#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zstd compression of the per-instance MITHrIL drug outputs
(<id>.perturbations.txt -> <id>.perturbations.txt.zst) of the LINCS cell-line
runs, for cell lines whose outputs have all been converted to CS inputs.

Design
------
- Cell lines: every conf.MITH_OUT/drug_signature_2025/<cell_id>_LM folder with
  an input matrix conf.MITH_IN_DRUG/LINCS_<cell_id>_LM.mi. A line is eligible
  only when every instance of its .mi has <id>.pkl and <id>_pathways.tsv in
  conf.CS_DIR/input/drug_signature_2025/<cell_id>_LM/ and a MITHrIL
  perturbation file (plain or .zst), and it is not the line MITHrIL is
  running (the last 'MITHrIL on' line of MITHRIL_LOG with no later
  'CS inputs' line for it).
- Per file: zstd -ZSTD_LEVEL to <file>.zst, integrity check (zstd -t), and
  only then removal of the plain file. A file that fails either step is kept
  and reported. N_WORKERS files are compressed at a time, one zstd thread each.
- Only *.perturbations.txt files are touched; *.output.txt and the CS inputs
  are left as they are. Folders without the '_LM' suffix are ignored.
- Loop: a pass over all eligible lines every INTERVAL_SEC seconds while the
  MITHrIL driver (run_drug_mithril_all_cell_lines.py) is running, then one
  last pass. Run under nice / ionice idle class (see usage).
- step7_drug_mith_out_to_cs_in and run_drug_mithril_all_cell_lines read the
  .zst files directly.

Inputs
------
- conf.MITH_OUT, conf.MITH_IN_DRUG, conf.CS_DIR, conf.LOGS_DIR
- MITHRIL_LOG: log of the running run_drug_mithril_all_cell_lines.py

Outputs
-------
- compressed files next to the originals (originals removed after the check)
- conf.LOGS_DIR/mithril_output_compression_<DD_MM_YYYY>.tsv: one row per
  (pass, cell line) with files compressed, bytes before / after, errors.

Usage
-----
    tmux new -d -s compress_mith "cd ~/netcos/src && PYTHONPATH=. nice -n 19 ionice -c3 python compress_mithril_drug_outputs.py"
"""

import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

from conf import CS_DIR, LOGS_DIR, MITH_IN_DRUG, MITH_OUT
from logger import append_run_metadata

ZSTD_LEVEL = 3
N_WORKERS = 2
INTERVAL_SEC = 1800
MITHRIL_LOG = LOGS_DIR / 'mithril_all_cell_lines_05_11_2026.log'
DRIVER = 'run_drug_mithril_all_cell_lines.py'

MITH_OUT_ROOT = MITH_OUT / 'drug_signature_2025'
CS_DRUG_ROOT = CS_DIR / 'input' / 'drug_signature_2025'
REGISTRY = LOGS_DIR / f'mithril_output_compression_{datetime.now():%d_%m_%Y}.tsv'
SUFFIX = '.perturbations.txt'


def running_line():
    """Run name MITHrIL is working on, from the driver log (None if none)."""
    if not MITHRIL_LOG.exists():
        return None
    current = None
    for line in open(MITHRIL_LOG, errors='replace'):
        m = re.search(r'\] (\S+_LM): MITHrIL on ', line)
        if m:
            current = m.group(1)
        m = re.search(r'\] (\S+_LM): CS inputs ', line)
        if m and m.group(1) == current:
            current = None
    return current


def driver_running():
    return subprocess.run(['pgrep', '-f', DRIVER], capture_output=True).returncode == 0


def eligible(run_name):
    mi = MITH_IN_DRUG / f'LINCS_{run_name}.mi'
    if not mi.exists():
        return False
    ids = open(mi).readline().rstrip('\n').split('\t')
    out_dir, cs_dir = MITH_OUT_ROOT / run_name, CS_DRUG_ROOT / run_name
    names_out = set(os.listdir(out_dir))
    names_cs = set(os.listdir(cs_dir)) if cs_dir.exists() else set()
    return all((f'{i}.pkl' in names_cs and f'{i}_pathways.tsv' in names_cs
                and (f'{i}{SUFFIX}' in names_out or f'{i}{SUFFIX}.zst' in names_out)) for i in ids)


def compress_one(path):
    zst = path.with_name(path.name + '.zst')
    before = path.stat().st_size
    r = subprocess.run(['zstd', f'-{ZSTD_LEVEL}', '-T1', '-q', '-f', str(path), '-o', str(zst)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return before, 0, f'{path.name}: compress failed: {r.stderr.strip()}'
    t = subprocess.run(['zstd', '-t', '-q', str(zst)], capture_output=True, text=True)
    if t.returncode != 0:
        return before, 0, f'{path.name}: test failed: {t.stderr.strip()}'
    after = zst.stat().st_size
    path.unlink()
    return before, after, None


def one_pass():
    skip = running_line()
    for out_dir in sorted(MITH_OUT_ROOT.glob('*_LM')):
        run_name = out_dir.name
        if run_name == skip or not eligible(run_name):
            continue
        files = sorted(out_dir.glob('*' + SUFFIX))
        if not files:
            continue
        start = time.time()
        with ThreadPoolExecutor(N_WORKERS) as ex:
            res = list(ex.map(compress_one, files))
        errors = [e for _, _, e in res if e]
        before = sum(b for b, a, e in res if not e)
        after = sum(a for b, a, e in res if not e)
        append_run_metadata(REGISTRY, {
            'timestamp': datetime.now().strftime('%d_%m_%Y_%H_%M_%S'), 'cell_line_run_name': run_name,
            'n_files_compressed': len(files) - len(errors), 'bytes_before': before, 'bytes_after': after,
            'ratio': round(before / after, 2) if after else None, 'elapsed_sec': round(time.time() - start, 1),
            'zstd_level': ZSTD_LEVEL, 'n_errors': len(errors), 'errors': ' | '.join(errors[:20])})
        print(f'[{datetime.now():%d/%m %H:%M}] {run_name}: {len(files) - len(errors)} files, '
              f'{before / 1e9:.1f} GB -> {after / 1e9:.1f} GB, {len(errors)} errors', flush=True)


if __name__ == '__main__':
    while True:
        alive = driver_running()
        one_pass()
        if not alive:
            print(f'[{datetime.now():%d/%m %H:%M}] MITHrIL driver not running: last pass done', flush=True)
            break
        time.sleep(INTERVAL_SEC)
