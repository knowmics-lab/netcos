#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Inspect the GSE152004 download (asthmatic vs healthy-control nasal airway
epithelium RNA-seq) and write the sample-group map used by the TSR disease
signature workflow.

Design
------
- Files: name and size of every file in P.data_dir.
- Count matrices: for the NCBI-recomputed counts (case_conf.NCBI_COUNTS,
  rows = NCBI GeneID, columns = GSM) and the submitter counts
  (case_conf.AUTHOR_COUNTS, rows = gene symbols, header without a name for
  the gene column): number of genes and samples, whether every value is a
  non-negative integer, number of all-zero genes, first gene ids.
- Sample metadata from the series matrix: '!Sample_geo_accession',
  '!Sample_title' and every '!Sample_characteristics_ch1' field, split on
  the first ':' into covariate name and level.
- Alignment: every GSM column of the NCBI counts must appear exactly once in
  the series matrix, and every submitter count column must match a
  '!Sample_title'; the script stops otherwise.
- Group map: one character per GSM column of the NCBI counts, in column
  order, '0' = asthmatic (disease) and '1' = healthy control, the coding of
  modules/disease_signature/DiseaseControlGroupsFactor.R.

Inputs (paths from case_conf.py)
--------------------------------
- P.data_dir (= DATA_DIR / asthma_GSE152004), NCBI_COUNTS, AUTHOR_COUNTS,
  SERIES_MATRIX.

Outputs
-------
- P.oout_dir/<slug>_download_inspection.tsv
    one row per check: section, item, value.
- P.oout_dir/<slug>_samples.tsv
    one row per GSM column of the NCBI counts: column index, GSM, title,
    every characteristic, group code.
- P.oout_dir/<slug>_samples_groups_map.txt
    the group-map string (single line), read by
    tsr/asthma_GSE152004/config/asthma_GSE152004_cfg.R.

Usage
-----
    python 00_inspect_geo_download.py   (from any directory)
"""

import sys
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P, NCBI_COUNTS, AUTHOR_COUNTS, SERIES_MATRIX


GROUP_CODE = {'asthmatic': '0', 'healthy control': '1'}
rows = []


def log(section, item, value):
    rows.append((section, item, value))
    print(f'[{section}] {item}: {value}')


def inspect_counts(path, label, index_col):
    m = pd.read_csv(path, sep='\t', index_col=index_col)
    v = m.to_numpy()
    log(label, 'file', path.name)
    log(label, 'genes x samples', f'{m.shape[0]} x {m.shape[1]}')
    log(label, 'all values non-negative integers',
        bool(np.issubdtype(v.dtype, np.integer) or (np.all(v >= 0) and np.all(v == np.round(v)))))
    log(label, 'genes with all-zero counts', int((v.sum(axis=1) == 0).sum()))
    log(label, 'duplicated gene ids', int(m.index.duplicated().sum()))
    log(label, 'first gene ids', ', '.join(map(str, m.index[:5])))
    log(label, 'first sample columns', ', '.join(m.columns[:3]))
    log(label, 'library size min / median / max',
        ' / '.join(f'{x:.0f}' for x in np.percentile(v.sum(axis=0), [0, 50, 100])))
    return m


def read_series_matrix():
    fields = {}
    chars = []
    with gzip.open(SERIES_MATRIX, 'rt') as f:
        for line in f:
            p = [x.strip().strip('"') for x in line.rstrip('\n').split('\t')]
            if p[0] in ('!Sample_geo_accession', '!Sample_title'):
                fields[p[0]] = p[1:]
            elif p[0] == '!Sample_characteristics_ch1':
                chars.append(p[1:])
    s = pd.DataFrame({'GSM': fields['!Sample_geo_accession'], 'title': fields['!Sample_title']})
    for c in chars:
        name = c[0].split(':', 1)[0].strip()
        s[name] = [x.split(':', 1)[1].strip() if ':' in x else x for x in c]
    log('series matrix', 'samples', len(s))
    for col in s.columns[2:]:
        lv = s[col].value_counts()
        log('series matrix', f'covariate "{col}"',
            '; '.join(f'{k} ({n})' for k, n in lv.head(10).items()) + (' ...' if len(lv) > 10 else ''))
    assert not s['GSM'].duplicated().any() and not s['title'].duplicated().any()
    return s


if __name__ == '__main__':
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    for p in sorted(P.data_dir.iterdir()):
        log('files', p.name, f'{p.stat().st_size / 2**20:.1f} MB')

    ncbi = inspect_counts(NCBI_COUNTS, 'NCBI counts', 'GeneID')
    author = inspect_counts(AUTHOR_COUNTS, 'submitter counts', 0)
    s = read_series_matrix()

    miss_gsm = sorted(set(ncbi.columns) - set(s['GSM']))
    miss_title = sorted(set(author.columns) - set(s['title']))
    log('alignment', 'NCBI GSM columns not in series matrix', len(miss_gsm))
    log('alignment', 'submitter columns not matching a title', len(miss_title))
    assert not miss_gsm and not miss_title

    s = s.set_index('GSM').loc[list(ncbi.columns)].rename_axis('GSM').reset_index()
    s.insert(0, 'column_index', range(1, len(s) + 1))
    assert set(s['asthma status']) == set(GROUP_CODE), set(s['asthma status'])
    s['group_code'] = s['asthma status'].map(GROUP_CODE)
    groups_map = ''.join(s['group_code'])
    log('group map', 'length', len(groups_map))
    log('group map', "'0' asthmatic / '1' healthy control",
        f"{groups_map.count('0')} / {groups_map.count('1')}")

    # same samples in both count files: library sizes, linked through the title
    lib_ncbi = ncbi.sum(axis=0).rename('lib_ncbi')
    lib_author = author.sum(axis=0).rename('lib_author')
    j = s[['GSM', 'title']].join(lib_ncbi, on='GSM').join(lib_author, on='title')
    log('alignment', 'Spearman of library sizes NCBI vs submitter (same sample)',
        round(j['lib_ncbi'].corr(j['lib_author'], method='spearman'), 4))

    pd.DataFrame(rows, columns=['section', 'item', 'value']).to_csv(
        P.oout_dir / f'{P.slug}_download_inspection.tsv', sep='\t', index=False)
    s.to_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', index=False)
    (P.oout_dir / f'{P.slug}_samples_groups_map.txt').write_text(groups_map + '\n')
