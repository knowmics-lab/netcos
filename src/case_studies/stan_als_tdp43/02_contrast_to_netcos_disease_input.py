#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Turn selected edgeR QLF contrasts of the TDP-43 iPSC motor-neuron RNA-seq
into the disease-side inputs of the pipeline (DEG CS input and MITHrIL
input), one disease run per contrast.

Design
------
- Source: the per-contrast logFC and FDR written by 01_edger_qlf_and_qc.py.
  Each contrast in conf.STAN_NETCOS_CONTRASTS becomes the disease run
  '<STAN_SLUG>_<contrast>'. logFC sign is numerator minus denominator, so a
  gene higher in the TDP-43 condition has DE_log2_FC > 0.
- Gene ids: human gene symbols (org.Hs.eg.db symbols of the counts file) are
  mapped to Entrez ids with conf.alias_2geneid, the dictionary used for the
  other disease signatures. Rows are dropped, with a status recorded, when:
  the row has no symbol (including the YFP transgene row); one ENSEMBL id
  carries several symbols (the counts file repeats the same counts under
  each symbol); the symbol is not in the dictionary; or several remaining
  rows map to the same Entrez id.
- The full tested gene list is written (no FDR filter); adj.p.value = edgeR
  BH FDR.

Inputs (paths from conf.py)
---------------------------
- STAN_OOUT_DIR/<STAN_SLUG>_edgeR_qlf_contrasts.tsv (from 01_).
- conf.alias_2geneid.

Outputs
-------
- CS_DIR/input/disease_signature_2025/<run>/<run>_signature_gene_id.csv
    ';'-separated, ',' decimal; columns gene, gene_id, DE_log2_FC,
    adj.p.value, ENSEMBL. DEG-level CS input, read by
    loader.load_disease_signature with DISEASE_SIGNATURE_SOURCE = "cs_in".
- MITH_IN_DISEASE/<run>_signature_gene_id.mi
    MITHrIL input, one 'gene_id<TAB>log2FC' line per gene, no header.
- STAN_OOUT_DIR/<STAN_SLUG>_gene_id_mapping_summary.tsv
    counts per mapping status (identical for every contrast), plus the
    number of genes written and of FDR genes (up/down) per disease run.

Usage
-----
    python 02_contrast_to_netcos_disease_input.py   (from any directory)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR))

from conf import (STAN_SLUG, STAN_OOUT_DIR, STAN_NETCOS_CONTRASTS, STAN_FDR,
                  CS_DIR, MITH_IN_DISEASE, alias_2geneid)


def map_rows(de):
    ens_multi = de.groupby('ENSEMBL')['SYMBOL'].transform('nunique') > 1
    gid = [str(alias_2geneid[s]) if isinstance(s, str) and s in alias_2geneid else None
           for s in de['SYMBOL']]
    de = de.assign(gene_id=pd.Series(gid, index=de.index, dtype=object))
    status = np.where(de['SYMBOL'].isna(), 'no_symbol',
             np.where(ens_multi, 'ensembl_with_several_symbols',
             np.where(de['gene_id'].isna(), 'symbol_not_in_alias_2geneid', 'mapped')))
    de['status'] = status
    m = de['status'] == 'mapped'
    dup = m & de['gene_id'].where(m).duplicated(keep=False) & de['gene_id'].notna()
    de.loc[dup, 'status'] = 'several_rows_same_gene_id'
    return de


if __name__ == '__main__':
    de = pd.read_csv(STAN_OOUT_DIR / f'{STAN_SLUG}_edgeR_qlf_contrasts.tsv', sep='\t')
    de = map_rows(de)
    summary = de['status'].value_counts().rename_axis('item').reset_index(name='n')
    ok = de[de['status'] == 'mapped']

    MITH_IN_DISEASE.mkdir(parents=True, exist_ok=True)
    extra = []
    for name in STAN_NETCOS_CONTRASTS:
        run = f'{STAN_SLUG}_{name}'
        out = pd.DataFrame({'gene': ok['SYMBOL'], 'gene_id': ok['gene_id'],
                            'DE_log2_FC': ok[f'{name}_logFC'], 'adj.p.value': ok[f'{name}_FDR'],
                            'ENSEMBL': ok['ENSEMBL']}).sort_values('gene_id')
        assert not out['gene_id'].duplicated().any()
        d = CS_DIR / 'input' / 'disease_signature_2025' / run
        d.mkdir(parents=True, exist_ok=True)
        out.to_csv(d / f'{run}_signature_gene_id.csv', sep=';', decimal=',', index=False)
        with open(MITH_IN_DISEASE / f'{run}_signature_gene_id.mi', 'w') as f:
            f.write('\n'.join(f'{g}\t{fc}' for g, fc in zip(out['gene_id'], out['DE_log2_FC'])))
        sig = out['adj.p.value'] <= STAN_FDR
        extra += [(f'{run}: genes written', len(out)),
                  (f'{run}: FDR<={STAN_FDR} up', int((sig & (out['DE_log2_FC'] > 0)).sum())),
                  (f'{run}: FDR<={STAN_FDR} down', int((sig & (out['DE_log2_FC'] < 0)).sum()))]
    summary = pd.concat([summary, pd.DataFrame(extra, columns=['item', 'n'])], ignore_index=True)
    summary.to_csv(STAN_OOUT_DIR / f'{STAN_SLUG}_gene_id_mapping_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))
