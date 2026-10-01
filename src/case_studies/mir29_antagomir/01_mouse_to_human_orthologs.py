#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map a mouse gene-level log2 fold-change signature to human Entrez gene ids
through 1:1 orthologs, and write the disease-side inputs of the pipeline.

Design
------
- Orthology: MGI homology report (HOM_MouseHumanSequence.rpt). A homology
  class (DB Class Key) is kept only if it contains exactly one mouse gene and
  exactly one human gene (1:1 orthologs). Many-to-many classes are dropped.
- Mouse genes are matched by MGI symbol; the human side is taken directly as
  the EntrezGene ID of the class (no symbol re-mapping).
- Mouse genes without a 1:1 human ortholog are dropped (e.g. most Gm*/Rik
  models, rodent-specific families).
- The input table carries no p-values: the adj.p.value column is written as
  NaN, and the full ranked list is used downstream.

Inputs (paths from conf.py)
---------------------------
- conf.MIR29_LFC_FILE : tab-separated, columns Gene (MGI symbol), LFC.
- conf.MGI_HOM_FILE   : MGI HOM_MouseHumanSequence.rpt.

Outputs
-------
- TSR_OUT_DISEASE/<DISEASE>_signature_gene_id.csv
    ';'-separated, ',' decimal; columns gene (human symbol), gene_id (human
    Entrez), DE_log2_FC, adj.p.value, mouse_gene. DEG-level CS input.
- MITH_IN_DISEASE/<disease_run_name>_signature_gene_id.mi
    MITHrIL input, one 'gene_id<TAB>log2FC' line per gene, no header.
- MIR29_OOUT_DIR/<DISEASE>_mouse_to_human_mapping.tsv
    one row per input mouse gene with its mapping status.
- MIR29_OOUT_DIR/<DISEASE>_mouse_to_human_mapping_summary.tsv
    counts per mapping status.

Usage
-----
    python 01_mouse_to_human_orthologs.py   (from any directory)
"""

import sys
from pathlib import Path

import pandas as pd
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR))

from conf import (DISEASE, disease_run_name, MIR29_LFC_FILE, MGI_HOM_FILE,
                  TSR_OUT_DISEASE, MITH_IN_DISEASE, MIR29_OOUT_DIR)

MOUSE = 'mouse, laboratory'
HUMAN = 'human'


def load_one_to_one_orthologs(path):
    hom = pd.read_csv(path, sep='\t', dtype=str)
    hom = hom[hom['Common Organism Name'].isin([MOUSE, HUMAN])]
    counts = hom.groupby(['DB Class Key', 'Common Organism Name']).size().unstack(fill_value=0)
    keep = counts[(counts.get(MOUSE, 0) == 1) & (counts.get(HUMAN, 0) == 1)].index
    hom = hom[hom['DB Class Key'].isin(keep)]
    mouse = hom[hom['Common Organism Name'] == MOUSE][['DB Class Key', 'Symbol']]
    human = hom[hom['Common Organism Name'] == HUMAN][['DB Class Key', 'Symbol', 'EntrezGene ID']]
    o2o = mouse.merge(human, on='DB Class Key', suffixes=('_mouse', '_human'))
    o2o = o2o.rename(columns={'Symbol_mouse': 'mouse_gene', 'Symbol_human': 'gene',
                              'EntrezGene ID': 'gene_id'})
    o2o = o2o.dropna(subset=['gene_id'])
    # a mouse symbol must identify one class, a human id one class
    o2o = o2o[~o2o['mouse_gene'].duplicated(keep=False) & ~o2o['gene_id'].duplicated(keep=False)]
    all_mouse = set(hom_all_mouse_symbols(path))
    return o2o[['mouse_gene', 'gene', 'gene_id']], all_mouse


def hom_all_mouse_symbols(path):
    hom = pd.read_csv(path, sep='\t', dtype=str, usecols=['Common Organism Name', 'Symbol'])
    return hom.loc[hom['Common Organism Name'] == MOUSE, 'Symbol']


def map_signature(sig, o2o, all_mouse):
    sig = sig.rename(columns={'Gene': 'mouse_gene', 'LFC': 'DE_log2_FC'}).copy()
    dup_mouse = sig['mouse_gene'].duplicated(keep=False)
    m = sig.merge(o2o, on='mouse_gene', how='left')
    status = np.where(dup_mouse, 'duplicated_mouse_symbol',
             np.where(m['gene_id'].notna(), 'mapped_1to1',
             np.where(m['mouse_gene'].isin(all_mouse), 'no_1to1_ortholog', 'symbol_not_in_MGI')))
    m['status'] = status
    return m


if __name__ == '__main__':
    sig = pd.read_csv(MIR29_LFC_FILE, sep='\t')
    o2o, all_mouse = load_one_to_one_orthologs(MGI_HOM_FILE)
    mapped = map_signature(sig, o2o, all_mouse)

    MIR29_OOUT_DIR.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(MIR29_OOUT_DIR / f'{DISEASE}_mouse_to_human_mapping.tsv', sep='\t', index=False)
    summary = mapped['status'].value_counts().rename_axis('status').reset_index(name='n_mouse_genes')
    summary.to_csv(MIR29_OOUT_DIR / f'{DISEASE}_mouse_to_human_mapping_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))

    out = mapped[mapped['status'] == 'mapped_1to1'].copy()
    out['adj.p.value'] = np.nan
    out = out[['gene', 'gene_id', 'DE_log2_FC', 'adj.p.value', 'mouse_gene']].sort_values('gene_id')
    assert not out['gene_id'].duplicated().any()

    TSR_OUT_DISEASE.mkdir(parents=True, exist_ok=True)
    out.to_csv(TSR_OUT_DISEASE / f'{DISEASE}_signature_gene_id.csv', sep=';', decimal=',', index=False)

    MITH_IN_DISEASE.mkdir(parents=True, exist_ok=True)
    with open(MITH_IN_DISEASE / f'{disease_run_name}_signature_gene_id.mi', 'w') as f:
        f.write('\n'.join(f'{g}\t{fc}' for g, fc in zip(out['gene_id'], out['DE_log2_FC'])))

    print(len(out), 'human genes written for', DISEASE)
