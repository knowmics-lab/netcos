#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transcriptomic similarity of the CCLE cell lines (LINCS cell lines among them)
to the GSE152004 nasal airway epithelium samples, i.e. the per-cell-line
weight of the Chen 2017 / OCTAD sRGES.

Design
------
- Samples: log2(TPM + 1) of the NCBI-recomputed TPM of GSE152004 (rows = NCBI
  GeneID), groups from the '<slug>_samples.tsv' table of 00_ (asthmatic,
  healthy control). TPM is used on both sides because CCLE is TPM and
  per-sample ranks across genes depend on gene-length normalisation.
- Cell lines: DepMap log2(TPM + 1) of every CCLE line with expression; LINCS
  lines are the 71 cell_id of the L1000 level-5 gold set
  (conf.LINCS_METADATA_PATH), matched to DepMap by cell_line_similarity.
- Marker genes: the three sets of cell_line_similarity.GENE_SETS, selected on
  the CCLE matrix only; 'binchen2017_iqr5000' is the primary set.
- Similarity of a line to a group: median over the group's samples of the
  per-sample Spearman correlations (Chen 2017 / OCTAD), for asthmatic,
  healthy control and all samples; Spearman correlation with the median
  profile of each group (the healthy-control one is the sRGES weight w).
- Among LINCS lines: OCTAD's one-sided paired Wilcoxon of the best line
  against every other line (all samples), BH FDR.
- Reference scale: the median tumour-cell line correlations of Chen 2017
  (BRCA, LIHC, COAD, ER+ BRCA; cell_line_<T>_tacle.csv, CCLE 2017 and TCGA,
  5,000 IQR genes) over all the CCLE lines they cover.

The data checks at the top print the sample alignment, the gene overlap of
the two matrices and the LINCS -> DepMap mapping.

Inputs
------
- case_conf.NCBI_TPM, P.oout_dir/<slug>_samples.tsv
- cell_line_similarity.CCLE_EXPRESSION, CCLE_MODEL, LINCS_LANDMARK
- conf.LINCS_METADATA_PATH, conf.BC_DISEASE_DATA/<T>/cell_line_<T>_tacle.csv

Outputs
-------
- P.oout_dir/<slug>_ccle_cell_line_similarity.tsv
    one row per (gene_set, CCLE line): ModelID, CellLineName, OncotreeLineage,
    lincs_cell_id (empty if not in LINCS), n_lincs_instances, n_genes,
    median_cor_all / _asthmatic / _control, cor_median_profile_asthmatic /
    _control (= w), rank_lincs (by median_cor_all, LINCS lines only),
    fdr_vs_best_lincs.
- P.oout_dir/<slug>_cell_line_weights.tsv
    sRGES weights: lincs_cell_id, cor = median_cor_asthmatic of the primary
    gene set (median over the asthmatic samples of the per-sample Spearman
    correlations, as in Chen 2017); LINCS lines with DepMap expression only,
    lines without a weight are dropped by sRGES.compute_sRGES.
- P.oout_dir/<slug>_lincs_ccle_mapping.tsv
    the 71 LINCS cell_id with their DepMap match.
- P.img_dir/<slug>_ccle_cell_line_similarity_panel.png
    A: per-sample correlations of the LINCS lines (primary set), asthmatic
    and healthy-control samples side by side, ordered by median over all
    samples; B: per-line median correlation, all CCLE vs LINCS lines vs the
    Chen 2017 tumour references; C: asthmatic vs healthy-control median per
    LINCS line; D: agreement of the LINCS-line ordering between gene sets.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/asthma_GSE152004/02_cell_line_weights.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P, NCBI_TPM
from conf import LINCS_METADATA_PATH, BC_DISEASE_DATA
import cell_line_similarity as cls

PRIMARY = 'binchen2017_iqr5000'
GROUPS = {'asthmatic': '0', 'control': '1'}
TACLE = ['BRCA', 'LIHC', 'COAD', 'ER']
BLUE, RED = '#1f4e9c', '#c0392b'

plt.rcParams.update({'axes.grid': True, 'grid.color': '0.7', 'grid.linestyle': '--',
                     'grid.linewidth': 0.6, 'font.size': 10})


def load_samples():
    s = pd.read_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', dtype={'group_code': str})
    x = pd.read_csv(NCBI_TPM, sep='\t', index_col=0)
    x.index = x.index.astype(int)
    print(f'TPM matrix: {x.shape[0]} genes x {x.shape[1]} samples')
    missing = set(s['GSM']) - set(x.columns)
    extra = set(x.columns) - set(s['GSM'])
    print(f'samples table {len(s)}; not in TPM {len(missing)}; TPM not in table {len(extra)}')
    assert not missing and not extra
    x = np.log2(x[s['GSM']] + 1)
    return s, x


def load_tacle():
    rows = []
    for t in TACLE:
        d = pd.read_csv(BC_DISEASE_DATA / t / f'cell_line_{t}_tacle.csv')
        rows.append(pd.DataFrame({'tumour': t, 'cor': d['cor'].values}))
    return pd.concat(rows)


def per_gene_set(x, s, ccle, gs, ref_models, lincs_models):
    genes = cls.marker_genes(ccle, gs, reference_models=ref_models)
    cor, n = cls.similarity(x, ccle, genes)
    out = pd.DataFrame(index=ccle.index)
    out['gene_set'] = gs
    out['n_genes'] = n
    out['median_cor_all'] = cor.median(axis=0)
    for g, code in GROUPS.items():
        idx = s.loc[s['group_code'] == code, 'GSM']
        out[f'median_cor_{g}'] = cor.loc[idx].median(axis=0)
        prof = x.loc[:, idx].median(axis=1).to_frame(f'median_{g}')
        out[f'cor_median_profile_{g}'] = cls.similarity(prof, ccle, genes)[0].iloc[0]
    lin = cor[lincs_models]
    best, fdr = cls.compare_with_best(lin)
    out['rank_lincs'] = lin.median(axis=0).rank(ascending=False)
    out['fdr_vs_best_lincs'] = fdr
    print(f'{gs}: {n} genes; best LINCS line {best}')
    return out, cor


def panel(tab, cor_primary, mapping, tacle, s):
    p = tab[tab['gene_set'] == PRIMARY]
    lin = p[p['lincs_cell_id'].notna()].sort_values('median_cor_all', ascending=False)
    fig, ax = plt.subplots(2, 2, figsize=(16, 12))

    a = ax[0, 0]
    pos = np.arange(1, len(lin) + 1)
    for g, col, dx in (('asthmatic', RED, -0.2), ('control', BLUE, 0.2)):
        idx = s.loc[s['group_code'] == GROUPS[g], 'GSM']
        data = [cor_primary.loc[idx, m].values for m in lin['ModelID']]
        a.boxplot(data, positions=pos + dx, widths=.35, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=col, alpha=.5, edgecolor=col),
                  whiskerprops=dict(color=col), capprops=dict(color=col),
                  medianprops=dict(color='k', linewidth=1), manage_ticks=False)
        a.plot([], [], 's', color=col, alpha=.5, label=f'{g} samples (n={len(idx)})')
    a.set_xticks(pos)
    a.set_xlim(0.3, len(lin) + 0.7)
    a.set_xticklabels([f'{c} ({l})' for c, l in zip(lin['lincs_cell_id'], lin['OncotreeLineage'])],
                      rotation=90, fontsize=7)
    a.set_ylabel('Spearman correlation with GSE152004 samples')
    a.set_title(f'A. LINCS lines, per-sample correlations ({PRIMARY}, {int(p.n_genes.iloc[0])} genes)')
    a.legend(loc='lower left')

    a = ax[0, 1]
    bins = np.linspace(min(p['median_cor_all'].min(), tacle['cor'].min()) - .02,
                       max(p['median_cor_all'].max(), tacle['cor'].max()) + .02, 50)
    a.hist(p['median_cor_all'], bins=bins, color='0.6', alpha=.7, label=f'asthma, all CCLE lines (n={len(p)})')
    a.hist(lin['median_cor_all'], bins=bins, color=RED, alpha=.8, label=f'asthma, LINCS lines (n={len(lin)})')
    for t, ls in zip(TACLE, ['-', '--', ':', '-.']):
        v = tacle.loc[tacle['tumour'] == t, 'cor']
        a.hist(v, bins=bins, histtype='step', color=BLUE, ls=ls, lw=1.4,
               label=f'Chen 2017 {t}, all CCLE lines (n={len(v)})')
    a.set_xlabel('median Spearman correlation per cell line')
    a.set_ylabel('cell lines')
    a.set_title('B. Similarity scale: asthma vs Chen 2017 tumour references')
    a.legend(fontsize=8)

    a = ax[1, 0]
    a.scatter(lin['median_cor_control'], lin['median_cor_asthmatic'], color=RED, s=18)
    for _, r in lin.head(8).iterrows():
        a.annotate(r['lincs_cell_id'], (r['median_cor_control'], r['median_cor_asthmatic']), fontsize=7)
    lo, hi = lin[['median_cor_control', 'median_cor_asthmatic']].values.min(), lin[['median_cor_control', 'median_cor_asthmatic']].values.max()
    a.plot([lo, hi], [lo, hi], color='k', lw=.6)
    a.set_xlabel('median correlation, healthy controls')
    a.set_ylabel('median correlation, asthmatics')
    a.set_title('C. LINCS lines: asthmatic vs healthy-control samples')

    a = ax[1, 1]
    piv = tab[tab['lincs_cell_id'].notna()].pivot_table(index='lincs_cell_id', columns='gene_set',
                                                         values='median_cor_all')
    for gs, c in zip([g for g in cls.GENE_SETS if g != PRIMARY], [BLUE, RED]):
        rho = spearmanr(piv[PRIMARY], piv[gs]).statistic
        a.scatter(piv[PRIMARY], piv[gs], s=16, color=c, label=f'{gs} (Spearman of orderings {rho:.2f})')
    a.set_xlabel(f'median correlation, {PRIMARY}')
    a.set_ylabel('median correlation, other gene set')
    a.set_title('D. LINCS lines: agreement between marker-gene sets')
    a.legend(fontsize=8)

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_ccle_cell_line_similarity_panel.png', dpi=130)


if __name__ == '__main__':
    s, x = load_samples()
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    print(f'CCLE matrix: {ccle.shape[0]} lines x {ccle.shape[1]} genes; '
          f'genes shared with GSE152004: {len(set(ccle.columns) & set(x.index))}')

    lincs = pd.read_csv(LINCS_METADATA_PATH)
    n_inst = lincs.groupby('cell_id').size()
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    mapping['n_lincs_instances'] = mapping['lincs_cell_id'].map(n_inst)
    print(f'LINCS lines {len(mapping)}; with DepMap expression {int(mapping.has_expression.sum())}; '
          f'without: {", ".join(mapping.loc[~mapping.has_expression, "lincs_cell_id"])}')
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t', index=False)

    lincs_models = mapping.loc[mapping.has_expression, 'ModelID'].tolist()
    meta = models.set_index('ModelID').loc[ccle.index, ['CellLineName', 'OncotreeLineage']]
    m2l = mapping[mapping.has_expression].set_index('ModelID')

    tabs = []
    for gs in cls.GENE_SETS:
        out, cor = per_gene_set(x, s, ccle, gs, lincs_models, lincs_models)
        if gs == PRIMARY:
            cor_primary = cor
        out = meta.join(out)
        out['lincs_cell_id'] = m2l['lincs_cell_id'].reindex(out.index)
        out['n_lincs_instances'] = m2l['n_lincs_instances'].reindex(out.index)
        tabs.append(out)
    tab = pd.concat(tabs).rename_axis('ModelID').reset_index()
    cols = ['gene_set', 'ModelID', 'CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'n_lincs_instances',
            'n_genes', 'median_cor_all', 'median_cor_asthmatic', 'median_cor_control',
            'cor_median_profile_asthmatic', 'cor_median_profile_control', 'rank_lincs', 'fdr_vs_best_lincs']
    tab = tab[cols].sort_values(['gene_set', 'median_cor_all'], ascending=[True, False])
    tab.to_csv(P.oout_dir / f'{P.slug}_ccle_cell_line_similarity.tsv', sep='\t', index=False)

    w = tab[(tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
    w = w[['lincs_cell_id', 'median_cor_asthmatic']].rename(columns={'median_cor_asthmatic': 'cor'})
    w.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t', index=False)
    print(f'sRGES weights written for {len(w)} LINCS lines')

    tacle = load_tacle()
    panel(tab, cor_primary, mapping, tacle, s)

    p = tab[tab.gene_set == PRIMARY]
    lin = p[p.lincs_cell_id.notna()]
    print('\nprimary set, median per-line correlation (all samples):')
    print(f'  all CCLE lines: median {p.median_cor_all.median():.3f}, max {p.median_cor_all.max():.3f} '
          f'({p.loc[p.median_cor_all.idxmax(), "CellLineName"]})')
    print(f'  LINCS lines:    median {lin.median_cor_all.median():.3f}, max {lin.median_cor_all.max():.3f}')
    for t in TACLE:
        v = tacle.loc[tacle.tumour == t, 'cor']
        print(f'  Chen 2017 {t}: median {v.median():.3f}, max {v.max():.3f}')
    print('\ntop 10 LINCS lines (primary set):')
    print(lin.head(10)[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'n_lincs_instances',
                        'median_cor_all', 'median_cor_asthmatic', 'median_cor_control',
                        'cor_median_profile_control', 'fdr_vs_best_lincs']].to_string(index=False))
    print('\ntop 10 CCLE lines overall (primary set):')
    print(p.head(10)[['CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'median_cor_all']].to_string(index=False))
