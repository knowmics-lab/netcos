#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transcriptomic similarity of the CCLE cell lines (LINCS cell lines among them)
to the GSE92592 lung tissue samples (IPF and controls), i.e. the per-cell-line
weight of the Chen 2017 / OCTAD sRGES.

Design
------
- Samples: log2(TPM + 1) of the NCBI-recomputed TPM of GSE92592 (rows = NCBI
  GeneID). Groups from the GEO series matrix ('disease state'): IPF and
  control. The data checks at the top verify that the series-matrix order
  equals the group coding of the disease signature
  (case_conf.SAMPLES_GROUPS_MAP, 0 = IPF, 1 = control) and that the TPM and
  raw-count columns are the same samples. TPM is used on both sides because
  CCLE is TPM and per-sample ranks across genes depend on gene-length
  normalisation.
- Cell lines: DepMap log2(TPM + 1) of every CCLE line with expression; LINCS
  lines are the 71 cell_id of the L1000 level-5 gold set
  (conf.LINCS_METADATA_PATH), matched to DepMap by cell_line_similarity.
- Marker genes: the three sets of cell_line_similarity.GENE_SETS, selected on
  the CCLE matrix only; 'binchen2017_iqr5000' is the primary set.
- Similarity of a line to a group: median over the group's samples of the
  per-sample Spearman correlations (Chen 2017 / OCTAD), for IPF, control and
  all samples; Spearman correlation with the median profile of each group.
- sRGES weight of a LINCS line: its median correlation with the IPF samples
  (median_cor_ipf) on the primary gene set.
- Among LINCS lines: OCTAD's one-sided paired Wilcoxon of the best line
  against every other line (all samples), BH FDR.
- Reference scale: the median tumour-cell line correlations of Chen 2017
  (BRCA, LIHC, COAD, ER+ BRCA; cell_line_<T>_tacle.csv, CCLE 2017 and TCGA,
  5,000 IQR genes) over all the CCLE lines they cover.

Inputs
------
- case_conf.NCBI_TPM, case_conf.SERIES_MATRIX, case_conf.NCBI_COUNTS (header only)
- cell_line_similarity.CCLE_EXPRESSION, CCLE_MODEL, LINCS_LANDMARK
- conf.LINCS_METADATA_PATH, conf.BC_DISEASE_DATA/<T>/cell_line_<T>_tacle.csv

Outputs
-------
- P.oout_dir/<slug>_samples.tsv
    GSM, title, disease_state, group (ipf / control), group_code (0 / 1).
- P.oout_dir/<slug>_ccle_cell_line_similarity.tsv
    one row per (gene_set, CCLE line): ModelID, CellLineName, OncotreeLineage,
    lincs_cell_id (empty if not in LINCS), n_lincs_instances, n_genes,
    median_cor_all / _ipf / _control, cor_median_profile_ipf / _control,
    rank_lincs (by median_cor_all, LINCS lines only), fdr_vs_best_lincs.
- P.oout_dir/<slug>_cell_line_weights.tsv
    sRGES weights: lincs_cell_id, cor = median_cor_ipf of the primary gene
    set; LINCS lines with DepMap expression only (lines without a weight are
    dropped by sRGES.compute_sRGES).
- P.oout_dir/<slug>_lincs_ccle_mapping.tsv
    the 71 LINCS cell_id with their DepMap match.
- P.img_dir/<slug>_ccle_cell_line_similarity_panel.png
    A: per-sample correlations of the LINCS lines (primary set), IPF and
    control samples side by side, ordered by median over all samples;
    B: per-line median correlation, all CCLE vs LINCS lines vs the Chen 2017
    tumour references; C: IPF vs control median per LINCS line; D: agreement
    of the LINCS-line ordering between gene sets.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/ipf/02_cell_line_weights.py
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P, NCBI_TPM, SERIES_MATRIX, NCBI_COUNTS, SAMPLES_GROUPS_MAP
from conf import LINCS_METADATA_PATH, BC_DISEASE_DATA
import cell_line_similarity as cls
import lfs_style

PRIMARY = 'binchen2017_iqr5000'
CASE = 'ipf'
GROUPS = {CASE: '0', 'control': '1'}
TACLE = ['BRCA', 'LIHC', 'COAD', 'ER']
BLUE, RED = lfs_style.BEFORE, lfs_style.AFTER


def read_series_matrix(path):
    """GSM, title and 'disease state' of every sample of a GEO series matrix."""
    rows = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_geo_accession'):
                rows['GSM'] = [v.strip('"') for v in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_title'):
                rows['title'] = [v.strip('"') for v in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1') and 'disease state' in line:
                rows['disease_state'] = [v.strip('"').split(': ', 1)[1]
                                         for v in line.rstrip('\n').split('\t')[1:]]
    return pd.DataFrame(rows)


def build_samples():
    s = read_series_matrix(SERIES_MATRIX)
    print('series matrix, disease state:', s['disease_state'].value_counts().to_dict())
    s['group'] = np.where(s['disease_state'] == 'Control', 'control', CASE)
    s['group_code'] = s['group'].map(GROUPS)
    coding = ''.join(s['group_code'])
    print(f'group coding from GEO == disease-signature coding: {coding == SAMPLES_GROUPS_MAP}')
    assert coding == SAMPLES_GROUPS_MAP
    counts_cols = pd.read_csv(NCBI_COUNTS, sep='\t', nrows=0).columns[1:].tolist()
    print(f'raw-count columns in the same order as the series matrix: {counts_cols == s["GSM"].tolist()}')
    assert counts_cols == s['GSM'].tolist()
    return s


def load_samples(s):
    x = pd.read_csv(NCBI_TPM, sep='\t', index_col=0)
    x.index = x.index.astype(int)
    print(f'TPM matrix: {x.shape[0]} genes x {x.shape[1]} samples; '
          f'column sums {x.sum().min():.0f}-{x.sum().max():.0f}')
    missing = set(s['GSM']) - set(x.columns)
    extra = set(x.columns) - set(s['GSM'])
    print(f'samples table {len(s)}; not in TPM {len(missing)}; TPM not in table {len(extra)}')
    assert not missing and not extra
    return np.log2(x[s['GSM']] + 1)


def load_tacle():
    rows = []
    for t in TACLE:
        d = pd.read_csv(BC_DISEASE_DATA / t / f'cell_line_{t}_tacle.csv')
        rows.append(pd.DataFrame({'tumour': t, 'cor': d['cor'].values}))
    return pd.concat(rows)


def per_gene_set(x, s, ccle, gs, ref_models, lincs_models, groups):
    genes = cls.marker_genes(ccle, gs, reference_models=ref_models)
    cor, n = cls.similarity(x, ccle, genes)
    out = pd.DataFrame(index=ccle.index)
    out['gene_set'] = gs
    out['n_genes'] = n
    out['median_cor_all'] = cor.median(axis=0)
    for g, code in groups.items():
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


def similarity_tables(x, s, groups, ccle, models, lincs):
    """Similarity table of every gene set, per-sample correlations of the primary set, LINCS mapping."""
    n_inst = lincs.groupby('cell_id').size()
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    mapping['n_lincs_instances'] = mapping['lincs_cell_id'].map(n_inst)
    lincs_models = mapping.loc[mapping.has_expression, 'ModelID'].tolist()
    meta = models.set_index('ModelID').loc[ccle.index, ['CellLineName', 'OncotreeLineage']]
    m2l = mapping[mapping.has_expression].set_index('ModelID')
    tabs, cor_primary = [], None
    for gs in cls.GENE_SETS:
        out, cor = per_gene_set(x, s, ccle, gs, lincs_models, lincs_models, groups)
        if gs == PRIMARY:
            cor_primary = cor
        out = meta.join(out)
        out['lincs_cell_id'] = m2l['lincs_cell_id'].reindex(out.index)
        out['n_lincs_instances'] = m2l['n_lincs_instances'].reindex(out.index)
        tabs.append(out)
    tab = pd.concat(tabs).rename_axis('ModelID').reset_index()
    g = list(groups)
    cols = (['gene_set', 'ModelID', 'CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'n_lincs_instances',
             'n_genes', 'median_cor_all'] + [f'median_cor_{k}' for k in g]
            + [f'cor_median_profile_{k}' for k in g] + ['rank_lincs', 'fdr_vs_best_lincs'])
    tab = tab[cols].sort_values(['gene_set', 'median_cor_all'], ascending=[True, False])
    return tab, cor_primary, mapping


def panel(tab, cor_primary, tacle, s):
    p = tab[tab['gene_set'] == PRIMARY]
    lin = p[p['lincs_cell_id'].notna()].sort_values('median_cor_all', ascending=False)
    fig, ax = plt.subplots(2, 2, figsize=(16, 12))

    a = ax[0, 0]
    pos = np.arange(1, len(lin) + 1)
    for g, col, dx in ((CASE, RED, -0.2), ('control', BLUE, 0.2)):
        idx = s.loc[s['group_code'] == GROUPS[g], 'GSM']
        data = [cor_primary.loc[idx, m].values for m in lin['ModelID']]
        a.boxplot(data, positions=pos + dx, widths=.35, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=col, alpha=.4, edgecolor=col),
                  whiskerprops=dict(color=col), capprops=dict(color=col),
                  medianprops=dict(color='k', linewidth=1), manage_ticks=False)
        a.plot([], [], 's', color=col, alpha=.4, label=f'{"IPF" if g == CASE else g} samples (n={len(idx)})')
    a.set_xticks(pos)
    a.set_xlim(0.3, len(lin) + 0.7)
    a.set_xticklabels([f'{c} ({l})' for c, l in zip(lin['lincs_cell_id'], lin['OncotreeLineage'])],
                      rotation=90, fontsize=7)
    a.set_ylabel('Spearman correlation with GSE92592 samples')
    a.set_title(f'A. LINCS lines, per-sample correlations ({PRIMARY}, {int(p.n_genes.iloc[0])} genes)')
    a.legend(loc='lower left')

    a = ax[0, 1]
    bins = np.linspace(min(p['median_cor_all'].min(), tacle['cor'].min()) - .02,
                       max(p['median_cor_all'].max(), tacle['cor'].max()) + .02, 50)
    a.hist(p['median_cor_all'], bins=bins, color='0.6', alpha=.7, label=f'IPF, all CCLE lines (n={len(p)})')
    a.hist(lin['median_cor_all'], bins=bins, color=RED, alpha=.6, label=f'IPF, LINCS lines (n={len(lin)})')
    for t, ls in zip(TACLE, ['-', '--', ':', '-.']):
        v = tacle.loc[tacle['tumour'] == t, 'cor']
        a.hist(v, bins=bins, histtype='step', color=BLUE, ls=ls, lw=1.4,
               label=f'Chen 2017 {t}, all CCLE lines (n={len(v)})')
    a.set_xlabel('median Spearman correlation per cell line')
    a.set_ylabel('cell lines')
    a.set_title('B. Similarity scale: IPF vs Chen 2017 tumour references')
    a.legend(fontsize=8)

    a = ax[1, 0]
    a.scatter(lin['median_cor_control'], lin[f'median_cor_{CASE}'], color=RED, s=18)
    for _, r in lin.head(8).iterrows():
        a.annotate(r['lincs_cell_id'], (r['median_cor_control'], r[f'median_cor_{CASE}']), fontsize=7)
    v = lin[['median_cor_control', f'median_cor_{CASE}']].values
    a.plot([v.min(), v.max()], [v.min(), v.max()], color='k', lw=.6)
    a.set_xlabel('median correlation, controls')
    a.set_ylabel('median correlation, IPF')
    a.set_title('C. LINCS lines: IPF vs control samples')

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
    lfs_style.apply()
    s = build_samples()
    x = load_samples(s)
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    print(f'CCLE matrix: {ccle.shape[0]} lines x {ccle.shape[1]} genes; '
          f'genes shared with GSE92592: {len(set(ccle.columns) & set(x.index))}')
    lincs = pd.read_csv(LINCS_METADATA_PATH)

    tab, cor_primary, mapping = similarity_tables(x, s, GROUPS, ccle, models, lincs)
    print(f'LINCS lines {len(mapping)}; with DepMap expression {int(mapping.has_expression.sum())}; '
          f'without: {", ".join(mapping.loc[~mapping.has_expression, "lincs_cell_id"])}')

    P.oout_dir.mkdir(parents=True, exist_ok=True)
    s.to_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', index=False)
    mapping.to_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t', index=False)
    tab.to_csv(P.oout_dir / f'{P.slug}_ccle_cell_line_similarity.tsv', sep='\t', index=False)
    w = tab[(tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
    w = w[['lincs_cell_id', f'median_cor_{CASE}']].rename(columns={f'median_cor_{CASE}': 'cor'})
    w.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t', index=False)
    print(f'sRGES weights written for {len(w)} LINCS lines (range {w.cor.min():.3f}-{w.cor.max():.3f})')

    tacle = load_tacle()
    panel(tab, cor_primary, tacle, s)

    p = tab[tab.gene_set == PRIMARY]
    lin = p[p.lincs_cell_id.notna()]
    d = lin[f'median_cor_{CASE}'] - lin['median_cor_control']
    print(f'\nIPF - control median correlation per LINCS line: min {d.min():+.3f}, max {d.max():+.3f}; '
          f'IPF higher in {(d > 0).sum()}/{len(d)} lines')
    print('\nprimary set, median per-line correlation (all samples):')
    print(f'  all CCLE lines: median {p.median_cor_all.median():.3f}, max {p.median_cor_all.max():.3f} '
          f'({p.loc[p.median_cor_all.idxmax(), "CellLineName"]})')
    print(f'  LINCS lines:    median {lin.median_cor_all.median():.3f}, max {lin.median_cor_all.max():.3f}')
    for t in TACLE:
        v = tacle.loc[tacle.tumour == t, 'cor']
        print(f'  Chen 2017 {t}: median {v.median():.3f}, max {v.max():.3f}')
    for gs in cls.GENE_SETS:
        q = tab[(tab.gene_set == gs) & tab.lincs_cell_id.notna()]
        print(f'  {gs}: LINCS median_cor_{CASE} range {q[f"median_cor_{CASE}"].min():.3f}-'
              f'{q[f"median_cor_{CASE}"].max():.3f}')
    print('\ntop 10 LINCS lines (primary set):')
    print(lin.head(10)[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'n_lincs_instances',
                        'median_cor_all', f'median_cor_{CASE}', 'median_cor_control',
                        'fdr_vs_best_lincs']].to_string(index=False))
    print('\nbottom 5 LINCS lines (primary set):')
    print(lin.tail(5)[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'median_cor_all']].to_string(index=False))
    print('\ntop 10 CCLE lines overall (primary set):')
    print(p.head(10)[['CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'median_cor_all']].to_string(index=False))
