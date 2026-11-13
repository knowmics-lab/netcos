#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transcriptomic similarity of the CCLE cell lines (LINCS cell lines among them)
to the GSE153960 postmortem CNS samples behind the als_NYGC signature, i.e.
the per-cell-line weight of the Chen 2017 / OCTAD sRGES.

Design
------
- Samples: log2(TPM + 1) of the NCBI-recomputed TPM of GSE153960 (rows = NCBI
  GeneID). Groups from the GSE153960 metadata column tissue_status, as in the
  als_NYGC disease signature: case_conf.CASE_STATUS ('ALS Spectrum MND',
  group 'als') and case_conf.CONTROL_STATUS ('Non-Neurological Control',
  group 'control'); every other status is left out. The
  samples cover 12 CNS regions (metadata column tissue), pooled. TPM is used
  on both sides because CCLE is TPM and per-sample ranks across genes depend
  on gene-length normalisation.
- Cell lines: DepMap log2(TPM + 1) of every CCLE line with expression; LINCS
  lines are the 71 cell_id of the L1000 level-5 gold set
  (conf.LINCS_METADATA_PATH), matched to DepMap by cell_line_similarity.
- Marker genes: the three sets of cell_line_similarity.GENE_SETS, selected on
  the CCLE matrix only; 'binchen2017_iqr5000' is the primary set.
- Similarity of a line to a group: median over the group's samples of the
  per-sample Spearman correlations (Chen 2017 / OCTAD), for als, control and
  all samples; Spearman correlation with the median profile of each group.
  On the primary set, also the median per CNS region class (cortex, spinal
  cord, cerebellum, hippocampus; als and control samples together).
- sRGES weight of a LINCS line: its median correlation with the als samples
  (median_cor_als) on the primary gene set. LINCS lines without DepMap
  expression have no weight.
- Among LINCS lines: OCTAD's one-sided paired Wilcoxon of the best line
  against every other line (all samples), BH FDR.
- Reference scale: the median tumour-cell line correlations of Chen 2017
  (BRCA, LIHC, COAD, ER+ BRCA; cell_line_<T>_tacle.csv, CCLE 2017 and TCGA,
  5,000 IQR genes) over all the CCLE lines they cover.

The data checks at the top print the sample alignment, the group and region
counts, the gene overlap of the two matrices and the LINCS -> DepMap mapping.

Inputs
------
- case_conf.NCBI_TPM, case_conf.METADATA
- cell_line_similarity.CCLE_EXPRESSION, CCLE_MODEL, LINCS_LANDMARK
- conf.LINCS_METADATA_PATH, conf.BC_DISEASE_DATA/<T>/cell_line_<T>_tacle.csv

Outputs
-------
- P.oout_dir/<slug>_samples.tsv
    the als and control samples used: GSM, tissue_status, group, tissue,
    region_class, subject_id.
- P.oout_dir/<slug>_ccle_cell_line_similarity.tsv
    one row per (gene_set, CCLE line): ModelID, CellLineName, OncotreeLineage,
    lincs_cell_id (empty if not in LINCS), n_lincs_instances, n_genes,
    median_cor_all / _als / _control, cor_median_profile_als / _control,
    rank_lincs (by median_cor_all, LINCS lines only), fdr_vs_best_lincs;
    primary set only: median_cor_<region_class>.
- P.oout_dir/<slug>_cell_line_weights.tsv
    sRGES weights: lincs_cell_id, cor = median_cor_als of the primary gene set;
    LINCS lines with DepMap expression only, lines without a weight are
    dropped by sRGES.compute_sRGES.
- P.oout_dir/<slug>_lincs_ccle_mapping.tsv
    the 71 LINCS cell_id with their DepMap match and n_lincs_instances.
- P.img_dir/<slug>_ccle_cell_line_similarity_panel.png
    A: per-sample correlations of the LINCS lines (primary set), als and
    control samples side by side, ordered by median over all samples;
    B: per-line median correlation, all CCLE vs LINCS lines vs the Chen 2017
    tumour references; C: als vs control median per LINCS line; D: agreement
    of the LINCS-line ordering between gene sets; E: median correlation of the
    LINCS lines per CNS region class.

Usage
-----
    cd src && PYTHONPATH=. python case_studies/als_NYGC/02_cell_line_weights.py
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

from case_conf import P, NCBI_TPM, METADATA, CASE_STATUS, CONTROL_STATUS
from conf import LINCS_METADATA_PATH, BC_DISEASE_DATA
import cell_line_similarity as cls
import lfs_style

PRIMARY = 'binchen2017_iqr5000'
GROUPS = {'als': CASE_STATUS, 'control': CONTROL_STATUS}
REGION_CLASSES = ['cortex', 'spinal_cord', 'cerebellum', 'hippocampus']
TACLE = ['BRCA', 'LIHC', 'COAD', 'ER']
BLUE, RED = lfs_style.BEFORE, lfs_style.AFTER

lfs_style.apply(**{'savefig.dpi': 130, 'font.size': 10, 'legend.shadow': False})


def region_class(tissue):
    t = str(tissue).lower()
    if t.startswith('cortex'):
        return 'cortex'
    if t.startswith('spinal cord'):
        return 'spinal_cord'
    return t


def load_samples():
    m = pd.read_csv(METADATA, sep=';', encoding='utf-8-sig')
    print(f'metadata: {len(m)} samples; tissue_status: {m.tissue_status.value_counts().to_dict()}')
    s = m[m['tissue_status'].isin(GROUPS.values())].copy()
    s['group'] = s['tissue_status'].map({v: k for k, v in GROUPS.items()})
    s['region_class'] = s['tissue'].map(region_class)
    s = s.rename(columns={'accession': 'GSM'})[['GSM', 'tissue_status', 'group', 'tissue', 'region_class', 'subject_id']]
    print('samples by group and region class:')
    print(pd.crosstab(s['region_class'], s['group'], margins=True).to_string())
    print(f'subjects: als {s[s.group == "als"].subject_id.nunique()}, '
          f'control {s[s.group == "control"].subject_id.nunique()}')
    x = pd.read_csv(NCBI_TPM, sep='\t', index_col=0)
    x.index = x.index.astype(int)
    print(f'TPM matrix: {x.shape[0]} genes x {x.shape[1]} samples')
    missing = set(s['GSM']) - set(x.columns)
    print(f'samples used {len(s)}; not in TPM {len(missing)}')
    assert not missing
    x = np.log2(x[s['GSM']] + 1)
    return s.reset_index(drop=True), x


def load_tacle():
    rows = []
    for t in TACLE:
        d = pd.read_csv(BC_DISEASE_DATA / t / f'cell_line_{t}_tacle.csv')
        rows.append(pd.DataFrame({'tumour': t, 'cor': d['cor'].values}))
    return pd.concat(rows)


def per_gene_set(x, s, ccle, gs, lincs_models):
    genes = cls.marker_genes(ccle, gs, reference_models=lincs_models)
    cor, n = cls.similarity(x, ccle, genes)
    out = pd.DataFrame(index=ccle.index)
    out['gene_set'] = gs
    out['n_genes'] = n
    out['median_cor_all'] = cor.median(axis=0)
    for g in GROUPS:
        idx = s.loc[s['group'] == g, 'GSM']
        out[f'median_cor_{g}'] = cor.loc[idx].median(axis=0)
        prof = x.loc[:, idx].median(axis=1).to_frame(f'median_{g}')
        out[f'cor_median_profile_{g}'] = cls.similarity(prof, ccle, genes)[0].iloc[0]
    if gs == PRIMARY:
        for r in REGION_CLASSES:
            out[f'median_cor_{r}'] = cor.loc[s.loc[s['region_class'] == r, 'GSM']].median(axis=0)
    lin = cor[lincs_models]
    best, fdr = cls.compare_with_best(lin)
    out['rank_lincs'] = lin.median(axis=0).rank(ascending=False)
    out['fdr_vs_best_lincs'] = fdr
    print(f'{gs}: {n} genes; best LINCS line {best}')
    return out, cor


def panel(tab, cor_primary, tacle, s):
    p = tab[tab['gene_set'] == PRIMARY]
    lin = p[p['lincs_cell_id'].notna()].sort_values('median_cor_all', ascending=False)
    fig = plt.figure(figsize=(16, 15))
    gsp = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.75])
    ax = {k: fig.add_subplot(v) for k, v in
          {'A': gsp[0, :], 'B': gsp[1, 0], 'C': gsp[1, 1], 'D': gsp[2, 0], 'E': gsp[2, 1]}.items()}

    a = ax['A']
    pos = np.arange(1, len(lin) + 1)
    for g, col, dx in (('als', RED, -0.2), ('control', BLUE, 0.2)):
        idx = s.loc[s['group'] == g, 'GSM']
        data = [cor_primary.loc[idx, m].values for m in lin['ModelID']]
        a.boxplot(data, positions=pos + dx, widths=.35, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=col, alpha=.4, edgecolor=col),
                  whiskerprops=dict(color=col), capprops=dict(color=col),
                  medianprops=dict(color='k', linewidth=1), manage_ticks=False)
        a.plot([], [], 's', color=col, alpha=.4, label=f'{g} samples (n={len(idx)})')
    a.set_xticks(pos)
    a.set_xlim(0.3, len(lin) + 0.7)
    a.set_xticklabels([f'{c} ({l})' for c, l in zip(lin['lincs_cell_id'], lin['OncotreeLineage'])],
                      rotation=90, fontsize=7)
    a.set_ylabel('Spearman correlation with GSE153960 samples')
    a.set_title(f'A. LINCS lines with CCLE expression, per-sample correlations '
                f'({PRIMARY}, {int(p.n_genes.iloc[0])} genes)')
    a.legend(loc='upper right')

    a = ax['B']
    bins = np.linspace(min(p['median_cor_all'].min(), tacle['cor'].min()) - .02,
                       max(p['median_cor_all'].max(), tacle['cor'].max()) + .02, 50)
    a.hist(p['median_cor_all'], bins=bins, color='0.6', alpha=.7, label=f'ALS CNS, all CCLE lines (n={len(p)})')
    a.hist(lin['median_cor_all'], bins=bins, color=RED, alpha=.6, label=f'ALS CNS, LINCS lines (n={len(lin)})')
    for t, ls in zip(TACLE, ['-', '--', ':', '-.']):
        v = tacle.loc[tacle['tumour'] == t, 'cor']
        a.hist(v, bins=bins, histtype='step', color=BLUE, ls=ls, lw=1.4,
               label=f'Chen 2017 {t}, all CCLE lines (n={len(v)})')
    a.set_xlabel('median Spearman correlation per cell line')
    a.set_ylabel('cell lines')
    a.set_title('B. Similarity scale: ALS CNS vs Chen 2017 tumour references')
    a.legend(fontsize=8)

    a = ax['C']
    a.scatter(lin['median_cor_control'], lin['median_cor_als'], color=RED, s=18)
    for _, r in lin.head(8).iterrows():
        a.annotate(r['lincs_cell_id'], (r['median_cor_control'], r['median_cor_als']), fontsize=7)
    v = lin[['median_cor_control', 'median_cor_als']].values
    a.plot([v.min(), v.max()], [v.min(), v.max()], color='k', lw=.6)
    a.set_xlabel('median correlation, control samples')
    a.set_ylabel('median correlation, als samples')
    a.set_title('C. LINCS lines: als vs control samples')

    a = ax['D']
    piv = tab[tab['lincs_cell_id'].notna()].pivot_table(index='lincs_cell_id', columns='gene_set',
                                                         values='median_cor_all')
    for gs, c in zip([g for g in cls.GENE_SETS if g != PRIMARY], [BLUE, RED]):
        rho = spearmanr(piv[PRIMARY], piv[gs]).statistic
        a.scatter(piv[PRIMARY], piv[gs], s=16, color=c, label=f'{gs} (Spearman of orderings {rho:.2f})')
    a.set_xlabel(f'median correlation, {PRIMARY}')
    a.set_ylabel('median correlation, other gene set')
    a.set_title('D. LINCS lines: agreement between marker-gene sets')
    a.legend(fontsize=8)

    a = ax['E']
    n_reg = s['region_class'].value_counts()
    for r, mk in zip(REGION_CLASSES, ['o', 's', '^', 'D']):
        a.scatter(lin['median_cor_all'], lin[f'median_cor_{r}'], s=16, marker=mk,
                  label=f'{r.replace("_", " ")} (n={n_reg.get(r, 0)})')
    v = lin['median_cor_all']
    a.plot([v.min(), v.max()], [v.min(), v.max()], color='k', lw=.6)
    a.set_xlabel('median correlation, all samples')
    a.set_ylabel('median correlation, samples of one region class')
    a.set_title('E. LINCS lines: median correlation per CNS region class')
    a.legend(fontsize=8)

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_ccle_cell_line_similarity_panel.png')


if __name__ == '__main__':
    s, x = load_samples()
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    print(f'CCLE matrix: {ccle.shape[0]} lines x {ccle.shape[1]} genes; '
          f'genes shared with GSE153960: {len(set(ccle.columns) & set(x.index))}')

    lincs = pd.read_csv(LINCS_METADATA_PATH)
    n_inst = lincs.groupby('cell_id').size()
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    mapping['n_lincs_instances'] = mapping['lincs_cell_id'].map(n_inst)
    print(f'LINCS lines {len(mapping)}; with DepMap expression {int(mapping.has_expression.sum())}; '
          f'without: {", ".join(mapping.loc[~mapping.has_expression, "lincs_cell_id"])} '
          f'({int(mapping.loc[~mapping.has_expression, "n_lincs_instances"].sum())} instances)')
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    s.to_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', index=False)
    mapping.to_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t', index=False)

    lincs_models = mapping.loc[mapping.has_expression, 'ModelID'].tolist()
    meta = models.set_index('ModelID').loc[ccle.index, ['CellLineName', 'OncotreeLineage']]
    m2l = mapping[mapping.has_expression].set_index('ModelID')

    tabs = []
    for gs in cls.GENE_SETS:
        out, cor = per_gene_set(x, s, ccle, gs, lincs_models)
        if gs == PRIMARY:
            cor_primary = cor
        out = meta.join(out)
        out['lincs_cell_id'] = m2l['lincs_cell_id'].reindex(out.index)
        out['n_lincs_instances'] = m2l['n_lincs_instances'].reindex(out.index)
        tabs.append(out)
    tab = pd.concat(tabs).rename_axis('ModelID').reset_index()
    cols = ['gene_set', 'ModelID', 'CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'n_lincs_instances',
            'n_genes', 'median_cor_all', 'median_cor_als', 'median_cor_control',
            'cor_median_profile_als', 'cor_median_profile_control', 'rank_lincs', 'fdr_vs_best_lincs'] + \
           [f'median_cor_{r}' for r in REGION_CLASSES]
    tab = tab[cols].sort_values(['gene_set', 'median_cor_all'], ascending=[True, False])
    tab.to_csv(P.oout_dir / f'{P.slug}_ccle_cell_line_similarity.tsv', sep='\t', index=False)

    w = tab[(tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
    w = w[['lincs_cell_id', 'median_cor_als']].rename(columns={'median_cor_als': 'cor'})
    w.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t', index=False)
    print(f'sRGES weights written for {len(w)} LINCS lines; min {w.cor.min():.3f}, max {w.cor.max():.3f}')

    tacle = load_tacle()
    panel(tab, cor_primary, tacle, s)

    p = tab[tab.gene_set == PRIMARY]
    lin = p[p.lincs_cell_id.notna()]
    print('\nprimary set, median per-line correlation (all samples):')
    print(f'  all CCLE lines: median {p.median_cor_all.median():.3f}, max {p.median_cor_all.max():.3f} '
          f'({p.loc[p.median_cor_all.idxmax(), "CellLineName"]})')
    print(f'  LINCS lines:    median {lin.median_cor_all.median():.3f}, max {lin.median_cor_all.max():.3f}')
    for t in TACLE:
        v = tacle.loc[tacle.tumour == t, 'cor']
        print(f'  Chen 2017 {t}: median {v.median():.3f}, max {v.max():.3f}')
    d = lin.median_cor_als - lin.median_cor_control
    print(f'\nals - control median correlation over LINCS lines: min {d.min():+.4f}, max {d.max():+.4f}; '
          f'als higher in {(d > 0).sum()}/{len(d)}')
    print('\nLINCS lines (primary set), all with a weight:')
    print(lin[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'n_lincs_instances',
               'median_cor_all', 'median_cor_als', 'median_cor_control'] +
              [f'median_cor_{r}' for r in REGION_CLASSES] + ['fdr_vs_best_lincs']].round(3).to_string(index=False))
    print('\ntop 15 CCLE lines overall (primary set):')
    print(p.head(15)[['CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'median_cor_all']].round(3).to_string(index=False))
