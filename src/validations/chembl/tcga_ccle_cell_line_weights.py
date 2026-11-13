#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Transcriptomic similarity of the CCLE cell lines (LINCS cell lines among them)
to the TCGA primary tumours of BRCA, LIHC and COAD, i.e. the per-cell-line
weight of the Chen 2017 sRGES for the ChEMBL IC50 validation.

Design
------
- Tumour samples: TCGA primary-tumour samples (barcode sample code 01) of each
  cancer, GDC STAR gene-level log2(TPM + 1) from the UCSC Xena GDC hub; rows
  are mapped Ensembl id -> gene symbol (GENCODE v36 probemap) -> NCBI GeneID
  (DepMap column labels 'SYMBOL (GeneID)'), symbols mapped to several Ensembl
  ids are averaged.
- Cell lines: DepMap log2(TPM + 1) of every CCLE line with expression; LINCS
  lines are the 71 cell_id of conf.LINCS_METADATA_PATH, matched to DepMap by
  cell_line_similarity.
- Marker genes: the three sets of cell_line_similarity.GENE_SETS, selected on
  the CCLE matrix only; 'binchen2017_iqr5000' is the primary set.
- Similarity of a line to a cancer: median over the tumour samples of the
  per-sample Spearman correlations (Chen 2017 / OCTAD).
- sRGES weight of a LINCS line for a cancer: its median correlation with the
  tumour samples on the primary gene set.
- Known answer: Chen 2017 per-line median tumour correlations
  (cell_line_<T>_tacle.csv: CCLE 2017 and TCGA, 5,000 IQR genes). Agreement
  is the Spearman correlation of the two orderings, over the LINCS lines both
  cover (LINCS -> CCLE names of raw/cell_line_lincs_ccle.csv) and over all
  CCLE lines both cover (names matched after upper-casing and removing
  non-alphanumerics). Absolute levels are not compared (different releases,
  normalisations and gene sets).

The data checks at the top print the samples kept, the gene mapping and the
LINCS -> DepMap mapping.

Inputs
------
- TCGA_DIR/TCGA-<T>.star_tpm.tsv.gz, TCGA_DIR/gencode.v36.annotation.gtf.gene.probemap
  (TCGA_DIR = conf.DATA_DIR / 'TCGA_xena_gdc')
- cell_line_similarity.CCLE_EXPRESSION, CCLE_MODEL
- conf.LINCS_METADATA_PATH
- conf.BC_DISEASE_DATA/<T>/cell_line_<T>_tacle.csv,
  conf.BC_DISEASE_DATA/raw/cell_line_lincs_ccle.csv

Outputs (OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR, IMG_OUT = conf.IMG_DIR / 'binchen2017_validation')
-------
- OUT_DIR/tcga_ccle_cell_line_similarity.tsv
    one row per (cancer, gene_set, CCLE line): ModelID, CellLineName,
    OncotreeLineage, lincs_cell_id (empty if not in LINCS), n_tumour_samples,
    n_genes, median_cor, rank_lincs (LINCS lines only, 1 = most similar),
    tacle_cor (Chen 2017, when the line is covered).
- OUT_DIR/tcga_cell_line_weights.tsv
    sRGES weights: cancer, lincs_cell_id, cor (median correlation with the
    tumour samples, primary gene set); LINCS lines with DepMap expression only.
- OUT_DIR/tcga_cell_line_weights_vs_chen2017.tsv
    per cancer and gene set: Spearman of the orderings vs Chen 2017, over the
    shared LINCS lines and over all shared CCLE lines, with the numbers of lines.
- IMG_OUT/tcga_ccle_cell_line_similarity_panel.png
    top: per-sample correlations of the LINCS lines per cancer (primary set),
    ordered by median; bottom: median correlation per LINCS line vs the Chen
    2017 value per cancer.

Usage
-----
    cd src && PYTHONPATH=. python validations/chembl/tcga_ccle_cell_line_weights.py
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

import conf
import cell_line_similarity as cls
import lfs_style

TCGA_DIR = conf.DATA_DIR / 'TCGA_xena_gdc'
PROBEMAP = TCGA_DIR / 'gencode.v36.annotation.gtf.gene.probemap'
OUT_DIR = conf.BINCHEN_VALIDATION_OUT_DIR
IMG_OUT = conf.IMG_DIR / 'binchen2017_validation'
CANCERS = {'BRCA': 'MCF7', 'LIHC': 'HEPG2', 'COAD': 'HT29'}  # cancer: Chen 2017 validation line
PRIMARY = 'binchen2017_iqr5000'
TUMOUR_CODE = '01'


def _strip(name):
    return re.sub(r'[^0-9A-Z]', '', str(name).upper())


def ccle_symbol_to_geneid():
    cols = pd.read_csv(cls.CCLE_EXPRESSION, nrows=0).columns[1:]
    m = cols.str.extract(r'^(.*) \((\d+)\)$').dropna()
    return dict(zip(m[0], m[1].astype(int)))


def load_tumours(cancer, sym2id):
    x = pd.read_csv(TCGA_DIR / f'TCGA-{cancer}.star_tpm.tsv.gz', sep='\t', index_col=0)
    codes = x.columns.str[13:15]
    print(f'{cancer}: {x.shape[0]} genes x {x.shape[1]} samples; sample codes '
          f'{pd.Series(codes).value_counts().to_dict()}')
    x = x.loc[:, codes == TUMOUR_CODE]
    pm = pd.read_csv(PROBEMAP, sep='\t', usecols=['id', 'gene']).set_index('id')['gene']
    sym = x.index.map(pm)
    gid = pd.Series(sym).map(sym2id).values
    keep = pd.notna(gid)
    x = x.loc[keep]
    x.index = gid[keep].astype(int)
    n_ens = int(keep.sum())
    x = x.groupby(level=0).mean()
    print(f'{cancer}: {x.shape[1]} primary tumours; {n_ens} Ensembl ids -> {x.shape[0]} DepMap GeneIDs')
    return x


def load_tacle(cancer, models, lincs_ccle):
    """Chen 2017 per-line cor, indexed by DepMap ModelID (CCLE name match) and by LINCS cell_id."""
    t = pd.read_csv(conf.BC_DISEASE_DATA / cancer / f'cell_line_{cancer}_tacle.csv')
    by_strip = models.drop_duplicates('stripped').set_index('stripped')['ModelID']
    t['ModelID'] = t['Cell.line.primary.name'].map(_strip).map(by_strip)
    by_model = t.dropna(subset=['ModelID']).drop_duplicates('ModelID').set_index('ModelID')['cor']
    by_lincs = t.merge(lincs_ccle, left_on='Cell.line.primary.name', right_on='ccle_cell_line_name')
    by_lincs = by_lincs.drop_duplicates('lincs_cell_id').set_index('lincs_cell_id')['cor']
    return by_model, by_lincs, len(t)


def panel(tab, cors):
    lfs_style.apply()
    fig, ax = plt.subplots(2, 3, figsize=(20, 11), gridspec_kw={'height_ratios': [1.3, 1]})
    for j, cancer in enumerate(CANCERS):
        p = tab[(tab.cancer == cancer) & (tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
        p = p.sort_values('median_cor', ascending=False)
        a = ax[0, j]
        data = [cors[cancer][m].values for m in p['ModelID']]
        colors = [lfs_style.AFTER if c == CANCERS[cancer] else lfs_style.BEFORE for c in p['lincs_cell_id']]
        bp = a.boxplot(data, showfliers=False, patch_artist=True, medianprops=dict(color='k', linewidth=1))
        for b, c in zip(bp['boxes'], colors):
            b.set_facecolor(c); b.set_alpha(.45)
        a.set_xticks(np.arange(1, len(p) + 1))
        a.set_xticklabels(p['lincs_cell_id'], rotation=90, fontsize=6)
        a.set_ylabel('Spearman correlation with tumour samples')
        a.set_title(f'{cancer}: {int(p.n_tumour_samples.iloc[0])} primary tumours, '
                    f'{len(p)} LINCS lines ({int(p.n_genes.iloc[0])} genes)', fontsize=10)
        a = ax[1, j]
        q = p.dropna(subset=['tacle_cor'])
        a.scatter(q['tacle_cor'], q['median_cor'], s=18, color=lfs_style.BEFORE)
        for _, r in q.iterrows():
            a.annotate(r['lincs_cell_id'], (r['tacle_cor'], r['median_cor']), fontsize=6,
                       color=lfs_style.AFTER if r['lincs_cell_id'] == CANCERS[cancer] else 'k')
        rho = spearmanr(q['tacle_cor'], q['median_cor']).statistic
        a.set_xlabel('Chen 2017 median correlation (CCLE 2017, TCGA)')
        a.set_ylabel('median correlation (DepMap 24Q4, TCGA GDC)')
        a.set_title(f'{cancer}: LINCS lines in both, n={len(q)}, Spearman {rho:.2f}', fontsize=10)
    fig.tight_layout()
    IMG_OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG_OUT / 'tcga_ccle_cell_line_similarity_panel.png', dpi=150)


if __name__ == '__main__':
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    sym2id = ccle_symbol_to_geneid()
    lincs = pd.read_csv(conf.LINCS_METADATA_PATH)
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    print(f'CCLE matrix: {ccle.shape[0]} lines x {ccle.shape[1]} genes; LINCS lines {len(mapping)}, '
          f'with DepMap expression {int(mapping.has_expression.sum())}')
    lincs_models = mapping.loc[mapping.has_expression, 'ModelID'].tolist()
    m2l = mapping[mapping.has_expression].set_index('ModelID')['lincs_cell_id']
    meta = models.set_index('ModelID').loc[ccle.index, ['CellLineName', 'OncotreeLineage']]
    lincs_ccle = pd.read_csv(conf.BC_DISEASE_DATA / 'raw' / 'cell_line_lincs_ccle.csv')
    genes = {gs: cls.marker_genes(ccle, gs, reference_models=lincs_models) for gs in cls.GENE_SETS}

    tabs, agree, cors = [], [], {}
    for cancer in CANCERS:
        x = load_tumours(cancer, sym2id)
        tac_model, tac_lincs, n_tac = load_tacle(cancer, models, lincs_ccle)
        print(f'{cancer}: Chen 2017 lines {n_tac}; matched to DepMap {len(tac_model)}; to LINCS {len(tac_lincs)}')
        for gs in cls.GENE_SETS:
            cor, n = cls.similarity(x, ccle, genes[gs])
            if gs == PRIMARY:
                cors[cancer] = cor
            out = meta.copy()
            out.insert(0, 'gene_set', gs)
            out.insert(0, 'cancer', cancer)
            out['lincs_cell_id'] = m2l.reindex(out.index)
            out['n_tumour_samples'] = cor.shape[0]
            out['n_genes'] = n
            out['median_cor'] = cor.median(axis=0)
            out['rank_lincs'] = out.loc[lincs_models, 'median_cor'].rank(ascending=False)
            out['tacle_cor'] = tac_model.reindex(out.index)
            tabs.append(out)
            lin = out[out.lincs_cell_id.notna()].set_index('lincs_cell_id')
            sl = lin['median_cor'].to_frame().join(tac_lincs.rename('tacle'), how='inner')
            sa = out[['median_cor', 'tacle_cor']].dropna()
            agree.append({'cancer': cancer, 'gene_set': gs,
                          'spearman_lincs_lines': spearmanr(sl['median_cor'], sl['tacle']).statistic,
                          'n_lincs_lines': len(sl),
                          'spearman_all_ccle_lines': spearmanr(sa['median_cor'], sa['tacle_cor']).statistic,
                          'n_ccle_lines': len(sa)})
            print(f'  {gs}: {n} genes; vs Chen 2017: LINCS lines rho={agree[-1]["spearman_lincs_lines"]:.3f} '
                  f'(n={len(sl)}), all CCLE lines rho={agree[-1]["spearman_all_ccle_lines"]:.3f} (n={len(sa)})')

    tab = pd.concat(tabs).rename_axis('ModelID').reset_index()
    tab = tab.sort_values(['cancer', 'gene_set', 'median_cor'], ascending=[True, True, False])
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT_DIR / 'tcga_ccle_cell_line_similarity.tsv', sep='\t', index=False)
    w = tab[(tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()][['cancer', 'lincs_cell_id', 'median_cor']]
    w.rename(columns={'median_cor': 'cor'}).to_csv(OUT_DIR / 'tcga_cell_line_weights.tsv', sep='\t', index=False)
    ag = pd.DataFrame(agree)
    ag.to_csv(OUT_DIR / 'tcga_cell_line_weights_vs_chen2017.tsv', sep='\t', index=False)
    print('\n' + ag.round(3).to_string(index=False))

    for cancer, line in CANCERS.items():
        p = tab[(tab.cancer == cancer) & (tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
        print(f'\n{cancer}: weights min {p.median_cor.min():.3f}, max {p.median_cor.max():.3f}; '
              f'{line} rank {int(p.loc[p.lincs_cell_id == line, "rank_lincs"].iloc[0])}/{len(p)}')
        print(p.head(8)[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'median_cor', 'tacle_cor']]
              .to_string(index=False))
    panel(tab, cors)
