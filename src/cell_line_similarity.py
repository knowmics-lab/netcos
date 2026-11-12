# -*- coding: utf-8 -*-
"""
Transcriptomic similarity between disease samples and CCLE cell lines.

Baseline expression of the cell lines (DepMap, log2(TPM + 1), protein-coding
genes) is compared with the expression of disease samples by Spearman
correlation over a set of marker genes chosen on the cell lines only.

Marker-gene sets (all selected on the CCLE matrix, then intersected with the
genes measured in the samples):
- 'binchen2017_iqr5000': the 5,000 genes with the largest interquartile range
  of log2(TPM + 1) across all CCLE cell lines (Chen et al., Nat. Commun. 2017,
  Supplementary Methods).
- 'octad_ranksd1000': OCTAD computeCellLine (Zeng et al., Nat. Protoc. 2021;
  R package octad): genes with median expression > 1 across the reference
  lines are kept, every cell line is rank-transformed across those genes, and
  the 1,000 genes with the largest standard deviation of the rank across all
  CCLE cell lines are kept.
- 'lincs_landmark978': the 978 L1000 landmark genes.

Per cell line, the similarity summary follows Chen 2017 / OCTAD: the median,
over the samples of a group, of the per-sample Spearman correlations. The
Spearman correlation with the median expression profile of a group is also
returned. OCTAD's comparison of every cell line with the best one (one-sided
paired Wilcoxon test on the per-sample correlations, Benjamini-Hochberg FDR)
is reproduced.

Cell-line identifiers: LINCS cell_id is matched to DepMap
StrippedCellLineName (upper case, alphanumerics only), with the aliases in
LINCS_TO_DEPMAP_ALIAS.
"""

import re

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon, false_discovery_control

from conf import CCLE_DIR, LINCS_BC_DATA

DEPMAP_RELEASE = 'DepMap_24Q4'
CCLE_EXPRESSION = CCLE_DIR / DEPMAP_RELEASE / 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'
CCLE_MODEL = CCLE_DIR / DEPMAP_RELEASE / 'Model.csv'
LINCS_LANDMARK = LINCS_BC_DATA / 'lincs_landmark.csv'

LINCS_TO_DEPMAP_ALIAS = {'H1299': 'NCIH1299', 'U266': 'U266B1'}

GENE_SETS = ['binchen2017_iqr5000', 'octad_ranksd1000', 'lincs_landmark978']


def _strip(name):
    return re.sub(r'[^0-9A-Z]', '', str(name).upper())


def load_ccle_expression():
    """CCLE log2(TPM + 1): rows = DepMap ModelID, columns = NCBI GeneID (int)."""
    x = pd.read_csv(CCLE_EXPRESSION, index_col=0)
    gene_ids = x.columns.str.extract(r'\((\d+)\)$')[0]
    keep = gene_ids.notna().values
    x = x.loc[:, keep]
    x.columns = gene_ids[keep].astype(int).values
    x = x.T.groupby(level=0).mean().T  # a GeneID appears once; guards against repeats
    return x


def load_ccle_models():
    m = pd.read_csv(CCLE_MODEL)
    m['stripped'] = m['StrippedCellLineName'].map(_strip)
    return m


def map_lincs_to_ccle(lincs_cell_ids, models, expressed_models):
    """One row per LINCS cell_id: ModelID, name, lineage, has_expression."""
    by_name = models.drop_duplicates('stripped').set_index('stripped')
    rows = []
    for cid in sorted(set(lincs_cell_ids)):
        key = _strip(LINCS_TO_DEPMAP_ALIAS.get(cid, cid))
        if key in by_name.index:
            r = by_name.loc[key]
            rows.append((cid, r['ModelID'], r['CellLineName'], r['OncotreeLineage'],
                         r['ModelID'] in expressed_models))
        else:
            rows.append((cid, None, None, None, False))
    return pd.DataFrame(rows, columns=['lincs_cell_id', 'ModelID', 'CellLineName',
                                       'OncotreeLineage', 'has_expression'])


def marker_genes(ccle, gene_set, reference_models=None, n=None):
    """Marker genes (NCBI GeneIDs) selected on the CCLE matrix (lines x genes)."""
    if gene_set == 'binchen2017_iqr5000':
        n = n or 5000
        q75, q25 = np.nanpercentile(ccle.values, [75, 25], axis=0)
        iqr = pd.Series(q75 - q25, index=ccle.columns)
        return iqr.sort_values(ascending=False).index[:n].tolist()
    if gene_set == 'octad_ranksd1000':
        n = n or 1000
        ref = ccle.loc[reference_models] if reference_models is not None else ccle
        expressed = ref.columns[ref.median(axis=0) > 1]
        ranks = pd.DataFrame(rankdata(ccle[expressed].values, axis=1), columns=expressed)
        return ranks.std(axis=0).sort_values(ascending=False).index[:n].tolist()
    if gene_set == 'lincs_landmark978':
        return pd.read_csv(LINCS_LANDMARK)['gene_id'].astype(int).tolist()
    raise ValueError(gene_set)


def spearman_rows(a, b):
    """Spearman correlation of every row of a (n x g) with every row of b (m x g)."""
    ra = rankdata(a, axis=1)
    rb = rankdata(b, axis=1)
    ra = (ra - ra.mean(axis=1, keepdims=True)) / ra.std(axis=1, keepdims=True)
    rb = (rb - rb.mean(axis=1, keepdims=True)) / rb.std(axis=1, keepdims=True)
    return ra @ rb.T / a.shape[1]


def similarity(samples, ccle, genes):
    """
    samples: genes x samples (log2 expression, NCBI GeneID index);
    ccle: lines x genes. Returns the samples x lines Spearman matrix over the
    marker genes present in both, and the number of genes used.
    """
    g = [x for x in genes if x in samples.index and x in ccle.columns]
    cor = spearman_rows(samples.loc[g].T.values, ccle[g].values)
    return pd.DataFrame(cor, index=samples.columns, columns=ccle.index), len(g)


def compare_with_best(cor):
    """OCTAD: one-sided paired Wilcoxon of the best line against every other line, BH FDR."""
    med = cor.median(axis=0)
    best = med.idxmax()
    others = [c for c in cor.columns if c != best]
    p = [wilcoxon(cor[best], cor[c], alternative='greater').pvalue for c in others]
    fdr = pd.Series(false_discovery_control(p, method='bh'), index=others)
    fdr[best] = np.nan
    return best, fdr
