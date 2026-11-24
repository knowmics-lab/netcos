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
LINCS_TO_DEPMAP_ALIAS. A LINCS line without DepMap expression can take the
profile of a Human Protein Atlas cell line (HPA_CELL_LINE, TPM; HPA_TISSUE
for tissue references), listed in EXTERNAL_PROFILES together with how close
the HPA line is to the LINCS line ('same line', 'parent line', 'same cell
type'). HPA profiles enter as log2(TPM + 1) over the genes whose symbol is
unique in HPA and in the DepMap header; they are scored by the same
marker genes and the same Spearman summary as the DepMap lines. A LINCS line
in neither DepMap nor HPA can take a GEO RNA-seq profile (GEO_PROFILES): raw
gene counts on Ensembl gene IDs, turned into TPM with the union-of-exons gene
length of ENSEMBL75_GTF, averaged over the listed samples and mapped to NCBI
GeneID through NCBI gene2ensembl (one-to-one pairs).
"""

import re

import numpy as np
import pandas as pd
from scipy.stats import rankdata, wilcoxon, false_discovery_control

import conf as _conf
from conf import CCLE_DIR, LINCS_BC_DATA, DATA_DIR

DEPMAP_RELEASE = 'DepMap_24Q4'
CCLE_EXPRESSION = CCLE_DIR / DEPMAP_RELEASE / 'OmicsExpressionProteinCodingGenesTPMLogp1.csv'
CCLE_MODEL = CCLE_DIR / DEPMAP_RELEASE / 'Model.csv'
LINCS_LANDMARK = LINCS_BC_DATA / 'lincs_landmark.csv'

LINCS_TO_DEPMAP_ALIAS = {'H1299': 'NCIH1299', 'U266': 'U266B1'}

# paths: conf.HPA_CELL_LINE, conf.HPA_TISSUE, conf.DEPMAP_PROFILES_DIR when defined, else under DATA_DIR
HPA_CELL_LINE = getattr(_conf, 'HPA_CELL_LINE', DATA_DIR / 'HPA' / 'rna_celline.tsv.zip')
HPA_TISSUE = getattr(_conf, 'HPA_TISSUE', DATA_DIR / 'HPA' / 'rna_tissue_consensus.tsv.zip')

# LINCS cell_id -> (HPA cell line, how close it is to the LINCS line, DepMap OncotreeLineage
# of the DepMap lines used to compare the two sources, or None)
EXTERNAL_PROFILES = {
    'MCF10A': ('MCF-10A', 'same line', 'Breast'),
    'SKLU1': ('SK-LU-1', 'same line', 'Lung'),
    'HEK293T': ('HEK293', 'parent line', 'Kidney'),
    'MCH58': ('fHDF/TERT166', 'same cell type', 'Fibroblast'),
}
# profiles that enter the candidate weights only when requested (external_optional)
EXTERNAL_PROFILES_OPTIONAL = {
    'ASC': ('ASC52telo', 'same cell type, hTERT-immortalised', None),
}
# RNA-seq from GEO of LINCS lines that are in neither DepMap nor HPA: raw gene counts on Ensembl
# GRCh37 release 75 gene IDs, turned into TPM with the union-of-exons length of each gene.
# paths: conf.GEO_GSE144279_COUNTS, conf.ENSEMBL75_GTF when defined, else under DATA_DIR
GEO_GSE144279_COUNTS = getattr(_conf, 'GEO_GSE144279_COUNTS',
                               DATA_DIR / 'GEO' / 'GSE144279' / 'GSE144279_all.cnts.txt.gz')
ENSEMBL75_GTF = getattr(_conf, 'ENSEMBL75_GTF',
                        DATA_DIR / 'Ensembl' / 'release75' / 'Homo_sapiens.GRCh37.75.gtf.gz')
GENE2ENSEMBL = getattr(_conf, 'GENE2ENSEMBL', DATA_DIR / 'NCBI_gene' / 'gene2ensembl.gz')
# LINCS cell_id -> (profile name, how close it is to the LINCS line, count file, wild-type samples
# averaged into the profile). LINCS NPC / NEU: aliases NPC-8330-8 NEURAL PROGENITORS / NEURONS
# (lincs_cell_info.csv), from the iPSC line 8330-iPS-8 (Cellosaurus CVCL_VK56, parent GM08330).
GEO_PROFILES = {
    'NPC': ('GSE144279 8330-8 NPC', 'same donor iPSC, external differentiation', GEO_GSE144279_COUNTS,
            ['NPC_GM_P32', 'NPC_GM_R_B1_P10', 'NPC_H8_R_B2_P10']),
    'NEU': ('GSE144279 8330-8 neurons', 'same donor iPSC, external differentiation', GEO_GSE144279_COUNTS,
            ['Neuron_GM_P37', 'Neuron_GM_R_B1_P14', 'Neuron_H8_R_B2_P14']),
}
# Cellosaurus accession of the iPSC line the LINCS line derives from
PARENT_CELLOSAURUS_ID = {'NPC': 'CVCL_VK56', 'NEU': 'CVCL_VK56'}

# tissue references of HPA_TISSUE (not used as weights)
TISSUE_REFERENCES = ['liver', 'skeletal muscle', 'adipose tissue', 'bone marrow', 'lung']

# Cellosaurus accession of the LINCS lines that are cell lines (RRID of DepMap Model.csv)
CELLOSAURUS_ID = {'MCF10A': 'CVCL_0598', 'SKLU1': 'CVCL_0629', 'HEK293T': 'CVCL_0063',
                  'HS27A': 'CVCL_3719', 'MCH58': 'CVCL_U486'}
DEPMAP_PROFILES_DIR = getattr(_conf, 'DEPMAP_PROFILES_DIR', CCLE_DIR / 'DepMap_OmicsProfiles')
DEPMAP_PROFILE_RELEASES = ['22Q4', '23Q2', '23Q4', '24Q2', '24Q4']
LINCS_CELL_INFO = LINCS_BC_DATA / 'lincs_cell_info.csv'

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


def ccle_gene_symbols():
    """DepMap gene symbol -> NCBI GeneID, for the symbols that appear once in the expression header."""
    head = pd.read_csv(CCLE_EXPRESSION, nrows=0, index_col=0).columns
    sym = head.str.extract(r'^(.*) \((\d+)\)$')
    sym = sym.dropna()
    sym = sym[~sym[0].duplicated(keep=False)]
    return dict(zip(sym[0], sym[1].astype(int)))


def _hpa_matrix(path, label, names, symbols, value):
    """log2(value + 1) of the HPA rows whose label is in names: rows = label, columns = NCBI GeneID."""
    parts = []
    for ch in pd.read_csv(path, sep='\t', chunksize=2_000_000, usecols=['Gene name', label, value]):
        parts.append(ch if names is None else ch[ch[label].isin(names)])
    d = pd.concat(parts)
    d = d[~d.duplicated(['Gene name', label], keep=False)]  # symbols with several Ensembl genes are ambiguous
    d['gene_id'] = d['Gene name'].map(symbols)
    d = d.dropna(subset=['gene_id'])
    d = d[~d.duplicated(['gene_id', label], keep=False)]
    x = d.pivot(index=label, columns='gene_id', values=value)
    x.columns = x.columns.astype(int)
    return np.log2(x + 1)


def load_hpa_cell_lines(names, symbols):
    """HPA cell-line log2(TPM + 1): rows = HPA cell line, columns = NCBI GeneID; names None = every line."""
    return _hpa_matrix(HPA_CELL_LINE, 'Cell line', None if names is None else set(names), symbols, 'TPM')


def ensembl_genes(gtf=None):
    """Ensembl gene ID -> (gene name, union-of-exons length in bp) from a GTF."""
    import gzip
    exons, names = {}, {}
    with gzip.open(gtf or ENSEMBL75_GTF, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.split('\t', 9)
            if p[2] != 'exon':
                continue
            gid = re.search(r'gene_id "([^"]+)"', p[8]).group(1)
            names[gid] = re.search(r'gene_name "([^"]+)"', p[8]).group(1)
            exons.setdefault(gid, []).append((int(p[3]), int(p[4])))
    length = {}
    for gid, iv in exons.items():
        iv.sort()
        total, (s, e) = 0, iv[0]
        for a, b in iv[1:]:
            if a > e:
                total += e - s + 1
                s, e = a, b
            else:
                e = max(e, b)
        length[gid] = total + e - s + 1
    return pd.DataFrame({'name': pd.Series(names), 'length': pd.Series(length)})


def ensembl_to_gene_id(path=None):
    """Human Ensembl gene ID -> NCBI GeneID from NCBI gene2ensembl, pairs that are one-to-one only."""
    parts = []
    for ch in pd.read_csv(path or GENE2ENSEMBL, sep='\t', usecols=[0, 1, 2], chunksize=5_000_000,
                          names=['tax_id', 'GeneID', 'ensembl'], header=0, dtype={'tax_id': int, 'GeneID': int}):
        parts.append(ch[ch['tax_id'] == 9606])
    d = pd.concat(parts)[['GeneID', 'ensembl']].drop_duplicates()
    d = d[~d['ensembl'].duplicated(keep=False) & ~d['GeneID'].duplicated(keep=False)]
    return d.set_index('ensembl')['GeneID']


def load_geo_profiles(profiles=None):
    """
    GEO profiles of GEO_PROFILES as log2(TPM + 1), rows = profile name, columns = NCBI GeneID; the TPM of
    a profile is the mean TPM of its samples. Ensembl gene IDs are mapped to NCBI GeneID through
    gene2ensembl (one-to-one pairs only).
    """
    profiles = profiles or GEO_PROFILES
    genes = ensembl_genes()
    to_id = ensembl_to_gene_id()
    mean = {}
    for path in {v[2] for v in profiles.values()}:
        c = pd.read_csv(path, sep='\t', index_col=0)
        c = c.loc[c.index.isin(genes.index)]
        rpk = c.div(genes.loc[c.index, 'length'] / 1e3, axis=0)
        tpm = rpk / rpk.sum() * 1e6
        for name, _, p, samples in profiles.values():
            if p == path:
                mean[name] = tpm[samples].mean(axis=1)
    t = pd.DataFrame(mean)
    t = t.loc[t.index.isin(to_id.index)]
    t.index = to_id.loc[t.index].astype(int).values
    return np.log2(t.T + 1)


def load_hpa_tissues(names, symbols):
    """HPA consensus tissue log2(nTPM + 1): rows = tissue, columns = NCBI GeneID."""
    return _hpa_matrix(HPA_TISSUE, 'Tissue', set(names), symbols, 'nTPM')


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
