"""
Step 'weights': transcriptomic similarity of the CCLE cell lines (LINCS cell
lines among them) to the disease samples of a case study, i.e. the
per-cell-line weight of the Chen 2017 / OCTAD sRGES.

Design
------
- Samples: log2(TPM + 1) of the NCBI-recomputed TPM of the study
  (cc.NCBI_TPM, rows = NCBI GeneID). Sample groups from cc.SAMPLES; the first
  group of SAMPLES['groups'] is the case group, the other the control group.
  SAMPLES['format']:
  * 'samples_table': the existing P.oout_dir/<slug>_samples.tsv, group from
    its column SAMPLES['group_column'];
  * 'series_matrix': the GEO series matrix SAMPLES['file'], group from the
    sample characteristic SAMPLES['characteristic'] (SAMPLES['control_value']
    = control, any other value = case); optional checks that the group coding
    equals SAMPLES['check_group_coding'] and that the columns of
    SAMPLES['check_sample_order'] (header only) are the same samples in the
    same order;
  * 'csv': a metadata table SAMPLES['file'] (read with SAMPLES['read_csv']),
    sample id in SAMPLES['sample_column'], group from SAMPLES['group_column'];
    samples of any other group are left out.
  Optional SAMPLES['strata']: a sample class derived from a metadata column
  (values lower-cased; a value starting with a key of 'prefixes' takes that
  class, any other value is its own class; 'label' / 'short_label' name the
  class in the figure); on the primary set the median correlation is also
  given per class in 'classes'.
  TPM is used on both sides because CCLE is TPM and per-sample ranks across
  genes depend on gene-length normalisation.
- Cell lines: DepMap log2(TPM + 1) of every CCLE line with expression; LINCS
  lines are the cell_id of the L1000 level-5 gold set
  (conf.LINCS_METADATA_PATH), matched to DepMap by cell_line_similarity.
- Marker genes: the three sets of cell_line_similarity.GENE_SETS, selected on
  the CCLE matrix only; 'binchen2017_iqr5000' is the primary set.
- Similarity of a line to a group: median over the group's samples of the
  per-sample Spearman correlations (Chen 2017 / OCTAD), for each group and
  all samples; Spearman correlation with the median profile of each group.
- sRGES weight of a LINCS line: its median correlation with the case samples
  (median_cor_<case>) on the primary gene set.
- Among LINCS lines: OCTAD's one-sided paired Wilcoxon of the best line
  against every other line (all samples), BH FDR.
- Cell line vs cell line reference (primary set): Spearman correlation of
  every pair between two disjoint sets of REF_N random CCLE lines (seed
  REF_SEED), and of every pair of LINCS lines; quartiles of each.
- Top cell lines: the LINCS lines whose sRGES weight is above a cutoff of the
  random-pair distribution, for each cutoff in TOP_CUTOFFS ('q3' = third
  quartile, 'median'); these lines feed the top-line sRGES of step 'srges'.
- Reference scale: the median tumour-cell line correlations of Chen 2017
  (BRCA, LIHC, COAD, ER+ BRCA; cell_line_<T>_tacle.csv, CCLE 2017 and TCGA,
  5,000 IQR genes) over all the CCLE lines they cover.

- Candidate weights (optional, run(cc, external_weights=True, external_optional=...) or
  run_candidate(cc, optional), which writes only the candidate files; the lines named in
  external_optional, among cell_line_similarity.EXTERNAL_PROFILES_OPTIONAL, are added to
  the candidate set, otherwise they are only reported in the effect table): LINCS lines
  without DepMap expression that have an entry in
  cell_line_similarity.EXTERNAL_PROFILES take the same weight computed from
  the Human Protein Atlas (HPA) cell-line profile (log2(TPM + 1), same marker
  genes, same median over the case samples); those with an entry in
  cell_line_similarity.GEO_PROFILES take it from GEO RNA-seq counts (TPM by
  union-of-exons gene length, mean of the listed samples). Comparability of the two
  sources: every DepMap line with expression that also has an HPA profile
  (match on the stripped cell-line name) is weighted from both sources; the
  gene-level Spearman between the two profiles of each paired line, the paired weights, their Spearman and Pearson correlation, the median and
  median absolute difference (HPA - DepMap) and the least-squares line
  HPA = a + b * DepMap are given over all shared lines and over the shared
  lines of the DepMap lineage of each external line, and the mean and standard
  deviation of the HPA weights of the paired lines and of the HPA lines not in DepMap. HPA consensus tissue
  profiles (cell_line_similarity.TISSUE_REFERENCES) are scored the same way
  and set against the distribution of the DepMap-line weights; they are
  references and never weights. The effect of the candidate weights is the
  number of L1000 instances (dose > 0, 6 h or 24 h, the instances sRGES
  uses) and of compounds that enter the weighted sRGES.

The data checks at the top print the sample alignment, the group counts, the
gene overlap of the two matrices and the LINCS -> DepMap mapping.

Inputs
------
- cc.NCBI_TPM, cc.SAMPLES (and the files it names)
- cell_line_similarity.CCLE_EXPRESSION, CCLE_MODEL, LINCS_LANDMARK
- conf.LINCS_METADATA_PATH, conf.BC_DISEASE_DATA/<T>/cell_line_<T>_tacle.csv

Outputs
-------
- P.oout_dir/<slug>_samples.tsv (formats 'series_matrix' and 'csv')
    the samples used, with their group.
- P.oout_dir/<slug>_ccle_cell_line_similarity.tsv
    one row per (gene_set, CCLE line): ModelID, CellLineName, OncotreeLineage,
    lincs_cell_id (empty if not in LINCS), n_lincs_instances, n_genes,
    median_cor_all / _<case> / _control, cor_median_profile_<case> /
    _control, rank_lincs (by median_cor_all, LINCS lines only),
    fdr_vs_best_lincs; with strata, primary set only: median_cor_<class>.
- P.oout_dir/<slug>_cell_line_weights.tsv
    sRGES weights: lincs_cell_id, cor = median_cor_<case> of the primary gene
    set, source; the LINCS lines with DepMap expression (source 'DepMap'),
    plus, with external_weights=True, the lines weighted from HPA or GEO
    (lines without a weight are dropped by sRGES.compute_sRGES). The top cell
    lines below are selected from this file.
- P.oout_dir/<slug>_external_cell_line_similarity.tsv   (external_weights=True)
    the LINCS lines weighted from HPA or GEO, primary gene set, in the columns
    of the similarity table: ModelID = CellLineName = profile name,
    OncotreeLineage = 'HPA' or 'GEO', lincs_cell_id, n_genes, median_cor_all,
    median_cor_<group>; they also enter panels A-C (and E) of the figure.
- P.oout_dir/<slug>_lincs_ccle_mapping.tsv
    the LINCS cell_id with their DepMap match and n_lincs_instances.
- P.oout_dir/<slug>_ccle_reference_scale.tsv
    one row per comparison (random CCLE pairs, LINCS pairs): n_pairs, q1,
    median, q3.
- P.oout_dir/<slug>_top_cell_lines.tsv
    one row per (cutoff, LINCS line above it): cutoff, threshold,
    lincs_cell_id, cor; a cutoff with no line above it has one row with
    empty lincs_cell_id.
- P.oout_dir/<slug>_cell_line_weights_candidate.tsv   (external_weights=True)
    lincs_cell_id, cor, source: the weights of the file above (source
    'DepMap') plus those of the lines weighted from HPA (source 'HPA <line>,
    <match>') or GEO (source 'GEO <profile>, <match>'); with external_weights=True
    it is also written as <slug>_cell_line_weights.tsv.
- P.oout_dir/<slug>_cell_line_weights_candidate_check.tsv   (external_weights=True)
    one row per paired line (kind 'shared_line': CellLineName,
    OncotreeLineage, n_genes of each source, cor_depmap, cor_hpa, diff =
    cor_hpa - cor_depmap) and per tissue reference (kind 'tissue_reference':
    cor_hpa, percentile of that value among the DepMap-line weights).
- P.oout_dir/<slug>_cell_line_weights_candidate_check_summary.tsv   (external_weights=True)
    one row per set of paired lines ('all' and one per lineage): n, Spearman,
    Pearson, median, median / 95th percentile / maximum absolute diff, slope and intercept of HPA on
    DepMap; and the recomputation check of the DepMap weights against the
    file above.
- P.oout_dir/<slug>_cell_line_weights_candidate_diagnosis.tsv   (external_weights=True)
    one row per LINCS line without DepMap expression: its LINCS annotation
    (cell type, histology, lineage), Cellosaurus accession where the line is
    a cell line, the DepMap model (name match or Model.RRID) and the number of
    DepMap RNA-seq profiles of that model in each release of
    cell_line_similarity.DEPMAP_PROFILE_RELEASES.
- P.oout_dir/<slug>_cell_line_weights_candidate_effect.tsv   (external_weights=True)
    kind 'line': one row per LINCS line without a DepMap weight (n_instances,
    candidate weight and its percentile among the DepMap-line weights);
    kind 'known_drug': per known drug, instances on lines weighted now / with
    the candidate weights / still unweighted, with the instances per unweighted
    line; kind 'all_compounds': the same
    counts over every compound (instances, and the number of compounds with
    at least one weighted instance).
- P.img_dir/<slug>_ccle_cell_line_similarity_panel.png
    A: per-sample correlations of the LINCS lines (primary set), case and
    control samples side by side, ordered by median over all samples;
    B: per-line median correlation, all CCLE vs LINCS lines vs the Chen 2017
    tumour references; C: case vs control median per LINCS line;
    D: agreement of the LINCS-line ordering between gene sets; with strata,
    E: median correlation of the LINCS lines per sample class.
    Figure text from cc.FIGURES; house style (case_study.common).
"""

import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, rankdata

from conf import LINCS_METADATA_PATH, BC_DISEASE_DATA
import cell_line_similarity as cls
from case_study.common import figures, colours

PRIMARY = 'binchen2017_iqr5000'
TACLE = ['BRCA', 'LIHC', 'COAD', 'ER']
REF_N = 200      # lines in each of the two random sets of the cell line vs cell line reference
REF_SEED = 0
TOP_CUTOFFS = {'q3': 'q3', 'median': 'median'}  # cutoff name: quantile of the random-pair distribution


# ---------------------------------------------------------------- samples

def read_series_matrix(path, characteristic):
    """GSM, title and one sample characteristic of every sample of a GEO series matrix
    (the !Sample_characteristics_ch1 row whose every value is '<characteristic>: <value>')."""
    col = characteristic.replace(' ', '_')
    rows = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('!Sample_geo_accession'):
                rows['GSM'] = [v.strip('"') for v in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_title'):
                rows['title'] = [v.strip('"') for v in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'):
                vals = [v.strip('"').split(': ', 1) for v in line.rstrip('\n').split('\t')[1:]]
                if vals and all(len(v) == 2 and v[0] == characteristic for v in vals):
                    rows[col] = [v[1] for v in vals]
    if col not in rows:
        raise ValueError(f'{path}: no !Sample_characteristics_ch1 row with the field {characteristic!r} '
                         f"in every sample (SAMPLES['characteristic'])")
    return pd.DataFrame(rows)


def samples_from_table(cc, spec):
    P = cc.P
    col = spec['group_column']
    s = pd.read_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', dtype={col: str})
    s['group'] = s[col].map({v: k for k, v in spec['groups'].items()})
    unknown = sorted(s.loc[s['group'].isna(), col].astype(str).unique())
    if unknown:
        print(f"WARNING: {col} values not in SAMPLES['groups'] (samples in no group): {unknown}")
    print(f'samples table: {s["group"].value_counts().to_dict()}')
    return s


def samples_from_series_matrix(cc, spec):
    groups = spec['groups']
    case = list(groups)[0]
    col = spec['characteristic'].replace(' ', '_')
    s = read_series_matrix(spec['file'], spec['characteristic'])
    print(f'series matrix, {spec["characteristic"]}:', s[col].value_counts().to_dict())
    s['group'] = np.where(s[col] == spec['control_value'], 'control', case)
    s['group_code'] = s['group'].map(groups)
    if 'check_group_coding' in spec:
        coding = ''.join(s['group_code'])
        print(f'group coding from GEO == disease-signature coding: {coding == spec["check_group_coding"]}')
        assert coding == spec['check_group_coding']
    if 'check_sample_order' in spec:
        counts_cols = pd.read_csv(spec['check_sample_order'], sep='\t', nrows=0).columns[1:].tolist()
        print(f'raw-count columns in the same order as the series matrix: {counts_cols == s["GSM"].tolist()}')
        assert counts_cols == s['GSM'].tolist()
    return s


def stratum_of(value, strata):
    t = str(value).lower()
    for prefix, name in strata['prefixes'].items():
        if t.startswith(prefix):
            return name
    return t


def samples_from_csv(cc, spec):
    groups = spec['groups']
    col = spec['group_column']
    m = pd.read_csv(spec['file'], **spec.get('read_csv', {}))
    print(f'metadata: {len(m)} samples; {col}: {m[col].value_counts().to_dict()}')
    s = m[m[col].isin(groups.values())].copy()
    s['group'] = s[col].map({v: k for k, v in groups.items()})
    st = spec.get('strata')
    if st:
        s[st['name']] = s[st['from_column']].map(lambda v: stratum_of(v, st))
    s = s.rename(columns={spec['sample_column']: 'GSM'})[spec['columns']]
    if st:
        print(f'samples by group and {st["name"]}:')
        print(pd.crosstab(s[st['name']], s['group'], margins=True).to_string())
    if 'subject_id' in s:
        print('subjects: ' + ', '.join(f'{g} {s[s.group == g].subject_id.nunique()}' for g in groups))
    return s


def load_samples(cc):
    """(samples table with a 'group' column, log2(TPM + 1) matrix of those samples)."""
    spec = cc.SAMPLES
    fmt = spec['format']
    if fmt == 'samples_table':
        s = samples_from_table(cc, spec)
    elif fmt == 'series_matrix':
        s = samples_from_series_matrix(cc, spec)
    elif fmt == 'csv':
        s = samples_from_csv(cc, spec)
    else:
        raise ValueError(f"SAMPLES['format'] = {fmt!r}: expected 'samples_table', 'series_matrix' or 'csv'")
    x = pd.read_csv(cc.NCBI_TPM, sep='\t', index_col=0)
    x.index = x.index.astype(int)
    print(f'TPM matrix: {x.shape[0]} genes x {x.shape[1]} samples; '
          f'column sums {x.sum().min():.0f}-{x.sum().max():.0f}')
    missing = set(s['GSM']) - set(x.columns)
    extra = set(x.columns) - set(s['GSM'])
    print(f'samples used {len(s)}; not in TPM {len(missing)}; TPM not in the samples used {len(extra)}')
    assert not missing
    if not spec.get('subset_of_tpm', False):
        assert not extra
    x = np.log2(x[s['GSM']] + 1)
    return s.reset_index(drop=True), x


# ------------------------------------------------------------- similarity

def load_tacle():
    rows = []
    for t in TACLE:
        d = pd.read_csv(BC_DISEASE_DATA / t / f'cell_line_{t}_tacle.csv')
        rows.append(pd.DataFrame({'tumour': t, 'cor': d['cor'].values}))
    return pd.concat(rows)


def per_gene_set(x, s, ccle, gs, lincs_models, groups, strata):
    genes = cls.marker_genes(ccle, gs, reference_models=lincs_models)
    cor, n = cls.similarity(x, ccle, genes)
    out = pd.DataFrame(index=ccle.index)
    out['gene_set'] = gs
    out['n_genes'] = n
    out['median_cor_all'] = cor.median(axis=0)
    for g in groups:
        idx = s.loc[s['group'] == g, 'GSM']
        out[f'median_cor_{g}'] = cor.loc[idx].median(axis=0)
        prof = x.loc[:, idx].median(axis=1).to_frame(f'median_{g}')
        out[f'cor_median_profile_{g}'] = cls.similarity(prof, ccle, genes)[0].iloc[0]
    if strata and gs == PRIMARY:
        for r in strata['classes']:
            out[f'median_cor_{r}'] = cor.loc[s.loc[s[strata['name']] == r, 'GSM']].median(axis=0)
    lin = cor[lincs_models]
    best, fdr = cls.compare_with_best(lin)
    out['rank_lincs'] = lin.median(axis=0).rank(ascending=False)
    out['fdr_vs_best_lincs'] = fdr
    print(f'{gs}: {n} genes; best LINCS line {best}')
    return out, cor


def similarity_tables(x, s, groups, strata, ccle, models, mapping):
    """Similarity table of every gene set and per-sample correlations of the primary set."""
    lincs_models = mapping.loc[mapping.has_expression, 'ModelID'].tolist()
    meta = models.set_index('ModelID').loc[ccle.index, ['CellLineName', 'OncotreeLineage']]
    m2l = mapping[mapping.has_expression].set_index('ModelID')
    tabs, cor_primary = [], None
    for gs in cls.GENE_SETS:
        out, cor = per_gene_set(x, s, ccle, gs, lincs_models, groups, strata)
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
    if strata:
        cols += [f'median_cor_{r}' for r in strata['classes']]
    tab = tab[cols].sort_values(['gene_set', 'median_cor_all'], ascending=[True, False])
    return tab, cor_primary


def external_similarity(x, s, groups, strata, ccle, w):
    """
    Per-sample correlations (primary gene set) of the LINCS lines weighted from an external profile
    (rows of w with source other than 'DepMap'), and their summary table in the columns of the
    similarity table (OncotreeLineage = 'HPA' or 'GEO'; ModelID = profile name).
    """
    ext = w[w['source'] != 'DepMap']
    if ext.empty:
        return None, None
    symbols = cls.ccle_gene_symbols()
    allp = {**cls.EXTERNAL_PROFILES, **cls.EXTERNAL_PROFILES_OPTIONAL}
    hpa_names = {lid: allp[lid][0] for lid in ext['lincs_cell_id'] if lid in allp}
    parts, meta = [], []
    if hpa_names:
        hpa = cls.load_hpa_cell_lines(set(hpa_names.values()), symbols)
        for lid, h in hpa_names.items():
            parts.append(hpa.loc[[h]])
            meta.append({'ModelID': h, 'lincs_cell_id': lid, 'OncotreeLineage': 'HPA'})
    geo_ids = {lid: cls.GEO_PROFILES[lid][0] for lid in ext['lincs_cell_id'] if lid in cls.GEO_PROFILES}
    if geo_ids:
        geo = cls.load_geo_profiles({lid: cls.GEO_PROFILES[lid] for lid in geo_ids})
        for lid, g in geo_ids.items():
            parts.append(geo.loc[[g]])
            meta.append({'ModelID': g, 'lincs_cell_id': lid, 'OncotreeLineage': 'GEO'})
    genes = cls.marker_genes(ccle, PRIMARY)
    cors, ns = [], {}
    for prof in parts:  # one source block at a time: each over the marker genes it has
        c, n = cls.similarity(x, prof, genes)
        cors.append(c)
        ns[prof.index[0]] = n
    cor = pd.concat(cors, axis=1)
    out = pd.DataFrame(meta)
    out['CellLineName'] = out['ModelID']
    out['gene_set'] = PRIMARY
    out['n_genes'] = out['ModelID'].map(ns)
    out['median_cor_all'] = cor.median(axis=0).reindex(out['ModelID']).values
    for g in groups:
        idx = s.loc[s['group'] == g, 'GSM']
        out[f'median_cor_{g}'] = cor.loc[idx].median(axis=0).reindex(out['ModelID']).values
    if strata:
        for c in strata['classes']:
            out[f'median_cor_{c}'] = cor.loc[s.loc[s[strata['name']] == c, 'GSM']].median(axis=0).reindex(out['ModelID']).values
    return out, cor


def reference_scale(ccle, lincs_models):
    """Quartiles of cell line vs cell line correlations on the primary set: random CCLE pairs, LINCS pairs."""
    x = ccle[cls.marker_genes(ccle, PRIMARY)].values
    idx = np.random.default_rng(REF_SEED).choice(len(x), 2 * REF_N, replace=False)
    rnd = cls.spearman_rows(x[idx[:REF_N]], x[idx[REF_N:]]).ravel()
    lin = cls.spearman_rows(x[ccle.index.get_indexer(lincs_models)], x[ccle.index.get_indexer(lincs_models)])
    lin = lin[np.triu_indices(len(lincs_models), 1)]
    rows = []
    for name, v in [(f'random CCLE pairs ({REF_N} x {REF_N})', rnd), ('LINCS line pairs', lin)]:
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        rows.append({'comparison': name, 'n_pairs': len(v), 'q1': q1, 'median': med, 'q3': q3})
    return pd.DataFrame(rows)


def top_cell_lines(w, ref):
    """LINCS lines with a weight above each cutoff of the random-pair distribution."""
    r = ref.iloc[0]
    rows = []
    for cut, q in TOP_CUTOFFS.items():
        top = w[w['cor'] > r[q]].sort_values('cor', ascending=False)
        if top.empty:
            rows.append({'cutoff': cut, 'threshold': r[q], 'lincs_cell_id': None, 'cor': None})
        rows += [{'cutoff': cut, 'threshold': r[q], 'lincs_cell_id': c, 'cor': v}
                 for c, v in zip(top['lincs_cell_id'], top['cor'])]
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- figure

def panel(cc, tab, cor_primary, tacle, s, groups, strata, ext_tab=None, cor_ext=None):
    P = cc.P
    F = figures(cc)
    BLUE, RED = colours()
    case = list(groups)[0]
    glabel = F.get('group_labels', {})
    box_alpha = .4
    p = tab[tab['gene_set'] == PRIMARY]
    lin = p[p['lincs_cell_id'].notna()]
    if ext_tab is not None:  # LINCS lines weighted from an external profile (HPA, GEO)
        lin = pd.concat([lin, ext_tab[[c for c in ext_tab.columns if c in lin.columns]]], ignore_index=True)
        cor_primary = pd.concat([cor_primary, cor_ext], axis=1)
    lin = lin.sort_values('median_cor_all', ascending=False)
    if strata:
        fig = plt.figure(figsize=(16, 15))
        gsp = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.75])
        ax = {k: fig.add_subplot(v) for k, v in
              {'A': gsp[0, :], 'B': gsp[1, 0], 'C': gsp[1, 1], 'D': gsp[2, 0], 'E': gsp[2, 1]}.items()}
    else:
        fig, a2 = plt.subplots(2, 2, figsize=(16, 12))
        ax = {'A': a2[0, 0], 'B': a2[0, 1], 'C': a2[1, 0], 'D': a2[1, 1]}

    a = ax['A']
    pos = np.arange(1, len(lin) + 1)
    for g, col, dx in ((case, RED, -0.2), ('control', BLUE, 0.2)):
        idx = s.loc[s['group'] == g, 'GSM']
        data = [cor_primary.loc[idx, m].values for m in lin['ModelID']]
        a.boxplot(data, positions=pos + dx, widths=.35, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=col, alpha=box_alpha, edgecolor=col),
                  whiskerprops=dict(color=col), capprops=dict(color=col),
                  medianprops=dict(color='k', linewidth=1), manage_ticks=False)
        a.plot([], [], 's', color=col, alpha=box_alpha, label=f'{glabel.get(g, g)} samples (n={len(idx)})')
    a.set_xticks(pos)
    a.set_xlim(0.3, len(lin) + 0.7)
    a.set_xticklabels([f'{c} ({l})' for c, l in zip(lin['lincs_cell_id'], lin['OncotreeLineage'])],
                      rotation=90, fontsize=7)
    a.set_ylabel(f'Spearman correlation with {F["dataset"]} samples')
    a.set_title(f'A. LINCS lines with a weight (CCLE, HPA, GEO), per-sample correlations\n'
                f'({PRIMARY}, {int(p.n_genes.iloc[0])} genes'
                + (f'; HPA/GEO profiles: the {", ".join(str(n) for n in sorted(ext_tab["n_genes"].unique(), reverse=True))} '
                   f'of them they contain' if ext_tab is not None else '') + ')')
    a.legend(loc='best')

    a = ax['B']
    bins = np.linspace(min(p['median_cor_all'].min(), tacle['cor'].min()) - .02,
                       max(p['median_cor_all'].max(), tacle['cor'].max()) + .02, 50)
    a.hist(p['median_cor_all'], bins=bins, color='0.6', alpha=.7,
           label=f'{F["disease"]}, all CCLE lines (n={len(p)})')
    a.hist(lin['median_cor_all'], bins=bins, color=RED, alpha=.6,
           label=f'{F["disease"]}, LINCS lines (n={len(lin)})')
    for t, ls in zip(TACLE, ['-', '--', ':', '-.']):
        v = tacle.loc[tacle['tumour'] == t, 'cor']
        a.hist(v, bins=bins, histtype='step', color=BLUE, ls=ls, lw=1.4,
               label=f'Chen 2017 {t}, all CCLE lines (n={len(v)})')
    a.set_xlabel('median Spearman correlation per cell line')
    a.set_ylabel('cell lines')
    a.set_title(f'B. Similarity scale: {F["disease"]} vs Chen 2017 tumour references')
    a.legend(fontsize=8)

    a = ax['C']
    a.scatter(lin['median_cor_control'], lin[f'median_cor_{case}'], color=RED, s=18)
    for _, r in lin.head(8).iterrows():
        a.annotate(r['lincs_cell_id'], (r['median_cor_control'], r[f'median_cor_{case}']), fontsize=7)
    v = lin[['median_cor_control', f'median_cor_{case}']].values
    a.plot([v.min(), v.max()], [v.min(), v.max()], color='k', lw=.6)
    a.set_xlabel(f'median correlation, {F["control_axis"]}')
    a.set_ylabel(f'median correlation, {F["case_axis"]}')
    a.set_title(f'C. LINCS lines: {F["panel_c_title"]}')

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

    if strata:
        a = ax['E']
        n_reg = s[strata['name']].value_counts()
        for r, mk in zip(strata['classes'], ['o', 's', '^', 'D']):
            a.scatter(lin['median_cor_all'], lin[f'median_cor_{r}'], s=16, marker=mk,
                      label=f'{r.replace("_", " ")} (n={n_reg.get(r, 0)})')
        v = lin['median_cor_all']
        a.plot([v.min(), v.max()], [v.min(), v.max()], color='k', lw=.6)
        regs = [f'median_cor_{r}' for r in strata['classes']]
        for _, r in lin.head(3).iterrows():  # the three lines most similar to the samples
            a.annotate(r['lincs_cell_id'], (r['median_cor_all'], r[regs].max()),
                       xytext=(3, 3), textcoords='offset points', fontsize=7)
        a.set_xlabel('median correlation, all samples')
        a.set_ylabel(f'median correlation, samples of one {strata["short_label"]}')
        a.set_title(f'E. LINCS lines: median correlation per {strata["label"]}')
        a.legend(fontsize=8)

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_ccle_cell_line_similarity_panel.png')


# ------------------------------------------------------- external weights

def instances_used(md):
    """Instances sRGES can use: dose > 0 and time 6 h or 24 h."""
    dose = pd.to_numeric(md['pert_dose'], errors='coerce')
    time = pd.to_numeric(md['pert_time'], errors='coerce')
    return md[(dose > 0) & time.isin([6, 24])]


def paired_summary(name, d):
    """Agreement of the HPA and DepMap weights of the same cell lines (rows of d)."""
    b, a = np.polyfit(d['cor_depmap'], d['cor_hpa'], 1) if len(d) > 1 else (np.nan, np.nan)
    diff = d['cor_hpa'] - d['cor_depmap']
    return {'set': name, 'n_lines': len(d),
            'spearman': spearmanr(d['cor_depmap'], d['cor_hpa']).statistic if len(d) > 2 else np.nan,
            'pearson': d['cor_depmap'].corr(d['cor_hpa']) if len(d) > 2 else np.nan,
            'median_diff': diff.median(), 'median_abs_diff': diff.abs().median(),
            'median_gene_spearman': d['gene_spearman'].median() if 'gene_spearman' in d else np.nan,
            'p95_abs_diff': diff.abs().quantile(.95), 'max_abs_diff': diff.abs().max(),
            'slope': b, 'intercept': a}


def diagnose_lines(models, mapping):
    """One row per LINCS line without DepMap expression: LINCS annotation, DepMap match by name and by
    Cellosaurus accession (Model.RRID), and the number of RNA-seq profiles of that model in each release."""
    info = pd.read_csv(cls.LINCS_CELL_INFO).set_index('cell_id')
    rna = {}
    for r in cls.DEPMAP_PROFILE_RELEASES:
        pr = pd.read_csv(cls.DEPMAP_PROFILES_DIR / f'OmicsProfiles_{r}.csv')
        rna[r] = pr[pr['Datatype'] == 'rna'].groupby('ModelID').size()
    by_rrid = models.dropna(subset=['RRID']).drop_duplicates('RRID').set_index('RRID')
    rows = []
    for lid in mapping.loc[~mapping['has_expression'], 'lincs_cell_id']:
        i = info.loc[lid] if lid in info.index else None
        acc = cls.CELLOSAURUS_ID.get(lid)
        mid = mapping.loc[mapping['lincs_cell_id'] == lid, 'ModelID'].iloc[0]
        if (not isinstance(mid, str)) and acc in by_rrid.index:
            mid = by_rrid.loc[acc, 'ModelID']
        row = {'lincs_cell_id': lid, 'lincs_cell_type': i['cell_type'] if i is not None else '',
               'lincs_histology': i['cell_histology'] if i is not None else '',
               'lincs_lineage': i['cell_lineage'] if i is not None else '',
               'lincs_aliases': i['cell_aliases'] if i is not None else '',
               'lincs_source': i['cell_source'] if i is not None else '',
               'parent_cellosaurus': cls.PARENT_CELLOSAURUS_ID.get(lid, ''),
               'cellosaurus': acc or '', 'depmap_model': mid if isinstance(mid, str) else ''}
        for r in cls.DEPMAP_PROFILE_RELEASES:
            row[f'rna_profiles_{r}'] = int(rna[r].get(mid, 0)) if isinstance(mid, str) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def candidate_weights(cc, x, s, ccle, models, mapping, optional=()):
    """
    Candidate weights of the LINCS lines without DepMap expression (docstring: 'Candidate
    weights'); writes the *_candidate* tables and returns the candidate weight table.
    x: log2(TPM + 1) of the samples (genes x samples), s: samples table.
    """
    from case_study.known_drugs import key
    from loader import read_lincs_metadata
    P = cc.P
    case = list(cc.SAMPLES['groups'])[0]
    case_idx = s.loc[s['group'] == case, 'GSM']
    genes = cls.marker_genes(ccle, PRIMARY)

    def weights_of(profiles):
        cor, n = cls.similarity(x[case_idx], profiles, genes)
        return cor.median(axis=0), n

    w_ccle, n_ccle = weights_of(ccle)
    symbols = cls.ccle_gene_symbols()
    hpa = cls.load_hpa_cell_lines(None, symbols)
    w_hpa, n_hpa = weights_of(hpa)
    print(f'external weights: DepMap {len(ccle)} lines over {n_ccle} marker genes; '
          f'HPA {len(hpa)} lines over {n_hpa} marker genes')
    geo = cls.load_geo_profiles()
    w_geo, n_geo = weights_of(geo)
    print(f'GEO {len(geo)} profiles over {n_geo} marker genes')

    # DepMap weights of the LINCS lines, checked against the weights file of the step
    m = mapping[mapping['has_expression']].set_index('ModelID')
    cur = pd.DataFrame({'lincs_cell_id': m['lincs_cell_id'], 'cor': w_ccle.reindex(m.index)})
    cur['source'] = 'DepMap'
    f = P.oout_dir / f'{P.slug}_cell_line_weights.tsv'
    recomputation = np.nan
    if f.exists():
        old = pd.read_csv(f, sep='\t').set_index('lincs_cell_id')['cor']
        recomputation = (cur.set_index('lincs_cell_id')['cor'] - old.reindex(cur['lincs_cell_id'].values)).abs().max()
        print(f'DepMap weights recomputed vs {f.name}: max |difference| {recomputation:.2e}')

    # paired lines: DepMap lines with expression whose stripped name is the stripped name of an HPA line
    hpa_by_name = {cls._strip(n): n for n in hpa.index}
    mod = models.set_index('ModelID')
    rows = []
    for mid in ccle.index:
        h = hpa_by_name.get(mod.loc[mid, 'stripped'])
        if h is not None:
            rows.append({'kind': 'shared_line', 'name': mod.loc[mid, 'CellLineName'], 'hpa_line': h,
                         'OncotreeLineage': mod.loc[mid, 'OncotreeLineage'], 'n_genes_depmap': n_ccle,
                         'n_genes_hpa': n_hpa, 'cor_depmap': w_ccle[mid], 'cor_hpa': w_hpa[h]})
    pair = pd.DataFrame(rows)
    g = [k for k in genes if k in ccle.columns and k in hpa.columns]
    mid_of = {n: m for m, n in zip(mod.index, mod['CellLineName'])}
    a = ccle.loc[[mid_of[n] for n in pair['name']], g].values
    b = hpa.loc[pair['hpa_line'], g].values
    ra, rb = rankdata(a, axis=1), rankdata(b, axis=1)
    ra = (ra - ra.mean(1, keepdims=True)) / ra.std(1, keepdims=True)
    rb = (rb - rb.mean(1, keepdims=True)) / rb.std(1, keepdims=True)
    pair['gene_spearman'] = (ra * rb).sum(1) / len(g)
    pair['diff'] = pair['cor_hpa'] - pair['cor_depmap']
    summ = [paired_summary('all', pair)]
    for lin in sorted({v[2] for v in cls.EXTERNAL_PROFILES.values() if v[2]}):
        summ.append(paired_summary(f'lineage {lin}', pair[pair['OncotreeLineage'] == lin]))
    only = w_hpa.drop(index=pair['hpa_line'])
    for name, v in (('scale: HPA weights of the paired lines', pair['cor_hpa']),
                    ('scale: HPA weights of the HPA lines not in DepMap', only)):
        summ.append({'set': name, 'n_lines': len(v), 'mean_cor': v.mean(), 'sd_cor': v.std()})
    summ = pd.DataFrame(summ)
    summ['recomputation_max_abs_diff_depmap_weights'] = np.nan
    summ.loc[0, 'recomputation_max_abs_diff_depmap_weights'] = recomputation
    print('\nHPA vs DepMap weights of the paired lines:')
    print(summ.round(3).to_string(index=False))

    diag = diagnose_lines(models, mapping)
    diag.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights_candidate_diagnosis.tsv', sep='\t', index=False)
    print('\nLINCS lines without DepMap expression:')
    print(diag.to_string(index=False))

    # tissue references against the distribution of DepMap-line weights
    tis = cls.load_hpa_tissues(cls.TISSUE_REFERENCES, symbols)
    w_tis, _ = weights_of(tis)
    pct = lambda v: 100 * (w_ccle < v).mean()
    trows = [{'kind': 'tissue_reference', 'name': t, 'cor_hpa': v, 'percentile_among_depmap_lines': pct(v)}
             for t, v in w_tis.items()]
    # candidate weights
    new, opt = [], []
    allp = {**cls.EXTERNAL_PROFILES, **cls.EXTERNAL_PROFILES_OPTIONAL}
    for lid, (h, match, _) in allp.items():
        r = mapping[mapping['lincs_cell_id'] == lid]
        if len(r) and not bool(r['has_expression'].iloc[0]) and h in w_hpa.index:
            row = {'lincs_cell_id': lid, 'cor': w_hpa[h], 'source': f'HPA {h}, {match}'}
            (new if lid in cls.EXTERNAL_PROFILES or lid in optional else opt).append(row)
    for lid, (name, match, _, _) in cls.GEO_PROFILES.items():
        r = mapping[mapping['lincs_cell_id'] == lid]
        if len(r) and not bool(r['has_expression'].iloc[0]):
            new.append({'lincs_cell_id': lid, 'cor': w_geo[name], 'source': f'GEO {name}, {match}'})
    cand = pd.concat([cur.reset_index(drop=True), pd.DataFrame(new)], ignore_index=True)
    cand.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights_candidate.tsv', sep='\t', index=False)
    print(f'\ncandidate weights: {len(cur)} DepMap + {len(new)} external (HPA, GEO) = {len(cand)} LINCS lines '
          f'(of {len(mapping)})')

    # effect
    md = read_lincs_metadata(dtype={'id': str})
    md = instances_used(md[['id', 'pert_iname', 'cell_id', 'pert_dose', 'pert_time']].copy())
    md['drug_key'] = key(md['pert_iname'])
    cur_set, cand_set = set(cur['lincs_cell_id']), set(cand['lincs_cell_id'])
    opt_set = cand_set | {o['lincs_cell_id'] for o in opt}
    optd = {o['lincs_cell_id']: o for o in opt}
    erows = []
    for lid in mapping.loc[~mapping['has_expression'], 'lincs_cell_id']:
        c = cand[cand['lincs_cell_id'] == lid]
        erows.append({'kind': 'line', 'name': lid, 'n_instances': int((md['cell_id'] == lid).sum()),
                      'candidate_cor': c['cor'].iloc[0] if len(c) else np.nan,
                      'candidate_source': c['source'].iloc[0] if len(c) else '',
                      'optional_cor': optd[lid]['cor'] if lid in optd else np.nan,
                      'optional_source': optd[lid]['source'] if lid in optd else '',
                      'optional_percentile_among_depmap_lines': pct(optd[lid]['cor']) if lid in optd else np.nan,
                      'percentile_among_depmap_lines': pct(c['cor'].iloc[0]) if len(c) else np.nan})
    kd = pd.DataFrame({'drug': list(cc.KNOWN_DRUGS)}) if hasattr(cc, 'KNOWN_DRUGS') else pd.DataFrame({'drug': []})
    kd['drug_key'] = key(kd['drug'])
    for _, r in kd.iterrows():
        d = md[md['drug_key'] == r['drug_key']]
        c_now, c_new = d['cell_id'].isin(cur_set), d['cell_id'].isin(cand_set)
        erows.append({'kind': 'known_drug', 'name': r['drug'], 'n_instances': len(d),
                      'instances_weighted_now': int(c_now.sum()), 'instances_weighted_candidate': int(c_new.sum()),
                      'instances_unweighted_candidate': int((~c_new).sum()),
                      'instances_weighted_with_optional': int(d['cell_id'].isin(opt_set).sum()),
                      'unweighted_now_by_line': ';'.join(f'{c}:{n}' for c, n in
                                                         d.loc[~c_now, 'cell_id'].value_counts().items()),
                      'unweighted_candidate_by_line': ';'.join(f'{c}:{n}' for c, n in
                                                               d.loc[~c_new, 'cell_id'].value_counts().items()),
                      'lines_weighted_now': d.loc[c_now, 'cell_id'].nunique(),
                      'lines_weighted_candidate': d.loc[c_new, 'cell_id'].nunique()})
    g = md.groupby('drug_key')['cell_id']
    now = g.apply(lambda v: v.isin(cur_set).sum())
    new_ = g.apply(lambda v: v.isin(cand_set).sum())
    erows.append({'kind': 'all_compounds', 'name': 'all compounds', 'n_instances': len(md),
                  'instances_weighted_now': int(md['cell_id'].isin(cur_set).sum()),
                  'instances_weighted_candidate': int(md['cell_id'].isin(cand_set).sum()),
                  'instances_unweighted_candidate': int((~md['cell_id'].isin(cand_set)).sum()),
                  'n_compounds': len(now), 'compounds_weighted_now': int((now > 0).sum()),
                  'compounds_weighted_candidate': int((new_ > 0).sum()),
                  'instances_weighted_with_optional': int(md['cell_id'].isin(opt_set).sum()),
                  'compounds_weighted_with_optional': int((g.apply(lambda v: v.isin(opt_set).sum()) > 0).sum())})
    eff = pd.DataFrame(erows)
    eff.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights_candidate_effect.tsv', sep='\t', index=False)
    print('\neffect of the candidate weights:')
    print(eff.round(3).to_string(index=False))

    pd.concat([pair, pd.DataFrame(trows)], ignore_index=True).to_csv(
        P.oout_dir / f'{P.slug}_cell_line_weights_candidate_check.tsv', sep='\t', index=False)
    summ.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights_candidate_check_summary.tsv', sep='\t', index=False)
    return cand


def run_candidate(cc, optional=()):
    """Only the candidate weights and their checks (no other output of the step is written)."""
    s, x = load_samples(cc)
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    lincs = pd.read_csv(LINCS_METADATA_PATH)
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    candidate_weights(cc, x, s, ccle, models, mapping, optional)


# -------------------------------------------------------------------- run

def run(cc, external_weights=False, external_optional=()):
    """Step 'weights'; external_weights=True also writes the candidate weights and their checks."""
    P = cc.P
    groups = cc.SAMPLES['groups']
    strata = cc.SAMPLES.get('strata')
    case = list(groups)[0]
    dataset = figures(cc)['dataset']

    s, x = load_samples(cc)
    ccle = cls.load_ccle_expression()
    models = cls.load_ccle_models()
    print(f'CCLE matrix: {ccle.shape[0]} lines x {ccle.shape[1]} genes; '
          f'genes shared with {dataset}: {len(set(ccle.columns) & set(x.index))}')

    lincs = pd.read_csv(LINCS_METADATA_PATH)
    n_inst = lincs.groupby('cell_id').size()
    mapping = cls.map_lincs_to_ccle(lincs['cell_id'], models, set(ccle.index))
    mapping['n_lincs_instances'] = mapping['lincs_cell_id'].map(n_inst)
    print(f'LINCS lines {len(mapping)}; with DepMap expression {int(mapping.has_expression.sum())}; '
          f'without: {", ".join(mapping.loc[~mapping.has_expression, "lincs_cell_id"])} '
          f'({int(mapping.loc[~mapping.has_expression, "n_lincs_instances"].sum())} instances)')
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(P.oout_dir / f'{P.slug}_lincs_ccle_mapping.tsv', sep='\t', index=False)
    if cc.SAMPLES['format'] != 'samples_table':
        s.to_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', index=False)

    tab, cor_primary = similarity_tables(x, s, groups, strata, ccle, models, mapping)
    tab.to_csv(P.oout_dir / f'{P.slug}_ccle_cell_line_similarity.tsv', sep='\t', index=False)

    w = tab[(tab.gene_set == PRIMARY) & tab.lincs_cell_id.notna()]
    w = w[['lincs_cell_id', f'median_cor_{case}']].rename(columns={f'median_cor_{case}': 'cor'})
    w['source'] = 'DepMap'
    if external_weights:
        print('\n=== candidate weights from HPA and GEO')
        w = candidate_weights(cc, x, s, ccle, models, mapping, external_optional)[['lincs_cell_id', 'cor', 'source']]
    w.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t', index=False)
    print(f'sRGES weights written for {len(w)} LINCS lines (range {w.cor.min():.3f}-{w.cor.max():.3f})')
    ext_tab, cor_ext = external_similarity(x, s, groups, strata, ccle, w)
    if ext_tab is not None:
        ext_tab.to_csv(P.oout_dir / f'{P.slug}_external_cell_line_similarity.tsv', sep='\t', index=False)

    ref = reference_scale(ccle, mapping.loc[mapping.has_expression, 'ModelID'].tolist())
    ref.to_csv(P.oout_dir / f'{P.slug}_ccle_reference_scale.tsv', sep='\t', index=False)
    print('\ncell line vs cell line reference (primary set):')
    print(ref.round(3).to_string(index=False))
    top = top_cell_lines(w, ref)
    top.to_csv(P.oout_dir / f'{P.slug}_top_cell_lines.tsv', sep='\t', index=False)
    for cut, g in top.groupby('cutoff', sort=False):
        lines = g.dropna(subset=['lincs_cell_id'])
        print(f'top cell lines, cutoff {cut} ({g.threshold.iloc[0]:.3f}): {len(lines)}'
              + (f': {", ".join(f"{c} {v:.3f}" for c, v in zip(lines.lincs_cell_id, lines.cor))}' if len(lines) else ''))

    tacle = load_tacle()
    panel(cc, tab, cor_primary, tacle, s, groups, strata, ext_tab, cor_ext)

    p = tab[tab.gene_set == PRIMARY]
    lin = p[p.lincs_cell_id.notna()]
    d = lin[f'median_cor_{case}'] - lin['median_cor_control']
    print(f'\n{case} - control median correlation per LINCS line: min {d.min():+.4f}, max {d.max():+.4f}; '
          f'{case} higher in {(d > 0).sum()}/{len(d)} lines')
    print('\nprimary set, median per-line correlation (all samples):')
    print(f'  all CCLE lines: median {p.median_cor_all.median():.3f}, max {p.median_cor_all.max():.3f} '
          f'({p.loc[p.median_cor_all.idxmax(), "CellLineName"]})')
    print(f'  LINCS lines:    median {lin.median_cor_all.median():.3f}, max {lin.median_cor_all.max():.3f}')
    for t in TACLE:
        v = tacle.loc[tacle.tumour == t, 'cor']
        print(f'  Chen 2017 {t}: median {v.median():.3f}, max {v.max():.3f}')
    for gs in cls.GENE_SETS:
        q = tab[(tab.gene_set == gs) & tab.lincs_cell_id.notna()]
        print(f'  {gs}: LINCS median_cor_{case} range {q[f"median_cor_{case}"].min():.3f}-'
              f'{q[f"median_cor_{case}"].max():.3f}')
    extra = [f'median_cor_{r}' for r in strata['classes']] if strata else []
    print('\nLINCS lines (primary set):')
    print(lin[['lincs_cell_id', 'CellLineName', 'OncotreeLineage', 'n_lincs_instances',
               'median_cor_all', f'median_cor_{case}', 'median_cor_control'] + extra +
              ['fdr_vs_best_lincs']].round(3).to_string(index=False))
    print('\ntop 15 CCLE lines overall (primary set):')
    print(p.head(15)[['CellLineName', 'OncotreeLineage', 'lincs_cell_id', 'median_cor_all']]
          .round(3).to_string(index=False))
