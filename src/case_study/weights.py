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
- Reference scale: the median tumour-cell line correlations of Chen 2017
  (BRCA, LIHC, COAD, ER+ BRCA; cell_line_<T>_tacle.csv, CCLE 2017 and TCGA,
  5,000 IQR genes) over all the CCLE lines they cover.

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
    set; LINCS lines with DepMap expression only (lines without a weight are
    dropped by sRGES.compute_sRGES).
- P.oout_dir/<slug>_lincs_ccle_mapping.tsv
    the LINCS cell_id with their DepMap match and n_lincs_instances.
- P.img_dir/<slug>_ccle_cell_line_similarity_panel.png
    A: per-sample correlations of the LINCS lines (primary set), case and
    control samples side by side, ordered by median over all samples;
    B: per-line median correlation, all CCLE vs LINCS lines vs the Chen 2017
    tumour references; C: case vs control median per LINCS line;
    D: agreement of the LINCS-line ordering between gene sets; with strata,
    E: median correlation of the LINCS lines per sample class.
    Text and style from cc.FIGURES.
"""

import gzip

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from conf import LINCS_METADATA_PATH, BC_DISEASE_DATA
import cell_line_similarity as cls
from case_study.common import figures, colours

PRIMARY = 'binchen2017_iqr5000'
TACLE = ['BRCA', 'LIHC', 'COAD', 'ER']


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


# ----------------------------------------------------------------- figure

def panel(cc, tab, cor_primary, tacle, s, groups, strata):
    P = cc.P
    F = figures(cc)
    BLUE, RED = colours(cc)
    case = list(groups)[0]
    glabel = F.get('group_labels', {})
    box_alpha = F.get('box_alpha', .4)
    p = tab[tab['gene_set'] == PRIMARY]
    lin = p[p['lincs_cell_id'].notna()].sort_values('median_cor_all', ascending=False)
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
    a.set_title(f'A. {F.get("panel_a_title", "LINCS lines")}, per-sample correlations '
                f'({PRIMARY}, {int(p.n_genes.iloc[0])} genes)')
    a.legend(loc=F.get('panel_a_legend', 'lower left'))

    a = ax['B']
    bins = np.linspace(min(p['median_cor_all'].min(), tacle['cor'].min()) - .02,
                       max(p['median_cor_all'].max(), tacle['cor'].max()) + .02, 50)
    a.hist(p['median_cor_all'], bins=bins, color='0.6', alpha=.7,
           label=f'{F["disease"]}, all CCLE lines (n={len(p)})')
    a.hist(lin['median_cor_all'], bins=bins, color=RED, alpha=F.get('lincs_hist_alpha', .6),
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
        a.set_xlabel('median correlation, all samples')
        a.set_ylabel(f'median correlation, samples of one {strata["short_label"]}')
        a.set_title(f'E. LINCS lines: median correlation per {strata["label"]}')
        a.legend(fontsize=8)

    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_ccle_cell_line_similarity_panel.png')


# -------------------------------------------------------------------- run

def run(cc):
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
    w.to_csv(P.oout_dir / f'{P.slug}_cell_line_weights.tsv', sep='\t', index=False)
    print(f'sRGES weights written for {len(w)} LINCS lines (range {w.cor.min():.3f}-{w.cor.max():.3f})')

    tacle = load_tacle()
    panel(cc, tab, cor_primary, tacle, s, groups, strata)

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
