#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Re-run the edgeR QLF differential expression of the TDP-43 iPSC motor-neuron
bulk RNA-seq from the raw counts, check it against the original edgeR tables,
and write a sample-level QC panel.

Design
------
- 32 samples, 8 groups x 4 replicates: EV, Y, T, N and EV_AKO, Y_AKO, T_AKO,
  N_AKO (EV = empty vector, Y = YFP, T = TDP-43:YFP, N = TDP-43dNLS:YFP;
  _AKO = ATXN2 knockout line). Group = sample name without the trailing
  replicate number.
- edgeR, as in the original analysis: keep genes with CPM > STAN_MIN_CPM in
  >= STAN_MIN_SAMPLES samples, TMM normalisation, estimateDisp on the
  8-group no-intercept design, glmQLFit (legacy = TRUE), glmQLFTest per
  contrast, Benjamini-Hochberg FDR. Rows are the rows of the counts file
  (one row per ENSEMBL x SYMBOL pair, as in the original tables).
- Contrasts: the nine of the original analysis (TvN, T_AKOvN_AKO, TvY,
  T_AKOvY_AKO, TvEV, T_AKOvEV_AKO, TvAKO, YvAKO, NvAKO) plus NvY, NvEV, YvEV.
  logFC is numerator minus denominator (e.g. TvY > 0 = higher in T).
- Known-answer check: for every contrast present as a sheet of the original
  workbook, genes are matched on ENSEMBL|SYMBOL (keys duplicated on either
  side are excluded) and logFC, PValue and the FDR <= STAN_FDR gene sets are
  compared.
- R is called through Rscript (edgeR required); intermediate files live in a
  temporary directory that is removed at the end.

Inputs (paths from conf.py)
---------------------------
- conf.STAN_COUNTS_FILE : raw counts, columns SYMBOL, ENSEMBL, one per sample.
- conf.STAN_EDGER_XLSX  : original edgeR QLF workbook, one sheet per contrast.

Outputs
-------
- STAN_OOUT_DIR/stan_als_tdp43_edgeR_qlf_contrasts.tsv
    one row per tested gene row; SYMBOL, ENSEMBL, logCPM and, per contrast,
    <contrast>_logFC, <contrast>_PValue, <contrast>_FDR.
- STAN_OOUT_DIR/stan_als_tdp43_edgeR_known_answer.tsv
    per contrast: genes compared, logFC correlation and max/median absolute
    difference, Spearman correlation of log10 PValue, number of FDR genes on
    each side and their Jaccard index.
- STAN_OOUT_DIR/stan_als_tdp43_sample_qc.tsv
    per sample: group, library size, CPM of YFP, TARDBP, ATXN2 and a few
    marker genes.
- STAN_IMG_DIR/stan_als_tdp43_qc_panel.png
    A: PCA of TMM log-CPM; B: marker-gene CPM per group; C: correlation of
    logFC between contrasts; D: number of FDR genes per contrast.

Usage
-----
    python 01_edger_qlf_and_qc.py   (from any directory; needs Rscript + edgeR)
"""

import sys
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SRC_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SRC_DIR))

from conf import (STAN_SLUG, STAN_COUNTS_FILE, STAN_EDGER_XLSX, STAN_OOUT_DIR,
                  STAN_IMG_DIR, STAN_MIN_CPM, STAN_MIN_SAMPLES, STAN_FDR)

GROUPS = ['EV', 'Y', 'T', 'N', 'EV_AKO', 'Y_AKO', 'T_AKO', 'N_AKO']
CONTRASTS = {'TvN': ('T', 'N'), 'T_AKOvN_AKO': ('T_AKO', 'N_AKO'),
             'TvY': ('T', 'Y'), 'T_AKOvY_AKO': ('T_AKO', 'Y_AKO'),
             'TvEV': ('T', 'EV'), 'T_AKOvEV_AKO': ('T_AKO', 'EV_AKO'),
             'TvAKO': ('T', 'T_AKO'), 'YvAKO': ('Y', 'Y_AKO'), 'NvAKO': ('N', 'N_AKO'),
             'NvY': ('N', 'Y'), 'NvEV': ('N', 'EV'), 'YvEV': ('Y', 'EV')}
MARKERS = ['YFP', 'TARDBP', 'ATXN2', 'STMN2', 'UNC13A', 'MNX1', 'ISL1', 'CHAT',
           'MAP2', 'MKI67', 'PEG3', 'CHCHD2']
GROUP_COLORS = dict(zip(GROUPS, ['#7f7f7f', '#e3a21a', '#2ca02c', '#1f77b4',
                                 '#c7c7c7', '#f5d58a', '#98df8a', '#aec7e8']))

R_CODE = r'''
suppressMessages(library(edgeR))
a <- commandArgs(trailingOnly=TRUE); counts_file <- a[1]; tmp <- a[2]
min_cpm <- as.numeric(a[3]); min_samples <- as.numeric(a[4])
x <- read.delim(counts_file, check.names=FALSE)
cnt <- as.matrix(x[, -(1:2)])
grp <- factor(sub("[0-9]+$", "", colnames(cnt)),
              levels=c("EV","Y","T","N","EV_AKO","Y_AKO","T_AKO","N_AKO"))
y <- DGEList(cnt, group=grp, genes=x[, 1:2])
write.table(data.frame(sample=colnames(cnt), group=grp, lib_size=colSums(cnt)),
            file.path(tmp, "samples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cpm_all <- cpm(y)
write.table(data.frame(x[, 1:2], cpm_all, check.names=FALSE), file.path(tmp, "cpm_all.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
keep <- rowSums(cpm_all > min_cpm) >= min_samples
y <- y[keep, , keep.lib.sizes=FALSE]; y <- calcNormFactors(y)
write.table(data.frame(y$genes, cpm(y, log=TRUE), check.names=FALSE), file.path(tmp, "logcpm.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
d <- model.matrix(~0 + grp); colnames(d) <- levels(grp)
y <- estimateDisp(y, d)
fit <- glmQLFit(y, d, legacy=TRUE)
cons <- read.delim(file.path(tmp, "contrasts.tsv"), stringsAsFactors=FALSE)
for (i in seq_len(nrow(cons))) {
  cv <- setNames(rep(0, ncol(d)), colnames(d)); cv[cons$num[i]] <- 1; cv[cons$den[i]] <- -1
  tt <- topTags(glmQLFTest(fit, contrast=cv), n=Inf, sort.by="none")$table
  write.table(tt, file.path(tmp, paste0("tt_", cons$name[i], ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
}
cat("genes tested:", nrow(y), " common BCV:", sqrt(y$common.dispersion), "\n")
'''


def run_edger(tmp):
    pd.DataFrame([(k, v[0], v[1]) for k, v in CONTRASTS.items()],
                 columns=['name', 'num', 'den']).to_csv(tmp / 'contrasts.tsv', sep='\t', index=False)
    (tmp / 'qlf.R').write_text(R_CODE)
    subprocess.run(['Rscript', str(tmp / 'qlf.R'), str(STAN_COUNTS_FILE), str(tmp),
                    str(STAN_MIN_CPM), str(STAN_MIN_SAMPLES)], check=True)
    tabs = {k: pd.read_csv(tmp / f'tt_{k}.tsv', sep='\t') for k in CONTRASTS}
    first = next(iter(tabs.values()))
    wide = first[['SYMBOL', 'ENSEMBL', 'logCPM']].copy()
    for k, t in tabs.items():
        for c in ['logFC', 'PValue', 'FDR']:
            wide[f'{k}_{c}'] = t[c].values
    return wide


def _key(df):
    return df['ENSEMBL'].astype(str) + '|' + df['SYMBOL'].astype(str)


def known_answer(wide):
    sheets = pd.read_excel(STAN_EDGER_XLSX, sheet_name=None, usecols=range(7))
    ours = wide.assign(key=_key(wide))
    ours = ours[~ours['key'].duplicated(keep=False)].set_index('key')
    rows = []
    for name, ref in sheets.items():
        if name not in CONTRASTS:
            continue
        ref = ref.dropna(subset=['ENSEMBL', 'logFC']).assign(key=lambda d: _key(d))
        ref = ref[~ref['key'].duplicated(keep=False)].set_index('key')
        common = ref.index.intersection(ours.index)
        r, o = ref.loc[common], ours.loc[common]
        dfc = (r['logFC'] - o[f'{name}_logFC']).abs()
        a = set(common[r['FDR'] <= STAN_FDR]); b = set(common[o[f'{name}_FDR'] <= STAN_FDR])
        rows.append({'contrast': name, 'n_rows_original': len(ref), 'n_rows_rerun': len(ours),
                     'n_common': len(common),
                     'pearson_logFC': np.corrcoef(r['logFC'], o[f'{name}_logFC'])[0, 1],
                     'max_abs_diff_logFC': dfc.max(), 'median_abs_diff_logFC': dfc.median(),
                     'spearman_log10P': pd.Series(np.log10(r['PValue'].values)).corr(
                         pd.Series(np.log10(o[f'{name}_PValue'].values)), method='spearman'),
                     'n_FDR_original': int((ref['FDR'] <= STAN_FDR).sum()),
                     'n_FDR_rerun': int((ours[f'{name}_FDR'] <= STAN_FDR).sum()),
                     'jaccard_FDR_common': len(a & b) / len(a | b)})
    return pd.DataFrame(rows)


def sample_qc(tmp):
    s = pd.read_csv(tmp / 'samples.tsv', sep='\t')
    cp = pd.read_csv(tmp / 'cpm_all.tsv', sep='\t')
    cp['label'] = np.where(cp['ENSEMBL'] == 'YFP', 'YFP', cp['SYMBOL'])
    cp = cp[cp['label'].isin(MARKERS)].drop_duplicates('label').set_index('label')
    for m in MARKERS:
        s[f'CPM_{m}'] = cp.loc[m, s['sample']].values
    return s


def qc_panel(tmp, samples, wide):
    lc = pd.read_csv(tmp / 'logcpm.tsv', sep='\t').iloc[:, 2:]
    X = lc.T.values - lc.T.values.mean(axis=0)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pcs = U * S; ve = 100 * S ** 2 / (S ** 2).sum()
    grp = samples.set_index('sample').loc[lc.columns, 'group']

    fig, ax = plt.subplots(2, 2, figsize=(15, 12))
    a = ax[0, 0]
    for g in GROUPS:
        m = (grp == g).values
        a.scatter(pcs[m, 0], pcs[m, 1], s=60, color=GROUP_COLORS[g], edgecolor='k', label=g)
    a.set_xlabel(f'PC1 ({ve[0]:.1f}%)'); a.set_ylabel(f'PC2 ({ve[1]:.1f}%)')
    a.set_title('A. PCA, TMM log-CPM (tested genes)'); a.legend(ncol=2, fontsize=9)

    a = ax[0, 1]
    mk = ['YFP', 'TARDBP', 'ATXN2', 'PEG3', 'CHCHD2', 'MKI67', 'STMN2', 'UNC13A']
    means = samples.groupby('group')[[f'CPM_{m}' for m in mk]].mean().loc[GROUPS]
    w = 0.1
    for i, g in enumerate(GROUPS):
        a.bar(np.arange(len(mk)) + (i - 3.5) * w, np.log2(means.loc[g].values + 1), w,
              color=GROUP_COLORS[g], edgecolor='k', linewidth=.3, label=g)
    a.set_xticks(range(len(mk))); a.set_xticklabels(mk)
    a.set_ylabel('log2(mean CPM + 1)'); a.set_title('B. Marker genes per group'); a.legend(ncol=4, fontsize=8)

    a = ax[1, 0]
    names = list(CONTRASTS)
    C = wide[[f'{k}_logFC' for k in names]].corr().values
    im = a.imshow(C, cmap='RdBu_r', vmin=-1, vmax=1)
    a.set_xticks(range(len(names))); a.set_xticklabels(names, rotation=90)
    a.set_yticks(range(len(names))); a.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            a.text(j, i, f'{C[i, j]:.2f}', ha='center', va='center', fontsize=7)
    fig.colorbar(im, ax=a, fraction=.046); a.set_title('C. Pearson correlation of logFC between contrasts')

    a = ax[1, 1]
    up = [int(((wide[f'{k}_FDR'] <= STAN_FDR) & (wide[f'{k}_logFC'] > 0)).sum()) for k in names]
    dn = [int(((wide[f'{k}_FDR'] <= STAN_FDR) & (wide[f'{k}_logFC'] < 0)).sum()) for k in names]
    a.barh(names, up, color='#d62728', label='up'); a.barh(names, [-d for d in dn], color='#1f77b4', label='down')
    a.axvline(0, color='k', lw=.5); a.invert_yaxis(); a.legend(loc='lower left')
    a.set_xlabel(f'genes at FDR <= {STAN_FDR} (down negative)'); a.set_title('D. DE genes per contrast (edgeR QLF)')
    fig.tight_layout()
    STAN_IMG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(STAN_IMG_DIR / f'{STAN_SLUG}_qc_panel.png', dpi=130)


if __name__ == '__main__':
    STAN_OOUT_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        wide = run_edger(tmp)
        samples = sample_qc(tmp)
        qc_panel(tmp, samples, wide)
    wide.to_csv(STAN_OOUT_DIR / f'{STAN_SLUG}_edgeR_qlf_contrasts.tsv', sep='\t', index=False)
    samples.to_csv(STAN_OOUT_DIR / f'{STAN_SLUG}_sample_qc.tsv', sep='\t', index=False)
    ka = known_answer(wide)
    ka.to_csv(STAN_OOUT_DIR / f'{STAN_SLUG}_edgeR_known_answer.tsv', sep='\t', index=False)
    print(ka.round(4).to_string(index=False))
