#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Disease signature of asthmatic vs healthy-control nasal airway epithelium
(GSE152004) with the R classes of modules/disease_signature/, written
directly as the DEG CS input and the MITHrIL input, plus a QC panel.

Design
------
- Samples and groups: the '<slug>_samples.tsv' table written by
  00_inspect_geo_download.py (one row per GSM column of the NCBI counts, in
  column order, group_code '0' = asthmatic, '1' = healthy control). The
  script stops if its GSM order differs from the header of the counts.
- R (Rscript, working directory = BASE_DIR, packages R6, data.table, edgeR,
  limma), with the classes used for the IPF signature:
  SamplesGroups (group-map string -> factor, levels Control, disease),
  LogarithmGeneFilter (keep genes with mean log2(count + 1) > 6),
  DiseaseSignature$compute(..., filter_by_proteing_coding = FALSE), i.e.
  DGEList(remove.zeros = TRUE), calcNormFactors(method = 'upperquartile'),
  voom, lmFit and eBayes on ~sample_type, topTable(coef = 2): logFC is
  asthmatic minus healthy control. Row names carry a 'GeneID_' prefix
  inside R (integer-like row names do not survive data.frame() and the
  topTable column order DiseaseSignature relies on); it is removed on
  output. The counts are read as in
  IlluminaHiSeqMapper (data.table::fread, rows = NCBI GeneID), but the rows
  keep their GeneID: there is no translation to gene symbols.
- QC: library sizes, PCA of the upper-quartile-normalised log-CPM of the
  filtered genes, volcano plot, and the airway type-2 epithelial genes in
  KNOWN_UP (expected higher in asthma), looked up by their Entrez id via
  conf.alias_2geneid.

Inputs
------
- case_conf.NCBI_COUNTS (rows = NCBI GeneID, columns = GSM).
- P.oout_dir/<slug>_samples.tsv (from 00_).

Outputs
-------
- P.cs_in_disease/<slug>_signature_gene_id.csv
    ';'-separated, ',' decimal; columns gene_id, DE_log2_FC, std.error,
    t.value, p.value, adj.p.value. DEG-level CS input, read by
    loader.load_disease_signature with DISEASE_SIGNATURE_SOURCE = "cs_in".
- conf.MITH_IN_DISEASE/<slug>_signature_gene_id.mi
    MITHrIL input, one 'gene_id<TAB>log2FC' line per gene, no header.
- P.oout_dir/<slug>_disease_signature_summary.tsv
    genes in the counts, genes after the filter, samples per group, genes at
    adj.p.value <= FDR up / down, and the KNOWN_UP genes with logFC and
    adj.p.value.
- P.img_dir/<slug>_disease_signature_qc_panel.png
    A: library sizes per group; B: PCA; C: volcano; D: KNOWN_UP genes,
    log-CPM per group.

Usage
-----
    python 01_disease_signature.py   (from any directory)
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

sys.path.insert(0, str(Path(__file__).resolve().parent))

from case_conf import P, NCBI_COUNTS
from conf import BASE_DIR, MITH_IN_DISEASE, alias_2geneid

FDR = 0.05
KNOWN_UP = ['POSTN', 'CLCA1', 'SERPINB2', 'CST1', 'CCL26']
GROUP_NAME = {'0': 'asthmatic', '1': 'healthy control'}
GROUP_COLORS = {'asthmatic': '#d62728', 'healthy control': '#1f77b4'}

R_CODE = r'''
suppressMessages({library(R6); library(data.table); library(edgeR); library(limma)})
source("modules/disease_signature/SamplesGroups.R")
source("modules/disease_signature/filter/LogarithmGeneFilter.R")
source("modules/disease_signature/DiseaseSignature.R")
a <- commandArgs(trailingOnly=TRUE)
counts_file <- a[1]; samples_groups_map <- readLines(a[2], n=1); tmp <- a[3]; disease_name <- a[4]
x <- as.matrix(data.table::fread(counts_file, header=TRUE, colClasses="integer"), rownames="GeneID")
sg <- SamplesGroups$new()$groups(samples_groups_map, disease_name)
stopifnot(length(sg$sample_groups_idxes) == ncol(x))
x <- x[, sg$sample_groups_idxes]
# non-numeric row names: data.frame() prefixes integer-like row names with "X", and a
# leading ID column then appears in topTable, where DiseaseSignature reads column 1 as logFC
rownames(x) <- paste0("GeneID_", rownames(x))
write.table(data.frame(GSM=colnames(x), lib_size=colSums(x)), file.path(tmp, "libs.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
n_all <- nrow(x)
x <- LogarithmGeneFilter$new()$filter(x, sg$disease_control_groups)
sig <- DiseaseSignature$new()$compute(
  list(gene_expressions=data.frame(x, check.names=FALSE), disease_control_groups=sg$disease_control_groups), FALSE)
write.table(sig, file.path(tmp, "signature.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
y <- calcNormFactors(DGEList(x, remove.zeros=TRUE), method="upperquartile")
write.table(data.frame(gene_id=sub("^GeneID_", "", rownames(y)), cpm(y, log=TRUE), check.names=FALSE), file.path(tmp, "logcpm.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
writeLines(as.character(c(n_all, nrow(x))), file.path(tmp, "n_genes.txt"))
cat("genes in counts:", n_all, " after LogarithmGeneFilter:", nrow(x), " in signature:", nrow(sig), "\n")
'''


def read_samples():
    s = pd.read_csv(P.oout_dir / f'{P.slug}_samples.tsv', sep='\t', dtype={'group_code': str})
    header = pd.read_csv(NCBI_COUNTS, sep='\t', nrows=0).columns.tolist()
    assert header[0] == 'GeneID' and header[1:] == s['GSM'].tolist(), \
        'GSM order of the samples table differs from the counts header'
    s['group'] = s['group_code'].map(GROUP_NAME)
    return s


def run_r(tmp, s):
    (tmp / 'map.txt').write_text(''.join(s['group_code']) + '\n')
    (tmp / 'signature.R').write_text(R_CODE)
    subprocess.run(['Rscript', str(tmp / 'signature.R'), str(NCBI_COUNTS), str(tmp / 'map.txt'),
                    str(tmp), P.slug], check=True, cwd=BASE_DIR)
    sig = pd.read_csv(tmp / 'signature.tsv', sep='\t', dtype={'gene': str})
    sig['gene'] = sig['gene'].str.removeprefix('GeneID_')
    assert sig['gene'].str.fullmatch(r'\d+').all(), 'signature gene ids are not plain NCBI GeneIDs'
    libs = pd.read_csv(tmp / 'libs.tsv', sep='\t')
    lc = pd.read_csv(tmp / 'logcpm.tsv', sep='\t', dtype={'gene_id': str}).set_index('gene_id')
    n_all, n_kept = map(int, (tmp / 'n_genes.txt').read_text().split())
    return sig, libs, lc, n_all, n_kept


def qc_panel(s, sig, lc, known):
    X = lc.T.values - lc.T.values.mean(axis=0)
    U, S, _ = np.linalg.svd(X, full_matrices=False)
    pcs = U * S; ve = 100 * S ** 2 / (S ** 2).sum()
    grp = s.set_index('GSM').loc[lc.columns, 'group'].values
    groups = ['healthy control', 'asthmatic']

    fig, ax = plt.subplots(2, 2, figsize=(14, 11))
    a = ax[0, 0]
    a.boxplot([s.loc[s['group'] == g, 'lib_size'] / 1e6 for g in groups])
    a.set_xticks([1, 2]); a.set_xticklabels(groups)
    a.set_ylabel('library size (millions of reads)'); a.set_title('A. Library sizes')

    a = ax[0, 1]
    for g in groups:
        m = grp == g
        a.scatter(pcs[m, 0], pcs[m, 1], s=10, alpha=.6, color=GROUP_COLORS[g], label=f'{g} (n={m.sum()})')
    a.set_xlabel(f'PC1 ({ve[0]:.1f}%)'); a.set_ylabel(f'PC2 ({ve[1]:.1f}%)')
    a.set_title('B. PCA, upper-quartile log-CPM (filtered genes)'); a.legend()

    a = ax[1, 0]
    sg = sig['adj.p.value'] <= FDR
    a.scatter(sig['DE_log2_FC'], -np.log10(sig['p.value'].clip(lower=1e-300)), s=3, color=np.where(sg, '#d62728', '#bbbbbb'))
    for sym, gid in known.items():
        r = sig[sig['gene_id'] == gid]
        if len(r):
            a.annotate(sym, (r['DE_log2_FC'].iloc[0], -np.log10(max(r['p.value'].iloc[0], 1e-300))), fontsize=9)
    a.axvline(0, color='k', lw=.5)
    a.set_xlabel('log2 FC (asthmatic - healthy control)'); a.set_ylabel('-log10 p')
    a.set_title(f'C. Volcano ({int(sg.sum())} genes at adj.p <= {FDR})')

    a = ax[1, 1]
    present = [(sym, gid) for sym, gid in known.items() if gid in lc.index]
    for i, (sym, gid) in enumerate(present):
        for j, g in enumerate(groups):
            a.boxplot([lc.loc[gid].values[grp == g]], positions=[i * 3 + j], widths=.8, patch_artist=True,
                      boxprops=dict(facecolor=GROUP_COLORS[g]), showfliers=False)
    a.set_xticks([i * 3 + .5 for i in range(len(present))]); a.set_xticklabels([p[0] for p in present])
    a.set_ylabel('log-CPM'); a.set_title('D. Known-answer genes (blue control, red asthma)')
    fig.tight_layout()
    P.img_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(P.img_dir / f'{P.slug}_disease_signature_qc_panel.png', dpi=130)


if __name__ == '__main__':
    s = read_samples()
    with tempfile.TemporaryDirectory() as t:
        sig, libs, lc, n_all, n_kept = run_r(Path(t), s)
    s = s.merge(libs, on='GSM')

    out = sig.rename(columns={'gene': 'gene_id'})[
        ['gene_id', 'DE_log2_FC', 'std.error', 't.value', 'p.value', 'adj.p.value']].sort_values('gene_id')
    assert not out['gene_id'].duplicated().any()
    P.cs_in_disease.mkdir(parents=True, exist_ok=True)
    out.to_csv(P.cs_in_disease / f'{P.slug}_signature_gene_id.csv', sep=';', decimal=',', index=False)
    MITH_IN_DISEASE.mkdir(parents=True, exist_ok=True)
    with open(MITH_IN_DISEASE / f'{P.slug}_signature_gene_id.mi', 'w') as f:
        f.write('\n'.join(f'{g}\t{fc}' for g, fc in zip(out['gene_id'], out['DE_log2_FC'])))

    known = {sym: str(alias_2geneid[sym]) for sym in KNOWN_UP if sym in alias_2geneid}
    qc_panel(s, out, lc, known)

    sg = out['adj.p.value'] <= FDR
    rows = [('genes in NCBI counts', n_all), ('genes after LogarithmGeneFilter', n_kept),
            ('genes in signature', len(out))]
    rows += [(f'samples {g}', int((s['group'] == g).sum())) for g in ['asthmatic', 'healthy control']]
    rows += [(f'genes adj.p.value <= {FDR} up', int((sg & (out['DE_log2_FC'] > 0)).sum())),
             (f'genes adj.p.value <= {FDR} down', int((sg & (out['DE_log2_FC'] < 0)).sum()))]
    for sym in KNOWN_UP:
        r = out[out['gene_id'] == known.get(sym)]
        rows.append((f'{sym} (gene_id {known.get(sym)}) log2FC / adj.p.value',
                     f"{r['DE_log2_FC'].iloc[0]:.3f} / {r['adj.p.value'].iloc[0]:.2e}" if len(r) else 'not in signature'))
    summary = pd.DataFrame(rows, columns=['item', 'value'])
    P.oout_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(P.oout_dir / f'{P.slug}_disease_signature_summary.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))
