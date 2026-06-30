#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map per-drug DEG meta-analysis files to NetCos connectivity-score input.

Source (one file per drug, already collapsed to one signature per drug via the
LMM, then gene-wise -> drug-wise mapped by step4 and de-duplicated on gene_id):

    TSR_OUT_DRUG / 'LINCS' / 'metanalysis_drug_wise_filtered' / '<drug>_metanalysis.pkl'

Each such file holds ALL three perturbation timepoints side by side, with the
column naming produced by the R meta-analysis pipeline:

    gene_id, gene, drug,
    DE_log2_FC_6h,     std.error_6h,     t.value_6h,     p.value_6h,     adj.p.value_6h,
    DE_log2_FC_24h,    std.error_24h,    t.value_24h,    p.value_24h,    adj.p.value_24h,
    DE_log2_FC_6h_24h, std.error_6h_24h, t.value_6h_24h, p.value_6h_24h

NOTE: the pooled 6h_24h timepoint has NO 'adj.p.value_6h_24h' column (this is a
known quirk of the meta-analysis output). See --pooled-pval below.

Output (mirrors validations/chembl Preprocess_binChen2017 `write_cs_drug_inputs`,
i.e. the DEG CS-input convention that the *current* cs_batch.py reads):

    CS_IN_DRUG / '<drug>_<timestep>_signature_gene_id.pkl'

with columns:  signature_id, gene_id (str), DE_log2_FC, adj.p.value

The signature id (= '<drug>_<timestep>') becomes the value of the 'pert_id'
column in the cs_batch output, because the IPF/ALS pipeline runs with
conf.drug_collapsed_before_cs = True.

cs_batch.py (DEG branch, rank_on='magnitude') only ranks on DE_log2_FC
(columns[1]); the adj.p.value column (columns[2]) must merely EXIST so the
`columns_of_interest = [gene_id, DE_log2_FC, adj.p.value]` selection does not
raise. It is used for ranking only when rank_on='p_value'.

Run on the server (paths come from the active conf.py / local.py), e.g.:

    cd src
    python step8_DEG_drug_metanalysis_to_cs_in.py --limit 1   # dry-run one drug
    python step8_DEG_drug_metanalysis_to_cs_in.py             # full run
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path
import pandas as pd
from conf import TSR_OUT_DRUG, CS_IN_DRUG

DEFAULT_SRC_DIR = Path(TSR_OUT_DRUG) / "LINCS" / "metanalysis_drug_wise_filtered"
DEFAULT_OUT_DIR = Path(CS_IN_DRUG)

TIMESTEPS = ["6h", "24h", "6h_24h"]
SRC_SUFFIX = "_metanalysis"  # source files are '<drug>_metanalysis.{pkl,csv}'
OUT_SUFFIX = "_signature_gene_id.pkl"  # cs_batch DEG-input convention

# cs_batch DEG branch expects exactly these output columns
ID_COL = "gene_id"
VALUE_COL = "DE_log2_FC"
PVAL_COL = "adj.p.value"


def read_source(path: Path) -> pd.DataFrame:
    """Load one '<drug>_metanalysis' file (pickle preferred, csv fallback)."""
    if path.suffix == ".pkl":
        with open(path, "rb") as fh:
            return pickle.load(fh)
    # tab-separated, matching loader.load_single_drug_signature(pkl=False)
    return pd.read_csv(path, sep="\t", dtype={ID_COL: "str"})


def resolve_pval(df: pd.DataFrame, timestep: str, pooled_pval: str):
    """
    Return the Series to write as 'adj.p.value' for this timestep.

    For 6h / 24h the real 'adj.p.value_<timestep>' column is used.
    For the pooled 6h_24h timestep (which has no adjusted p-value) the
    behaviour is controlled by `pooled_pval`:
        'pvalue'  -> fall back to 'p.value_6h_24h'
        'one'     -> constant 1.0 (matches the BinChen DEG preprocessing)
    """
    adj_col = f"adj.p.value_{timestep}"
    if adj_col in df.columns:
        return df[adj_col]

    if pooled_pval == "pvalue":
        raw_col = f"p.value_{timestep}"
        if raw_col in df.columns:
            print(f"    [{timestep}] no '{adj_col}'; falling back to '{raw_col}'")
            return df[raw_col]
        print(f"    [{timestep}] no '{adj_col}' nor '{f'p.value_{timestep}'}'; "
              f"using constant 1.0")
        return 1.0

    # pooled_pval == 'one'
    print(f"    [{timestep}] no '{adj_col}'; using constant 1.0")
    return 1.0


def split_one_drug(df: pd.DataFrame, drug: str, out_dir: Path,
                   timesteps, pooled_pval: str, overwrite: bool) -> int:
    """Write one CS-input pickle per timestep for a single drug. Returns n written."""
    if ID_COL not in df.columns:
        raise KeyError(
            f"{drug}: source has no '{ID_COL}' column (found: {list(df.columns)})")

    written = 0
    for ts in timesteps:
        fc_col = f"{VALUE_COL}_{ts}"
        if fc_col not in df.columns:
            print(f"    [{ts}] missing '{fc_col}', skipping this timestep")
            continue

        signature_id = f"{drug}_{ts}"
        out_pkl = out_dir / f"{signature_id}{OUT_SUFFIX}"
        if out_pkl.exists() and not overwrite:
            print(f"    skipping existing {out_pkl.name}")
            continue

        out = pd.DataFrame({
            "signature_id": signature_id,
            ID_COL: df[ID_COL].astype(str),
            VALUE_COL: pd.to_numeric(df[fc_col], errors="coerce"),
            PVAL_COL: resolve_pval(df, ts, pooled_pval),
        }).dropna(subset=[ID_COL, VALUE_COL]).reset_index(drop=True)

        with open(out_pkl, "wb") as fh:
            pickle.dump(out, fh)
        print(f"    wrote {out_pkl.name}  ({out.shape[0]} genes)")
        written += 1
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src-dir", type=Path, default=DEFAULT_SRC_DIR,
                    help=f"metanalysis_drug_wise_filtered dir (default: {DEFAULT_SRC_DIR})")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                    help=f"CS drug-input dir (default: {DEFAULT_OUT_DIR})")
    ap.add_argument("--timesteps", nargs="+", default=TIMESTEPS, choices=TIMESTEPS,
                    help="timesteps to emit (default: all three)")
    ap.add_argument("--pooled-pval", choices=["pvalue", "one"], default="pvalue",
                    help="adj.p.value source for the pooled 6h_24h timestep that "
                         "lacks one: 'pvalue' = use p.value_6h_24h (default), "
                         "'one' = constant 1.0")
    ap.add_argument("--overwrite", action="store_true",
                    help="overwrite existing output pickles")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N drugs (dry-run / smoke test)")
    args = ap.parse_args(argv)

    src_dir = args.src_dir
    out_dir = args.out_dir
    if not src_dir.is_dir():
        sys.exit(f"source dir not found: {src_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    # prefer .pkl, fall back to .csv only for drugs without a .pkl
    pkls = sorted(p for p in src_dir.iterdir() if p.name.endswith(SRC_SUFFIX + ".pkl"))
    pkl_stems = {p.name[:-len(SRC_SUFFIX + ".pkl")] for p in pkls}
    csvs = sorted(p for p in src_dir.iterdir()
                  if p.name.endswith(SRC_SUFFIX + ".csv")
                  and p.name[:-len(SRC_SUFFIX + ".csv")] not in pkl_stems)
    src_files = pkls + csvs
    if args.limit is not None:
        src_files = src_files[:args.limit]

    print(f"source : {src_dir}")
    print(f"output : {out_dir}")
    print(f"drugs  : {len(src_files)}  timesteps: {args.timesteps}  "
          f"pooled_pval: {args.pooled_pval}")

    total_written = 0
    for n, path in enumerate(src_files, 1):
        if path.name.endswith(SRC_SUFFIX + ".pkl"):
            drug = path.name[:-len(SRC_SUFFIX + ".pkl")]
        else:
            drug = path.name[:-len(SRC_SUFFIX + ".csv")]
        print(f"[{n}/{len(src_files)}] {drug}")
        df = read_source(path)
        total_written += split_one_drug(df, drug, out_dir, args.timesteps,
                                        args.pooled_pval, args.overwrite)

    print(f"\ndone: wrote {total_written} CS-input files for {len(src_files)} drugs "
          f"into {out_dir}")


if __name__ == "__main__":
    main()
