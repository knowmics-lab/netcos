#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map MITHrIL batch output files to NetCos connectivity-score input, 1:1,
using analysis-neutral column names.

Compared with the older NetCos step7 script, this version:
- keeps a 1:1 mapping between each `*.perturbations.txt` file and one output pickle;
- preserves the old cleaning rule that removes special gene names;
- does NOT encode timepoint labels in perturbation column names;
- can optionally append selected metadata columns from `lincs_sig_info_new.csv`.

This makes each output file represent one signature-item, while the choice to run
6h-only, 24h-only, or pooled analyses is deferred to downstream analysis.
"""

from __future__ import annotations

import argparse
import pickle
import re
from pathlib import Path
from typing import Iterable, Optional
from conf import MITH_OUT_DRUG, CS_IN_DRUG
import pandas as pd
from time import time


MITH_REQUIRED_COLUMNS = {
    "gene_id": 2,
    "gene": 3,
    "perturbation": 4,
    "p_value": 6,
    "adj_p_value": 7,
}

SPECIAL_GENE_CHARS = ';/:*?"|'
DEFAULT_METADATA_COLUMNS = [
    "pert_iname",
    "pert_desc",
    "pert_id",
    "sig_id",
    "cell_id",
    "pert_time",
    "pert_dose",
    "pert_type",
    "is_gold",
    "DrugBank.ID",
]


def remove_special_gene_names(df: pd.DataFrame, gene_col: str = "gene") -> pd.DataFrame:
    """Reproduce the old NetCos cleaning rule for special gene names."""
    out = df.copy()
    for char in SPECIAL_GENE_CHARS:
        mask = out[gene_col].astype(str).str.contains(re.escape(char), na=False)
        out = out.loc[~mask]
    return out


def infer_signature_id(path: Path) -> str:
    suffix = ".perturbations.txt"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected MITHrIL output filename: {path.name}")
    return path.name[: -len(suffix)]


def load_metadata(metadata_csv: Optional[Path], id_column: str) -> Optional[pd.DataFrame]:
    if metadata_csv is None:
        return None
    meta = pd.read_csv(metadata_csv)
    if id_column not in meta.columns:
        raise ValueError(f"Column '{id_column}' not found in metadata file {metadata_csv}")
    meta = meta.copy()
    meta[id_column] = meta[id_column].astype(str)
    if meta[id_column].duplicated().any():
        dup = meta.loc[meta[id_column].duplicated(), id_column].astype(str).tolist()[:10]
        raise ValueError(
            f"Metadata column '{id_column}' contains duplicated IDs. Examples: {dup}"
        )
    return meta


def build_output_frame(
    mith_path: Path,
    metadata_row: Optional[pd.Series] = None,
    metadata_columns: Optional[list[str]] = None,
    drop_special_gene_names: bool = True,
) -> pd.DataFrame:
    raw = pd.read_csv(mith_path, sep="\t", header=None, skiprows=1)
    max_required_col = max(MITH_REQUIRED_COLUMNS.values())
    if raw.shape[1] <= max_required_col:
        raise ValueError(
            f"{mith_path.name} has only {raw.shape[1]} columns, but step7 needs at least "
            f"{max_required_col + 1}."
        )

    signature_id = infer_signature_id(mith_path)

    frame = raw[[
        MITH_REQUIRED_COLUMNS["gene_id"],
        MITH_REQUIRED_COLUMNS["gene"],
        MITH_REQUIRED_COLUMNS["perturbation"],
        MITH_REQUIRED_COLUMNS["p_value"],
        MITH_REQUIRED_COLUMNS["adj_p_value"],
    ]].copy()

    # Old step7 behaviour: keep the first occurrence of each gene_id across pathways.
    frame.drop_duplicates(subset=MITH_REQUIRED_COLUMNS["gene_id"], keep="first", inplace=True)

    frame.columns = [
        "gene_id",
        "gene",
        "Perturbation",
        "p.value",
        "adj.p.value",
    ]
    frame.insert(0, "signature_id", signature_id)

    if drop_special_gene_names:
        frame = remove_special_gene_names(frame, gene_col="gene")

    frame["gene_id"] = frame["gene_id"].astype(str)

    if metadata_row is not None:
        columns_to_add = metadata_columns or DEFAULT_METADATA_COLUMNS
        for col in columns_to_add:
            if col in metadata_row.index:
                frame[col] = metadata_row[col]

    return frame.reset_index(drop=True)


def iter_input_files(input_dir: Path) -> list[Path]:
    files = sorted(input_dir.glob("*.perturbations.txt"))
    if not files:
        raise FileNotFoundError(f"No '*.perturbations.txt' files found in {input_dir}")
    return files


def filter_files_by_metadata(
    files: Iterable[Path],
    metadata_df: Optional[pd.DataFrame],
    id_column: str,
) -> list[Path]:
    files = list(files)
    if metadata_df is None:
        return files

    allowed = set(metadata_df[id_column].astype(str))
    filtered = [p for p in files if infer_signature_id(p) in allowed]
    missing = sorted(allowed.difference({infer_signature_id(p) for p in files}))
    if missing:
        print(f"Warning: {len(missing)} metadata IDs do not have a matching MITHrIL output file.")
        print("First missing IDs:", missing[:10])
    return filtered


def convert_directory(
    input_dir: Path,
    output_dir: Path,
    metadata_csv: Optional[Path] = None,
    id_column: str = "id",
    metadata_columns: Optional[list[str]] = None,
    save_tsv: bool = False,
    overwrite: bool = False,
    drop_special_gene_names: bool = True,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    files = iter_input_files(input_dir)
    metadata_df = load_metadata(metadata_csv, id_column=id_column)
    files = filter_files_by_metadata(files, metadata_df=metadata_df, id_column=id_column)

    print(f"Input directory : {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Files to convert: {len(files)}")
    start = time()
    if metadata_df is not None:
        cols = metadata_columns or DEFAULT_METADATA_COLUMNS
        existing = [c for c in cols if c in metadata_df.columns]
        print(f"Metadata columns appended: {existing}")

    for idx, mith_path in enumerate(files, start=1):
        signature_id = infer_signature_id(mith_path)
        out_pkl = output_dir / f"{signature_id}.pkl"
        out_tsv = output_dir / f"{signature_id}.tsv"

        if out_pkl.exists() and not overwrite:
            print(f"[{idx}/{len(files)}] skipping existing {out_pkl.name}")
            continue

        metadata_row = None
        if metadata_df is not None:
            match = metadata_df.loc[metadata_df[id_column] == signature_id]
            if not match.empty:
                metadata_row = match.iloc[0]

        df = build_output_frame(
            mith_path,
            metadata_row=metadata_row,
            metadata_columns=metadata_columns,
            drop_special_gene_names=drop_special_gene_names,
        )

        with open(out_pkl, "wb") as fh:
            pickle.dump(df, fh)
        if save_tsv:
            df.to_csv(out_tsv, sep="\t", index=False)

        print(f"[{idx}/{len(files)}] wrote {out_pkl.name} ({len(df)} rows)")
    print('DONE. time passed:', time()-start)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Map MITHrIL perturbation outputs to NetCos CS-input files with 1:1 filename preservation "
            "and neutral perturbation column names."
        )
    )
    
    parser.add_argument("--input_dir", type=Path, default=MITH_OUT_DRUG, help="Directory containing *.perturbations.txt files")
    parser.add_argument("--output_dir", type=Path, default=CS_IN_DRUG, help="Directory where .pkl outputs will be written")

    parser.add_argument(
        "--metadata-csv",
        type=Path,
        default=None,
        help="Optional metadata table (e.g. lincs_sig_info_new.csv)",
    )
    parser.add_argument(
        "--id-column",
        default="id",
        help="Column in --metadata-csv containing the signature IDs. Default: id",
    )
    parser.add_argument(
        "--metadata-columns",
        nargs="+",
        default=None,
        help=(
            "Optional list of metadata columns to append to every output dataframe. "
            "Default: a standard LINCS/DrugBank subset."
        ),
    )
    parser.add_argument(
        "--save-tsv",
        action="store_true",
        help="Also save a human-readable TSV next to each pickle",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files",
    )
    parser.add_argument(
        "--keep-special-gene-names",
        action="store_true",
        help="Do not apply the old NetCos special-gene-name cleaning rule",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert_directory(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        metadata_csv=args.metadata_csv,
        id_column=args.id_column,
        metadata_columns=args.metadata_columns,
        save_tsv=args.save_tsv,
        overwrite=args.overwrite,
        drop_special_gene_names=not args.keep_special_gene_names,
    )
