#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the distribution of connectivity scores for the CS run selected by conf.py.

The script resolves the correct connectivity-score run from logs/cs_runs.tsv using
current conf.py parameters, opens the corresponding file under
connectivity_score/output/<disease_run_name>/, and saves a simple histogram.

Typical use, from the NetCos repo root:
    python src/plot_drug_ranking.py

Optional overrides:
    python src/plot_drug_ranking.py --cs-run-id 31_03_2026_16_25_mith_connectivity_score
    python src/plot_drug_ranking.py --score-col connectivity_score --bins 100
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent if (HERE.parent / "src").exists() else HERE
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))
else:
    sys.path.insert(0, str(REPO_ROOT))

from conf import (  # noqa: E402
    CS_METHOD,
    CS_ON_PATHWAYS,
    CS_OUT,
    DISEASE,
    IMG_DIR,
    LOGS_DIR,
    cell_line_run_name,
    cs_log_filename,
    cs_mith,
    cs_on_LM,
    disease_run_name,
    selected_cs_run_id,
)


def _to_int_flag(value) -> int:
    """Convert bool / 0-1 / string bool-like values to 0 or 1."""
    if isinstance(value, str):
        value = value.strip().lower()
        if value in {"true", "t", "yes", "y"}:
            return 1
        if value in {"false", "f", "no", "n"}:
            return 0
    return int(value)


def resolve_cs_run(
    cs_runs_tsv: Path,
    disease_run_id: str,
    drug_run_id: str,
    mith: int,
    cs_on_lm: int,
    cs_on_pathways: int,
    cs_method: str,
    explicit_cs_run_id: Optional[str] = None,
) -> pd.Series:
    """
    Resolve the CS run row from cs_runs.tsv.

    If explicit_cs_run_id is supplied, it takes precedence. Otherwise, match the
    current conf.py hyperparameters. If more than one run matches, use the most
    recent timestamp, mirroring the downstream validation scripts.
    """
    cs_runs_tsv = Path(cs_runs_tsv)
    if not cs_runs_tsv.exists():
        raise FileNotFoundError(f"CS run log not found: {cs_runs_tsv}")

    runs = pd.read_csv(cs_runs_tsv, sep="\t")
    required = {
        "cs_run_id",
        "timestamp",
        "disease_run_id",
        "drug_run_id",
        "mith",
        "cs_on_LM",
        "CS_ON_PATHWAYS",
        "CS_METHOD",
        "output_file",
    }
    missing = sorted(required - set(runs.columns))
    if missing:
        raise ValueError(f"Missing required columns in {cs_runs_tsv}: {missing}")

    runs["mith"] = runs["mith"].map(_to_int_flag)
    runs["cs_on_LM"] = runs["cs_on_LM"].map(_to_int_flag)
    runs["CS_ON_PATHWAYS"] = runs["CS_ON_PATHWAYS"].map(_to_int_flag)

    if explicit_cs_run_id:
        hit = runs[runs["cs_run_id"].astype(str) == str(explicit_cs_run_id)].copy()
        if hit.empty:
            raise ValueError(f"cs_run_id '{explicit_cs_run_id}' not found in {cs_runs_tsv}")
    else:
        hit = runs[
            (runs["disease_run_id"].astype(str) == str(disease_run_id))
            & (runs["drug_run_id"].astype(str) == str(drug_run_id))
            & (runs["mith"] == int(mith))
            & (runs["cs_on_LM"] == int(cs_on_lm))
            & (runs["CS_ON_PATHWAYS"] == int(cs_on_pathways))
            & (runs["CS_METHOD"].astype(str) == str(cs_method))
        ].copy()
        if hit.empty:
            raise ValueError(
                "No matching CS run found in cs_runs.tsv for current conf.py parameters:\n"
                f"  disease_run_id={disease_run_id}\n"
                f"  drug_run_id={drug_run_id}\n"
                f"  mith={int(mith)}\n"
                f"  cs_on_LM={int(cs_on_lm)}\n"
                f"  CS_ON_PATHWAYS={int(cs_on_pathways)}\n"
                f"  CS_METHOD={cs_method}"
            )

    # Prefer the most recent matching run. dayfirst handles your dd/mm/YYYY log.
    hit["_timestamp"] = pd.to_datetime(hit["timestamp"], errors="coerce", dayfirst=True)
    hit = hit.sort_values(["_timestamp", "timestamp"], ascending=False, na_position="last")
    return hit.iloc[0]


def resolve_cs_file(run_row: pd.Series, cs_out: Path) -> Path:
    """
    Resolve the CS output file in the current checkout.

    The log stores absolute paths from the machine that ran the CS job. To keep
    this portable across machines, first try current conf.CS_OUT/<filename>, then
    fall back to the absolute path recorded in the log.
    """
    logged_output = Path(str(run_row["output_file"]))
    filename = logged_output.name if logged_output.name else f"{run_row['cs_run_id']}.tsv"

    candidates = [
        Path(cs_out) / filename,
        Path(cs_out) / f"{run_row['cs_run_id']}.tsv",
        logged_output,
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "Could not find the resolved CS output file. Tried:\n"
        + "\n".join(f"  - {candidate}" for candidate in candidates)
    )


def load_cs_table(cs_file: Path, score_col: str) -> pd.DataFrame:
    """Load a CS table and validate the score column."""
    cs_file = Path(cs_file)
    df = pd.read_csv(cs_file, sep="\t")
    if score_col not in df.columns:
        raise ValueError(
            f"Score column '{score_col}' not found in {cs_file}. "
            f"Available columns: {list(df.columns)}"
        )
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col]).copy()
    if df.empty:
        raise ValueError(f"No numeric values found in column '{score_col}' of {cs_file}")
    return df


def plot_score_distribution(
    df: pd.DataFrame,
    score_col: str,
    title: str,
    output_file: Path,
    bins: int = 100,
    xlim: Optional[tuple[float, float]] = None,
    split_by_perturbation_time: bool = True,
) -> None:
    """Plot the connectivity-score distribution, optionally split by perturbation_time."""
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))

    if split_by_perturbation_time and "perturbation_time" in df.columns:
        for pert_time, sub_df in df.groupby("perturbation_time", sort=True):
            ax.hist(sub_df[score_col], bins=bins, alpha=0.45, label=str(pert_time))
    else:
        ax.hist(df[score_col], bins=bins, alpha=0.75, label=score_col)

    ax.axvline(0, linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("Connectivity score")
    ax.set_ylabel("Number of drug profiles")
    ax.set_title(title)
    ax.grid(linestyle="--", alpha=0.35)
    ax.legend(frameon=False)

    if xlim is not None:
        ax.set_xlim(*xlim)

    # Lightweight summary annotation.
    scores = df[score_col].to_numpy(dtype=float)
    txt = (
        f"n={len(scores)}\n"
        f"mean={np.mean(scores):.3g}\n"
        f"median={np.median(scores):.3g}\n"
        f"min={np.min(scores):.3g}\n"
        f"max={np.max(scores):.3g}"
    )
    ax.text(
        0.98,
        0.97,
        txt,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the connectivity-score distribution for the CS run selected by conf.py."
    )
    parser.add_argument(
        "--cs-run-id",
        default=selected_cs_run_id,
        help="Optional explicit cs_run_id. Overrides automatic resolution from conf.py.",
    )
    parser.add_argument(
        "--cs-runs-tsv",
        type=Path,
        default=Path(cs_log_filename) if "cs_log_filename" in globals() else Path(LOGS_DIR) / "cs_runs.tsv",
        help="Path to cs_runs.tsv.",
    )
    parser.add_argument("--score-col", default="connectivity_score", help="Column to plot.")
    parser.add_argument("--bins", type=int, default=100, help="Number of histogram bins.")
    parser.add_argument(
        "--xlim",
        nargs=2,
        type=float,
        default=None,
        metavar=("XMIN", "XMAX"),
        help="Optional x-axis limits, e.g. --xlim -2 2.",
    )
    parser.add_argument(
        "--no-split-pert-time",
        action="store_true",
        help="Do not split histograms by perturbation_time even if the column exists.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(IMG_DIR) / "drug_ranking_distributions",
        help="Directory where the plot will be saved.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    run_row = resolve_cs_run(
        cs_runs_tsv=args.cs_runs_tsv,
        disease_run_id=disease_run_name,
        drug_run_id=cell_line_run_name,
        mith=cs_mith,
        cs_on_lm=cs_on_LM,
        cs_on_pathways=CS_ON_PATHWAYS,
        cs_method=CS_METHOD,
        explicit_cs_run_id=args.cs_run_id,
    )

    cs_file = resolve_cs_file(run_row, CS_OUT)
    df = load_cs_table(cs_file, args.score_col)

    title = (
        f"{DISEASE} drug ranking CS distribution\n"
        f"{run_row['cs_run_id']} | {disease_run_name} vs {cell_line_run_name}"
    )
    output_file = args.out_dir / f"{run_row['cs_run_id']}_{args.score_col}_distribution.png"

    plot_score_distribution(
        df=df,
        score_col=args.score_col,
        title=title,
        output_file=output_file,
        bins=args.bins,
        xlim=tuple(args.xlim) if args.xlim else None,
        split_by_perturbation_time=not args.no_split_pert_time,
    )

    print("Resolved cs_run_id:", run_row["cs_run_id"])
    print("Using CS file:", cs_file)
    print("Rows plotted:", len(df))
    print("Saved plot to:", output_file)


if __name__ == "__main__":
    main()
