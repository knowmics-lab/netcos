#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the distribution of connectivity scores for the CS run selected by conf.py.

The script delegates CS-run discovery to `loader.load_cs_run`, which uses the
hyperparameters set in src/conf.py (DISEASE, landmark_disease, cell_line, mith,
cs_on_LM, CS_ON_PATHWAYS, CS_METHOD, selected_cs_run_id) and returns the latest
matching run from logs/cs_runs.tsv whose .tsv still exists on disk.

Typical use (run from anywhere; the script self-locates src/):
    python src/plot_drug_ranking.py                              # use conf.py defaults
    python src/plot_drug_ranking.py --cs-run-id 13_05_2026_14_57_..._connectivity_score
    python src/plot_drug_ranking.py --cs-file path/to/run.tsv
    python src/plot_drug_ranking.py --cs-method bin_chen_disease_sorted
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

from conf import DISEASE, IMG_DIR, cell_line_run_name, disease_run_name  # noqa: E402
from loader import load_cs_run  # noqa: E402


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
        description="Plot the connectivity-score distribution for the CS run "
                    "selected by conf.py (or a CLI override)."
    )
    # CS-run selection overrides (all optional; default = conf.py)
    parser.add_argument(
        "--cs-run-id",
        default=None,
        help="Optional explicit cs_run_id. Overrides conf-based resolution.",
    )
    parser.add_argument(
        "--cs-file",
        type=Path,
        default=None,
        help="Optional explicit path to a CS .tsv. Overrides cs_run_id.",
    )
    parser.add_argument("--disease", default=None,
                        help="Override conf.DISEASE.")
    parser.add_argument("--cs-method", default=None,
                        help="Override conf.CS_METHOD (e.g. 'bin_chen', 'bin_chen_disease_sorted').")
    parser.add_argument("--mith", type=int, default=None,
                        help="Override conf.cs_mith (0 or 1).")
    parser.add_argument("--cs-on-lm", type=int, default=None,
                        help="Override conf.cs_on_LM (0 or 1).")
    parser.add_argument("--cs-on-pathways", type=int, default=None,
                        help="Override conf.CS_ON_PATHWAYS (0 or 1).")
    parser.add_argument("--landmark-disease", type=int, default=None,
                        help="Override conf.landmark_disease (0 or 1).")
    # Plot knobs
    parser.add_argument("--score-col", default="connectivity_score",
                        help="Column to plot.")
    parser.add_argument("--bins", type=int, default=100,
                        help="Number of histogram bins.")
    parser.add_argument("--xlim", nargs=2, type=float, default=None,
                        metavar=("XMIN", "XMAX"),
                        help="Optional x-axis limits, e.g. --xlim -2 2.")
    parser.add_argument("--no-split-pert-time", action="store_true",
                        help="Do not split histograms by perturbation_time "
                             "even if the column exists.")
    parser.add_argument("--out-dir", type=Path,
                        default=Path(IMG_DIR) / "drug_ranking_distributions",
                        help="Directory where the plot will be saved.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    df, cs_file, cs_run_id = load_cs_run(
        # CLI overrides win; None falls back to conf.py
        disease=args.disease,
        landmark_disease=(None if args.landmark_disease is None
                          else bool(args.landmark_disease)),
        mith=args.mith,
        cs_on_LM=args.cs_on_lm,
        CS_ON_PATHWAYS=args.cs_on_pathways,
        CS_METHOD=args.cs_method,
        selected_cs_run_id=args.cs_run_id,
        cs_file_path=args.cs_file,
        score_col=args.score_col,
    )

    title = (
        f"{args.disease or DISEASE} drug ranking CS distribution\n"
        f"{cs_run_id} | {disease_run_name} vs {cell_line_run_name}"
    )
    output_file = args.out_dir / f"{cs_run_id}_{args.score_col}_distribution.png"

    plot_score_distribution(
        df=df,
        score_col=args.score_col,
        title=title,
        output_file=output_file,
        bins=args.bins,
        xlim=tuple(args.xlim) if args.xlim else None,
        split_by_perturbation_time=not args.no_split_pert_time,
    )

    print("Resolved cs_run_id:", cs_run_id)
    print("Using CS file:", cs_file)
    print("Rows plotted:", len(df))
    print("Saved plot to:", output_file)


if __name__ == "__main__":
    main()
