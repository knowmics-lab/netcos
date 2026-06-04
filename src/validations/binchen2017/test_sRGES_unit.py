#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for src/sRGES.py.

These exercise the algorithm on a small hand-crafted synthetic dataset where
the expected outputs can be hand-computed, so that we can be confident the
Python port matches the R semantics of Bin Chen 2017's getsRGES.

Run from anywhere in the repo:

    cd src && python -m validations.binchen2017.test_sRGES_unit

or:

    cd src && python validations/binchen2017/test_sRGES_unit.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent.parent  # .../src
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


def hand_diff(df, score_col="RGES"):
    """Independent reference implementation of build_reference_diff, for cross-check."""
    sub = df[(df["pert_dose"] > 0) & (df["pert_time"].isin([6, 24]))].copy()
    merged = sub.merge(sub, on=["pert_iname", "cell_id"], suffixes=("_x", "_y"))
    merged = merged[(merged["pert_time_x"] == 24) & (merged["pert_dose_x"] == 10)]
    merged = merged[merged["id_x"] != merged["id_y"]]
    cmap_diff = merged[f"{score_col}_x"] - merged[f"{score_col}_y"]
    dose_bin = np.where(merged["pert_dose_y"] < 10, "low", "high")
    bucket = [f"{d} {int(t)}" for d, t in zip(dose_bin, merged["pert_time_y"])]
    return pd.Series(cmap_diff.values).groupby(bucket).mean().sort_index()


def build_synthetic():
    rows = []
    drug_levels = {
        "DrugA": [
            ("HEPG2", 10, 24, -0.50),  # reference
            ("HEPG2", 10,  6, -0.30),
            ("HEPG2",  1, 24, -0.20),
            ("HEPG2",  1,  6, -0.10),
            ("MCF7",  10, 24, -0.60),  # reference
            ("MCF7",  10,  6, -0.40),
            ("MCF7",   1, 24, -0.25),
            ("MCF7",   1,  6, -0.15),
        ],
        "DrugB": [
            ("HEPG2", 10, 24, +0.10),
            ("HEPG2", 10,  6, +0.20),
            ("HEPG2",  1, 24, +0.30),
            ("HEPG2",  1,  6, +0.40),
            ("MCF7",  10, 24, +0.00),
            ("MCF7",  10,  6, +0.05),
            ("MCF7",   1, 24, +0.10),
            ("MCF7",   1,  6, +0.15),
        ],
    }
    next_id = 0
    for drug, instances in drug_levels.items():
        for (cell, dose, time, rges) in instances:
            rows.append({
                "id": next_id, "pert_iname": drug, "cell_id": cell,
                "pert_dose": dose, "pert_time": time, "RGES": float(rges),
                "is_gold": "1",
            })
            next_id += 1
    return pd.DataFrame(rows)


def test_diff_matches_hand_computed():
    from sRGES import build_reference_diff
    df = build_synthetic()
    mine = build_reference_diff(df, score_col="RGES")
    expected = hand_diff(df, score_col="RGES")
    print("\n[test_diff] diff vector mine vs hand-computed:")
    for k in sorted(set(mine.index) | set(expected.index)):
        print(f"  bucket={k!r:>12}  mine={mine.get(k, np.nan):+.6f}  "
              f"expected={expected.get(k, np.nan):+.6f}")
    for k in expected.index:
        assert abs(mine[k] - expected[k]) < 1e-12, \
            f"bucket {k}: mine={mine[k]} vs expected={expected[k]}"
    print("[test_diff] PASS")


def test_full_pipeline_no_corweight():
    from sRGES import compute_sRGES, build_reference_diff
    df = build_synthetic()
    diff = build_reference_diff(df, score_col="RGES")

    out = compute_sRGES(
        df, score_col="RGES", drug_col="pert_iname",
        cell_col="cell_id", dose_col="pert_dose", time_col="pert_time",
        cell_lines=None, cell_line_weights=None,
        apply_cor_weight=False, srges_diff_mode="binchen_r", 
        filter_is_gold=False,
    )

    def expected_srges(drug):
        sub = df[df["pert_iname"] == drug]
        diff_d = dict(diff)
        adj = []
        for _, r in sub.iterrows():
            v = r["RGES"]
            bt = "short" if r["pert_time"] < 24 else "long"
            bd = "low" if r["pert_dose"] < 10 else "high"
            if bt == "short" and bd == "low":   v += diff_d.get("low 6", 0.0)
            if bd == "low" and bt == "long":    v += diff_d.get("high 6", 0.0)
            if bd == "high" and bt == "short":  v += diff_d.get("high 24", 0.0)
            adj.append(v)
        return float(np.mean(adj))

    for drug in df["pert_iname"].unique():
        my = float(out.loc[out["pert_iname"] == drug, "sRGES"].iloc[0])
        ex = expected_srges(drug)
        assert abs(my - ex) < 1e-12, f"{drug}: mine={my} expected={ex}"
    print("[test_full] PASS")


def test_per_cell_line_restriction():
    from sRGES import compute_sRGES
    df = build_synthetic()
    only_hep = compute_sRGES(df, score_col="RGES", cell_lines=["HEPG2"],
                              apply_cor_weight=False, filter_is_gold=False)
    full = compute_sRGES(df, score_col="RGES", cell_lines=None,
                          apply_cor_weight=False, filter_is_gold=False)
    # per-cell-line and full diffs are typically different -> sRGES differs.
    same = (only_hep.sort_values("pert_iname").reset_index(drop=True)["sRGES"]
            == full.sort_values("pert_iname").reset_index(drop=True)["sRGES"]).all()
    assert not same, "per-cell-line and full sRGES should generally differ"
    print("[test_per_cell_line] PASS")


def test_cor_weighting_smoke():
    from sRGES import compute_sRGES
    df = build_synthetic()
    out = compute_sRGES(df, score_col="RGES",
                        cell_line_weights={"HEPG2": 0.5, "MCF7": 0.3},
                        apply_cor_weight=True, filter_is_gold=False)
    out_noweights = compute_sRGES(df, score_col="RGES",
                                   apply_cor_weight=False,
                                   filter_is_gold=False)
    # cor weighting shrinks magnitudes (since cor/max_cor < 1 for some cells).
    assert abs(out["sRGES"].abs().sum()) < abs(out_noweights["sRGES"].abs().sum()), \
        "cor weighting should shrink overall sRGES magnitudes"
    print("[test_cor_weighting] PASS")


if __name__ == "__main__":
    test_diff_matches_hand_computed()
    test_full_pipeline_no_corweight()
    test_per_cell_line_restriction()
    test_cor_weighting_smoke()
    print("\nALL SYNTHETIC UNIT TESTS PASSED")
