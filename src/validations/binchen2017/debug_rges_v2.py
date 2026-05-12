#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Debug & repair of connectivity_score.bin_chen_connectivity against the
Bin Chen 2017 reference RGES, using 3 HEPG2 LIHC drug signatures as test
cases:

    LINCS_id   drug             ref RGES   (from all_lincs_score.csv)
    8726       AG-957           -0.0390461497984907
    8728       albendazole      +0.187899879214054
    8729       alclometasone    +0.103086619069829

The current bin_chen_connectivity returns +0.0390, -0.1879, -0.1031 for
these three — a perfect sign flip to floating-point precision.

HYPOTHESIS: rank_genes() sorts drug_signature by drug_col_name with the
pandas default (ascending), but Lamb 2006 / Bin Chen 2017 / CMap convention
is descending (rank 1 = most upregulated reference gene). The np.sort
inside compute_KS sorts the resulting V values but does not invert their
semantic meaning — V[i] is still the ASCENDING position of disease gene
i in the drug profile, half-way through a Lamb mirror.

v2 here flips rank_genes to descending and leaves everything else
identical. If the hypothesis is right, v2 RGES should equal the reference
exactly (within float noise from the underlying merges).

Run from the project root:
    python src/validations/binchen2017/debug_rges_v2.py

If v2 matches, we promote the change into connectivity_score.py.
Otherwise the per-signature dump at the bottom shows which intermediate
(a_up / b_up / a_down / b_down) diverges, and we iterate.
"""
import sys
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent.parent
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import conf  # noqa: E402
from connectivity_score import (  # noqa: E402
    bin_chen_connectivity,
    rank_genes as rank_genes_v1,
    compute_KS,
    calculate_RGES,
    get_common_genes,
)
from loader import load_disease_signature, load_landmark_gene_ids  # noqa: E402


# ============================================================================
# Reference: Bin Chen's RGES for LINCS signatures 8726 / 8728 / 8729, scored
# against the LIHC disease signature. Values taken verbatim from
# data/BinChen2017/data/data/LIHC/all_lincs_score.csv (column `cmap_score`
# of that file == column `RGES` of lincs_score_0/1 for any overlapping id).
# ============================================================================
REFERENCE = {
    "8726": {"pert_iname": "AG-957",        "rges": -0.0390461497984907},
    "8728": {"pert_iname": "albendazole",   "rges":  0.187899879214054},
    "8729": {"pert_iname": "alclometasone", "rges":  0.103086619069829},
}


# ============================================================================
# v2 candidate
# ============================================================================
def rank_genes_v2(drug_signature, disease_signature,
                  drug_col_name, disease_col_name, id_col="gene_id"):
    """
    Lamb 2006 / Bin Chen 2017 convention.

    Reference (drug) profile is ordered DESCENDING by drug_col_name so that
    position 0 corresponds to the most upregulated gene in the drug
    signature. The disease genes' V values are then "position in
    descending-FC drug profile" — small V means up-by-drug, large V means
    down-by-drug, which is what Lamb's KS formula assumes.

    Only one character changed vs v1: ascending=False on the drug sort.
    """
    sorted_drug = (
        drug_signature
        .sort_values(by=drug_col_name, ascending=False)
        .reset_index(drop=True)
    )
    sorted_disease = disease_signature.sort_values(by=disease_col_name)
    merged = pd.merge(sorted_disease, sorted_drug.reset_index(), on=id_col)
    V = merged["index"].to_numpy()
    return V


def bin_chen_connectivity_v2(disease_signature, drug_signature,
                             rank_on="magnitude", id_col="gene_id",
                             drug_ranking_col_name=None,
                             disease_ranking_col_name=None,
                             return_intermediates=False):
    """
    v2 of bin_chen_connectivity. Identical to v1 except:
        - uses rank_genes_v2 (drug profile sorted DESCENDING, Lamb/Bin
          Chen convention)
        - returns RGES only (no Monte Carlo p-value) to keep the
          benchmark deterministic
        - optionally returns the KS intermediates for debugging
    """
    if drug_ranking_col_name is None or disease_ranking_col_name is None:
        if rank_on == "p_value":
            drug_ranking_col_name = drug_signature.columns[2]
            disease_ranking_col_name = disease_signature.columns[2]
        if rank_on == "magnitude":
            drug_ranking_col_name = drug_signature.columns[1]
            disease_ranking_col_name = disease_signature.columns[1]

    d_up   = disease_signature[disease_signature[disease_ranking_col_name] > 0]
    d_down = disease_signature[disease_signature[disease_ranking_col_name] < 0]

    V_up   = rank_genes_v2(drug_signature, d_up,   drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)
    V_down = rank_genes_v2(drug_signature, d_down, drug_ranking_col_name, disease_ranking_col_name, id_col=id_col)

    a_up,   b_up,   s_up,   r = compute_KS(d_up[id_col],   V_up,   drug_signature[id_col])
    a_down, b_down, s_down, _ = compute_KS(d_down[id_col], V_down, drug_signature[id_col])

    rges = calculate_RGES(a_up, a_down, b_up, b_down)

    if return_intermediates:
        return rges, {
            "s_up": s_up, "s_down": s_down, "r": r,
            "a_up": a_up, "b_up": b_up,
            "a_down": a_down, "b_down": b_down,
        }
    return rges


# ============================================================================
# v1 dump (current connectivity_score.py, with intermediates exposed)
# ============================================================================
def bin_chen_connectivity_v1_with_intermediates(disease_signature, drug_signature,
                                                rank_on="magnitude",
                                                id_col="gene_id"):
    """
    Re-runs the v1 logic but exposes a/b/s/r alongside the score. Mirrors
    bin_chen_connectivity() in connectivity_score.py without the Monte
    Carlo step.
    """
    if rank_on == "magnitude":
        drug_col = drug_signature.columns[1]
        dis_col  = disease_signature.columns[1]
    else:
        drug_col = drug_signature.columns[2]
        dis_col  = disease_signature.columns[2]

    d_up   = disease_signature[disease_signature[dis_col] > 0]
    d_down = disease_signature[disease_signature[dis_col] < 0]

    V_up   = rank_genes_v1(drug_signature, d_up,   drug_col, dis_col, id_col=id_col)
    V_down = rank_genes_v1(drug_signature, d_down, drug_col, dis_col, id_col=id_col)

    a_up,   b_up,   s_up,   r = compute_KS(d_up[id_col],   V_up,   drug_signature[id_col])
    a_down, b_down, s_down, _ = compute_KS(d_down[id_col], V_down, drug_signature[id_col])

    rges = calculate_RGES(a_up, a_down, b_up, b_down)
    return rges, {
        "s_up": s_up, "s_down": s_down, "r": r,
        "a_up": a_up, "b_up": b_up,
        "a_down": a_down, "b_down": b_down,
    }


# ============================================================================
# Main
# ============================================================================
def main():
    # Mirror the canonical Bin Chen RGES replication config:
    #   cs_mith=0, cs_on_LM=1, CS_METHOD='bin_chen', landmark_drug=True,
    #   landmark_disease=False (the 08_05_2026_10_26 run config)
    cell_line         = "HEPG2"
    landmark_drug     = True
    landmark_disease  = False
    cs_on_LM          = True

    DISEASE             = "LIHC"
    disease_run_name    = DISEASE + ("_LM" if landmark_disease else "")
    cell_line_run_name  = cell_line + ("_LM" if landmark_drug else "")

    CS_IN_DRUG    = conf.CS_DIR / "input" / "drug_signature_2025"    / cell_line_run_name
    CS_IN_DISEASE = conf.CS_DIR / "input" / "disease_signature_2025" / disease_run_name
    LINCS_BC_DATA = conf.LINCS_BC_DATA

    # --- load disease signature ---
    disease_sig = load_disease_signature(
        disease_run_name, cs_input_dir=CS_IN_DISEASE, mith=False
    )
    print(f"[load] disease signature: {len(disease_sig):,} rows from {CS_IN_DISEASE}")

    # filter disease to landmark (matches cs_on_LM=1 step in cs_batch)
    if cs_on_LM:
        lm_ids = [str(x) for x in load_landmark_gene_ids(LINCS_BC_DATA)]
        disease_sig = disease_sig[disease_sig["gene_id"].isin(lm_ids)].reset_index(drop=True)
        print(f"       after landmark filter: {len(disease_sig):,} rows")

    cols = ["gene_id", "DE_log2_FC", "adj.p.value"]
    disease_sig = disease_sig[cols]

    # --- iterate over the three test signatures ---
    print()
    header = f"{'sig':<6} {'drug':<15} {'v1':>22} {'v2':>22} {'ref':>22}   {'v2 - ref':>12}   {'v2 == -v1'}"
    print(header)
    print("-" * len(header))

    rows = []
    for sig_id, info in REFERENCE.items():
        drug_path = CS_IN_DRUG / f"{sig_id}_signature_gene_id.pkl"
        with open(drug_path, "rb") as f:
            drug_sig = pickle.load(f)

        # filter drug to landmark + landmark intersect already applied
        if cs_on_LM:
            drug_sig = drug_sig[drug_sig["gene_id"].isin(lm_ids)].reset_index(drop=True)

        # common-genes filtering, mirrors cs_batch.run_connectivity_score_drugs_batch
        d_idx, _x_idx = get_common_genes(disease_sig, drug_sig, id_col="gene_id")
        d_filt = disease_sig.iloc[d_idx].reset_index(drop=True)

        # v1
        rges_v1, intern_v1 = bin_chen_connectivity_v1_with_intermediates(
            d_filt, drug_sig[cols], rank_on="magnitude", id_col="gene_id"
        )
        # v2
        rges_v2, intern_v2 = bin_chen_connectivity_v2(
            d_filt, drug_sig[cols], rank_on="magnitude", id_col="gene_id",
            return_intermediates=True,
        )
        ref = info["rges"]
        delta_v2 = rges_v2 - ref
        is_anti  = abs(rges_v1 + rges_v2) < 1e-12

        print(f"{sig_id:<6} {info['pert_iname']:<15} "
              f"{rges_v1:>22.16f} {rges_v2:>22.16f} {ref:>22.16f}   "
              f"{delta_v2:>12.2e}   {'YES' if is_anti else 'NO'}")

        rows.append({
            "sig_id": sig_id, "drug": info["pert_iname"],
            "rges_v1": rges_v1, "rges_v2": rges_v2, "rges_ref": ref,
            "delta_v2": delta_v2,
            **{f"v1_{k}": v for k, v in intern_v1.items()},
            **{f"v2_{k}": v for k, v in intern_v2.items()},
        })

    # --- diagnostic dump ---
    print()
    print("---- per-signature intermediates ----")
    df = pd.DataFrame(rows)
    # Pretty-print the KS intermediates so deviations are visible.
    show_cols = ["sig_id", "drug",
                 "v1_s_up", "v1_s_down", "v1_r",
                 "v1_a_up", "v1_b_up", "v1_a_down", "v1_b_down",
                 "v2_a_up", "v2_b_up", "v2_a_down", "v2_b_down",
                 "rges_v1", "rges_v2", "rges_ref"]
    show_cols = [c for c in show_cols if c in df.columns]
    with pd.option_context("display.max_columns", None,
                           "display.width", 220,
                           "display.float_format", "{:.10f}".format):
        print(df[show_cols].to_string(index=False))

    # --- verdict ---
    print()
    max_abs_err_v2 = df["delta_v2"].abs().max()
    if max_abs_err_v2 < 1e-10:
        print(f"VERDICT: v2 reproduces Bin Chen reference RGES "
              f"(max |Δ| = {max_abs_err_v2:.2e}). Promote to "
              f"connectivity_score.py.")
    else:
        print(f"VERDICT: v2 does NOT match reference yet. "
              f"max |Δ| = {max_abs_err_v2:.2e}. See intermediates above.")


if __name__ == "__main__":
    main()
