#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Intersect the repoDB drug-indication pairs with the disease-signature inventory
(DiSignAtlas + CREEDS) on the UMLS CUI: which labelled pairs have a usable
disease signature, at each Psamples threshold.

Design
------
Stage 2, part 2 of the drug-disease-pairs validation (disease leg only). Reads
the repoDB table (as characterised by 02_) and the per-disease signature summary
written by 01_ (`creeds_disignatlas_merged_summary_by_disease_ge<x>Psamples.csv`)
for every x in PSAMPLES_THRESHOLDS, and joins them on the disease key.

Join key. repoDB indications are UMLS CUIs (`ind_id`). A signature disease
matches when one of its `umls_cuis` equals the repoDB CUI (match_method = cui).
Signature diseases with no CUI (mostly CREEDS) fall back to an exact match of the
normalised name against repoDB `ind_name` (match_method = name, same conservative
normaliser as 01_); these are flagged, never silently merged. No fuzzy or
ontology-based mapping is attempted: a CUI that differs only by granularity
(e.g. a subtype vs its parent) is a miss here, so every count is a lower bound.

Labels (distinct drug-disease pairs, as in 02_):
- approved      : >=1 repoDB record with status Approved.
- not_approved  : >=1 record Terminated / Withdrawn / Suspended and NO approved
                  record for the same pair.
- conflict      : the pair is both approved and carries a failed trial (e.g. a
                  failed trial on an already-approved indication). Reported, and
                  excluded from both classes.
A not-approved label is a clinical-trial outcome, NOT a biological negative: it
mixes efficacy/safety stops with administrative/recruitment stops. Separating
them (AACT stop reason) is a later stage.

The drug side is NOT joined here (DrugBank -> LINCS is a separate mapping, the
second bottleneck). Pair counts are therefore an upper bound on what NetCoS can
score: a pair still needs its drug in the LINCS drug universe.

Two views, as asked by Pulvirenti (22/09): the pooled multi-disease benchmark
(all surviving pairs) and the single-disease view (diseases with enough approved
AND not-approved pairs to validate within one disease, which avoids pooling CS
across diseases with different per-disease scale).

Inputs
------
- repoDB table from conf (REPODB_FILE), v2 `full.csv`.
- 01_ outputs `creeds_disignatlas_merged_summary_by_disease_ge<x>Psamples.csv`
  for each x in PSAMPLES_THRESHOLDS (missing thresholds are skipped with a
  warning; run 01_ with that MIN_REPLICATES to produce them).

Outputs (flat inside the slug; one file serves every x through the column
`max_psamples` = the largest threshold at which the disease still has a
signature; filter `max_psamples >= x` to get the ge<x>Psamples subset)
------
- repodb_signature_intersection_summary.tsv -- one row per x: signature
  diseases, matched repoDB diseases, pairs per label, not-approved pairs per
  trial status, distinct drugs, diseases with both labels, single-disease
  candidates.
- repodb_signature_pairs.tsv -- one row per distinct repoDB pair whose disease
  has a signature at x >= min(PSAMPLES_THRESHOLDS): drug, disease, label,
  not-approved statuses / trial accessions, match method, max_psamples, and the
  signature evidence of the disease at the smallest threshold.
- repodb_signature_by_disease.tsv -- one row per matched disease: pair counts
  per label, max_psamples, per-database signature evidence, and the flag
  single_disease_candidate (>= SINGLE_DISEASE_MIN_PAIRS in both classes).
- repodb_signature_intersection_panel_ge<x>Psamples.png -- panel A-D, x =
  SELECTED_PSAMPLES (A/B show every threshold; C/D the selected one).

Usage
-----
Run from the source directory (so `import conf` resolves), e.g. in Spyder:
    runfile('validations/drug_disease_pairs_validation/03_repodb_signature_intersection.py')

Metadata only: no connectivity score is computed, no run ID is produced.
"""

import re
import textwrap
import unicodedata
import warnings

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

from conf import REPODB_FILE, OOUT_DIR, IMG_DIR

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
IMG_OUT_DIR = IMG_DIR / "validations" / SLUG


def signature_summary_path(x):
    """Per-disease signature summary written by 01_ at threshold x."""
    return OUT_DIR / f"creeds_disignatlas_merged_summary_by_disease_ge{x}Psamples.csv"


SUMMARY_TSV = OUT_DIR / "repodb_signature_intersection_summary.tsv"
PAIRS_TSV = OUT_DIR / "repodb_signature_pairs.tsv"
BY_DISEASE_TSV = OUT_DIR / "repodb_signature_by_disease.tsv"

# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
PSAMPLES_THRESHOLDS = (3, 10, 20, 50)   # x values read from 01_ outputs
SELECTED_PSAMPLES = 10                  # x shown in panels C/D (file-name tag)
SINGLE_DISEASE_MIN_PAIRS = 10           # min approved AND not-approved pairs
LABEL_TOP_N = 20                        # diseases shown in panel D

APPROVED_STATUS = "Approved"
NOT_APPROVED_STATUSES = ("Terminated", "Withdrawn", "Suspended")

COL_DRUG_NAME = "drug_name"
COL_DRUG_ID = "drugbank_id"
COL_IND_NAME = "ind_name"
COL_IND_ID = "ind_id"
COL_TRIAL = "NCT"
COL_STATUS = "status"
COL_PHASE = "phase"

PANEL_TAG = f"ge{SELECTED_PSAMPLES}Psamples"
PANEL_PNG = IMG_OUT_DIR / f"repodb_signature_intersection_panel_{PANEL_TAG}.png"

C_APPROVED = "#2b6cb0"
C_NOT_APPROVED = "#c05621"
C_CONFLICT = "#718096"
C_NEUTRAL = "#4a5568"
GRID_COLOUR = "#d9dde3"


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _normalise_disease_name(value):
    """Same conservative normaliser as 01_ (case, accents, punctuation only)."""
    text = unicodedata.normalize("NFKD", value or "")
    text = text.encode("ascii", "ignore").decode("ascii").lower()
    text = text.replace("&", " and ").replace("'s", "s")
    return " ".join(re.findall(r"[a-z0-9]+", text))


def _join(values):
    return "|".join(sorted({str(v) for v in values if isinstance(v, str) and v}))


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_repodb(path=REPODB_FILE):
    """repoDB with 'NA' / '' as missing (same reader as 02_)."""
    last_err = None
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding=enc)
            break
        except (UnicodeDecodeError, LookupError) as err:  # pragma: no cover
            last_err = err
    else:  # pragma: no cover
        raise last_err
    df.columns = [c.strip() for c in df.columns]
    return df.replace({"NA": pd.NA, "": pd.NA})


def build_labelled_pairs(df):
    """
    One row per distinct (drug, disease) pair with its label
    (approved / not_approved / conflict) and the failed-trial details.
    """
    n_missing = int(df[[COL_DRUG_ID, COL_IND_ID]].isna().any(axis=1).sum())
    print(f"repoDB records dropped for a missing drug/disease id: {n_missing}")
    unknown = set(df[COL_STATUS].dropna()) - {APPROVED_STATUS, *NOT_APPROVED_STATUSES}
    if unknown:
        warnings.warn(f"repoDB statuses ignored (not in any class): {sorted(unknown)}")
    df = df.dropna(subset=[COL_DRUG_ID, COL_IND_ID])
    keys = [COL_DRUG_ID, COL_IND_ID]
    appr = df[df[COL_STATUS] == APPROVED_STATUS]
    fail = df[df[COL_STATUS].isin(NOT_APPROVED_STATUSES)]

    appr_keys = set(map(tuple, appr[keys].drop_duplicates().to_numpy()))
    fail_agg = (fail.groupby(keys)
                    .agg(n_failed_records=(COL_STATUS, "size"),
                         failed_statuses=(COL_STATUS, _join),
                         failed_trials=(COL_TRIAL, _join),
                         failed_phases=(COL_PHASE, _join))
                    .reset_index())
    fail_keys = set(map(tuple, fail_agg[keys].to_numpy()))

    names = df.groupby(keys)[[COL_DRUG_NAME, COL_IND_NAME]].first().reset_index()
    pairs = names[names.apply(lambda r: (r[COL_DRUG_ID], r[COL_IND_ID]) in
                              (appr_keys | fail_keys), axis=1)]
    pairs = pairs.merge(fail_agg, on=keys, how="left")
    pairs["n_failed_records"] = pairs["n_failed_records"].fillna(0).astype(int)

    def _label(r):
        k = (r[COL_DRUG_ID], r[COL_IND_ID])
        a, f = k in appr_keys, k in fail_keys
        return "conflict" if (a and f) else ("approved" if a else "not_approved")

    pairs["label"] = pairs.apply(_label, axis=1)
    return pairs.fillna("")


def load_signature_diseases(repodb_names):
    """
    For each available threshold, the set of repoDB CUIs with a signature and
    how they matched. Returns (per_x dict of DataFrames, list of thresholds read).

    Each per-x frame has one row per (repoDB CUI, source_database) with the
    signature evidence of that database's disease entry.
    """
    name_to_cui = {}
    for cui, name in repodb_names.items():
        name_to_cui.setdefault(_normalise_disease_name(name), set()).add(cui)

    per_x, read, audit = {}, [], []
    for x in PSAMPLES_THRESHOLDS:
        path = signature_summary_path(x)
        if not path.exists():
            warnings.warn(f"missing {path.name}: run 01_ with MIN_REPLICATES={x}")
            continue
        s = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
        rows = []
        for _, r in s.iterrows():
            cuis = [c for c in r["umls_cuis"].split("|") if c]
            hits = [(c, "cui") for c in cuis if c in repodb_names]
            name_cuis = name_to_cui.get(_normalise_disease_name(r["disease"]), set())
            if not cuis:  # name fallback only when the signature has no CUI
                hits += [(c, "name") for c in name_cuis]
            elif not hits and name_cuis:
                # CUI misses repoDB but the name matches under another CUI:
                # not merged automatically, printed for manual review.
                audit.append((x, r["source_database"], r["disease"],
                              "|".join(cuis), "|".join(sorted(name_cuis))))
            for cui, method in hits:
                rows.append({
                    COL_IND_ID: cui,
                    "source_database": r["source_database"],
                    "signature_disease": r["disease"],
                    "match_method": method,
                    "n_datasets": int(r["n_datasets"] or 0),
                    "n_gse": int(r["n_gse"] or 0),
                    "biggest_ctrl_case": r["biggest_ctrl_case"],
                    "gse_accessions": r["gse_accessions"],
                })
        per_x[x] = (pd.DataFrame(rows), len(s))
        read.append(x)
    if not read:
        raise FileNotFoundError("no 01_ signature summaries found")
    # max_psamples assumes nested sets (a disease at x=50 is also at x=10).
    for lo, hi in zip(read, read[1:]):
        extra = set(per_x[hi][0][COL_IND_ID]) - set(per_x[lo][0][COL_IND_ID])
        assert not extra, f"CUIs at >= {hi} but not at >= {lo}: {sorted(extra)}"
    if audit:
        print("\n--- CUI mismatch, name matches repoDB (review in 01_; NOT merged) ---")
        for row in sorted(set(audit)):
            print("  ge{}Psamples  {:<11} {:<40} sig CUI {}  repoDB CUI {}".format(*row))
    return per_x, read


# --------------------------------------------------------------------------- #
# Intersection
# --------------------------------------------------------------------------- #
def _signature_evidence(hits, x0):
    """Collapse the (CUI, database) hits into one evidence row per CUI."""
    hits = hits.drop_duplicates()   # guard: same entry hitting a CUI twice

    def per_cui(g):
        out = {"match_method": _join(g["match_method"]),
               "signature_databases": _join(g["source_database"]),
               "signature_disease_names": _join(g["signature_disease"])}
        for db, tag in (("DiSignAtlas", "dsa"), ("CREEDS", "creeds")):
            sub = g[g["source_database"] == db]
            out[f"{tag}_n_datasets"] = int(sub["n_datasets"].sum())
            out[f"{tag}_biggest_ctrl_case"] = _join(sub["biggest_ctrl_case"])
        out["gse_accessions"] = _join(
            a for v in g["gse_accessions"] for a in v.split("|"))
        return pd.Series(out)
    ev = hits.groupby(COL_IND_ID).apply(per_cui, include_groups=False)
    # Evidence is described at the smallest threshold read (x0), say so.
    tag = f"_ge{x0}Psamples"
    return ev.rename(columns={c: c + tag for c in ev.columns
                              if c.startswith(("dsa_", "creeds_", "gse_"))})


def intersect(pairs, per_x, thresholds):
    """Pairs and diseases with a signature, annotated with max_psamples."""
    max_ps = {}
    for x in thresholds:
        for cui in set(per_x[x][0][COL_IND_ID]):
            max_ps[cui] = max(max_ps.get(cui, 0), x)

    x0 = min(thresholds)
    evidence = _signature_evidence(per_x[x0][0], x0)
    sig_pairs = pairs[pairs[COL_IND_ID].isin(max_ps)].copy()
    sig_pairs["max_psamples"] = sig_pairs[COL_IND_ID].map(max_ps)
    sig_pairs = sig_pairs.merge(evidence, left_on=COL_IND_ID,
                                right_index=True, how="left")

    counts = (sig_pairs.pivot_table(index=COL_IND_ID, columns="label",
                                    values=COL_DRUG_ID, aggfunc="count",
                                    fill_value=0)
                       .reindex(columns=["approved", "not_approved", "conflict"],
                                fill_value=0)
                       .add_prefix("n_").add_suffix("_pairs"))
    by_dis = (sig_pairs.groupby(COL_IND_ID)
                       .agg(ind_name=(COL_IND_NAME, "first"),
                            max_psamples=("max_psamples", "first"))
                       .join(counts).join(evidence))
    by_dis["n_min_class"] = by_dis[["n_approved_pairs",
                                    "n_not_approved_pairs"]].min(axis=1)
    by_dis["single_disease_candidate"] = (
        by_dis["n_min_class"] >= SINGLE_DISEASE_MIN_PAIRS).astype(int)
    by_dis = by_dis.sort_values(["n_min_class", "n_approved_pairs"],
                                ascending=False).reset_index()
    return sig_pairs, by_dis


def summarise(sig_pairs, by_dis, per_x, thresholds, pairs):
    rows = []
    for x in thresholds:
        hits, n_sig_entries = per_x[x]
        p = sig_pairs[sig_pairs["max_psamples"] >= x]
        d = by_dis[by_dis["max_psamples"] >= x]
        napp = p[p["label"] == "not_approved"]
        row = {
            "min_psamples": x,
            "signature_disease_entries": n_sig_entries,
            "repodb_diseases_matched": d.shape[0],
            "matched_by_cui": int(hits.loc[hits.match_method == "cui", COL_IND_ID].nunique()),
            "matched_by_name_only": int(len(
                set(hits.loc[hits.match_method == "name", COL_IND_ID])
                - set(hits.loc[hits.match_method == "cui", COL_IND_ID]))),
            "pairs_total": p.shape[0],
            "pairs_approved": int((p["label"] == "approved").sum()),
            "pairs_not_approved": napp.shape[0],
            "pairs_conflict": int((p["label"] == "conflict").sum()),
        }
        for st in NOT_APPROVED_STATUSES:  # a pair may carry several statuses
            row[f"not_approved_with_{st.lower()}"] = int(
                napp["failed_statuses"].str.split("|").apply(lambda t: st in t).sum())
        row.update({
            "distinct_drugs": p[COL_DRUG_ID].nunique(),
            "distinct_drugs_approved": p.loc[p.label == "approved", COL_DRUG_ID].nunique(),
            "distinct_drugs_not_approved": napp[COL_DRUG_ID].nunique(),
            "diseases_with_approved": int((d.n_approved_pairs > 0).sum()),
            "diseases_with_not_approved": int((d.n_not_approved_pairs > 0).sum()),
            # conflict pairs excluded from both classes ...
            "diseases_with_both": int(((d.n_approved_pairs > 0) &
                                       (d.n_not_approved_pairs > 0)).sum()),
            # ... sensitivity: conflict counted as approved
            "diseases_with_both_conflict_as_approved": int((
                ((d.n_approved_pairs + d.n_conflict_pairs) > 0) &
                (d.n_not_approved_pairs > 0)).sum()),
            f"single_disease_candidates_ge{SINGLE_DISEASE_MIN_PAIRS}": int(
                d.single_disease_candidate.sum()),
            "share_of_repodb_pairs": round(p.shape[0] / len(pairs), 4),
        })
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def _style_axis(ax, grid_axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID_COLOUR)
    ax.grid(axis=grid_axis, color=GRID_COLOUR, linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    fmt = StrMethodFormatter("{x:,.0f}")
    (ax.yaxis if grid_axis == "y" else ax.xaxis).set_major_formatter(fmt)


def _panel_letter(ax, letter):
    ax.text(-0.10, 1.06, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="left")


def _short(name, width=40):
    return textwrap.shorten(" ".join(str(name).split()), width=width, placeholder="…")


def plot_panel(summary, by_dis, total_pairs):
    x = SELECTED_PSAMPLES
    d = by_dis[by_dis["max_psamples"] >= x]
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    labels = [f"≥{t}" for t in summary.min_psamples]
    pos = range(len(labels))

    # A -- labelled pairs with a disease signature, per threshold
    ax = axes[0, 0]
    w = 0.27
    series = (("pairs_approved", "approved", C_APPROVED, -w),
              ("pairs_not_approved", "not approved", C_NOT_APPROVED, 0),
              ("pairs_conflict", "conflict (excluded)", C_CONFLICT, w))
    for col, lab, colour, off in series:
        bars = ax.bar([p + off for p in pos], summary[col], w, color=colour, label=lab)
        for b, v in zip(bars, summary[col]):
            ax.annotate(f"{v:,}", (b.get_x() + b.get_width() / 2, b.get_height()),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8)
    ax.set_xticks(list(pos), labels)
    ax.set_xlabel("Psamples threshold (paired case/control samples)")
    ax.set_ylabel("distinct repoDB drug–disease pairs")
    ax.set_title(f"Pairs whose disease has a signature\n(of {total_pairs:,} "
                 "distinct repoDB pairs)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    _style_axis(ax)
    _panel_letter(ax, "A")

    # B -- diseases, per threshold
    ax = axes[0, 1]
    series = (("repodb_diseases_matched", "with a signature", C_NEUTRAL, -w),
              ("diseases_with_both", "both classes", C_NOT_APPROVED, 0),
              (f"single_disease_candidates_ge{SINGLE_DISEASE_MIN_PAIRS}",
               f"≥{SINGLE_DISEASE_MIN_PAIRS} pairs in both classes", C_APPROVED, w))
    for col, lab, colour, off in series:
        bars = ax.bar([p + off for p in pos], summary[col], w, color=colour, label=lab)
        for b, v in zip(bars, summary[col]):
            ax.annotate(f"{v:,}", (b.get_x() + b.get_width() / 2, b.get_height()),
                        xytext=(0, 3), textcoords="offset points", ha="center",
                        fontsize=8)
    ax.set_xticks(list(pos), labels)
    ax.set_xlabel("Psamples threshold (paired case/control samples)")
    ax.set_ylabel("repoDB diseases (UMLS CUI)")
    ax.set_title("repoDB diseases with a signature", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    _style_axis(ax)
    _panel_letter(ax, "B")

    # C -- approved vs not approved per disease, selected threshold
    ax = axes[1, 0]
    cand = d["single_disease_candidate"] == 1
    ax.scatter(d.loc[~cand, "n_approved_pairs"], d.loc[~cand, "n_not_approved_pairs"],
               s=18, color=C_CONFLICT, alpha=0.6, label="other diseases")
    ax.scatter(d.loc[cand, "n_approved_pairs"], d.loc[cand, "n_not_approved_pairs"],
               s=30, color=C_NOT_APPROVED, label=f"≥{SINGLE_DISEASE_MIN_PAIRS} in both")
    for _, r in d[cand].head(8).iterrows():
        ax.annotate(_short(r["ind_name"], 26),
                    (r["n_approved_pairs"], r["n_not_approved_pairs"]),
                    xytext=(4, 3), textcoords="offset points", fontsize=7)
    ax.axvline(SINGLE_DISEASE_MIN_PAIRS, color=GRID_COLOUR, ls="--", lw=1)
    ax.axhline(SINGLE_DISEASE_MIN_PAIRS, color=GRID_COLOUR, ls="--", lw=1)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlim(left=-0.3)
    ax.set_ylim(bottom=-0.3)
    ax.set_xlabel("approved pairs per disease")
    ax.set_ylabel("not-approved pairs per disease")
    ax.set_title(f"Per-disease label balance (≥{x} Psamples, {d.shape[0]} diseases)",
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    _style_axis(ax, "both")
    _panel_letter(ax, "C")

    # D -- top single-disease candidates, selected threshold
    ax = axes[1, 1]
    top = d.head(LABEL_TOP_N).iloc[::-1]
    yy = range(len(top))
    ax.barh([y + 0.2 for y in yy], top["n_approved_pairs"], 0.4,
            color=C_APPROVED, label="approved")
    ax.barh([y - 0.2 for y in yy], top["n_not_approved_pairs"], 0.4,
            color=C_NOT_APPROVED, label="not approved")
    ax.set_yticks(list(yy), [f"{_short(n, 34)}  (≥{m})" for n, m in
                             zip(top["ind_name"], top["max_psamples"])], fontsize=8)
    ax.set_xlabel("distinct pairs")
    ax.axvline(SINGLE_DISEASE_MIN_PAIRS, color=GRID_COLOUR, ls="--", lw=1)
    ax.set_title(f"Top {LABEL_TOP_N} diseases by the smaller class (≥{x} Psamples)\n"
                 f"dashed = {SINGLE_DISEASE_MIN_PAIRS} pairs; in brackets the highest "
                 "threshold reached", fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    _style_axis(ax, "x")
    _panel_letter(ax, "D")

    fig.suptitle("repoDB × disease signatures (DiSignAtlas + CREEDS), join on UMLS CUI"
                 " — disease leg only, drug side (LINCS) not yet joined",
                 fontsize=13, y=1.0)
    fig.tight_layout()
    PANEL_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PANEL_PNG, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
PAIR_COLUMNS = [
    COL_DRUG_ID, COL_DRUG_NAME, COL_IND_ID, COL_IND_NAME, "label",
    "max_psamples", "match_method", "failed_statuses", "n_failed_records",
    "failed_trials", "failed_phases", "signature_databases",
    "signature_disease_names",
]  # + the per-database evidence columns (suffix _ge<x0>Psamples), added in main


def main():
    repodb = load_repodb()
    pairs = build_labelled_pairs(repodb)
    names = (repodb.dropna(subset=[COL_IND_ID])
                   .groupby(COL_IND_ID)[COL_IND_NAME].first().to_dict())
    per_x, thresholds = load_signature_diseases(names)
    sig_pairs, by_dis = intersect(pairs, per_x, thresholds)
    summary = summarise(sig_pairs, by_dis, per_x, thresholds, pairs)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    evid_cols = [c for c in sig_pairs.columns
                 if c.startswith(("dsa_", "creeds_"))]
    (sig_pairs[PAIR_COLUMNS + evid_cols]
        .sort_values(["max_psamples", COL_IND_NAME, "label", COL_DRUG_NAME],
                     ascending=[False, True, True, True])
        .to_csv(PAIRS_TSV, sep="\t", index=False))
    by_dis.to_csv(BY_DISEASE_TSV, sep="\t", index=False)
    if SELECTED_PSAMPLES in thresholds:
        plot_panel(summary, by_dis, len(pairs))
    else:
        warnings.warn(f"SELECTED_PSAMPLES={SELECTED_PSAMPLES} not available; no panel")

    print(f"repoDB distinct pairs: {len(pairs):,} "
          f"({pairs.label.value_counts().to_dict()})")
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print(summary.T.to_string(header=False))
        cols = [COL_IND_ID, "ind_name", "max_psamples", "n_approved_pairs",
                "n_not_approved_pairs", "n_conflict_pairs", "signature_databases"]
        print("\n--- single-disease candidates (all thresholds) ---")
        print(by_dis[by_dis.single_disease_candidate == 1][cols].to_string(index=False))


if __name__ == "__main__":
    main()
