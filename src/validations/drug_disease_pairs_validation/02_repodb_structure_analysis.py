#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Characterise the repoDB drug-indication table: how the drug-disease pairs are
built, on which identifiers, and how the outcome labels are distributed.

Design
------
Reads the repoDB gold-standard table (one row per drug-indication-trial record)
and describes its structure without joining any external resource. A row is NOT
a distinct drug-disease pair: an approved indication appears once, but a failed
indication appears once per clinical trial, so the same drug-disease pair can
carry several not-approved rows (different trial accessions). Every count here is
therefore reported at two levels -- rows (records) and distinct drug-disease
pairs -- and the pair is the unit that matters for a benchmark.

It answers three questions about the table itself: (1) what a pair is keyed on --
a DrugBank drug id and a UMLS-CUI indication id -- and how complete and clean
those keys are; (2) how the outcome labels are distributed (Approved vs the three
not-approved clinical-trial statuses Terminated / Withdrawn / Suspended), at both
row and distinct-pair level, and how many not-approved rows carry a trial
accession that can later be linked to a stop reason; (3) the shape of the table
by disease and by drug (distinct approved / not-approved pairs per disease and per
drug), and how many diseases carry both an approved pair and at least one
not-approved pair within repoDB alone.

This is a structural characterisation only. It deliberately does not intersect
repoDB with disease-signature availability and does not decide which pairs are
usable positives or negatives. In particular a not-approved status is NOT treated
as a biological negative here: the not-approved class mixes efficacy/safety stops
with administrative and recruitment stops, and separating them is a downstream
stage (trial stop-reason via AACT + a released classifier), not this script.
The disease key emitted here (UMLS CUI) is the join key the next stage uses
against the disease-signature inventory.

Inputs
------
- The repoDB table, resolved from conf (REPODB_FILE). Downloaded once from the
  repoDB distribution (the "Download the full repoDB Dataset" export of the repoDB
  web application) into that location; it is not regenerated here. repoDB is the
  gold standard from Brown & Patel, Sci. Data 4:170029 (2017), maintained with
  later DrugCentral / AACT / UMLS releases.
- Analysis parameters below (the set of statuses treated as not-approved).

Outputs (written to the results dir resolved from conf, flat inside the slug)
------
- repodb_schema_summary.tsv        -- one row per column: dtype, non-null count,
                                      distinct-value count, and an example value.
- repodb_status_distribution.tsv   -- pair counts per outcome status, with the
                                      approved / not-approved grouping and the
                                      count of pairs carrying a trial accession.
- repodb_pairs_by_disease.tsv      -- one row per disease (UMLS CUI): approved,
                                      not-approved and total pair counts.
- repodb_pairs_by_drug.tsv         -- one row per drug (DrugBank id): approved,
                                      not-approved and total pair counts.
- repodb_both_labels_diseases.tsv  -- diseases carrying >=1 approved AND >=1
                                      not-approved pair within repoDB (signatures
                                      not yet intersected -- an upper structural
                                      bound, not the benchmark ceiling).
- repodb_overview.tsv              -- headline counts (pairs, distinct drugs,
                                      distinct diseases, key completeness).
- Four annotated PNG figures: status distribution, pairs-per-disease and
  pairs-per-drug distributions, and approved-vs-not-approved pairs per disease.

Usage
-----
Run from the source directory (so `import conf` resolves), e.g. in Spyder:
    runfile('validations/drug_disease_pairs_validation/02_repodb_structure_analysis.py')

Convention
----------
ind_id is a UMLS CUI (e.g. C0272275); the drug id column (COL_DRUG_ID) is a
DrugBank id (e.g. DB00001). A row is a drug-indication-trial record, not a
distinct pair: pair-level counts
deduplicate on (drug_id, ind_id) within a label. The not-approved class =
{Terminated, Withdrawn, Suspended}; it is a clinical-trial outcome, not an
efficacy verdict.
"""

import csv as _csv
import unicodedata

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

from conf import REPODB_FILE, OOUT_DIR, IMG_DIR

# --------------------------------------------------------------------------- #
# Paths (flat inside the canonical slug, mirrored across other_outputs / imgs)
# --------------------------------------------------------------------------- #
SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
IMG_OUT_DIR = IMG_DIR / "validations" / SLUG

SCHEMA_TSV = OUT_DIR / "repodb_schema_summary.tsv"
STATUS_TSV = OUT_DIR / "repodb_status_distribution.tsv"
BY_DISEASE_TSV = OUT_DIR / "repodb_pairs_by_disease.tsv"
BY_DRUG_TSV = OUT_DIR / "repodb_pairs_by_drug.tsv"
BOTH_LABELS_TSV = OUT_DIR / "repodb_both_labels_diseases.tsv"
OVERVIEW_TSV = OUT_DIR / "repodb_overview.tsv"

STATUS_PNG = IMG_OUT_DIR / "repodb_status_distribution.png"
DISEASE_BUCKETS_PNG = IMG_OUT_DIR / "repodb_pairs_per_disease.png"
DRUG_BUCKETS_PNG = IMG_OUT_DIR / "repodb_pairs_per_drug.png"
BOTH_LABELS_PNG = IMG_OUT_DIR / "repodb_approved_vs_notapproved_per_disease.png"

# --------------------------------------------------------------------------- #
# Analysis parameters
# --------------------------------------------------------------------------- #
APPROVED_STATUS = "Approved"
NOT_APPROVED_STATUSES = ("Terminated", "Withdrawn", "Suspended")

# In the approved-vs-not-approved scatter, name diseases carrying strictly more
# than this many pairs in BOTH classes.
LABEL_MIN_PAIRS = 10

# Column roles in the repoDB table.
COL_DRUG_NAME = "drug_name"
COL_DRUG_ID = "drugbank_id" #    # DrugBank id,  column name  for : RepoDB v1 "drug_id"
COL_IND_NAME = "ind_name"
COL_IND_ID = "ind_id"         # UMLS CUI  <-- the disease join key
COL_TRIAL = "NCT"             # clinical-trial accession (blank for approved rows)
COL_STATUS = "status"
COL_PHASE = "phase"
COL_DETAIL = "DetailedStatus"


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_repodb(path):
    """
    Load the repoDB table as a DataFrame.

    The distributed file is comma-separated with a header and encodes empty
    trial accessions / phases as the string 'NA'. We read with an encoding
    fallback (some drug/indication names carry non-UTF-8 bytes) and normalise
    'NA' to missing so key-completeness counts are honest.
    """
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
    # Treat the literal 'NA' and empty strings as missing.
    df = df.replace({"NA": pd.NA, "": pd.NA})
    return df


def _norm_name(name):
    """Lower-case, strip accents and collapse whitespace for a display fallback."""
    if not isinstance(name, str):
        return ""
    n = unicodedata.normalize("NFKD", name)
    n = "".join(c for c in n if not unicodedata.combining(c))
    return " ".join(n.lower().split())


# --------------------------------------------------------------------------- #
# Structural summaries
# --------------------------------------------------------------------------- #
def schema_summary(df):
    """One row per column: non-null count, distinct-value count, example value."""
    rows = [("column", "n_nonnull", "n_distinct", "example")]
    for col in df.columns:
        nonnull = df[col].notna().sum()
        ndist = df[col].dropna().nunique()
        ex = df[col].dropna().iloc[0] if nonnull else ""
        rows.append((col, int(nonnull), int(ndist), str(ex)))
    return rows


def _distinct_pairs(sub):
    """The distinct (drug_id, ind_id) pairs present in a row subset."""
    return sub[[COL_DRUG_ID, COL_IND_ID]].drop_duplicates()


def status_distribution(df):
    """
    Record counts per outcome status, plus the approved / not-approved grouping
    at both row and distinct-pair level, and how many rows in each carry a
    linkable trial accession.
    """
    counts = df[COL_STATUS].value_counts(dropna=False)
    rows = [("status", "n_rows", "n_distinct_pairs", "n_rows_with_trial_accession")]
    for status, n in counts.items():
        sub = df[df[COL_STATUS] == status]
        rows.append((str(status), int(n), int(_distinct_pairs(sub).shape[0]),
                     int(sub[COL_TRIAL].notna().sum())))

    approved = df[df[COL_STATUS] == APPROVED_STATUS]
    not_appr = df[df[COL_STATUS].isin(NOT_APPROVED_STATUSES)]
    rows.append(("__APPROVED_TOTAL__", int(len(approved)),
                 int(_distinct_pairs(approved).shape[0]),
                 int(approved[COL_TRIAL].notna().sum())))
    rows.append(("__NOT_APPROVED_TOTAL__", int(len(not_appr)),
                 int(_distinct_pairs(not_appr).shape[0]),
                 int(not_appr[COL_TRIAL].notna().sum())))
    return rows, approved, not_appr


def _per_key_counts(df, key_id, key_name, approved, not_appr):
    """
    One row per distinct key (disease CUI or drug id) with the number of DISTINCT
    approved / not-approved drug-disease pairs and a display name, sorted by total
    descending. Pairs are deduplicated within each label, so a pair failing in
    several trials is counted once (n_not_approved is pairs, not trials).
    """
    appr_c = _distinct_pairs(approved).groupby(key_id).size()
    napr_c = _distinct_pairs(not_appr).groupby(key_id).size()
    names = (df.dropna(subset=[key_id])
               .groupby(key_id)[key_name].first())
    keys = names.index

    header = (key_id, key_name, "n_approved_pairs", "n_not_approved_pairs",
              "n_total_pairs")
    out = []
    for k in keys:
        na = int(appr_c.get(k, 0))
        nn = int(napr_c.get(k, 0))
        out.append((k, str(names.get(k, "")), na, nn, na + nn))
    out.sort(key=lambda r: (-r[4], r[0]))
    return [header] + out


def both_labels_diseases(by_disease_rows):
    """
    Diseases carrying >=1 approved AND >=1 not-approved pair within repoDB.

    This is a repoDB-internal, signature-free count: an upper structural bound on
    the eventual benchmark ceiling, which additionally requires an admissible
    disease signature (intersection with the inventory -- the next stage).
    """
    header = by_disease_rows[0]
    body = [r for r in by_disease_rows[1:] if r[2] >= 1 and r[3] >= 1]
    body.sort(key=lambda r: (-r[4], r[0]))
    return [header] + body


def overview(df, approved, not_appr, both_rows):
    """Headline counts + key completeness, one metric per row."""
    rows = [
        ("metric", "value"),
        ("n_rows", int(len(df))),
        ("n_distinct_drug_disease_pairs",
         int(_distinct_pairs(df).shape[0])),
        ("n_distinct_drugs_drugbank", int(df[COL_DRUG_ID].dropna().nunique())),
        ("n_distinct_diseases_umls_cui", int(df[COL_IND_ID].dropna().nunique())),
        ("n_approved_rows", int(len(approved))),
        ("n_approved_pairs", int(_distinct_pairs(approved).shape[0])),
        ("n_not_approved_rows", int(len(not_appr))),
        ("n_not_approved_pairs", int(_distinct_pairs(not_appr).shape[0])),
        ("n_diseases_with_approved", int(approved[COL_IND_ID].dropna().nunique())),
        ("n_diseases_with_not_approved", int(not_appr[COL_IND_ID].dropna().nunique())),
        ("n_diseases_with_both_repodb_internal", len(both_rows) - 1),
        ("drug_id_completeness", f"{df[COL_DRUG_ID].notna().mean():.4f}"),
        ("ind_id_completeness", f"{df[COL_IND_ID].notna().mean():.4f}"),
        ("n_not_approved_rows_with_trial_accession",
         int(not_appr[COL_TRIAL].notna().sum())),
    ]
    return rows


# --------------------------------------------------------------------------- #
# Writing
# --------------------------------------------------------------------------- #
def _write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerows(rows)


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _savefig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _bucket(counts, edges):
    """Bucket a Series of integer counts into labelled ranges."""
    labels, values = [], []
    for lo, hi in edges:
        if hi is None:
            labels.append(f">={lo}")
            values.append(int((counts >= lo).sum()))
        elif lo == hi:
            labels.append(str(lo))
            values.append(int((counts == lo).sum()))
        else:
            labels.append(f"{lo}-{hi}")
            values.append(int(((counts >= lo) & (counts <= hi)).sum()))
    return labels, values


def fig_status(status_rows):
    # Plot DISTINCT drug-disease pairs per status (column n_distinct_pairs),
    # not records: a pair can recur across records (a failed pair over several
    # trials; an approved indication under several drug-name synonyms sharing one
    # DrugBank id), so record counts are higher. The pair is the benchmark unit.
    real = [r for r in status_rows[1:] if not r[0].startswith("__")]
    labels = [r[0] for r in real]
    values = [r[2] for r in real]
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, values, color=["#3d6b9c", "#c05a3e", "#c88a2e", "#7a7a7a"][:len(labels)])
    ax.bar_label(bars, padding=2, fontsize=9)
    ax.set_ylabel("distinct drug-disease pairs")
    ax.set_title("repoDB outcome-status distribution")
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    fig.tight_layout()
    _savefig(fig, STATUS_PNG)


def fig_buckets(counts, title, path):
    edges = [(1, 1), (2, 4), (5, 9), (10, 19), (20, None)]
    labels, values = _bucket(counts, edges)
    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, values, color="#3d6b9c")
    ax.bar_label(bars, padding=2, fontsize=9)
    ax.set_xlabel("pairs per key")
    ax.set_ylabel("number of keys")
    ax.set_title(title)
    fig.tight_layout()
    _savefig(fig, path)


def fig_both(by_disease_rows):
    body = by_disease_rows[1:]
    appr = [r[2] for r in body]
    napr = [r[3] for r in body]
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(appr, napr, s=14, alpha=0.4, color="#3d6b9c", edgecolor="none")
    lim = max(max(appr, default=1), max(napr, default=1)) + 1
    ax.plot([0, lim], [0, lim], color="#bbbbbb", lw=0.8, ls="--")

    # Name the diseases that are dense in BOTH classes (> LABEL_MIN pairs each):
    # these are the well-populated both-labels diseases worth calling out.
    labelled = [r for r in body
                if r[2] > LABEL_MIN_PAIRS and r[3] > LABEL_MIN_PAIRS]
    ax.scatter([r[2] for r in labelled], [r[3] for r in labelled],
               s=32, color="#c05a3e", edgecolor="white", linewidth=0.5, zorder=3)
    for r in labelled:
        ax.annotate(r[1], (r[2], r[3]), xytext=(6, 4), textcoords="offset points",
                    fontsize=8, color="#2b2b2b", ha="left", va="bottom", zorder=4)

    ax.set_xlabel("approved pairs (per disease)")
    ax.set_ylabel("not-approved pairs (per disease)")
    ax.set_title("repoDB: approved vs not-approved pairs per disease")
    ax.set_xlim(-0.5, lim)
    ax.set_ylim(-0.5, lim)
    fig.tight_layout()
    _savefig(fig, BOTH_LABELS_PNG)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def main():
    df = load_repodb(REPODB_FILE)

    schema_rows = schema_summary(df)
    status_rows, approved, not_appr = status_distribution(df)
    by_disease = _per_key_counts(df, COL_IND_ID, COL_IND_NAME, approved, not_appr)
    by_drug = _per_key_counts(df, COL_DRUG_ID, COL_DRUG_NAME, approved, not_appr)
    both_rows = both_labels_diseases(by_disease)
    overview_rows = overview(df, approved, not_appr, both_rows)

    _write_tsv(SCHEMA_TSV, schema_rows)
    _write_tsv(STATUS_TSV, status_rows)
    _write_tsv(BY_DISEASE_TSV, by_disease)
    _write_tsv(BY_DRUG_TSV, by_drug)
    _write_tsv(BOTH_LABELS_TSV, both_rows)
    _write_tsv(OVERVIEW_TSV, overview_rows)

    disease_counts = pd.Series([r[4] for r in by_disease[1:]])
    drug_counts = pd.Series([r[4] for r in by_drug[1:]])

    fig_status(status_rows)
    fig_buckets(disease_counts, "repoDB: pairs per disease", DISEASE_BUCKETS_PNG)
    fig_buckets(drug_counts, "repoDB: pairs per drug", DRUG_BUCKETS_PNG)
    fig_both(by_disease)

    # Console recap (kept out of the written artifacts).
    print(f"repoDB: {len(df):,} rows | "
          f"{_distinct_pairs(df).shape[0]:,} distinct drug-disease pairs | "
          f"{df[COL_DRUG_ID].dropna().nunique():,} drugs | "
          f"{df[COL_IND_ID].dropna().nunique():,} diseases (UMLS CUI)")
    print(f"  approved pairs={_distinct_pairs(approved).shape[0]:,}  "
          f"not-approved pairs={_distinct_pairs(not_appr).shape[0]:,} "
          f"({'/'.join(NOT_APPROVED_STATUSES)}; from {len(not_appr):,} trial rows)")
    print(f"  diseases with both labels (repoDB-internal, no signatures): "
          f"{len(both_rows) - 1:,}")


if __name__ == "__main__":
    main()
