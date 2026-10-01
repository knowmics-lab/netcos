#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build and compare human disease-signature inventories from DiSignAtlas and CREEDS.

Design
------
Reads the DiSignAtlas "Disease information of All Datasets" table (one row per
dataset) and produces a filtered inventory of case/control transcriptomic
signatures usable as disease signatures, plus a summary of how many distinct
diseases survive each filter (the signature-availability ceiling). Nothing is
computed on expression values here; this operates on the catalogue metadata only.

The catalogue gives, per dataset, the GEO/ArrayExpress/TCGA accession, the
disease and its DisGeNET concept id, the assay type, the organism, and the
control/case sample counts (as "n_control|n_case"). Per-sample accessions
(the GSM split) are NOT in this table; they live on each dataset's record page
and are fetched separately when a signature is actually built.

Filters, in order: organism == Homo sapiens; library_strategy in the requested
set (default bulk = RNA-Seq + Microarray); n_control >= MIN_REPLICATES and
n_case >= MIN_REPLICATES. Distinct diseases are counted on the DisGeNET/UMLS
id, with a normalised-name fallback when that id is missing.

MIN_REPLICATES (x) is the single Psamples threshold that defines the filtered
inventory used by every downstream table and figure. Every output built on the
filtered inventory carries the tag "ge<x>Psamples" in its file name and states
the threshold in its figure title, so outputs produced with different values
of x coexist side by side. The signature-ceiling table and figure are the
exception: they report every threshold in CEILING_THRESHOLDS at once (plus x)
and therefore carry no tag.

Inputs
------
- The catalogue CSV, resolved from conf (DATA_DIR / 'DiSignAtlas' /
  'dis_info_datasets.csv'). Download it once from the DiSignAtlas Download page
  (file "Disease Information of All DiSignAtlas Datasets", CSV).
- Filter parameters below (assay set, minimum paired samples, Psamples).
- CREEDS uses organism == "human" and counts samples from ctrl_ids / pert_ids;
  it has no assay/library-strategy field, so no RNA-Seq/Microarray filter is inferred.

Outputs (written to the results dir resolved from conf, flat inside the slug;
<tag> = "ge<x>Psamples")
------
- disignatlas_signature_ceiling.tsv — dataset/disease counts per filter stage,
  for every threshold in CEILING_THRESHOLDS (untagged).
- disignatlas_inventory_human_<tag>.tsv  — one row per surviving dataset.
- disignatlas_inventory_by_disease_<tag>.tsv — one row per disease, aggregated across
  its surviving datasets (dataset count, distinct-accession count, RNA-Seq count,
  the case/control counts of the largest single dataset, and a name-string neoplasm flag).
- creeds_signature_ceiling.tsv — dataset/disease counts per CREEDS filter stage,
  for every threshold in CEILING_THRESHOLDS (untagged).
- creeds_inventory_human_<tag>.tsv — CREEDS signatures surviving the human / replicate filter.
- creeds_inventory_by_disease_<tag>.tsv — one row per CREEDS disease, aggregated across
  surviving signatures (dataset count, distinct-GEO count, platform count,
  largest-signature case/control counts, and a name-string neoplasm flag).
- creeds_disignatlas_overlap_summary_<tag>.tsv — headline coverage and shared-GSE counts.
- creeds_disignatlas_disease_matches_<tag>.tsv — one row per canonical CREEDS disease,
  classified as UMLS-CUI match, exact-name-only match, or not found.
- creeds_disignatlas_common_gse_<tag>.tsv — one row per shared disease/GSE pair.
- creeds_not_found_in_disignatlas_<tag>.tsv — the CREEDS diseases not found
  conservatively in DiSignAtlas, with the closest name candidate for manual review.
- creeds_disignatlas_merged_summary_by_disease_<tag>.csv — one row per disease and
  source database, with record/accession counts, assay counts, largest study,
  and the complete list of contributing accessions. Built from an in-memory
  row-level union of both inventories, which is not written to disk.

Figures (written to the image dir resolved from conf, same slug)
-------
- <db>_signature_ceiling.png (db = disignatlas, creeds) — records and distinct
  diseases retained after each filter stage, with one bar per threshold in
  CEILING_THRESHOLDS; the bar for x is highlighted (untagged).
- <db>_inventory_panel_<tag>.png — one multi-panel figure per database for the
  filtered inventory: (A) datasets/signatures per disease, (B) top diseases by
  data volume, (C) technical composition (DiSignAtlas: source repository and
  assay; CREEDS: GPL platforms), (D) cross-database disease coverage
  (DiSignAtlas diseases matched by CREEDS; CREEDS diseases found in DiSignAtlas).

Usage
-----
Run from the source directory (so `import conf` resolves), e.g. in Spyder:
    runfile('validations/drug_disease_pairs_validation/01_disignatlas_creeds_data_analysis.py')

Convention
----------
control_case_sample_count is "n_control|n_case" — control first, case second.
"""

import csv as _csv
import difflib
import re
import textwrap
import unicodedata

import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

from conf import DSA_CATALOGUE, OOUT_DIR, CREED_DIR, IMG_DIR

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
ASSAYS = {"RNA-Seq", "Microarray"}   # bulk; drop to {"RNA-Seq"} for sequencing only
MIN_REPLICATES = 10                  # x: minimum samples required in EACH of case and control
CEILING_THRESHOLDS = (3, 10, 20, 50)  # Psamples thresholds shown in the ceiling figure
FILTER_TAG = f"ge{MIN_REPLICATES}Psamples"   # file-name tag of every filtered output
UMLS_CUI_RE = re.compile(r"C\d{7}", flags=re.IGNORECASE)

# ---------------------------------------------------------------------------
# Paths (resolved from conf)
# ---------------------------------------------------------------------------
DSA_CSV = DSA_CATALOGUE
OUT_DIR = OOUT_DIR / "validations" / "drug_disease_pairs_validation"
CEILING_TSV = OUT_DIR / "disignatlas_signature_ceiling.tsv"
INVENTORY_TSV = OUT_DIR / f"disignatlas_inventory_human_{FILTER_TAG}.tsv"
BY_DISEASE_TSV = OUT_DIR / f"disignatlas_inventory_by_disease_{FILTER_TAG}.tsv"

CREED_METADATA = CREED_DIR / "disease_signatures_metadata.csv"
CREED_CEILING_TSV = OUT_DIR / "creeds_signature_ceiling.tsv"
CREED_INVENTORY_TSV = OUT_DIR / f"creeds_inventory_human_{FILTER_TAG}.tsv"
CREED_BY_DISEASE_TSV = OUT_DIR / f"creeds_inventory_by_disease_{FILTER_TAG}.tsv"

OVERLAP_SUMMARY_TSV = OUT_DIR / f"creeds_disignatlas_overlap_summary_{FILTER_TAG}.tsv"
OVERLAP_BY_DISEASE_TSV = OUT_DIR / f"creeds_disignatlas_disease_matches_{FILTER_TAG}.tsv"
OVERLAP_GSE_TSV = OUT_DIR / f"creeds_disignatlas_common_gse_{FILTER_TAG}.tsv"
CREEDS_ONLY_TSV = OUT_DIR / f"creeds_not_found_in_disignatlas_{FILTER_TAG}.tsv"
MERGED_SUMMARY_CSV = OUT_DIR / f"creeds_disignatlas_merged_summary_by_disease_{FILTER_TAG}.csv"

IMG_OUT_DIR = IMG_DIR / "validations" / "drug_disease_pairs_validation"
DSA_STAGE_PNG = IMG_OUT_DIR / "disignatlas_signature_ceiling.png"
DSA_PANEL_PNG = IMG_OUT_DIR / f"disignatlas_inventory_panel_{FILTER_TAG}.png"
CREED_STAGE_PNG = IMG_OUT_DIR / "creeds_signature_ceiling.png"
CREED_PANEL_PNG = IMG_OUT_DIR / f"creeds_inventory_panel_{FILTER_TAG}.png"

# A restrained, colour-blind-friendly palette shared by all figures.  The
# database-specific colour is used for totals; lighter/darker shades encode
# subdivisions without making the charts visually busy.
DSA_COLOUR = "#2F6690"
DSA_LIGHT = "#8ECAE6"
CREEDS_COLOUR = "#7A5195"
CREEDS_LIGHT = "#C6A0D5"
NEUTRAL = "#A7A9AC"
GRID_COLOUR = "#D9DEE3"

# Name-string keywords flagging a likely neoplastic disease. Heuristic on the
# disease NAME only (no ontology walk); a flagged row is a candidate, not a
# verdict. Used to separate the proliferation-artefact-prone cancer datasets.
NEOPLASM_KEYWORDS = (
    "cancer", "carcinoma", "tumor", "tumour", "neoplasm", "neoplastic",
    "leukemia", "leukaemia", "lymphoma", "sarcoma", "melanoma", "myeloma",
    "glioma", "glioblastoma", "adenocarcinoma", "blastoma", "metasta",
    "malignan", "oma of",
)

COLUMNS = [
    "dsaid", "accession", "platform", "deg_count", "disease", "diseaseid",
    "tissue", "data_source", "library_strategy", "organism",
    "control_case_sample_count", "definition",
]


def _split_counts(value):
    """Return (n_control, n_case) from a 'n_control|n_case' string, or (None, None)."""
    try:
        control, case = value.split("|")
        return int(control), int(case)
    except (ValueError, AttributeError):
        return None, None


def _normalise_disease_name(value):
    """Return a conservative comparison key for a disease label.

    Only case, accents, possessive apostrophes, punctuation, and repeated
    whitespace are normalised.  No fuzzy or semantic synonym is accepted as
    an automatic match: those candidates are reported separately for review.
    """
    text = unicodedata.normalize("NFKD", value or "")
    text = text.encode("ascii", "ignore").decode("ascii").lower()
    text = text.replace("&", " and ").replace("'s", "s")
    return " ".join(re.findall(r"[a-z0-9]+", text))


def _representative_name(rows, field):
    """Choose the most frequent non-empty label with deterministic tie-breaking."""
    counts = {}
    for row in rows:
        name = (row.get(field) or "").strip()
        if name:
            counts[name] = counts.get(name, 0) + 1
    if not counts:
        return ""
    return min(counts, key=lambda name: (-counts[name], name.casefold(), name))


def _is_gse(value):
    """Return whether value is a GEO Series accession."""
    return bool(re.fullmatch(r"GSE\d+", (value or "").strip(), flags=re.IGNORECASE))


def _dsa_disease_key(row):
    """Use the DisGeNET/UMLS id, falling back to the normalised disease name.

    The catalogue contains many blank diseaseid values.  Treating the empty
    string as a concept collapses unrelated diseases into one and materially
    undercounts the inventory; the name fallback prevents that loss.
    """
    diseaseid = (row.get("diseaseid") or "").strip()
    if diseaseid:
        return diseaseid
    return f"NAME:{_normalise_disease_name(row.get('disease'))}"


def load_catalogue(path=DSA_CSV):
    """Read the DiSignAtlas catalogue CSV into a list of dict rows.

    The file is not valid UTF-8 (disease definitions carry Latin-1/Windows-1252
    bytes), so the encoding is tried in order and the first that decodes wins.
    """
    for encoding in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            with open(path, newline="", encoding=encoding) as fh:
                rows = list(_csv.DictReader(fh))
            break
        except UnicodeDecodeError:
            continue
    else:
        raise UnicodeDecodeError(
            "utf-8", b"", 0, 1, f"could not decode {path} with utf-8/cp1252/latin-1")
    if rows and set(COLUMNS) - set(rows[0].keys()):
        missing = set(COLUMNS) - set(rows[0].keys())
        raise ValueError(f"catalogue missing expected columns: {sorted(missing)}")
    for r in rows:
        r["n_control"], r["n_case"] = _split_counts(r["control_case_sample_count"])
    return rows


def _distinct_diseases(rows):
    return len({_dsa_disease_key(r) for r in rows})


def _ceiling_thresholds():
    """Thresholds reported in the ceiling: CEILING_THRESHOLDS plus x, ascending."""
    return sorted(set(CEILING_THRESHOLDS) | {MIN_REPLICATES})


def _passes_psamples_filter(row, threshold):
    return (row["n_control"] is not None
            and row["n_control"] >= threshold
            and row["n_case"] >= threshold)


def build_inventory(rows):
    """Apply the human / assay / replicate filters and return the surviving rows.

    The returned inventory uses MIN_REPLICATES. The ceiling additionally counts
    what survives each threshold in CEILING_THRESHOLDS; each row carries the
    threshold (blank for the non-replicate stages) and whether it is the one
    used for the inventory.
    """
    human = [r for r in rows if r["organism"] == "Homo sapiens"]
    assay = [r for r in human if r["library_strategy"] in ASSAYS]
    stages = [
        ("all_datasets", "", rows),
        ("human", "", human),
        (f"human_{'+'.join(sorted(ASSAYS))}", "", assay),
    ]
    for t in _ceiling_thresholds():
        stages.append((f"human_assay_ge{t}_each_psamples", t,
                       [r for r in assay if _passes_psamples_filter(r, t)]))
    ceiling = [
        {"stage": name, "min_psamples": t, "selected": int(t == MIN_REPLICATES),
         "n_datasets": len(sub), "n_distinct_diseases": _distinct_diseases(sub)}
        for name, t, sub in stages
    ]
    adequate = [r for r in assay if _passes_psamples_filter(r, MIN_REPLICATES)]
    return adequate, ceiling


def _is_neoplasm(name):
    low = name.lower()
    return any(kw in low for kw in NEOPLASM_KEYWORDS)


def aggregate_by_disease(inventory):
    """Collapse the surviving datasets to one row per disease (DisGeNET id).

    Per disease we keep counts and, as a size proxy, the case/control counts of the single
    largest dataset (the one maximising the smaller of its case and control groups, i.e. the
    most balanced-and-large single study) — both counts therefore come from the
    SAME dataset, unlike a per-column maximum. n_accessions counts all distinct
    source accessions; n_gse is restricted to valid GEO Series accessions.
    """
    by_id = {}
    for r in inventory:
        disease_key = _dsa_disease_key(r)
        d = by_id.setdefault(disease_key, {
            "diseaseid": r["diseaseid"], "disease": r["disease"],
            "n_datasets": 0, "n_rnaseq": 0, "n_microarray": 0,
            "accessions": set(), "gse_accessions": set(),
            "best_control": 0, "best_case": 0, "_best_min": -1,
        })
        d["n_datasets"] += 1
        d["n_rnaseq"] += r["library_strategy"] == "RNA-Seq"
        d["n_microarray"] += r["library_strategy"] == "Microarray"
        d["accessions"].add(r["accession"])
        if _is_gse(r["accession"]):
            d["gse_accessions"].add(r["accession"].strip().upper())
        arm_min = min(r["n_control"], r["n_case"])
        if arm_min > d["_best_min"]:
            d["_best_min"] = arm_min
            d["best_control"], d["best_case"] = r["n_control"], r["n_case"]
    rows = []
    for d in by_id.values():
        d["n_accessions"] = len(d.pop("accessions"))
        d["n_gse"] = len(d.pop("gse_accessions"))
        d.pop("_best_min")
        d["neoplasm_kw"] = int(_is_neoplasm(d["disease"]))
        rows.append(d)
    return sorted(rows, key=lambda x: (-x["n_datasets"], x["disease"]))


def inspect_inventory(inventory, by_disease):
    """Print a descriptive summary of the filtered inventory and return it."""
    n_dis = len(by_disease)
    # datasets-per-disease buckets
    buckets = {"1": 0, "2-4": 0, "5-9": 0, "10+": 0}
    for d in by_disease:
        n = d["n_datasets"]
        key = "1" if n == 1 else "2-4" if n <= 4 else "5-9" if n <= 9 else "10+"
        buckets[key] += 1
    # breakdowns over datasets
    by_source, by_assay = {}, {}
    for r in inventory:
        by_source[r["data_source"]] = by_source.get(r["data_source"], 0) + 1
        by_assay[r["library_strategy"]] = by_assay.get(r["library_strategy"], 0) + 1
    neo_dis = sum(d["neoplasm_kw"] for d in by_disease)
    neo_ds = sum(1 for r in inventory if _is_neoplasm(r["disease"]))

    print(f"\n--- inventory inspection (human, assay, >= {MIN_REPLICATES} case and control) ---")
    print(f"datasets: {len(inventory)}   distinct diseases: {n_dis}")
    print("datasets per disease:  " +
          "  ".join(f"{k}:{v}" for k, v in buckets.items()))
    print("by source:  " + "  ".join(f"{k}:{v}" for k, v in sorted(by_source.items())))
    print("by assay:   " + "  ".join(f"{k}:{v}" for k, v in sorted(by_assay.items())))
    print(f"neoplasm (name-kw): {neo_dis}/{n_dis} diseases, {neo_ds}/{len(inventory)} datasets"
          f"  -> non-neoplastic diseases: {n_dis - neo_dis}")
    print("top 15 diseases by dataset count")
    print(f"  {'#ds':>4} {'#acc':>4} {'#RNAseq':>7}  {'biggest ctrl|case':>17}  disease")
    print(f"  {'':->4} {'':->4} {'':->7}  {'(same dataset)':>17}  {'':->30}")
    for d in by_disease[:15]:
        tag = " [neo]" if d["neoplasm_kw"] else ""
        biggest = f"{d['best_control']}|{d['best_case']}"
        print(f"  {d['n_datasets']:>4} {d['n_accessions']:>4} {d['n_rnaseq']:>7}  {biggest:>17}"
              f"  {d['disease'][:40]}{tag}")
    return buckets


def write_outputs(inventory, ceiling, by_disease=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = ["dsaid", "accession", "data_source", "disease", "diseaseid",
              "tissue", "library_strategy", "platform", "n_control", "n_case",
              "deg_count"]
    with open(INVENTORY_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(inventory, key=lambda x: (x["disease"], x["dsaid"])):
            w.writerow(r)
    with open(CEILING_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=["stage", "min_psamples", "selected", "n_datasets",
                                        "n_distinct_diseases"],
                            delimiter="\t")
        w.writeheader()
        w.writerows(ceiling)
    if by_disease is not None:
        fields = ["diseaseid", "disease", "n_datasets", "n_accessions", "n_gse", "n_rnaseq",
                  "n_microarray", "best_control", "best_case", "neoplasm_kw"]
        with open(BY_DISEASE_TSV, "w", newline="", encoding="utf-8") as fh:
            w = _csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(by_disease)


# ---------------------------------------------------------------------------
# CREEDS
# ---------------------------------------------------------------------------
def _count_gsm_ids(value):
    """Count distinct pipe-separated sample accessions in a CREEDS case or control group."""
    if value is None:
        return 0
    value = str(value).strip()
    if not value or value.lower() == "nan":
        return 0
    return len({x.strip().upper() for x in value.split("|") if x.strip()})


def load_creeds_metadata(path=CREED_METADATA):
    """Read CREEDS disease-signature metadata.

    CREEDS stores the control and disease-group sample accessions directly in
    ctrl_ids and pert_ids, so replicate counts are derived from those fields.
    There is no assay/library-strategy column analogous to DiSignAtlas.
    """
    encodings = ("utf-8-sig", "cp1252", "latin-1")
    rows = None
    for encoding in encodings:
        try:
            with open(path, newline="", encoding=encoding) as fh:
                rows = list(_csv.DictReader(fh))
            break
        except UnicodeDecodeError:
            continue
    if rows is None:
        raise UnicodeDecodeError(
            "utf-8", b"", 0, 1, f"could not decode {path} with utf-8/cp1252/latin-1")

    expected = {
        "id", "cell_type", "ctrl_ids", "disease_name", "do_id", "geo_id",
        "organism", "pert_ids", "platform"
    }
    if rows and expected - set(rows[0].keys()):
        raise ValueError(
            f"CREEDS metadata missing expected columns: {sorted(expected - set(rows[0].keys()))}"
        )

    seen_ids = set()
    for r in rows:
        sid = (r.get("id") or "").strip()
        if not sid:
            raise ValueError("CREEDS metadata contains a row with an empty signature id")
        if sid in seen_ids:
            raise ValueError(f"duplicate CREEDS signature id: {sid}")
        seen_ids.add(sid)
        r["n_control"] = _count_gsm_ids(r.get("ctrl_ids"))
        r["n_case"] = _count_gsm_ids(r.get("pert_ids"))
        # CREEDS has no assay/library-strategy field; keep this explicit
        # rather than guessing it from GPL/platform identifiers.
        r["data_source"] = "GEO" if r.get("geo_id") else ""
    return rows


def _creeds_disease_key(row):
    """Use DOID when present; otherwise fall back to the disease name."""
    doid = (row.get("do_id") or "").strip()
    return doid if doid else f"NAME:{_normalise_disease_name(row['disease_name'])}"


def _group_creeds_by_disease(rows):
    """Group CREEDS rows into canonical disease concepts.

    A few labels occur once with a DOID and once without it.  When a normalised
    name maps to exactly one observed DOID, the id-less rows are folded into
    that DOID instead of being counted as a second disease.
    """
    name_to_doids = {}
    for row in rows:
        name = _normalise_disease_name(row.get("disease_name"))
        doid = (row.get("do_id") or "").strip()
        if doid:
            name_to_doids.setdefault(name, set()).add(doid)

    groups = {}
    for row in rows:
        name = _normalise_disease_name(row.get("disease_name"))
        observed_doids = name_to_doids.get(name, set())
        if len(observed_doids) == 1:
            key = next(iter(observed_doids))
        else:
            key = _creeds_disease_key(row)
        groups.setdefault(key, []).append(row)
    return groups


def build_creeds_inventory(rows):
    """Apply the human and minimum-Psamples filters to CREEDS."""
    human = [r for r in rows if (r.get("organism") or "").strip().lower() == "human"]
    stages = [("all_datasets", "", rows), ("human", "", human)]
    for t in _ceiling_thresholds():
        stages.append((f"human_ge{t}_each_psamples", t,
                       [r for r in human if _passes_psamples_filter(r, t)]))
    ceiling = [
        {
            "stage": name,
            "min_psamples": t,
            "selected": int(t == MIN_REPLICATES),
            "n_datasets": len(sub),
            "n_distinct_diseases": len(_group_creeds_by_disease(sub)),
        }
        for name, t, sub in stages
    ]
    adequate = [r for r in human if _passes_psamples_filter(r, MIN_REPLICATES)]
    return adequate, ceiling


def aggregate_creeds_by_disease(inventory):
    """Collapse CREEDS signatures to one row per disease concept.

    The disease key is DOID where available, otherwise the disease name.
    n_gse counts distinct GEO accessions.  best_control/best_case are taken
    from the SAME signature, choosing the signature with the largest smaller
    group, matching the DiSignAtlas aggregation logic.
    """
    by_id = {}
    for key, disease_rows in _group_creeds_by_disease(inventory).items():
        # Prefer the most frequent spelling and a non-empty DOID for display.
        disease_name = _representative_name(disease_rows, "disease_name")
        diseaseid = next(
            ((r.get("do_id") or "").strip() for r in disease_rows
             if (r.get("do_id") or "").strip()),
            "",
        )
        d = by_id.setdefault(key, {
            "diseaseid": diseaseid,
            "disease": disease_name,
            "n_datasets": 0,
            "accessions": set(),
            "platforms": set(),
            "best_control": 0,
            "best_case": 0,
            "_best_min": -1,
        })
        for r in disease_rows:
            d["n_datasets"] += 1
            if r.get("geo_id"):
                d["accessions"].add(r["geo_id"])
            if r.get("platform"):
                d["platforms"].add(r["platform"])

            arm_min = min(r["n_control"], r["n_case"])
            if arm_min > d["_best_min"]:
                d["_best_min"] = arm_min
                d["best_control"], d["best_case"] = r["n_control"], r["n_case"]

    rows = []
    for d in by_id.values():
        d["n_gse"] = len(d.pop("accessions"))
        d["n_platforms"] = len(d.pop("platforms"))
        d.pop("_best_min")
        d["neoplasm_kw"] = int(_is_neoplasm(d["disease"]))
        rows.append(d)
    return sorted(rows, key=lambda x: (-x["n_datasets"], x["disease"].lower()))


def inspect_creeds_inventory(inventory, by_disease):
    """Print the CREEDS analogue of the DiSignAtlas inventory inspection."""
    n_dis = len(by_disease)

    buckets = {"1": 0, "2-4": 0, "5-9": 0, "10+": 0}
    for d in by_disease:
        n = d["n_datasets"]
        key = "1" if n == 1 else "2-4" if n <= 4 else "5-9" if n <= 9 else "10+"
        buckets[key] += 1

    by_source, by_platform = {}, {}
    for r in inventory:
        by_source[r["data_source"]] = by_source.get(r["data_source"], 0) + 1
        by_platform[r["platform"]] = by_platform.get(r["platform"], 0) + 1

    neo_dis = sum(d["neoplasm_kw"] for d in by_disease)
    neo_ds = sum(1 for r in inventory if _is_neoplasm(r["disease_name"]))

    print(f"\n--- CREEDS inventory inspection (human, >= {MIN_REPLICATES} case and control) ---")
    print(f"datasets: {len(inventory)}   distinct diseases: {n_dis}")
    print("datasets per disease:  " +
          "  ".join(f"{k}:{v}" for k, v in buckets.items()))
    print("by source:  " +
          "  ".join(f"{k}:{v}" for k, v in sorted(by_source.items()) if k))
    print("by platform: " +
          "  ".join(f"{k}:{v}" for k, v in sorted(by_platform.items())))
    print(f"neoplasm (name-kw): {neo_dis}/{n_dis} diseases, {neo_ds}/{len(inventory)} datasets"
          f"  -> non-neoplastic diseases: {n_dis - neo_dis}")
    print("top 15 diseases by dataset count")
    print(f"  {'#ds':>4} {'#GSE':>4} {'#platforms':>10}  {'biggest ctrl|case':>17}  disease")
    print(f"  {'':->4} {'':->4} {'':->10}  {'(same dataset)':>17}  {'':->30}")
    for d in by_disease[:15]:
        tag = " [neo]" if d["neoplasm_kw"] else ""
        biggest = f"{d['best_control']}|{d['best_case']}"
        print(f"  {d['n_datasets']:>4} {d['n_gse']:>4} {d['n_platforms']:>10}  {biggest:>17}"
              f"  {d['disease'][:40]}{tag}")

    return buckets


def write_creeds_outputs(inventory, ceiling, by_disease):
    """Write the three CREEDS inventory tables."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fields = [
        "id", "geo_id", "data_source", "disease_name", "do_id", "cell_type",
        "organism", "platform", "n_control", "n_case",
    ]
    with open(CREED_INVENTORY_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in sorted(inventory, key=lambda x: (x["disease_name"].lower(), x["id"])):
            w.writerow(r)

    with open(CREED_CEILING_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(
            fh,
            fieldnames=["stage", "min_psamples", "selected", "n_datasets",
                        "n_distinct_diseases"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(ceiling)

    fields = [
        "diseaseid", "disease", "n_datasets", "n_gse", "n_platforms",
        "best_control", "best_case", "neoplasm_kw",
    ]
    with open(CREED_BY_DISEASE_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        w.writerows(by_disease)


# ---------------------------------------------------------------------------
# CREEDS versus DiSignAtlas coverage
# ---------------------------------------------------------------------------
def _gsm_ids(row):
    """Return the distinct CREEDS sample accessions used by one signature."""
    ids = set()
    for field in ("ctrl_ids", "pert_ids"):
        ids.update(
            value.strip().upper()
            for value in (row.get(field) or "").split("|")
            if value.strip()
        )
    return ids


def build_creeds_disignatlas_overlap(dsa_inventory, creeds_inventory):
    """Compare the usable human disease inventories conservatively.

    Matching is hierarchical and auditable:
      1. exact UMLS CUI (CREEDS umls_cui vs DiSignAtlas diseaseid);
      2. if no CUI match exists, exact normalised disease label;
      3. otherwise ``not_found``.  Fuzzy similarity is reported only as a
         review candidate and never promoted to an automatic match.

    A GSE is shared only when it occurs in both databases for the same matched
    CREEDS disease.  True GSM overlap cannot be computed because the supplied
    DiSignAtlas catalogue contains case/control sample counts but not GSM accessions.
    """
    dsa_groups = {}
    for row in dsa_inventory:
        dsa_groups.setdefault(_dsa_disease_key(row), []).append(row)

    dsa_by_cui, dsa_by_name = {}, {}
    dsa_display_by_name = {}
    for key, rows in dsa_groups.items():
        for row in rows:
            diseaseid = (row.get("diseaseid") or "").strip()
            name_key = _normalise_disease_name(row.get("disease"))
            if UMLS_CUI_RE.fullmatch(diseaseid):
                dsa_by_cui.setdefault(diseaseid.upper(), set()).add(key)
            dsa_by_name.setdefault(name_key, set()).add(key)
            dsa_display_by_name.setdefault(name_key, row.get("disease") or "")

    creeds_groups = _group_creeds_by_disease(creeds_inventory)
    match_rows, shared_gse_rows = [], []
    all_shared_gse, all_candidate_gsm = set(), set()

    for creeds_key, rows in creeds_groups.items():
        disease_name = _representative_name(rows, "disease_name")
        doids = sorted({(r.get("do_id") or "").strip() for r in rows
                        if (r.get("do_id") or "").strip()})
        cuis = sorted({(r.get("umls_cui") or "").strip().upper() for r in rows
                       if UMLS_CUI_RE.fullmatch((r.get("umls_cui") or "").strip())})
        name_keys = sorted({_normalise_disease_name(r.get("disease_name")) for r in rows})

        cui_matches = set()
        for cui in cuis:
            cui_matches.update(dsa_by_cui.get(cui, set()))
        name_matches = set()
        for name_key in name_keys:
            name_matches.update(dsa_by_name.get(name_key, set()))
        # CUI is the stronger identifier.  Name matching is a fallback only;
        # unioning discordant name matches here would inflate downstream GSEs.
        matched_dsa_keys = cui_matches if cui_matches else name_matches

        if cui_matches and name_matches:
            match_status = "cui_and_exact_name"
        elif cui_matches:
            match_status = "cui_only"
        elif name_matches:
            match_status = "exact_name_only"
        else:
            match_status = "not_found"

        matched_dsa_rows = [
            row for key in matched_dsa_keys for row in dsa_groups[key]
        ]
        dsa_gse = {
            row["accession"].strip().upper()
            for row in matched_dsa_rows if _is_gse(row.get("accession"))
        }
        creeds_gse = {
            row["geo_id"].strip().upper()
            for row in rows if _is_gse(row.get("geo_id"))
        }
        shared_gse = sorted(dsa_gse & creeds_gse)
        all_shared_gse.update(shared_gse)

        candidate_gsm = set()
        for row in rows:
            if (row.get("geo_id") or "").strip().upper() in shared_gse:
                candidate_gsm.update(_gsm_ids(row))
        all_candidate_gsm.update(candidate_gsm)

        closest_name, closest_score = "", ""
        if match_status == "not_found" and dsa_display_by_name:
            query_name = _normalise_disease_name(disease_name)
            closest_key = min(
                dsa_display_by_name,
                key=lambda key: (
                    -difflib.SequenceMatcher(None, query_name, key).ratio(), key
                ),
            )
            closest_name = dsa_display_by_name[closest_key]
            closest_score = round(
                difflib.SequenceMatcher(None, query_name, closest_key).ratio(), 3
            )

        matched_names = sorted({
            row.get("disease") or "" for row in matched_dsa_rows
            if row.get("disease")
        })
        match_rows.append({
            "creeds_key": creeds_key,
            "creeds_disease": disease_name,
            "creeds_do_ids": "|".join(doids),
            "creeds_umls_cuis": "|".join(cuis),
            "match_status": match_status,
            "matched_dsa_keys": "|".join(sorted(matched_dsa_keys)),
            "matched_dsa_diseases": "|".join(matched_names),
            "n_creeds_signatures": len(rows),
            "n_dsa_datasets": len(matched_dsa_rows),
            "n_creeds_gse": len(creeds_gse),
            "n_dsa_gse": len(dsa_gse),
            "n_shared_gse": len(shared_gse),
            "shared_gse": "|".join(shared_gse),
            "n_creeds_candidate_gsm_in_shared_gse": len(candidate_gsm),
            "closest_dsa_name_if_unmatched": closest_name,
            "closest_name_similarity": closest_score,
        })

        for gse in shared_gse:
            creeds_rows_in_gse = [
                row for row in rows
                if (row.get("geo_id") or "").strip().upper() == gse
            ]
            dsa_rows_in_gse = [
                row for row in matched_dsa_rows
                if (row.get("accession") or "").strip().upper() == gse
            ]
            gse_candidate_gsm = set().union(
                *(_gsm_ids(row) for row in creeds_rows_in_gse)
            ) if creeds_rows_in_gse else set()
            shared_gse_rows.append({
                "creeds_key": creeds_key,
                "creeds_disease": disease_name,
                "gse": gse,
                "n_creeds_signatures": len(creeds_rows_in_gse),
                "n_dsa_datasets": len(dsa_rows_in_gse),
                "n_creeds_candidate_gsm": len(gse_candidate_gsm),
                "shared_gsm": "NA",
                "shared_gsm_status": "not_computable_without_disignatlas_gsm_ids",
            })

    n_cui = sum(row["match_status"] in {"cui_only", "cui_and_exact_name"}
                for row in match_rows)
    n_not_found = sum(row["match_status"] == "not_found" for row in match_rows)
    n_matched = len(match_rows) - n_not_found
    n_with_shared_gse = sum(row["n_shared_gse"] > 0 for row in match_rows)

    summary = [
        {"metric": "disignatlas_diseases", "value": len(dsa_groups),
         "definition": "UMLS/DisGeNET id; normalised-name fallback when id is blank"},
        {"metric": "creeds_diseases", "value": len(creeds_groups),
         "definition": "DOID; id-less duplicate labels folded into their unique DOID"},
        {"metric": "common_diseases_strict_cui", "value": n_cui,
         "definition": "at least one exact shared UMLS CUI"},
        {"metric": "common_diseases_cui_or_exact_name", "value": n_matched,
         "definition": "exact UMLS CUI or exact normalised disease label"},
        {"metric": "creeds_diseases_not_found", "value": n_not_found,
         "definition": "no exact CUI or exact normalised-name match; requires manual review"},
        {"metric": "creeds_fraction_found",
         "value": round(n_matched / len(match_rows), 6) if match_rows else 0.0,
         "definition": "common_diseases_cui_or_exact_name / creeds_diseases"},
        {"metric": "unique_shared_gse", "value": len(all_shared_gse),
         "definition": "distinct GSE accessions shared within a matched disease"},
        {"metric": "shared_disease_gse_pairs", "value": len(shared_gse_rows),
         "definition": "one row per matched CREEDS disease and shared GSE"},
        {"metric": "common_diseases_with_shared_gse", "value": n_with_shared_gse,
         "definition": "matched diseases with at least one same GSE"},
        {"metric": "creeds_candidate_gsm_in_shared_gse", "value": len(all_candidate_gsm),
         "definition": "CREEDS GSMs in shared disease/GSE pairs; not a GSM intersection"},
        {"metric": "shared_gsm", "value": "NA",
         "definition": "not computable: DiSignAtlas catalogue has counts but no GSM ids"},
    ]

    match_rows.sort(key=lambda row: (
        row["match_status"] != "not_found",
        -row["n_creeds_signatures"],
        row["creeds_disease"].lower(),
    ))
    shared_gse_rows.sort(key=lambda row: (row["creeds_disease"].lower(), row["gse"]))
    return match_rows, shared_gse_rows, summary


def write_creeds_disignatlas_overlap(match_rows, shared_gse_rows, summary):
    """Write the comparison summary, detailed matches, shared GSEs, and misses."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with open(OVERLAP_SUMMARY_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=["metric", "value", "definition"],
                            delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    match_fields = [
        "creeds_key", "creeds_disease", "creeds_do_ids", "creeds_umls_cuis",
        "match_status", "matched_dsa_keys", "matched_dsa_diseases",
        "n_creeds_signatures", "n_dsa_datasets", "n_creeds_gse", "n_dsa_gse",
        "n_shared_gse", "shared_gse", "n_creeds_candidate_gsm_in_shared_gse",
        "closest_dsa_name_if_unmatched", "closest_name_similarity",
    ]
    with open(OVERLAP_BY_DISEASE_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
        w.writeheader()
        w.writerows(match_rows)

    gse_fields = [
        "creeds_key", "creeds_disease", "gse", "n_creeds_signatures",
        "n_dsa_datasets", "n_creeds_candidate_gsm", "shared_gsm",
        "shared_gsm_status",
    ]
    with open(OVERLAP_GSE_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=gse_fields, delimiter="\t")
        w.writeheader()
        w.writerows(shared_gse_rows)

    creeds_only = [row for row in match_rows if row["match_status"] == "not_found"]
    with open(CREEDS_ONLY_TSV, "w", newline="", encoding="utf-8") as fh:
        w = _csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
        w.writeheader()
        w.writerows(creeds_only)


def build_merged_inventory(dsa_inventory, creeds_inventory, match_rows):
    """Return one harmonised row for every usable record in either database.

    The union is held in memory only: it feeds the disease-level summary and
    is not written to disk.

    The function intentionally retains records documented by both sources.
    ``shared_dataset_across_databases`` identifies same-disease GSEs present in
    both inventories so a downstream analysis can deduplicate them explicitly.

    A CREEDS disease is assigned to a DiSignAtlas disease key only when the
    conservative overlap analysis yields exactly one key.  Multiple candidates
    are preserved as ``AMBIGUOUS_DSA:...`` rather than resolved arbitrarily;
    unmatched CREEDS diseases retain a namespaced CREEDS key.
    """
    dsa_groups = {}
    for row in dsa_inventory:
        dsa_groups.setdefault(_dsa_disease_key(row), []).append(row)
    dsa_names = {
        key: _representative_name(rows, "disease")
        for key, rows in dsa_groups.items()
    }

    match_by_creeds_key = {row["creeds_key"]: row for row in match_rows}
    creeds_key_by_record_id = {
        row["id"]: key
        for key, rows in _group_creeds_by_disease(creeds_inventory).items()
        for row in rows
    }

    shared_dsa_gse = set()
    for match in match_rows:
        dsa_keys = [key for key in match["matched_dsa_keys"].split("|") if key]
        shared_gse = [gse for gse in match["shared_gse"].split("|") if gse]
        shared_dsa_gse.update((key, gse) for key in dsa_keys for gse in shared_gse)

    merged = []
    for row in dsa_inventory:
        disease_key = _dsa_disease_key(row)
        accession = (row.get("accession") or "").strip()
        merged.append({
            "source_database": "DiSignAtlas",
            "source_record_id": row.get("dsaid") or "",
            "source_disease_key": disease_key,
            "dataset_accession": accession,
            "merged_disease_key": f"DSA:{disease_key}",
            "merged_disease_name": dsa_names[disease_key],
            "merged_disease_name_normalized": _normalise_disease_name(
                dsa_names[disease_key]
            ),
            "source_disease_name": row.get("disease") or "",
            "source_disease_id": row.get("diseaseid") or "",
            "umls_cui": ((row.get("diseaseid") or "").strip().upper()
                         if UMLS_CUI_RE.fullmatch((row.get("diseaseid") or "").strip())
                         else ""),
            "doid": "",
            "disease_match_status": "native_disignatlas",
            "matched_dsa_keys": disease_key,
            "matched_dsa_diseases": dsa_names[disease_key],
            "shared_dataset_across_databases": int(
                (disease_key, accession.upper()) in shared_dsa_gse
            ),
            "data_source": row.get("data_source") or "",
            "assay": row.get("library_strategy") or "",
            "platform": row.get("platform") or "",
            "tissue_or_cell_type": row.get("tissue") or "",
            "organism": row.get("organism") or "",
            "n_control": row["n_control"],
            "n_case": row["n_case"],
            "control_sample_ids": "",
            "case_sample_ids": "",
            "deg_count": row.get("deg_count") or "",
            "definition": row.get("definition") or "",
        })

    for row in creeds_inventory:
        creeds_key = creeds_key_by_record_id[row["id"]]
        match = match_by_creeds_key[creeds_key]
        dsa_keys = [key for key in match["matched_dsa_keys"].split("|") if key]
        if len(dsa_keys) == 1:
            merged_key = f"DSA:{dsa_keys[0]}"
            merged_name = dsa_names[dsa_keys[0]]
        elif len(dsa_keys) > 1:
            merged_key = "AMBIGUOUS_DSA:" + "|".join(dsa_keys)
            merged_name = match["creeds_disease"]
        else:
            merged_key = f"CREEDS:{creeds_key}"
            merged_name = match["creeds_disease"]

        accession = (row.get("geo_id") or "").strip()
        shared_gse = {gse for gse in match["shared_gse"].split("|") if gse}
        merged.append({
            "source_database": "CREEDS",
            "source_record_id": row.get("id") or "",
            "source_disease_key": creeds_key,
            "dataset_accession": accession,
            "merged_disease_key": merged_key,
            "merged_disease_name": merged_name,
            "merged_disease_name_normalized": _normalise_disease_name(merged_name),
            "source_disease_name": row.get("disease_name") or "",
            "source_disease_id": row.get("do_id") or "",
            "umls_cui": row.get("umls_cui") or "",
            "doid": row.get("do_id") or "",
            "disease_match_status": match["match_status"],
            "matched_dsa_keys": match["matched_dsa_keys"],
            "matched_dsa_diseases": match["matched_dsa_diseases"],
            "shared_dataset_across_databases": int(accession.upper() in shared_gse),
            "data_source": row.get("data_source") or "",
            "assay": "",
            "platform": row.get("platform") or "",
            "tissue_or_cell_type": row.get("cell_type") or "",
            "organism": row.get("organism") or "",
            "n_control": row["n_control"],
            "n_case": row["n_case"],
            "control_sample_ids": row.get("ctrl_ids") or "",
            "case_sample_ids": row.get("pert_ids") or "",
            "deg_count": "",
            "definition": "",
        })

    return sorted(merged, key=lambda row: (
        _normalise_disease_name(row["merged_disease_name"]),
        row["source_database"],
        row["dataset_accession"],
        row["source_record_id"],
    ))


def build_merged_disease_summary(merged_rows):
    """Aggregate the merged inventory by disease and source database."""
    groups = {}
    for row in merged_rows:
        key = (row["source_database"], row["source_disease_key"])
        groups.setdefault(key, []).append(row)

    summary_rows = []
    for (source_database, source_disease_key), rows in groups.items():
        accessions = sorted({
            row["dataset_accession"].strip()
            for row in rows if row["dataset_accession"].strip()
        })
        gse_accessions = sorted({
            accession.upper() for accession in accessions if _is_gse(accession)
        })
        best = max(rows, key=lambda row: (
            min(int(row["n_control"]), int(row["n_case"])),
            int(row["n_control"]) + int(row["n_case"]),
            int(row["n_control"]),
            int(row["n_case"]),
            row["source_record_id"],
        ))
        is_dsa = source_database == "DiSignAtlas"
        summary_rows.append({
            "source_database": source_database,
            "source_disease_key": source_disease_key,
            "merged_disease_key": _representative_name(rows, "merged_disease_key"),
            "disease": _representative_name(rows, "merged_disease_name"),
            "disease_name_normalized": _representative_name(
                rows, "merged_disease_name_normalized"
            ),
            "source_disease_names": "|".join(sorted({
                row["source_disease_name"] for row in rows
                if row["source_disease_name"]
            })),
            "source_disease_ids": "|".join(sorted({
                row["source_disease_id"] for row in rows
                if row["source_disease_id"]
            })),
            "umls_cuis": "|".join(sorted({
                row["umls_cui"] for row in rows if row["umls_cui"]
            })),
            "doids": "|".join(sorted({
                row["doid"] for row in rows if row["doid"]
            })),
            "n_datasets": len(rows),
            "n_accessions": len(accessions),
            "n_gse": len(gse_accessions),
            # CREEDS has no assay field: blanks mean unknown, not zero.
            "n_rnaseq": sum(row["assay"] == "RNA-Seq" for row in rows)
            if is_dsa else "",
            "n_microarray": sum(row["assay"] == "Microarray" for row in rows)
            if is_dsa else "",
            "n_platforms": len({row["platform"] for row in rows if row["platform"]}),
            "biggest_ctrl_case": f"{best['n_control']}|{best['n_case']}",
            "biggest_n_control": best["n_control"],
            "biggest_n_case": best["n_case"],
            "dataset_accessions": "|".join(accessions),
            "gse_accessions": "|".join(gse_accessions),
            "n_records_shared_across_databases": sum(
                int(row["shared_dataset_across_databases"]) for row in rows
            ),
        })

    return sorted(summary_rows, key=lambda row: (
        _normalise_disease_name(row["disease"]),
        row["source_database"],
        row["source_disease_key"],
    ))


def write_merged_disease_summary(rows, path=MERGED_SUMMARY_CSV):
    """Write the disease-level summary as a UTF-8 CSV."""
    fields = [
        "source_database", "source_disease_key", "merged_disease_key", "disease",
        "disease_name_normalized", "source_disease_names", "source_disease_ids",
        "umls_cuis", "doids", "n_datasets", "n_accessions", "n_gse",
        "n_rnaseq", "n_microarray", "n_platforms", "biggest_ctrl_case",
        "biggest_n_control", "biggest_n_case", "dataset_accessions",
        "gse_accessions", "n_records_shared_across_databases",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = _csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Plots (shared visual language for DiSignAtlas and CREEDS)
# ---------------------------------------------------------------------------
def _savefig(fig, path):
    """Save one reproducible figure and close only that figure."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _style_axis(ax, grid_axis="y"):
    """Apply a quiet publication-style baseline without external packages."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID_COLOUR)
    ax.spines["bottom"].set_color(GRID_COLOUR)
    ax.grid(axis=grid_axis, color=GRID_COLOUR, linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    if grid_axis == "y":
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    else:
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))


def _annotate_vertical_bars(ax, bars, values, denominator=None, fontsize=9):
    """Label vertical bars with counts and, when requested, percentages."""
    for bar, value in zip(bars, values):
        label = f"{value:,}"
        if denominator:
            label += f"\n({value / denominator:.1%})"
        ax.annotate(
            label,
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            linespacing=1.15,
        )
    ax.margins(y=0.22)


def _short_disease_name(name, width=46):
    """Keep long disease names readable on horizontal-bar charts."""
    return textwrap.shorten(" ".join(name.split()), width=width, placeholder="…")


def _panel_letter(ax, letter):
    ax.text(-0.10, 1.08, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom", ha="left")


def _filter_phrase():
    return f"filtered inventory: ≥{MIN_REPLICATES} paired case/control samples (Psamples)"


# --- signature ceiling (one standalone figure per database, untagged) ------
def _plot_signature_ceiling(ceiling, base_labels, database, colour, light, path):
    """Records and distinct diseases retained after each filter stage.

    The non-replicate stages are drawn in the database colour; one bar per
    Psamples threshold follows in the light shade, and the threshold used for
    the filtered inventory (MIN_REPLICATES) is drawn in the full colour with
    a dark outline.
    """
    n_base = len(base_labels)
    labels = list(base_labels) + [f"≥{c['min_psamples']} Psamples" for c in ceiling[n_base:]]
    colours, edges = [], []
    for i, c in enumerate(ceiling):
        if i < n_base:
            colours.append(colour); edges.append(colour)
        elif c["selected"]:
            colours.append(colour); edges.append("#222222")
        else:
            colours.append(light); edges.append(colour)

    metrics = (
        ("n_datasets", "Datasets / signatures"),
        ("n_distinct_diseases", "Distinct diseases"),
    )
    x = list(range(len(labels)))
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.6))
    for ax, (field, panel_title) in zip(axes, metrics):
        values = [c[field] for c in ceiling]
        bars = ax.bar(x, values, color=colours, edgecolor=edges, linewidth=1.2,
                      width=0.70)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8.5)
        for tick, c in zip(ax.get_xticklabels(), ceiling):
            if c["selected"]:
                tick.set_fontweight("bold")
        if len(labels) > n_base:
            ax.axvline(n_base - 0.5, color=GRID_COLOUR, linewidth=1.0, linestyle="--")
        ax.set_title(panel_title, loc="left", fontweight="bold")
        ax.set_ylabel("Count")
        _style_axis(ax, "y")
        _annotate_vertical_bars(ax, bars, values, denominator=values[0], fontsize=8)

    thresholds = ", ".join(f"≥{c['min_psamples']}" for c in ceiling[n_base:])
    fig.suptitle(
        f"{database}: metadata retained after each filter "
        f"(Psamples thresholds {thresholds})",
        x=0.06, ha="left", fontsize=14, fontweight="bold", y=1.02,
    )
    fig.text(
        0.06, 0.965,
        "Right of the dashed line: records passing each Psamples threshold "
        "(paired case/control samples), each applied to the preceding stage independently. "
        f"Outlined bar = ≥{MIN_REPLICATES} Psamples, the threshold used for the "
        "filtered inventory. Percentages are relative to all metadata records.",
        ha="left", va="top", fontsize=9, color="#555555", wrap=True,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    _savefig(fig, path)


def plot_disignatlas_signature_ceiling(ceiling):
    labels = ["All\nrecords", "Human", "Human +\nbulk assay"]
    _plot_signature_ceiling(ceiling, labels, "DiSignAtlas", DSA_COLOUR,
                            DSA_LIGHT, DSA_STAGE_PNG)


def plot_creeds_signature_ceiling(ceiling):
    labels = ["All\nsignatures", "Human"]
    _plot_signature_ceiling(ceiling, labels, "CREEDS", CREEDS_COLOUR,
                            CREEDS_LIGHT, CREED_STAGE_PNG)


# --- panel building blocks (each draws on a given axis) --------------------
def _ax_datasets_per_disease(ax, by_disease, database, item_name, colour):
    labels = ["1", "2–4", "5–9", "10+"]
    values = [0, 0, 0, 0]
    for disease in by_disease:
        n = disease["n_datasets"]
        idx = 0 if n == 1 else 1 if n <= 4 else 2 if n <= 9 else 3
        values[idx] += 1
    total = sum(values)
    bars = ax.bar(labels, values, color=colour, width=0.68)
    ax.set_xlabel(f"Usable {item_name.lower()} per disease")
    ax.set_ylabel("Number of diseases")
    ax.set_title("Depth of evidence per disease", loc="left",
                 fontsize=12, fontweight="bold", pad=10)
    _style_axis(ax, "y")
    _annotate_vertical_bars(ax, bars, values, denominator=total)
    ax.text(1.0, 1.02, f"n = {total:,} diseases", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9, color="#555555")


def _ax_top_diseases(ax, by_disease, item_name, colour, light=None,
                     split=None, top_n=15):
    """Top diseases by volume; with split=(field_a, label_a, field_b, label_b)
    each bar is stacked into two components."""
    top = list(reversed(by_disease[:top_n]))
    labels = [_short_disease_name(d["disease"]) for d in top]
    totals = [d["n_datasets"] for d in top]
    y = list(range(len(top)))
    if split:
        fa, la, fb, lb = split
        a = [d[fa] for d in top]
        ax.barh(y, a, color=light, label=la)
        ax.barh(y, [d[fb] for d in top], left=a, color=colour, label=lb)
        ax.legend(frameon=False, ncol=2, loc="lower right", fontsize=9)
    else:
        ax.barh(y, totals, color=colour)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(f"Usable {item_name}")
    ax.set_title(f"Top {len(top)} diseases by {item_name[:-1]} count", loc="left",
                 fontsize=12, fontweight="bold", pad=10)
    _style_axis(ax, "x")
    for yi, total in zip(y, totals):
        ax.annotate(f"{total:,}", (total, yi), xytext=(5, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8.5)
    ax.margins(x=0.12)


def _ax_vertical_counts(ax, counts, title, ylabel, total, light, colour):
    items = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    labels = [item[0] for item in items]
    values = [item[1] for item in items]
    bars = ax.bar(labels, values, color=light, edgecolor=colour,
                  linewidth=0.8, width=0.68)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=10)
    ax.set_ylabel(ylabel)
    _style_axis(ax, "y")
    _annotate_vertical_bars(ax, bars, values, denominator=total)


def _ax_creeds_platforms(ax, inventory, top_n=10):
    counts = {}
    for row in inventory:
        platform = row.get("platform") or "Unknown"
        counts[platform] = counts.get(platform, 0) + 1
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    shown = ranked[:top_n]
    remainder = ranked[top_n:]
    if remainder:
        shown.append((f"Other ({len(remainder)} platforms)",
                      sum(value for _, value in remainder)))
    shown = list(reversed(shown))  # largest named platform on top, Other last
    labels = [label for label, _ in shown]
    values = [value for _, value in shown]
    colours = [NEUTRAL if label.startswith("Other") else CREEDS_LIGHT
               for label in labels]
    total = len(inventory)
    y = list(range(len(labels)))
    bars = ax.barh(y, values, color=colours, edgecolor=CREEDS_COLOUR, linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Usable signatures")
    ax.set_title(f"Most frequent platforms (top {min(top_n, len(ranked))})",
                 loc="left", fontsize=12, fontweight="bold", pad=10)
    _style_axis(ax, "x")
    for bar, value in zip(bars, values):
        ax.annotate(f"{value:,} ({value / total:.1%})",
                    (value, bar.get_y() + bar.get_height() / 2),
                    xytext=(5, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8.5)
    ax.margins(x=0.25)


def _ax_coverage(ax, labels, values, colours, title, ylabel, n_label):
    total = sum(values)
    bars = ax.bar(labels, values, color=colours, width=0.66)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=10)
    _style_axis(ax, "y")
    _annotate_vertical_bars(ax, bars, values, denominator=total)
    ax.text(1.0, 1.02, f"n = {total:,} {n_label}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9, color="#555555")


def _creeds_coverage_counts(match_rows):
    return [
        sum(r["match_status"] in {"cui_only", "cui_and_exact_name"} for r in match_rows),
        sum(r["match_status"] == "exact_name_only" for r in match_rows),
        sum(r["match_status"] == "not_found" for r in match_rows),
    ]


def _dsa_coverage_counts(dsa_by_disease_keys, match_rows):
    """DiSignAtlas disease keys matched by >=1 CREEDS disease, split by match type.

    A DiSignAtlas disease reached through a CUI match by any CREEDS disease is
    counted as a CUI match; otherwise as exact-name-only if reached by name.
    """
    via_cui, via_name = set(), set()
    for r in match_rows:
        keys = {k for k in r["matched_dsa_keys"].split("|") if k}
        if r["match_status"] in {"cui_only", "cui_and_exact_name"}:
            via_cui |= keys
        elif r["match_status"] == "exact_name_only":
            via_name |= keys
    via_cui &= dsa_by_disease_keys
    via_name = (via_name & dsa_by_disease_keys) - via_cui
    return [len(via_cui), len(via_name),
            len(dsa_by_disease_keys) - len(via_cui) - len(via_name)]


def _panel_figure(database):
    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], hspace=0.38, wspace=0.42)
    return fig, gs


def _finish_panel(fig, database, n_records, item_name, n_diseases, path):
    fig.suptitle(
        f"{database}: usable human inventory — {_filter_phrase()}",
        x=0.06, ha="left", fontsize=16, fontweight="bold", y=0.975,
    )
    fig.text(
        0.06, 0.945,
        f"{n_records:,} {item_name} on {n_diseases:,} distinct diseases survive the "
        f"filter (organism = human; case and control ≥{MIN_REPLICATES} samples). "
        "Every panel describes this filtered set.",
        ha="left", va="top", fontsize=10, color="#555555",
    )
    _savefig(fig, path)


def plot_disignatlas_panel(inventory, by_disease, match_rows):
    """(A) depth per disease, (B) top diseases, (C) source+assay, (D) CREEDS coverage."""
    fig, gs = _panel_figure("DiSignAtlas")
    ax_a = fig.add_subplot(gs[0, 0])
    _ax_datasets_per_disease(ax_a, by_disease, "DiSignAtlas", "datasets", DSA_COLOUR)
    _panel_letter(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    _ax_top_diseases(ax_b, by_disease, "datasets", DSA_COLOUR, DSA_LIGHT,
                     split=("n_microarray", "Microarray", "n_rnaseq", "RNA-Seq"))
    _panel_letter(ax_b, "B")

    sub = gs[1, 0].subgridspec(1, 2, wspace=0.45)
    source_counts, assay_counts = {}, {}
    for row in inventory:
        s = row.get("data_source") or "Unknown"
        a = row.get("library_strategy") or "Unknown"
        source_counts[s] = source_counts.get(s, 0) + 1
        assay_counts[a] = assay_counts.get(a, 0) + 1
    ax_c1 = fig.add_subplot(sub[0, 0])
    _ax_vertical_counts(ax_c1, source_counts, "Source repository", "Usable datasets",
                        len(inventory), DSA_LIGHT, DSA_COLOUR)
    _panel_letter(ax_c1, "C")
    ax_c2 = fig.add_subplot(sub[0, 1])
    _ax_vertical_counts(ax_c2, assay_counts, "Assay", "Usable datasets",
                        len(inventory), DSA_LIGHT, DSA_COLOUR)

    ax_d = fig.add_subplot(gs[1, 1])
    keys = {_dsa_disease_key(r) for r in inventory}
    _ax_coverage(ax_d, ["Matched by CREEDS\n(UMLS CUI)",
                        "Matched by CREEDS\n(exact name only)", "DiSignAtlas only"],
                 _dsa_coverage_counts(keys, match_rows),
                 [DSA_COLOUR, DSA_LIGHT, NEUTRAL],
                 "DiSignAtlas diseases also present in CREEDS",
                 "DiSignAtlas diseases", "DiSignAtlas diseases")
    _panel_letter(ax_d, "D")
    _finish_panel(fig, "DiSignAtlas", len(inventory), "datasets",
                  len(by_disease), DSA_PANEL_PNG)


def plot_creeds_panel(inventory, by_disease, match_rows):
    """(A) depth per disease, (B) top diseases, (C) platforms, (D) coverage in DiSignAtlas."""
    fig, gs = _panel_figure("CREEDS")
    ax_a = fig.add_subplot(gs[0, 0])
    _ax_datasets_per_disease(ax_a, by_disease, "CREEDS", "signatures", CREEDS_COLOUR)
    _panel_letter(ax_a, "A")
    ax_b = fig.add_subplot(gs[0, 1])
    _ax_top_diseases(ax_b, by_disease, "signatures", CREEDS_COLOUR)
    _panel_letter(ax_b, "B")
    ax_c = fig.add_subplot(gs[1, 0])
    _ax_creeds_platforms(ax_c, inventory)
    _panel_letter(ax_c, "C")
    ax_d = fig.add_subplot(gs[1, 1])
    _ax_coverage(ax_d, ["UMLS CUI match", "Exact name only", "Not found"],
                 _creeds_coverage_counts(match_rows),
                 [CREEDS_COLOUR, CREEDS_LIGHT, "#D46A6A"],
                 "CREEDS diseases found in DiSignAtlas",
                 "Canonical CREEDS diseases", "CREEDS diseases")
    _panel_letter(ax_d, "D")
    _finish_panel(fig, "CREEDS", len(inventory), "signatures",
                  len(by_disease), CREED_PANEL_PNG)


# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------
def _print_ceiling(ceiling):
    for c in ceiling:
        mark = "  <- inventory" if c["selected"] else ""
        print(f"  {c['stage']:<40} datasets={c['n_datasets']:>6}  "
              f"diseases={c['n_distinct_diseases']:>5}{mark}")


def run_disignatlas_analysis():
    rows = load_catalogue()
    inventory, ceiling = build_inventory(rows)
    by_disease = aggregate_by_disease(inventory)
    write_outputs(inventory, ceiling, by_disease)
    print(f"Loaded {len(rows)} datasets from {DSA_CSV.name}")
    _print_ceiling(ceiling)
    inspect_inventory(inventory, by_disease)
    plot_disignatlas_signature_ceiling(ceiling)
    print(f"\nWrote {CEILING_TSV}")
    print(f"Wrote {INVENTORY_TSV}")
    print(f"Wrote {BY_DISEASE_TSV}")
    print(f"Wrote {DSA_STAGE_PNG}")
    return inventory, ceiling, by_disease


def run_creeds_analysis():
    rows = load_creeds_metadata()
    inventory, ceiling = build_creeds_inventory(rows)
    by_disease = aggregate_creeds_by_disease(inventory)
    write_creeds_outputs(inventory, ceiling, by_disease)
    print(f"\nLoaded {len(rows)} datasets from {CREED_METADATA.name}")
    _print_ceiling(ceiling)
    inspect_creeds_inventory(inventory, by_disease)
    plot_creeds_signature_ceiling(ceiling)
    print(f"\nWrote {CREED_CEILING_TSV}")
    print(f"Wrote {CREED_INVENTORY_TSV}")
    print(f"Wrote {CREED_BY_DISEASE_TSV}")
    print(f"Wrote {CREED_STAGE_PNG}")
    return inventory, ceiling, by_disease


def run_creeds_disignatlas_overlap(dsa_inventory, creeds_inventory):
    """Run and report the conservative cross-database coverage analysis."""
    match_rows, shared_gse_rows, summary = build_creeds_disignatlas_overlap(
        dsa_inventory, creeds_inventory
    )
    merged_rows = build_merged_inventory(dsa_inventory, creeds_inventory, match_rows)
    merged_summary_rows = build_merged_disease_summary(merged_rows)
    write_creeds_disignatlas_overlap(match_rows, shared_gse_rows, summary)
    write_merged_disease_summary(merged_summary_rows)

    metrics = {row["metric"]: row["value"] for row in summary}
    print(f"\n--- CREEDS coverage in DiSignAtlas ({_filter_phrase()}) ---")
    print(f"DiSignAtlas diseases: {metrics['disignatlas_diseases']}")
    print(f"CREEDS diseases:      {metrics['creeds_diseases']}")
    print(f"common, strict CUI:   {metrics['common_diseases_strict_cui']}")
    print(f"common, CUI or exact normalised name: "
          f"{metrics['common_diseases_cui_or_exact_name']}")
    print(f"CREEDS not found conservatively: {metrics['creeds_diseases_not_found']}")
    print(f"shared GSEs: {metrics['unique_shared_gse']} unique accessions; "
          f"{metrics['shared_disease_gse_pairs']} disease-GSE pairs; "
          f"{metrics['common_diseases_with_shared_gse']} diseases with >=1 shared GSE")
    print("shared GSMs: not computable from the supplied DiSignAtlas catalogue "
          "(GSM ids are absent)")
    print(f"CREEDS candidate GSMs inside shared disease-GSE pairs: "
          f"{metrics['creeds_candidate_gsm_in_shared_gse']} (not an intersection)")

    creeds_only = [row for row in match_rows if row["match_status"] == "not_found"]
    print("top CREEDS diseases not found (automatic conservative match)")
    for row in creeds_only[:15]:
        print(f"  {row['n_creeds_signatures']:>3} signatures  "
              f"{row['creeds_disease'][:48]}")

    print(f"\nWrote {OVERLAP_SUMMARY_TSV}")
    print(f"Wrote {OVERLAP_BY_DISEASE_TSV}")
    print(f"Wrote {OVERLAP_GSE_TSV}")
    print(f"Wrote {CREEDS_ONLY_TSV}")
    print(f"Wrote merged disease summary ({len(merged_summary_rows)} rows, from "
          f"{len(merged_rows)} records) to {MERGED_SUMMARY_CSV}")
    return match_rows, shared_gse_rows, summary


def main():
    dsa_inventory, _, dsa_by_disease = run_disignatlas_analysis()
    creeds_inventory, _, creeds_by_disease = run_creeds_analysis()
    match_rows, _, _ = run_creeds_disignatlas_overlap(dsa_inventory, creeds_inventory)
    plot_disignatlas_panel(dsa_inventory, dsa_by_disease, match_rows)
    plot_creeds_panel(creeds_inventory, creeds_by_disease, match_rows)
    print(f"Wrote {DSA_PANEL_PNG}")
    print(f"Wrote {CREED_PANEL_PNG}")


if __name__ == "__main__":
    main()
