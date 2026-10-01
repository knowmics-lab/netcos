#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map the repoDB drugs (DrugBank ids) onto the LINCS L1000 perturbagens and keep
the repoDB x disease-signature pairs whose drug has an L1000 signature.

Design
------
Stage 2, part 3 of the drug-disease-pairs validation (drug leg). Input is the
pair table written by 03_ (repoDB pairs whose disease has a DiSignAtlas/CREEDS
signature, with `max_psamples`). Each pair's drug is mapped to LINCS L1000
compound perturbagens (`pert_id`, BRD-...) through the LINCS signature metadata
used by the pipeline (conf LINCS_METADATA_PATH, `lincs_sig_info_new.csv`, the
Bin Chen 2017 gold-signature table, trt_cp only).

Is the id the same? LINCS' native id is the BRD `pert_id`; the `DrugBank.ID`
column of `lincs_sig_info_new.csv` is an annotation added by the Bin Chen lab
(2016), so an id match is a match between two DrugBank annotations of
different ages, not a native identifier. It is therefore cross-checked by name.

Mapping, in order of trust (every repoDB name of the id is tried):
1. id          : repoDB `drugbank_id` == LINCS `DrugBank.ID`.
2. name        : no id match; a repoDB name of the id equals (normalised) a
                 LINCS `pert_iname` whose perturbagens carry NO DrugBank id.
3. parent_name : the name matches a perturbagen annotated with a DIFFERENT
                 DrugBank id (salt / ester / prodrug id vs parent). Same active
                 moiety, different entity: used only if ALLOW_PARENT_MATCH.
MANUAL_EXCLUDE removes adjudicated wrong id matches.

Pair-level check. repoDB v2 attaches indications of combination products and
prodrugs to their ingredient ids, so an id match can be right for the molecule
and wrong for the pair. Each pair's own repoDB drug names are compared with
the LINCS pert_iname (equal, equal after stripping salt/ester words, or in the
hand-adjudicated SYNONYMS list): agree / no_name / mismatch; mismatches are
removed when EXCLUDE_NAME_MISMATCH. Drugs that repoDB lists under two ids
(same perturbagens, same disease) are counted once.

Per mapped drug, the L1000 coverage is reported: number of gold signatures,
cell lines and the cell-line list, because NetCoS scores a drug in one cell
line at a time: a pair is scorable in cell line c only if the drug has a
signature in c. REPORT_CELL_LINES gives per-pair flags for the cell lines of
interest.

Metadata only: no connectivity score, no run ID. Coverage here is the table
the pipeline reads (`lincs_sig_info_new.csv`); a drug present in LINCS but not
in this table (or not in a given cell-line run) is not scorable by the
current pipeline.

Inputs
------
- `repodb_signature_pairs.tsv` from 03_ (run 03_ first).
- LINCS metadata from conf (LINCS_METADATA_PATH).

Outputs (flat inside the slug; filter `max_psamples >= x` for ge<x>Psamples)
------
- repodb_lincs_drug_map.tsv -- one row per repoDB drug in the 03_ pairs:
  match method (id / name / parent_name / manual_exclude / none), all repoDB
  names, name agreement, LINCS perturbagens,
  number of signatures and cell lines.
- repodb_signature_lincs_pairs.tsv -- the 03_ pairs whose drug maps to LINCS,
  with the LINCS columns and one 0/1 flag per REPORT_CELL_LINES.
- repodb_signature_lincs_summary.tsv -- one row per x: pairs / drugs /
  diseases per label before and after the LINCS filter, diseases with both
  labels, single-disease candidates.
- repodb_signature_lincs_removed_pairs.tsv -- audit: pairs removed for a name
  mismatch (or agree_approx when EXCLUDE_APPROX) and duplicate repoDB ids.
- repodb_signature_lincs_by_disease.tsv -- one row per disease after the
  filter: pair counts per label, max_psamples, single_disease_candidate.
- repodb_signature_lincs_panel_ge<x>Psamples.png -- panel A-D.

Usage
-----
    runfile('validations/drug_disease_pairs_validation/04_repodb_lincs_drug_mapping.py')
"""

import re
import textwrap

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import StrMethodFormatter

from conf import LINCS_METADATA_PATH, REPODB_FILE, OOUT_DIR, IMG_DIR

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
IMG_OUT_DIR = IMG_DIR / "validations" / SLUG

PAIRS_IN_TSV = OUT_DIR / "repodb_signature_pairs.tsv"          # from 03_
DRUG_MAP_TSV = OUT_DIR / "repodb_lincs_drug_map.tsv"
PAIRS_OUT_TSV = OUT_DIR / "repodb_signature_lincs_pairs.tsv"
SUMMARY_TSV = OUT_DIR / "repodb_signature_lincs_summary.tsv"
BY_DISEASE_TSV = OUT_DIR / "repodb_signature_lincs_by_disease.tsv"
AUDIT_TSV = OUT_DIR / "repodb_signature_lincs_removed_pairs.tsv"   # mismatches + duplicates

# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
PSAMPLES_THRESHOLDS = (3, 10, 20, 50)
SELECTED_PSAMPLES = 10
SINGLE_DISEASE_MIN_PAIRS = 10       # same definition as 03_
SINGLE_DISEASE_MIN_PAIRS_LOW = 5    # extra, looser count (few survive at 10)
ALLOW_NAME_MATCH = True
# Salt / ester / prodrug ids whose name matches a LINCS perturbagen annotated
# with the PARENT id (e.g. ramipril salt DB14208 vs ramipril DB00178). Same
# active moiety, but not the same DrugBank entity: off by default (rigorous).
ALLOW_PARENT_MATCH = False
# Pairs whose repoDB drug name disagrees with the LINCS perturbagen name
# (repoDB explodes combination products / prodrugs onto ingredient ids) are
# flagged; with this on they are removed from the scorable set.
EXCLUDE_NAME_MISMATCH = True

# Manual exclusions of wrong id matches (id -> reason).
MANUAL_EXCLUDE = {
    "DB05332": "repoDB romiplostim (peptibody) vs LINCS semaxanib: wrong annotation",
}

# Accepted synonyms / salts / INN variants between a repoDB name and a LINCS
# pert_iname (normalised keys). Every entry was adjudicated by hand (Critic,
# 24/09/2026); anything not listed and not equal after salt stripping is a
# mismatch.
SYNONYMS = [
    ("vitamine", "alphatocopherol"), ("acetylsalicylicacid", "aspirin"),
    ("indomethacin", "indometacin"), ("calcipotriene", "calcipotriol"),
    ("torasemide", "torsemide"), ("thioguanine", "tioguanine"),
    ("cyclacillin", "ciclacillin"), ("dyphylline", "diprophylline"),
    ("omacetaxinemepesuccinate", "homoharringtonine"),
    ("omacetaxine", "homoharringtonine"),
    ("mesalamine", "mesalazine"), ("epirubicin", "pidorubicine"),
    ("norethindrone", "norethisterone"), ("hydrocortisone", "cortisol"),
    ("flurandrenolide", "fludroxycortide"), ("colforsin", "forskolin"),
    ("salbutamol", "albuterol"), ("phytonadione", "phytomenadione"),
    ("vitamink1", "phytomenadione"), ("pyridostigmine", "mestinon"),
    ("sonidegib", "erismodegib"), ("flucloxacillin", "floxacillin"),
    ("tenofovirdisoproxil", "tenofovir"),
    ("mitomycin", "mitomycinc"), ("erythromycinethylsuccinate", "erythromycin"), ("paracetamol", "acetaminophen"),
    ("phytomenadione", "phylloquinone"), ("bacampicillin", "bacampicilline"),
    ("nicotinicacid", "niacin"),
]
# Near-equivalents, NOT the same molecule (enantiomer/racemate, prodrug/active
# moiety, ester/parent). Pairs agreeing only through these are tagged
# `agree_approx`: kept by default, dropped when EXCLUDE_APPROX = True.
APPROX_SYNONYMS = [
    ("esomeprazole", "omeprazole"), ("levonorgestrel", "norgestrel"),
    ("levalbuterol", "salbutamol"), ("enalapril", "enalaprilat"),
    ("etretinate", "acitretin"), ("prednicarbate", "prednisolone"),
    ("sulfasalazine", "mesalazine"), ("hydrocortamate", "hydrocortisone"),
    ("hydrocortamate", "cortisol"), ("penamecillin", "benzylpenicillin"),
]
EXCLUDE_APPROX = False
SALT_WORDS = (
    "hydrochloride", "dihydrochloride", "hydrobromide", "sodium", "disodium",
    "potassium", "calcium", "magnesium", "mesylate", "mesilate", "maleate",
    "fumarate", "succinate", "tartrate", "bitartrate", "citrate", "sulfate",
    "sulphate", "phosphate", "acetate", "propionate", "dipropionate",
    "furoate", "valerate", "butyrate", "besylate", "besilate", "tosylate",
    "bromide", "chloride", "nitrate", "lactate", "gluconate", "hyclate",
    "monohydrate", "dihydrate", "trihydrate", "anhydrous", "medoxomil",
    "cilexetil", "mofetil", "pivoxil", "axetil", "embonate", "pamoate",
    "hemifumarate", "oxalate", "napsylate", "xinafoate", "decanoate",
    "enanthate", "cypionate", "undecanoate", "palmitate", "hcl",
    "pivalate", "diacetate", "proxetil", "acetonide", "etabonate",
    "tocoferil",
)
REPORT_CELL_LINES = ("MCF7", "PC3", "A549", "HEPG2", "HT29", "A375", "VCAP", "HA1E")
LABEL_TOP_N = 20

LABELS = ("approved", "not_approved", "conflict")
PANEL_PNG = IMG_OUT_DIR / f"repodb_signature_lincs_panel_ge{SELECTED_PSAMPLES}Psamples.png"

C_APPROVED = "#2b6cb0"
C_NOT_APPROVED = "#c05621"
C_CONFLICT = "#718096"
C_LIGHT = "#cbd5e0"
GRID_COLOUR = "#d9dde3"


def _norm_drug(name):
    """Conservative drug-name key: lower case, alphanumerics only."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def _syn_map(pairs):
    m = {}
    for a, b in pairs:
        m.setdefault(a, set()).add(b)
        m.setdefault(b, set()).add(a)
    return m


_SYN, _APPROX = _syn_map(SYNONYMS), _syn_map(APPROX_SYNONYMS)


def _strip_salt(name):
    """Name key with salt / ester words removed (word-level, then normalised)."""
    words = [w for w in re.split(r"[^a-z0-9]+", str(name).lower())
             if w and w not in SALT_WORDS]
    return "".join(words)


def name_relation(repodb_name, lincs_iname):
    """
    'exact' | 'approx' | None between a repoDB name and a LINCS pert_iname.

    Salt/ester words are stripped from the repoDB side only (LINCS inames are
    parent compounds; stripping them would equate e.g. adenosine with
    adenosine-phosphate). An empty stripped key never matches.
    """
    a, b = _norm_drug(repodb_name), _norm_drug(lincs_iname)
    sa = _strip_salt(repodb_name)
    if not a or not b:
        return None
    keys = {a} | ({sa} if sa else set())
    if b in keys or any(b in _SYN.get(k, ()) for k in keys):
        return "exact"
    if any(b in _APPROX.get(k, ()) for k in keys):
        return "approx"
    return None


def names_agree(repodb_name, lincs_iname):
    return name_relation(repodb_name, lincs_iname) is not None


def load_repodb_names(path=REPODB_FILE):
    """All non-missing repoDB names per drug id and per (drug, disease) pair."""
    r = pd.read_csv(path, dtype=str, keep_default_na=False)
    r = r[~r["drug_name"].isin(["", "NA"])]
    by_id = r.groupby("drugbank_id")["drug_name"].agg(lambda v: sorted(set(v)))
    by_pair = r.groupby(["drugbank_id", "ind_id"])["drug_name"].agg(
        lambda v: sorted(set(v)))
    return by_id.to_dict(), by_pair.to_dict()


def _join(values):
    return "|".join(sorted({str(v) for v in values if isinstance(v, str) and v}))


# --------------------------------------------------------------------------- #
# LINCS
# --------------------------------------------------------------------------- #
def load_lincs(path=LINCS_METADATA_PATH):
    """LINCS gold compound signatures; 'NULL' DrugBank ids as missing."""
    lin = pd.read_csv(path, dtype=str, keep_default_na=False)
    lin = lin[(lin["pert_type"] == "trt_cp") & (lin["is_gold"] == "1")].copy()
    lin["DrugBank.ID"] = lin["DrugBank.ID"].replace({"NULL": "", "NA": ""}).str.strip()
    lin["name_key"] = lin["pert_iname"].map(_norm_drug)
    return lin


def _coverage(sub):
    return pd.Series({
        "lincs_pert_ids": _join(sub["pert_id"]),
        "lincs_pert_inames": _join(sub["pert_iname"]),
        "lincs_n_signatures": int(sub.shape[0]),
        "lincs_n_cell_lines": int(sub["cell_id"].nunique()),
        "lincs_cell_lines": _join(sub["cell_id"]),
    })


def build_drug_map(pairs, lin, names_by_id):
    """
    One row per repoDB drug: how (and whether) it maps to LINCS.

    Every repoDB name of the id is tried (repoDB v2 carries several names per
    id, and 'NA' names), not an arbitrary first one.
    """
    ids = sorted(set(pairs["drugbank_id"]))
    by_id = {k: g for k, g in lin[lin["DrugBank.ID"] != ""].groupby("DrugBank.ID")}
    by_name = {k: g for k, g in lin.groupby("name_key")}

    rows = []
    for db in ids:
        names = names_by_id.get(db, [])
        row = {"drugbank_id": db, "repodb_drug_names": "|".join(names)}
        sub, method = by_id.get(db), "id"
        if db in MANUAL_EXCLUDE:
            sub, method = None, "manual_exclude"
            row["note"] = MANUAL_EXCLUDE[db]
        elif sub is None:
            cands = [by_name[k] for k in {_norm_drug(n) for n in names} if k in by_name]
            if not cands:
                sub, method = None, "none"
            else:
                cand = pd.concat(cands).drop_duplicates("sig_id")
                other_ids = set(cand["DrugBank.ID"]) - {""}
                if other_ids:
                    row["conflicting_lincs_drugbank_ids"] = _join(other_ids)
                    method = "parent_name"
                    sub = cand if ALLOW_PARENT_MATCH else None
                    if sub is None:
                        row.update({k: v for k, v in _coverage(cand).items()
                                    if k in ("lincs_pert_ids", "lincs_pert_inames")})
                else:
                    sub, method = cand, "name"
        row["lincs_match_method"] = method
        if sub is not None:
            row.update(_coverage(sub).to_dict())
            row["name_agrees"] = int(any(names_agree(n, i) for n in names
                                         for i in set(sub["pert_iname"])))
        rows.append(row)

    dm = pd.DataFrame(rows)
    for c in ("lincs_n_signatures", "lincs_n_cell_lines", "name_agrees"):
        dm[c] = dm[c].fillna(0).astype(int) if c in dm else 0
    for c in ("lincs_pert_ids", "lincs_pert_inames", "lincs_cell_lines",
              "conflicting_lincs_drugbank_ids", "note"):
        dm[c] = dm[c].fillna("") if c in dm else ""
    usable = {"id"} | ({"name"} if ALLOW_NAME_MATCH else set()) | (
        {"parent_name"} if ALLOW_PARENT_MATCH else set())
    dm["mapped"] = dm["lincs_match_method"].isin(usable).astype(int)
    return dm.sort_values(["mapped", "lincs_n_signatures"], ascending=False)


# --------------------------------------------------------------------------- #
# Filter and summaries
# --------------------------------------------------------------------------- #
def filter_pairs(pairs, dm, names_by_pair):
    """
    Pairs whose drug maps to LINCS, with a pair-level name check and a
    de-duplication of drugs that repoDB lists under two ids.

    name_check (per pair, on the repoDB names of THAT pair's records):
      agree    -- a repoDB name of the pair matches a LINCS pert_iname;
      agree_approx -- only through APPROX_SYNONYMS (dropped if EXCLUDE_APPROX);
      no_name  -- the pair's records carry only 'NA' names (not checkable);
      mismatch -- names exist and none matches: repoDB attached this
                  indication to the id via another drug (combination product,
                  prodrug). Removed when EXCLUDE_NAME_MISMATCH.
    """
    cols = ["drugbank_id", "lincs_match_method", "lincs_pert_ids",
            "lincs_pert_inames", "lincs_n_signatures", "lincs_n_cell_lines",
            "lincs_cell_lines"]
    out = pairs.merge(dm[dm["mapped"] == 1][cols], on="drugbank_id", how="inner")

    def _check(r):
        names = names_by_pair.get((r["drugbank_id"], r["ind_id"]), [])
        if not names:
            return "no_name", ""
        rel = {name_relation(n, i) for n in names
               for i in r["lincs_pert_inames"].split("|")}
        tag = ("agree" if "exact" in rel else
               "agree_approx" if "approx" in rel else "mismatch")
        return tag, "|".join(names)

    chk = out.apply(_check, axis=1, result_type="expand")
    out["name_check"], out["repodb_pair_drug_names"] = chk[0], chk[1]
    mismatches = out[out["name_check"] == "mismatch"].copy()
    if EXCLUDE_NAME_MISMATCH:
        out = out[out["name_check"] != "mismatch"]
    if EXCLUDE_APPROX:
        mismatches = pd.concat([mismatches, out[out["name_check"] == "agree_approx"]])
        out = out[out["name_check"] != "agree_approx"]

    # Same perturbagen set + same disease under two repoDB ids = one pair.
    # Keep the id match first, then the lowest drugbank_id; a label conflict
    # between the two records is reported (label_conflict_with_duplicate).
    rank = {"id": 0, "name": 1, "parent_name": 2}
    out = (out.assign(_m=out["lincs_match_method"].map(rank))
              .sort_values(["_m", "drugbank_id"]))
    dup = out.duplicated(["lincs_pert_ids", "ind_id"], keep="first")
    duplicates = out[dup].drop(columns="_m")
    out = out[~dup].drop(columns="_m")
    key_labels = duplicates.groupby(["lincs_pert_ids", "ind_id"])["label"].agg(set)
    out["label_conflict_with_duplicate"] = [
        int(bool(key_labels.get((p, d), set()) - {l}))
        for p, d, l in zip(out["lincs_pert_ids"], out["ind_id"], out["label"])]

    cells = out["lincs_cell_lines"].str.split("|")
    for c in REPORT_CELL_LINES:
        out[f"in_{c}"] = cells.apply(lambda t, c=c: int(c in t))
    return out.reset_index(drop=True), mismatches, duplicates


def by_disease(lp):
    counts = (lp.pivot_table(index="ind_id", columns="label", values="drugbank_id",
                             aggfunc="count", fill_value=0)
                .reindex(columns=list(LABELS), fill_value=0)
                .add_prefix("n_").add_suffix("_pairs"))
    bd = (lp.groupby("ind_id").agg(ind_name=("ind_name", "first"),
                                   max_psamples=("max_psamples", "first"))
            .join(counts))
    bd["n_min_class"] = bd[["n_approved_pairs", "n_not_approved_pairs"]].min(axis=1)
    bd["single_disease_candidate"] = (bd["n_min_class"] >= SINGLE_DISEASE_MIN_PAIRS).astype(int)
    return bd.sort_values(["n_min_class", "n_approved_pairs"],
                          ascending=False).reset_index()


def summarise(pairs, lp, bd):
    rows = []
    for x in PSAMPLES_THRESHOLDS:
        p = pairs[pairs["max_psamples"] >= x]
        q = lp[lp["max_psamples"] >= x]
        d = bd[bd["max_psamples"] >= x]
        row = {"min_psamples": x}
        for lab in LABELS:
            row[f"pairs_{lab}_before"] = int((p["label"] == lab).sum())
            row[f"pairs_{lab}_lincs"] = int((q["label"] == lab).sum())
        row["pairs_lincs_by_name_only"] = int((q["lincs_match_method"] == "name").sum())
        row["drugs_before"] = p["drugbank_id"].nunique()
        row["drugs_lincs"] = q["drugbank_id"].nunique()
        row["diseases_lincs"] = d.shape[0]
        row["diseases_with_both"] = int(((d.n_approved_pairs > 0) &
                                         (d.n_not_approved_pairs > 0)).sum())
        row[f"single_disease_candidates_ge{SINGLE_DISEASE_MIN_PAIRS}"] = int(
            (d.n_min_class >= SINGLE_DISEASE_MIN_PAIRS).sum())
        row[f"single_disease_candidates_ge{SINGLE_DISEASE_MIN_PAIRS_LOW}"] = int(
            (d.n_min_class >= SINGLE_DISEASE_MIN_PAIRS_LOW).sum())
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def _style(ax, axis="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID_COLOUR)
    ax.grid(axis=axis, color=GRID_COLOUR, lw=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(StrMethodFormatter("{x:,.0f}"))


def _letter(ax, l):
    ax.text(-0.10, 1.06, l, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="bottom")


def _short(s, w=34):
    return textwrap.shorten(" ".join(str(s).split()), width=w, placeholder="…")


def _annot(ax, bars):
    for b in bars:
        ax.annotate(f"{int(b.get_height()):,}", (b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7)


def plot_panel(summary, dm, pairs, lp, bd):
    x = SELECTED_PSAMPLES
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    pos = list(range(len(summary)))
    ticks = [f"≥{t}" for t in summary.min_psamples]

    # A: pairs before/after LINCS, per label and threshold
    ax = axes[0, 0]
    w = 0.2
    for i, (lab, col) in enumerate((("approved", C_APPROVED),
                                    ("not_approved", C_NOT_APPROVED))):
        off = (i - 0.5) * 2 * w
        ax.bar([p + off - w / 2 for p in pos], summary[f"pairs_{lab}_before"], w,
               color=C_LIGHT, edgecolor=col, label=f"{lab.replace('_', ' ')}: disease signature")
        _annot(ax, ax.bar([p + off + w / 2 for p in pos], summary[f"pairs_{lab}_lincs"],
                          w, color=col, label=f"{lab.replace('_', ' ')}: + drug in L1000"))
    ax.set_xticks(pos, ticks)
    ax.set_xlabel("Psamples threshold (paired case/control samples)")
    ax.set_ylabel("distinct repoDB drug–disease pairs")
    ax.set_title("Pairs before and after requiring an L1000 drug signature", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    _style(ax)
    _letter(ax, "A")

    # B: drug mapping outcome (drugs of the >= x pairs)
    ax = axes[0, 1]
    ids = set(pairs.loc[pairs.max_psamples >= x, "drugbank_id"])
    sub = dm[dm.drugbank_id.isin(ids)]
    order = ["id", "name", "parent_name", "manual_exclude", "none"]
    vals = [int((sub.lincs_match_method == m).sum()) for m in order]
    bars = ax.bar(order, vals, color=[C_APPROVED, "#63b3ed", C_CONFLICT, C_CONFLICT, C_LIGHT])
    _annot(ax, bars)
    ax.set_ylabel("repoDB drugs")
    ax.set_title("How repoDB drugs map to L1000 perturbagens\n"
                 f"(id = DrugBank id; name = exact name, no id; drugs of the ≥{x} Psamples pairs)",
                 fontsize=11)
    _style(ax)
    _letter(ax, "B")

    # C: cell-line coverage of the >= x LINCS pairs
    ax = axes[1, 0]
    q = lp[lp.max_psamples >= x]
    cov = pd.DataFrame({lab: [int(q.loc[q.label == lab, f"in_{c}"].sum())
                              for c in REPORT_CELL_LINES]
                        for lab in ("approved", "not_approved")}, index=REPORT_CELL_LINES)
    yy = range(len(cov))
    ax.barh([y + 0.2 for y in yy], cov["approved"], 0.4, color=C_APPROVED, label="approved")
    ax.barh([y - 0.2 for y in yy], cov["not_approved"], 0.4, color=C_NOT_APPROVED,
            label="not approved")
    ax.set_yticks(list(yy), cov.index)
    ax.invert_yaxis()
    ax.set_xlabel("pairs whose drug has a signature in that cell line")
    ax.set_title(f"Scorable pairs per L1000 cell line (≥{x} Psamples)", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    _style(ax, "x")
    _letter(ax, "C")

    # D: top single-disease candidates after the filter
    ax = axes[1, 1]
    d = bd[bd.max_psamples >= x].head(LABEL_TOP_N).iloc[::-1]
    yy = range(len(d))
    ax.barh([y + 0.2 for y in yy], d["n_approved_pairs"], 0.4, color=C_APPROVED, label="approved")
    ax.barh([y - 0.2 for y in yy], d["n_not_approved_pairs"], 0.4, color=C_NOT_APPROVED,
            label="not approved")
    ax.axvline(SINGLE_DISEASE_MIN_PAIRS, color=GRID_COLOUR, ls="--", lw=1)
    ax.set_yticks(list(yy), [f"{_short(n)}  (≥{m})" for n, m in
                             zip(d["ind_name"], d["max_psamples"])], fontsize=8)
    ax.set_xlabel("distinct pairs with drug in L1000")
    ax.set_title(f"Top {LABEL_TOP_N} diseases by the smaller class, after the L1000 filter\n"
                 f"(≥{x} Psamples; dashed = {SINGLE_DISEASE_MIN_PAIRS} pairs)", fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    _style(ax, "x")
    _letter(ax, "D")

    fig.suptitle("repoDB × disease signatures × L1000 drug signatures "
                 "(lincs_sig_info_new.csv, gold trt_cp)", fontsize=13, y=1.0)
    fig.tight_layout()
    PANEL_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PANEL_PNG, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    pairs = pd.read_csv(PAIRS_IN_TSV, sep="\t", dtype=str, keep_default_na=False)
    pairs["max_psamples"] = pairs["max_psamples"].astype(int)
    lin = load_lincs()
    names_by_id, names_by_pair = load_repodb_names()

    dm = build_drug_map(pairs, lin, names_by_id)
    lp, mismatches, duplicates = filter_pairs(pairs, dm, names_by_pair)
    bd = by_disease(lp)
    summary = summarise(pairs, lp, bd)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dm.to_csv(DRUG_MAP_TSV, sep="\t", index=False)
    lp.to_csv(PAIRS_OUT_TSV, sep="\t", index=False)
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    bd.to_csv(BY_DISEASE_TSV, sep="\t", index=False)
    pd.concat([mismatches.assign(removed_because="name_" + mismatches["name_check"]),
               duplicates.assign(removed_because="duplicate_repodb_id")]
              ).to_csv(AUDIT_TSV, sep="\t", index=False)
    plot_panel(summary, dm, pairs, lp, bd)

    print(f"LINCS: {lin.pert_id.nunique():,} perturbagens, "
          f"{(lin['DrugBank.ID'] != '').groupby(lin.pert_id).any().sum():,} with a DrugBank id")
    print("repoDB drugs by match method:", dm.lincs_match_method.value_counts().to_dict())
    print(f"id/name matches whose names never agree (review): "
          f"{int(((dm.mapped == 1) & (dm.name_agrees == 0) & (dm.repodb_drug_names != '')).sum())}"
          f" (+{int(((dm.mapped == 1) & (dm.repodb_drug_names == '')).sum())} without any repoDB name)")
    print(f"pairs: name_check {lp.name_check.value_counts().to_dict()}; "
          f"mismatches {'removed' if EXCLUDE_NAME_MISMATCH else 'kept'}: {len(mismatches)}; "
          f"duplicate ids collapsed: {len(duplicates)}")
    if len(mismatches):
        print("\n--- name mismatches (repoDB pair names vs LINCS pert_iname) ---")
        print(mismatches[["drugbank_id", "repodb_pair_drug_names", "lincs_pert_inames",
                          "ind_name", "label"]].drop_duplicates(
                              ["drugbank_id", "repodb_pair_drug_names"]).to_string(index=False))
    with pd.option_context("display.width", 200):
        print(summary.T.to_string(header=False))
        cols = ["ind_name", "max_psamples", "n_approved_pairs", "n_not_approved_pairs",
                "n_conflict_pairs"]
        print(f"\n--- diseases with >= {SINGLE_DISEASE_MIN_PAIRS_LOW} pairs in both classes ---")
        print(bd[bd.n_min_class >= SINGLE_DISEASE_MIN_PAIRS_LOW][cols].to_string(index=False))


if __name__ == "__main__":
    main()
