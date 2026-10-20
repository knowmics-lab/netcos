#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recover repoDB drug-indication pairs lost by the exact UMLS CUI join because
the indication is coded at a finer granularity than the disease signature
(e.g. a disease subtype vs its parent), by climbing a disease hierarchy
(UMLS source vocabularies NCIt and MeSH, or MONDO), and count how many
recovered pairs survive the LINCS L1000 filter.

Design
------
A UMLS CUI is a flat concept identifier: the hierarchy lives in the source
vocabularies, not in the CUI. ROLLUP_HIERARCHY (conf) selects the hierarchy:

- "umls"  : the source-vocabulary hierarchies distributed in UMLS, read through
            the UTS REST API (version UMLS_VERSION). Each CUI is mapped to its
            NCIt and MeSH codes (atoms of the CUI); each code is climbed
            through its `parents` in the same source. Sources are tried in
            UMLS_SOURCES order: the first source in which the CUI reaches a
            signature CUI gives the target; the others are not used for that
            CUI. API responses are cached on disk (UMLS_CACHE_FILE).
- "mondo" : MONDO is_a edges, reached through MONDO's UMLS cross-references
            (`xrefs` UMLS:Cxxxxxxx, deprecated MONDO classes ignored).

A hierarchy node is a source code ("NCI:C2926", "MSH:D002289") or a MONDO id.
repoDB CUIs without a node in the selected hierarchy cannot be rolled up.

The roll-up is threshold-aware. At each Psamples threshold x, every repoDB
indication CUI WITHOUT its own signature at >= x (including CUIs whose own
signature exists only below x, flagged own_signature_below_x) is climbed from
each node carrying it, up to MAX_LEVELS parent steps; the nearest node that
carries a CUI with a DiSignAtlas/CREEDS signature at >= x is the
roll-up target (level 0 = the same node, i.e. two CUIs for one concept;
level 1 = parent; level 2 = grandparent). The nearest level wins even when a
farther ancestor has a larger signature. Ties at the same level keep the
target with the largest max_psamples; all tied targets are listed.

Every level-0/level-1 candidate carries a curated decision (manual
LLM-assisted inspection of child vs target). With "mondo" the decisions are
ROLLUP_DECISIONS, keyed by the child CUI. With "umls" they are
ROLLUP_DECISIONS_UMLS, keyed by (child CUI, target CUI); the exclude_infection
entries of ROLLUP_DECISIONS apply to the child whatever its target:
- keep              : the child is a genuine subtype / clinical form of the
                      target and the target signature is a fair proxy.
- exclude_infection : an infectious indication; the drug acts on the pathogen,
                      not by reversing the host transcriptional response.
- exclude_wrong     : the target is a different disease or a class too
                      generic to stand for the child.
Level-2 candidates are excluded by rule (EXCLUDE_LEVEL2); candidates without a
curated decision are excluded and listed as `uncurated`. CUIs whose target
changes across thresholds are printed for checking the curated note.

repoDB `conflict` pairs and the records in MANUAL_PAIR_EXCLUDE are not rolled
up. Kept pairs take the TARGET as their effective disease (its CUI,
max_psamples and signature evidence) and go through the drug leg of 04_
(DrugBank -> LINCS map, pair-level name check on the child's repoDB names,
duplicate-id collapse) with the 04_ parameters. Each recovered pair is then
compared with the exact pairs on the same target and perturbagen set,
conflict pairs included: conflict_with_exact, label_disagrees_with_exact and
redundant_with_exact pairs are removed, as are recovered pairs whose children
give one perturbagen set different labels (label_conflict); remaining
duplicates on one target are counted once.

The exact baseline in the summary is the 04_ set (conflict pairs excluded);
infectious exact indications are not removed from it.

Metadata only: no connectivity score, no run ID. Counts depend on the
hierarchy release (UMLS version or MONDO version IRI, printed and written to
the summary).

Inputs
------
- repoDB table from conf (REPODB_FILE).
- 01_ outputs `creeds_disignatlas_merged_summary_by_disease_ge<x>Psamples.csv`
  (read through 03_).
- LINCS metadata from conf (LINCS_METADATA_PATH, read through 04_).
- "umls": UTS REST API (UMLS_API_KEY from local, UMLS_VERSION, UMLS_SOURCES,
  UMLS_CACHE_FILE from conf).
- "mondo": MONDO JSON release from conf (MONDO_FILE; obographs `mondo.json`
  from the monarch-initiative/mondo GitHub releases).

Outputs (flat inside the slug; with "umls" every file name ends in `_umls`)
------
- umls_rollup_candidates.tsv -- one row per (threshold, repoDB CUI) with a
  roll-up target: level, target CUI(s)/name, hierarchy nodes, decision, pair
  counts per label.
- umls_rollup_lincs_pairs.tsv -- recovered pairs kept AND scorable in LINCS,
  one row per (threshold, pair), with the child indication, the target and
  the 04_ LINCS columns.
- umls_rollup_removed_pairs.tsv -- audit: recovered pairs removed in the drug
  leg, with removed_because and the threshold.
- umls_rollup_summary.tsv -- one row per x: exact baseline after LINCS, the
  recovered pairs per decision and label before and after LINCS, removals
  against the exact pairs, target diseases, new disease x label cells and
  the pairs carried by the three largest targets, the split between
  children with no own signature and children whose own signature is below
  x, and the gain without the child contributing most pairs.
- umls_rollup_by_target.tsv -- kept recovered pairs per (x, target disease),
  with the children and the exact pair counts of the same target.
- umls_rollup_panel_ge<x>Psamples.png (image dir) -- recovered pairs per
  decision class and label, before and after LINCS, vs the exact baseline.

Usage
-----
    runfile('validations/drug_disease_pairs_validation/05_umls_hierarchy_rollup.py')
"""

import importlib.util
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from conf import (IMG_DIR, MONDO_FILE, OOUT_DIR, ROLLUP_HIERARCHY,
                  UMLS_CACHE_FILE, UMLS_SOURCES, UMLS_VERSION)

# --------------------------------------------------------------------------- #
# Sibling stages (03_ disease leg, 04_ drug leg) reused as modules
# --------------------------------------------------------------------------- #
_HERE = Path(__file__).resolve().parent


def _load_stage(filename, name):
    spec = importlib.util.spec_from_file_location(name, _HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S03 = _load_stage("03_repodb_signature_intersection.py", "stage03")
S04 = _load_stage("04_repodb_lincs_drug_mapping.py", "stage04")

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
SLUG = "drug_disease_pairs_validation"
OUT_DIR = OOUT_DIR / "validations" / SLUG
IMG_OUT_DIR = IMG_DIR / "validations" / SLUG

_SFX = "_umls" if ROLLUP_HIERARCHY == "umls" else ""
CANDIDATES_TSV = OUT_DIR / f"umls_rollup_candidates{_SFX}.tsv"
PAIRS_OUT_TSV = OUT_DIR / f"umls_rollup_lincs_pairs{_SFX}.tsv"
AUDIT_TSV = OUT_DIR / f"umls_rollup_removed_pairs{_SFX}.tsv"
SUMMARY_TSV = OUT_DIR / f"umls_rollup_summary{_SFX}.tsv"
BY_TARGET_TSV = OUT_DIR / f"umls_rollup_by_target{_SFX}.tsv"

# --------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------- #
MAX_LEVELS = 2                 # parent steps climbed (candidates reported)
EXCLUDE_LEVEL2 = True          # grandparent targets are never used
SELECTED_PSAMPLES = S04.SELECTED_PSAMPLES
PANEL_PNG = IMG_OUT_DIR / f"umls_rollup_panel_ge{SELECTED_PSAMPLES}Psamples{_SFX}.png"
UTS_BASE = "https://uts-ws.nlm.nih.gov/rest"
UTS_SLEEP = 0.06               # seconds between uncached UTS calls (< 20 req/s)

LEVEL_NAMES = {0: "same_node", 1: "parent", 2: "grandparent"}
DECISION_CLASSES = ("keep", "exclude_infection", "exclude_wrong",
                    "exclude_level2", "uncurated")
LABELS = ("approved", "not_approved")

# Single repoDB records excluded from the roll-up ((DrugBank id, CUI) -> reason).
MANUAL_PAIR_EXCLUDE = {
    ("DB01229", "C0017612"): "paclitaxel - open-angle glaucoma: implausible repoDB record",
    ("DB00368", "C0034067"): "norepinephrine - pulmonary emphysema: implausible repoDB record",
}

# Curated decisions for level-0/1 candidates (repoDB CUI -> (decision, note)),
# from manual LLM-assisted inspection of each child indication vs its target.
ROLLUP_DECISIONS = {
    "C0004352": ("keep", "autistic disorder is within autism spectrum disorder"),
    "C0005283": ("keep", "beta-thalassaemia is a thalassaemia"),
    "C0005587": ("keep", "bipolar depression is a bipolar disorder phase"),
    "C0007102": ("keep", "malignant colon tumour within colon tumour"),
    "C0007115": ("keep", "malignant thyroid neoplasm within thyroid tumour"),
    "C0008513": ("keep", "chorioretinitis is posterior uveitis"),
    "C0008677": ("keep", "chronic bronchitis is a COPD phenotype"),
    "C0010308": ("keep", "congenital form of hypothyroidism"),
    "C0014547": ("keep", "focal epilepsy is an epilepsy"),
    "C0014869": ("keep", "peptic (reflux) oesophagitis is an oesophagitis"),
    "C0017605": ("keep", "angle-closure glaucoma is a glaucoma"),
    "C0017612": ("keep", "open-angle glaucoma is a glaucoma"),
    "C0018621": ("keep", "hay fever is seasonal allergic rhinitis"),
    "C0018802": ("keep", "congestive heart failure is heart failure"),
    "C0019147": ("keep", "hepatic coma is the severe stage of hepatic encephalopathy"),
    "C0022081": ("keep", "iritis is anterior uveitis"),
    "C0022575": ("keep", "keratoconjunctivitis sicca is dry eye"),
    "C0023465": ("keep", "acute monocytic leukaemia is an AML subtype"),
    "C0023493": ("keep", "adult T-cell leukaemia/lymphoma is a T-cell leukaemia"),
    "C0024305": ("keep", "non-Hodgkin lymphoma is a lymphoma"),
    "C0024623": ("keep", "malignant gastric neoplasm within gastric tumour"),
    "C0027145": ("keep", "myxoedema is severe hypothyroidism"),
    "C0030809": ("keep", "pemphigus vulgaris is the main pemphigus form"),
    "C0033246": ("keep", "proctitis (ulcerative proctitis) within IBD"),
    "C0034067": ("keep", "emphysema is a COPD phenotype"),
    "C0035457": ("keep", "perennial allergic rhinitis is allergic rhinitis"),
    "C0036421": ("keep", "systemic sclerosis is the systemic form of scleroderma"),
    "C0039520": ("keep", "tenosynovitis is tendon-sheath inflammation, close to tendinitis"),
    "C0042165": ("keep", "anterior uveitis is a uveitis"),
    "C0079731": ("keep", "B-cell lymphoma is a lymphoma"),
    "C0079773": ("keep", "cutaneous T-cell lymphoma is a T-cell lymphoma"),
    "C0085159": ("keep", "seasonal affective disorder is a depressive disorder"),
    "C0085580": ("keep", "essential hypertension is hypertension"),
    "C0085695": ("keep", "chronic gastritis is a gastritis"),
    "C0149782": ("keep", "lung squamous cell carcinoma is a lung cancer"),
    "C0152081": ("keep", "pustular psoriasis is a psoriasis"),
    "C0153633": ("keep", "malignant brain neoplasm within brain neoplasms"),
    "C0155877": ("keep", "allergic asthma is asthma"),
    "C0155880": ("keep", "intrinsic asthma is asthma"),
    "C0221074": ("keep", "postpartum depression is a depressive disorder"),
    "C0235480": ("keep", "paroxysmal atrial fibrillation is atrial fibrillation"),
    "C0238298": ("keep", "myxoedema coma is decompensated hypothyroidism"),
    "C0264115": ("keep", "drug-induced osteoporosis is osteoporosis"),
    "C0267390": ("keep", "ulcerative rectosigmoiditis is distal ulcerative colitis"),
    "C0269658": ("keep", "mild pre-eclampsia is pre-eclampsia"),
    "C0278701": ("keep", "gastric adenocarcinoma is gastric cancer"),
    "C0278996": ("keep", "malignant head and neck neoplasm within head and neck neoplasms"),
    "C0280324": ("keep", "laryngeal SCC is a head and neck SCC"),
    "C0302592": ("keep", "cervix carcinoma is cervical cancer"),
    "C0334579": ("keep", "anaplastic astrocytoma is an astrocytoma"),
    "C0338908": ("keep", "mixed anxiety-depressive disorder, depressive component"),
    "C0341950": ("keep", "severe pre-eclampsia is pre-eclampsia"),
    "C0342753": ("keep", "late-onset Pompe disease is Pompe disease"),
    "C0342881": ("keep", "homozygous FH is familial hypercholesterolaemia"),
    "C0345905": ("keep", "intrahepatic cholangiocarcinoma is a cholangiocarcinoma"),
    "C0346161": ("keep", "malignant epithelial ovarian tumour is ovarian cancer"),
    "C0349543": ("keep", "brain glioblastoma is glioblastoma"),
    "C0349636": ("keep", "pre-B ALL is an ALL"),
    "C0349649": ("keep", "pulmonary LAM is LAM"),
    "C0349790": ("keep", "asthma exacerbation is asthma"),
    "C0409818": ("keep", "CINCA/NOMID is the severe CAPS phenotype"),
    "C0555198": ("keep", "malignant glioma is a glioma"),
    "C0585442": ("keep", "osteosarcoma of bone is osteosarcoma"),
    "C0730285": ("keep", "diabetic macular oedema is a complication within diabetic retinopathy"),
    "C0745138": ("keep", "hypertensive urgency is hypertension"),
    "C0751362": ("keep", "narcolepsy with cataplexy is narcolepsy"),
    "C0751967": ("keep", "relapsing-remitting MS is MS"),
    "C0853193": ("keep", "bipolar I is a bipolar disorder"),
    "C0855112": ("keep", "refractory DLBCL is DLBCL"),
    "C0917796": ("keep", "two CUIs for Leber hereditary optic neuropathy"),
    "C1282500": ("keep", "colon cancer metastasis, same primary tumour"),
    "C1319317": ("keep", "pharyngeal SCC is a head and neck SCC"),
    "C1328504": ("keep", "hormone-refractory prostate cancer is prostate cancer"),
    "C1961099": ("keep", "precursor T-lymphoblastic leukaemia is a T-cell leukaemia"),
    "C1961835": ("keep", "Gaucher type 1 is Gaucher disease"),
    "C4509816": ("keep", "squamous NSCLC is NSCLC"),
    "C0002351": ("exclude_wrong", "altitude sickness is not a respiratory disease class proxy"),
    "C0003492": ("exclude_wrong", "congenital malformation vs generic aortic disease"),
    "C0003872": ("exclude_wrong", "psoriatic arthritis is not rheumatoid arthritis"),
    "C0006267": ("exclude_wrong", "bronchiectasis is not COPD"),
    "C0008526": ("exclude_wrong", "choroiditis is not a vasculitis"),
    "C0008728": ("exclude_wrong", "EGPA is not granulomatosis with polyangiitis"),
    "C0011606": ("exclude_wrong", "target too generic (inflammatory dermatosis)"),
    "C0011616": ("exclude_wrong", "target too generic (inflammatory dermatosis)"),
    "C0018671": ("exclude_wrong", "target too generic (neoplasm)"),
    "C0023443": ("exclude_wrong", "hairy cell leukaemia is not CLL"),
    "C0028840": ("exclude_wrong", "ocular hypertension is not systemic hypertension"),
    "C0029408": ("exclude_wrong", "osteoarthritis vs generic arthritis"),
    "C0032019": ("exclude_wrong", "target too generic (endocrine tumour)"),
    "C0035460": ("exclude_wrong", "vasomotor rhinitis is non-allergic"),
    "C0036508": ("exclude_wrong", "target too generic (inflammatory dermatosis)"),
    "C0038013": ("exclude_wrong", "ankylosing spondylitis is not rheumatoid arthritis"),
    "C0039103": ("exclude_wrong", "target too generic (arthritis)"),
    "C0042109": ("exclude_wrong", "target too generic (inflammatory dermatosis)"),
    "C0042373": ("exclude_wrong", "target too generic (cardiovascular diseases)"),
    "C0085697": ("exclude_wrong", "chronic pyelonephritis is not a CKD proxy"),
    "C0149526": ("exclude_wrong", "urticaria is not atopic dermatitis"),
    "C0206754": ("exclude_wrong", "neuroendocrine tumours are mostly gastro-entero-pancreatic, not endocrine-gland tumours"),
    "C0271093": ("exclude_wrong", "Stargardt disease is not age-related macular degeneration"),
    "C0281361": ("exclude_wrong", "target too generic (adenocarcinoma)"),
    "C0423775": ("exclude_wrong", "target too generic (exanthema)"),
    "C0750952": ("exclude_wrong", "biliary tract cancer is not liver cancer"),
    "C0752347": ("exclude_wrong", "target too generic (dementia)"),
    "C0917798": ("exclude_wrong", "target too generic (ischaemia)"),
    "C1292758": ("exclude_wrong", "precursor T-lymphoblastic lymphoma is a T-ALL-type disease, not a mature T-cell lymphoma"),
    "C0149925": ("exclude_wrong", "small cell lung carcinoma is neuroendocrine; the lung cancer signature is NSCLC-dominated"),
    "C0220756": ("exclude_wrong", "Niemann-Pick type C (NPC1/2) is distinct from types A/B; target dataset subtype not verified"),
    "C1301363": ("exclude_wrong", "BPDCN is myeloid-derived, not a T-cell lymphoma"),
    "C1318543": ("exclude_wrong", "tenosynovial giant cell tumour is not giant cell tumour of bone"),
    "C1527407": ("exclude_wrong", "eosinophilic pneumonia is non-infectious; pneumonia signature is not a proxy"),
    "C1561643": ("exclude_wrong", "target too generic (kidney disease)"),
    "C3495559": ("exclude_wrong", "juvenile idiopathic arthritis vs generic arthritis"),
    "C0006846": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0006849": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0009443": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0009768": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0013882": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0014121": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0014836": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0016436": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0023241": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0023290": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0024110": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0024535": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0024536": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0024537": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0025294": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0025306": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0030636": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0031157": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0032290": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0032302": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0032308": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0034686": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0035585": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0038160": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0040262": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0041321": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0041327": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0043167": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0085423": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0152068": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0152072": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0153252": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0154639": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0155862": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0206178": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0206525": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0239295": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0276075": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0276648": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0338395": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0341503": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0519030": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0524688": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0729524": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0854215": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0854248": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C0854328": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C1328252": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C1562585": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
    "C1692886": ("exclude_infection", "infectious indication: the drug acts on the pathogen"),
}

# Curated decisions for "umls" level-0/1 candidates
# ((repoDB CUI, target CUI) -> (decision, note)).
ROLLUP_DECISIONS_UMLS = {
}

C_APPROVED = S04.C_APPROVED
C_NOT_APPROVED = S04.C_NOT_APPROVED
GRID_COLOUR = S04.GRID_COLOUR


# --------------------------------------------------------------------------- #
# MONDO
# --------------------------------------------------------------------------- #
def _mondo_id(iri):
    return iri.rsplit("/", 1)[1] if "/obo/MONDO_" in iri else None


def load_mondo(path=MONDO_FILE):
    """
    (cui -> set of MONDO ids, MONDO id -> label, child -> set of parents,
    version IRI). Deprecated classes are dropped; only is_a edges between
    MONDO classes are kept.
    """
    with open(path, encoding="utf-8") as fh:
        graph = json.load(fh)["graphs"][0]
    cui2m, label, deprecated = {}, {}, set()
    for node in graph["nodes"]:
        m = _mondo_id(node["id"])
        if m is None:
            continue
        meta = node.get("meta", {})
        if meta.get("deprecated"):
            deprecated.add(m)
            continue
        label[m] = node.get("lbl", "")
        for x in meta.get("xrefs", []):
            if x["val"].startswith("UMLS:"):
                cui2m.setdefault(x["val"][5:], set()).add(m)
    parents = {}
    for e in graph["edges"]:
        if e["pred"] != "is_a":
            continue
        c, p = _mondo_id(e["sub"]), _mondo_id(e["obj"])
        if c and p and c not in deprecated and p not in deprecated:
            parents.setdefault(c, set()).add(p)
    version = graph.get("meta", {}).get("version", "")
    return cui2m, label, parents, version


# --------------------------------------------------------------------------- #
# UMLS source vocabularies (UTS REST API, cached)
# --------------------------------------------------------------------------- #
def _uts_get(session, path, api_key):
    """All result items of one UTS endpoint (paged); [] when UTS answers 404."""
    items, page = [], 1
    while True:
        r = session.get(f"{UTS_BASE}{path}",
                        params={"apiKey": api_key, "pageSize": 200, "pageNumber": page},
                        timeout=60)
        time.sleep(UTS_SLEEP)
        if r.status_code == 404:
            return items
        r.raise_for_status()
        body = r.json()
        res = body.get("result", [])
        items.extend(res if isinstance(res, list) else [res])
        if page >= int(body.get("pageCount", 1) or 1):
            return items
        page += 1


def load_umls(cuis, max_levels=MAX_LEVELS):
    """
    (cui -> set of nodes, node -> label, child node -> set of parent nodes,
    version) for the given CUIs, nodes being "SAB:code" in UMLS_SOURCES.
    Parents are fetched up to max_levels steps above the CUIs' own codes.
    """
    import requests
    from local import UMLS_API_KEY

    cache = {"version": UMLS_VERSION, "atoms": {}, "parents": {}}
    if UMLS_CACHE_FILE.exists():
        with open(UMLS_CACHE_FILE, encoding="utf-8") as fh:
            cache = json.load(fh)
    sabs = ",".join(UMLS_SOURCES)
    session = requests.Session()
    n_calls = 0

    todo = sorted(c for c in cuis if c not in cache["atoms"])
    for i, cui in enumerate(todo, 1):
        atoms = _uts_get_atoms(session, cui, sabs, UMLS_API_KEY)
        cache["atoms"][cui] = sorted({(a["rootSource"], a["code"].rsplit("/", 1)[1], a["name"])
                                      for a in atoms if a.get("rootSource") in UMLS_SOURCES})
        n_calls += 1
        if i % 200 == 0:
            print(f"UTS atoms: {i}/{len(todo)}")
            _save_cache(cache)

    cui2n, label = {}, {}
    for cui in cuis:
        for sab, code, name in cache["atoms"].get(cui, []):
            n = f"{sab}:{code}"
            cui2n.setdefault(cui, set()).add(n)
            label.setdefault(n, name)

    parents, frontier = {}, set().union(*cui2n.values()) if cui2n else set()
    for _ in range(max_levels):
        nxt = set()
        for n in sorted(frontier):
            if n not in cache["parents"]:
                sab, code = n.split(":", 1)
                res = _uts_get(session, f"/content/{UMLS_VERSION}/source/{sab}/{code}/parents",
                               UMLS_API_KEY)
                cache["parents"][n] = sorted({(sab, p["ui"], p.get("name", "")) for p in res})
                n_calls += 1
                if n_calls % 200 == 0:
                    print(f"UTS calls: {n_calls}")
                    _save_cache(cache)
            for sab, ui, name in cache["parents"][n]:
                p = f"{sab}:{ui}"
                parents.setdefault(n, set()).add(p)
                label.setdefault(p, name)
                nxt.add(p)
        frontier = nxt
    if n_calls:
        _save_cache(cache)
    return cui2n, label, parents, f"UMLS {cache['version']} ({','.join(UMLS_SOURCES)})"


def _uts_get_atoms(session, cui, sabs, api_key):
    """Atoms of one CUI restricted to the given sources ([] when none)."""
    items, page = [], 1
    while True:
        r = session.get(f"{UTS_BASE}/content/{UMLS_VERSION}/CUI/{cui}/atoms",
                        params={"apiKey": api_key, "sabs": sabs, "pageSize": 200,
                                "pageNumber": page}, timeout=60)
        time.sleep(UTS_SLEEP)
        if r.status_code == 404:
            return items
        r.raise_for_status()
        body = r.json()
        items.extend(body.get("result", []))
        if page >= int(body.get("pageCount", 1) or 1):
            return items
        page += 1


def _save_cache(cache):
    UMLS_CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = UMLS_CACHE_FILE.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(cache, fh)
    tmp.replace(UMLS_CACHE_FILE)


def _source(node):
    """Source of a hierarchy node ("NCI", "MSH"; "MONDO" for MONDO ids)."""
    return node.split(":", 1)[0] if ":" in node else "MONDO"


def _decision(cui, target):
    """Curated decision of a level-0/1 candidate for the selected hierarchy."""
    if ROLLUP_HIERARCHY == "mondo":
        return ROLLUP_DECISIONS.get(cui, ("uncurated", ""))
    if cui in ROLLUP_DECISIONS and ROLLUP_DECISIONS[cui][0] == "exclude_infection":
        return ROLLUP_DECISIONS[cui]
    return ROLLUP_DECISIONS_UMLS.get((cui, target), ("uncurated", ""))


def ancestors(m, parents, max_levels):
    """{ancestor: shortest is_a distance} up to max_levels steps."""
    out, frontier = {}, {m}
    for d in range(1, max_levels + 1):
        frontier = {p for x in frontier for p in parents.get(x, ())}
        for p in frontier:
            out.setdefault(p, d)
    return out


# --------------------------------------------------------------------------- #
# Candidates (threshold-aware)
# --------------------------------------------------------------------------- #
def find_candidates(pairs, max_ps, x, ind_names, hier):
    """
    Roll-up candidates at threshold x: one row per repoDB CUI WITHOUT its own
    signature at >= x whose nearest ancestor node (<= MAX_LEVELS) carries a
    CUI with a signature at >= x. A CUI that has its own signature only below
    x is a candidate at x (own_signature_below_x = 1). Nodes are searched one
    source at a time, in UMLS_SOURCES order for "umls"; the first source with
    a hit gives the target.
    """
    cui2m, label, parents, _ = hier
    order = list(UMLS_SOURCES) if ROLLUP_HIERARCHY == "umls" else ["MONDO"]
    sig = {c for c, v in max_ps.items() if v >= x}
    sig_by_m = {}
    for cui in sig:
        for m in cui2m.get(cui, ()):
            sig_by_m.setdefault(m, set()).add(cui)

    rows = []
    for cui in sorted(set(pairs["ind_id"]) - sig):
        hits = {}                                   # target CUI -> (level, node)
        for src in order:
            for m in (n for n in cui2m.get(cui, ()) if _source(n) == src):
                reach = {m: 0, **ancestors(m, parents, MAX_LEVELS)}
                for a, d in reach.items():
                    for t in sig_by_m.get(a, ()):
                        if t != cui and (t not in hits or d < hits[t][0]):
                            hits[t] = (d, a)
            if hits:
                break
        if not hits:
            continue
        level = min(d for d, _ in hits.values())
        tied = sorted((t for t, (d, _) in hits.items() if d == level),
                      key=lambda t: (-max_ps[t], t))
        target = tied[0]
        dec, note = _decision(cui, target)
        if level == 2 and EXCLUDE_LEVEL2:
            dec, note = "exclude_level2", "grandparent target"
        cnt = pairs[pairs["ind_id"] == cui]["label"].value_counts()
        rows.append({
            "min_psamples": x,
            "ind_id": cui, "ind_name": ind_names.get(cui, ""),
            "own_signature_below_x": int(cui in max_ps),
            "rollup_level": LEVEL_NAMES[level],
            "target_cui": target, "target_name": ind_names.get(target, ""),
            "target_source": _source(hits[target][1]),
            "target_node": hits[target][1],
            "target_node_label": label.get(hits[target][1], ""),
            "tied_target_cuis": "|".join(tied),
            "target_max_psamples": max_ps[target],
            "child_nodes": "|".join(sorted(cui2m[cui])),
            "decision": dec, "decision_note": note,
            "n_approved_pairs": int(cnt.get("approved", 0)),
            "n_not_approved_pairs": int(cnt.get("not_approved", 0)),
            "n_conflict_pairs": int(cnt.get("conflict", 0)),
        })
    return pd.DataFrame(rows)


def rolled_pairs(pairs, cand, evidence, max_ps):
    """
    repoDB pairs (conflict pairs and MANUAL_PAIR_EXCLUDE removed) of the
    candidate CUIs at one threshold, re-keyed to the target disease (effective
    ind_id = target CUI), with the target's signature evidence and max_psamples.
    """
    drop = pairs.apply(lambda r: (r["drugbank_id"], r["ind_id"]) in MANUAL_PAIR_EXCLUDE,
                       axis=1)
    pairs = pairs[(pairs["label"] != "conflict") & ~drop]
    rp = pairs.merge(cand[["ind_id", "rollup_level", "target_cui", "target_name",
                           "decision", "own_signature_below_x"]],
                     on="ind_id", how="inner")
    rp = rp.rename(columns={"ind_id": "child_ind_id", "ind_name": "child_ind_name"})
    rp["ind_id"], rp["ind_name"] = rp["target_cui"], rp["target_name"]
    rp["max_psamples"] = rp["ind_id"].map(max_ps)
    rp = rp.join(evidence, on="ind_id")
    rp["match_method"] = "rollup_" + rp["rollup_level"]
    return rp


# --------------------------------------------------------------------------- #
# Drug leg (04_) on the recovered pairs
# --------------------------------------------------------------------------- #
def drug_leg(rp_keep, exact_lp, dm, names_by_pair):
    """
    04_ filter on the kept recovered pairs of one threshold, then the roll-up
    specific removals against the exact pairs (conflict pairs included):
      conflict_with_exact        -- the exact pair on the target is a conflict;
      label_disagrees_with_exact -- the exact pair has the other label;
      redundant_with_exact       -- the exact pair has the same label;
      label_conflict             -- two children give one perturbagen set
                                    different labels on the same target;
      duplicate_on_target        -- same perturbagen set, target and label.
    The pair-level name check uses the CHILD indication's repoDB names (the
    records the label comes from).
    """
    empty = pd.DataFrame(columns=["removed_because"])
    if rp_keep.empty:
        return rp_keep.assign(lincs_pert_ids=""), empty
    tmp = rp_keep.assign(ind_id=rp_keep["child_ind_id"])
    lp, mism, dups = S04.filter_pairs(tmp, dm, names_by_pair)
    lp["ind_id"] = lp["target_cui"]
    removed = [mism.assign(removed_because="name_" + mism["name_check"]),
               dups.assign(removed_because="duplicate_repodb_id")]

    exact_label = (exact_lp.groupby(["lincs_pert_ids", "ind_id"])["label"]
                           .agg(lambda v: "conflict" if "conflict" in set(v) or len(set(v)) > 1
                                else next(iter(v))).to_dict())
    ex = [exact_label.get((p, d)) for p, d in zip(lp["lincs_pert_ids"], lp["ind_id"])]
    lp = lp.assign(exact_label=ex)
    why = lp.apply(lambda r: None if pd.isna(r["exact_label"]) else
                   ("conflict_with_exact" if r["exact_label"] == "conflict" else
                    "redundant_with_exact" if r["exact_label"] == r["label"] else
                    "label_disagrees_with_exact"), axis=1)
    hit = why.notna()
    removed.append(lp[hit].assign(removed_because=why[hit]))
    lp = lp[~hit].drop(columns="exact_label")

    labs = lp.groupby(["lincs_pert_ids", "ind_id"])["label"].agg(set)
    bad = {k for k, v in labs.items() if len(v) > 1}
    conf = pd.Series([(p, d) in bad for p, d in zip(lp["lincs_pert_ids"], lp["ind_id"])],
                     index=lp.index)
    removed.append(lp[conf].assign(removed_because="label_conflict"))
    lp = lp[~conf]
    dup = lp.duplicated(["lincs_pert_ids", "ind_id"], keep="first")
    removed.append(lp[dup].assign(removed_because="duplicate_on_target"))
    lp = lp[~dup]
    return lp.reset_index(drop=True), pd.concat(removed, ignore_index=True)


# --------------------------------------------------------------------------- #
# Summary and figure
# --------------------------------------------------------------------------- #
def summarise(exact_lp, rp, lp, removed, version):
    """One row per threshold x; every column refers to that x."""
    rows = []
    for x in sorted(rp["min_psamples"].unique()):
        row = {"min_psamples": x, "hierarchy": ROLLUP_HIERARCHY,
               "hierarchy_version": version}
        e = exact_lp[(exact_lp["max_psamples"] >= x) & (exact_lp["label"] != "conflict")]
        r = rp[rp["min_psamples"] == x]
        k = lp[lp["min_psamples"] == x] if len(lp) else lp.reindex(
            columns=["label", "child_ind_name", "child_ind_id", "ind_id",
                     "own_signature_below_x", "min_psamples"])
        d = (removed[removed["min_psamples"] == x] if len(removed)
             else pd.DataFrame(columns=["removed_because"]))
        for lab in LABELS:
            row[f"exact_{lab}_lincs"] = int((e["label"] == lab).sum())
            for dec in DECISION_CLASSES:
                row[f"rollup_{dec}_{lab}_before"] = int(
                    ((r["decision"] == dec) & (r["label"] == lab)).sum())
            row[f"rollup_keep_{lab}_lincs"] = int((k["label"] == lab).sum())
            base = row[f"exact_{lab}_lincs"]
            row[f"rollup_keep_{lab}_gain_pct"] = round(
                100 * row[f"rollup_keep_{lab}_lincs"] / base, 1) if base else None
        for why in ("conflict_with_exact", "label_disagrees_with_exact",
                    "redundant_with_exact", "label_conflict"):
            row[f"removed_{why}"] = int((d["removed_because"] == why).sum())
        for lab in LABELS:
            kl = k[k["label"] == lab]
            row[f"rollup_keep_{lab}_no_own_signature"] = int((kl["own_signature_below_x"] == 0).sum())
            row[f"rollup_keep_{lab}_own_signature_below_x"] = int((kl["own_signature_below_x"] == 1).sum())
        top_child = k["child_ind_name"].value_counts()
        row["rollup_keep_largest_child"] = top_child.index[0] if len(top_child) else ""
        row["rollup_keep_largest_child_pairs"] = int(top_child.iloc[0]) if len(top_child) else 0
        for lab in LABELS:
            row[f"rollup_keep_{lab}_excl_largest_child"] = int(
                ((k["label"] == lab) & (k["child_ind_name"] != row["rollup_keep_largest_child"])).sum())
        row["rollup_keep_child_diseases_lincs"] = k["child_ind_id"].nunique()
        row["rollup_keep_target_diseases_lincs"] = k["ind_id"].nunique()
        exact_cells = set(zip(e["ind_id"], e["label"]))
        row["rollup_keep_new_disease_label_cells"] = len(
            set(zip(k["ind_id"], k["label"])) - exact_cells)
        row["rollup_keep_pairs_in_top3_targets"] = int(
            k["ind_id"].value_counts().head(3).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def by_target(lp, exact_lp):
    """Kept recovered pairs per target disease and threshold, next to the exact ones."""
    if lp.empty:
        return pd.DataFrame()
    out = []
    for x, k in lp.groupby("min_psamples"):
        e = exact_lp[(exact_lp["max_psamples"] >= x) & (exact_lp["label"] != "conflict")]
        for t, g in k.groupby("ind_id"):
            ee = e[e["ind_id"] == t]
            out.append({"min_psamples": x, "target_cui": t,
                        "target_name": g["ind_name"].iloc[0],
                        "children": "|".join(sorted(set(g["child_ind_name"]))),
                        "rollup_approved": int((g["label"] == "approved").sum()),
                        "rollup_not_approved": int((g["label"] == "not_approved").sum()),
                        "exact_approved": int((ee["label"] == "approved").sum()),
                        "exact_not_approved": int((ee["label"] == "not_approved").sum())})
    return pd.DataFrame(out).sort_values(["min_psamples", "rollup_approved"],
                                         ascending=[True, False])


def plot_panel(summary):
    s = summary[summary["min_psamples"] == SELECTED_PSAMPLES].iloc[0]
    cats = ["keep", "exclude_infection", "exclude_wrong", "exclude_level2"]
    names = ["kept\n(subtype)", "excluded\n(infection)",
             "excluded\n(wrong/generic)", "excluded\n(grandparent)"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    ax = axes[0]
    xs = range(len(cats))
    w = 0.38
    for i, (lab, col) in enumerate(((LABELS[0], C_APPROVED), (LABELS[1], C_NOT_APPROVED))):
        vals = [s[f"rollup_{c}_{lab}_before"] for c in cats]
        bars = ax.bar([x + (i - 0.5) * w for x in xs], vals, w, color=col,
                      label=lab.replace("_", " "))
        S04._annot(ax, bars)
    ax.set_xticks(list(xs), names)
    ax.set_ylabel("recovered repoDB pairs (before LINCS)")
    ax.set_title("A  Roll-up candidates by curated decision", loc="left")
    ax.legend(frameon=False)
    S04._style(ax)

    ax = axes[1]
    groups = ["exact CUI join", "kept roll-up\n(after LINCS)"]
    for i, (lab, col) in enumerate(((LABELS[0], C_APPROVED), (LABELS[1], C_NOT_APPROVED))):
        vals = [s[f"exact_{lab}_lincs"], s[f"rollup_keep_{lab}_lincs"]]
        bars = ax.bar([x + (i - 0.5) * w for x in range(2)], vals, w, color=col,
                      label=lab.replace("_", " "))
        S04._annot(ax, bars)
    ax.set_xticks(range(2), groups)
    ax.set_ylabel("scorable pairs (LINCS)")
    ax.set_title("B  Gain over the exact CUI join", loc="left")
    S04._style(ax)
    via = "NCIt/MeSH (UMLS)" if ROLLUP_HIERARCHY == "umls" else "MONDO"
    fig.suptitle(f"UMLS granularity roll-up through {via}, "
                 f"≥{SELECTED_PSAMPLES} Psamples", x=0.01, ha="left")
    fig.tight_layout()
    IMG_OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(PANEL_PNG, dpi=200)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    repodb = S03.load_repodb()
    pairs = S03.build_labelled_pairs(repodb)
    ind_names = (repodb.dropna(subset=["ind_id"])
                       .groupby("ind_id")["ind_name"].first().to_dict())
    # Every signature CUI is a potential roll-up target, also when repoDB has
    # no indication coded with it: extend the CUI -> name table before the join.
    names_all = dict(ind_names)
    for x in S03.PSAMPLES_THRESHOLDS:
        path = S03.signature_summary_path(x)
        if path.exists():
            s = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
            for cuis, nm in zip(s["umls_cuis"], s["disease"]):
                for c in filter(None, cuis.split("|")):
                    names_all.setdefault(c, nm)
    per_x, thresholds = S03.load_signature_diseases(names_all)
    sig_pairs, _ = S03.intersect(pairs, per_x, thresholds)
    max_ps = {}
    for x in thresholds:
        for cui in set(per_x[x][0]["ind_id"]):
            max_ps[cui] = max(max_ps.get(cui, 0), x)
    x0 = min(thresholds)
    evidence = S03._signature_evidence(per_x[x0][0], x0)
    ind_names = names_all

    if ROLLUP_HIERARCHY == "umls":
        hier = load_umls(set(pairs["ind_id"]) | set(max_ps))
    else:
        hier = load_mondo()

    lin = S04.load_lincs()
    names_by_id, names_by_pair = S04.load_repodb_names()

    cands, rps = [], []
    for x in thresholds:
        c = find_candidates(pairs, max_ps, x, ind_names, hier)
        cands.append(c)
        rps.append(rolled_pairs(pairs, c, evidence, max_ps).assign(min_psamples=x))
    cand, rp = pd.concat(cands, ignore_index=True), pd.concat(rps, ignore_index=True)

    # Drug leg: one drug map over exact + recovered pairs (04_ logic); the
    # exact pairs keep their conflict rows for the redundancy check.
    dm = S04.build_drug_map(pd.concat([sig_pairs, rp]), lin, names_by_id)
    exact_lp, _, _ = S04.filter_pairs(sig_pairs, dm, names_by_pair)
    lps, rems = [], []
    for x in thresholds:
        keep = rp[(rp["min_psamples"] == x) & (rp["decision"] == "keep")]
        e = exact_lp[exact_lp["max_psamples"] >= x]
        lp_x, rem_x = drug_leg(keep, e, dm, names_by_pair)
        lps.append(lp_x.assign(min_psamples=x))
        rems.append(rem_x.assign(min_psamples=x))
    lp, removed = pd.concat(lps, ignore_index=True), pd.concat(rems, ignore_index=True)
    summary = summarise(exact_lp, rp, lp, removed, hier[3])
    bt = by_target(lp, exact_lp)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cand.sort_values(["min_psamples", "decision", "n_approved_pairs"],
                     ascending=[True, True, False]).to_csv(CANDIDATES_TSV, sep="\t", index=False)
    lp.to_csv(PAIRS_OUT_TSV, sep="\t", index=False)
    removed.to_csv(AUDIT_TSV, sep="\t", index=False)
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    bt.to_csv(BY_TARGET_TSV, sep="\t", index=False)
    plot_panel(summary)

    cui2m = hier[0]
    print(f"{hier[3]}: {len(cui2m):,} UMLS CUIs with a hierarchy node")
    print(f"repoDB indication CUIs with a node: "
          f"{sum(c in cui2m for c in set(pairs.ind_id)):,}/{pairs.ind_id.nunique():,}; "
          f"signature CUIs with a node: {sum(c in cui2m for c in max_ps)}/{len(max_ps)}")
    print("candidates by threshold x source:")
    print(cand.groupby(["min_psamples", "target_source"]).size().to_string())
    print("candidates by threshold x level x decision:")
    print(cand.groupby(["min_psamples", "rollup_level", "decision"]).size().to_string())
    multi = cand.groupby("ind_id")["target_cui"].nunique()
    if (multi > 1).any():
        print("\nCUIs whose roll-up target changes with x (check the curated note):",
              sorted(multi[multi > 1].index))
    unc = cand[cand["decision"] == "uncurated"].drop_duplicates(["ind_id", "target_cui"])
    if len(unc):
        print(f"\n--- {len(unc)} level-0/1 candidates WITHOUT a curated decision "
              "(excluded) ---")
        print(unc[["ind_id", "ind_name", "rollup_level", "target_cui", "target_name",
                   "target_source", "target_node_label",
                   "n_approved_pairs", "n_not_approved_pairs"]].to_string(index=False))
    print("\nremoved in the drug leg (all x):",
          removed.groupby(["min_psamples", "removed_because"]).size().to_dict()
          if len(removed) else {})
    with pd.option_context("display.width", 200):
        print(summary.T.to_string(header=False))


if __name__ == "__main__":
    main()
