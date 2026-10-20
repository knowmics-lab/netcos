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
            through its `parents` in the same source. A CUI marks a node as a
            target only when it is the node's preferred term (UMLS_PREFERRED_TTY,
            e.g. MeSH MH/NM, NCIt PT): a MeSH entry term carries the code of
            the broader descriptor that contains it. A child reached through
            a non-preferred atom starts one level up (the containing node is
            level 1). Sources are tried in
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
Level-2 candidates are excluded by rule for the sources in
EXCLUDE_LEVEL2_SOURCES; level-2 candidates of the other sources (NCIt, whose
stage and marker subtypes sit two steps below their histological type) are
curated like level 0/1. Candidates without a
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
                  UMLS_CACHE_FILE, UMLS_PREFERRED_TTY, UMLS_SOURCES,
                  UMLS_VERSION)

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
EXCLUDE_LEVEL2_SOURCES = ("MSH", "MONDO")  # grandparent targets never used for these
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
# ((repoDB CUI, target CUI) -> (decision, note)), from manual LLM-assisted
# inspection of each child indication vs its target; comment: names, source, level,
# approved/not-approved pairs at the lowest threshold where the candidate occurs.
ROLLUP_DECISIONS_UMLS = {
    ("C0276535", "C0036220"): ("keep", "genuine subtype / clinical form of the target"),  # AIDS with Kaposi's sarcoma -> Kaposi Sarcoma [NCI p] A3N0
    ("C0014553", "C0014544"): ("exclude_wrong", "target is a different disease or too generic"),  # Absence Epilepsy -> Epilepsy [NCI g] A10N0
    ("C0000880", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Acanthamoeba Keratitis -> Protozoal Infection [NCI p] A1N0
    ("C0001144", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Acne Vulgaris -> Dermatitis [NCI p] A34N1
    ("C0271693", "C0023787"): ("keep", "NCIt level 2"),  # Acquired generalized lipodystrophy -> Lipodystrophy [NCI g] A1N0
    ("C0220989", "C0023787"): ("keep", "NCIt level 2"),  # Acquired partial lipodystrophy -> Lipodystrophy [NCI g] A1N0
    ("C0001261", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Actinomycosis -> Bacterial Infection [NCI g] A10N0
    ("C0023440", "C0023467"): ("keep", "genuine subtype / clinical form of the target"),  # Acute Erythroblastic Leukemia -> Acute Myeloid Leukemia [NCI p] A1N0
    ("C0023487", "C0023467"): ("keep", "APL is a distinct AML subtype (ATRA/arsenic); AML signature a weak proxy"),  # Acute Promyelocytic Leukemia -> Acute Myeloid Leukemia [NCI g] A7N0
    ("C0023465", "C0023467"): ("keep", "genuine subtype / clinical form of the target"),  # Acute monocytic leukemia -> Acute Myeloid Leukemia [NCI g] A1N0
    ("C0023479", "C0023467"): ("keep", "genuine subtype / clinical form of the target"),  # Acute myelomonocytic leukemia -> Acute Myeloid Leukemia [NCI p] A1N0
    ("C0281361", "C0001418"): ("exclude_wrong", "target is a different disease or too generic"),  # Adenocarcinoma of pancreas -> Adenocarcinoma [NCI g] A4N2
    ("C0023493", "C0023492"): ("keep", "genuine subtype / clinical form of the target"),  # Adult T-Cell Lymphoma/Leukemia -> Leukemia, T-Cell [MSH p] A1N0
    ("C0023493", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # Adult T-Cell Lymphoma/Leukemia -> Lymphoma [NCI g] A1N0
    ("C0041228", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # African Trypanosomiasis -> Protozoal Infection [NCI p] A4N0
    ("C0085280", "C0008372"): ("exclude_wrong", "target is a different disease or too generic"),  # Alagille Syndrome -> Cholestasis, Intrahepatic [MSH p] A1N0
    ("C0002351", "C0035204"): ("exclude_wrong", "target is a different disease or too generic"),  # Altitude Sickness -> Respiration Disorders [MSH p] A1N0
    ("C0002438", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Amebiasis -> Protozoal Infection [NCI p] A5N0
    ("C0013370", "C0009319"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Amebic colitis -> Colitis [NCI g] A2N0
    ("C0206247", "C0002726"): ("keep", "genuine subtype / clinical form of the target"),  # Amyloid Neuropathies -> Amyloidosis [MSH p] A2N0
    ("C0334579", "C0004114"): ("keep", "genuine subtype / clinical form of the target"),  # Anaplastic astrocytoma -> Astrocytic Tumor [NCI p] A1N0
    ("C0158996", "C0002871"): ("exclude_wrong", "target is a different disease or too generic"),  # Anemia of prematurity -> Anemia [NCI p] A5N0
    ("C0002888", "C0002871"): ("exclude_wrong", "target is a different disease or too generic"),  # Anemia, Megaloblastic -> Anemia [NCI g] A2N0
    ("C0002962", "C0151744"): ("keep", "genuine subtype / clinical form of the target"),  # Angina Pectoris -> Ischemic Heart Disease [MSH p] A24N1
    ("C0002963", "C0002965"): ("exclude_wrong", "target is a different disease or too generic"),  # Angina Pectoris, Variant -> Angina, Unstable [MSH p] A5N0
    ("C0241961", "C0022665"): ("exclude_wrong", "target is a different disease or too generic"),  # Angiomyolipoma of kidney -> Kidney Neoplasm [NCI p] A1N0
    ("C0241961", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Angiomyolipoma of kidney -> Kidney Diseases [NCI g] A1N0
    ("C0017605", "C0017601"): ("keep", "genuine subtype / clinical form of the target"),  # Angle Closure Glaucoma -> Glaucoma [NCI p] A4N0
    ("C0038013", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Ankylosing spondylitis -> Arthritis [NCI g] A27N0
    ("C0042165", "C0042164"): ("keep", "genuine subtype / clinical form of the target"),  # Anterior uveitis -> Uveitis [NCI p] A4N0
    ("C0003175", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Anthrax disease -> Bacterial Infection [NCI p] A13N0
    ("C0003872", "C0033860"): ("exclude_wrong", "target is a different disease or too generic"),  # Arthritis, Psoriatic -> Psoriasis [MSH p] A19N0
    ("C0003872", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Arthritis, Psoriatic -> Arthritis [NCI p] A19N0
    ("C0004031", "C0004030"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Aspergillosis, Allergic Bronchopulmonary -> Aspergillosis [NCI p] A1N0
    ("C4038730", "C0004096"): ("keep", "genuine subtype / clinical form of the target"),  # Asthma-Chronic Obstructive Pulmonary Disea -> Asthma [MSH p] A3N0
    ("C0270846", "C0036572"): ("exclude_wrong", "target is a different disease or too generic"),  # Atonic seizure -> Seizure [MSH p] A2N0
    ("C0004352", "C1510586"): ("keep", "genuine subtype / clinical form of the target"),  # Autistic Disorder -> Autism Spectrum Disorder [NCI p] A2N0
    ("C0241910", "C0019189"): ("exclude_wrong", "target is a different disease or too generic"),  # Autoimmune Chronic Hepatitis -> Chronic Hepatitis [MSH p] A0N4
    ("C0079731", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # B-Cell Lymphomas -> Lymphoma [NCI g] A1N0
    ("C0085166", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Bacterial Vaginosis -> Bacterial Infection [MSH p] A3N1
    ("C0004771", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Bartonella Infections -> Bacterial Infection [NCI p] A3N0
    ("C0005138", "C0206062"): ("exclude_wrong", "target is a different disease or too generic"),  # Berylliosis -> Interstitial Lung Disease [NCI g] A10N1
    ("C0005699", "C0007621"): ("exclude_wrong", "target is a different disease or too generic"),  # Blast Phase -> Cell Transformation, Neoplastic [MSH p] A3N0
    ("C1301363", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # Blastic plasmacytoid dendritic cell neopla -> Lymphoma [MSH p] A1N0
    ("C0006057", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Botulism -> Bacterial Infection [NCI p] A1N0
    ("C0349543", "C1621958"): ("keep", "genuine subtype / clinical form of the target"),  # Brain Glioblastoma -> Glioblastoma [NCI p] A2N0
    ("C0006267", "C0024117"): ("exclude_wrong", "target is a different disease or too generic"),  # Bronchiectasis -> Chronic Obstructive Pulmonary Disease [NCI p] A6N2
    ("C0006277", "C0024115"): ("exclude_wrong", "target is a different disease or too generic"),  # Bronchitis -> Lung Diseases [NCI g] A30N0
    ("C0008677", "C0024117"): ("keep", "genuine subtype / clinical form of the target"),  # Bronchitis, Chronic -> Chronic Obstructive Pulmonary Disease [NCI p] A9N1
    ("C0006266", "C0024115"): ("exclude_wrong", "target is a different disease or too generic"),  # Bronchospasm -> Lung Diseases [NCI g] A6N0
    ("C0006413", "C0149678"): ("exclude_wrong", "MeSH places Burkitt under EBV infections; target is an infection"),  # Burkitt Lymphoma -> EBV Infection [MSH p] A4N4
    ("C5203670", "C0042769"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # COVID19 (disease) -> Viral Infection [NCI g] A11N51
    ("C0153254", "C1541923"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Candidal endocarditis -> Infective Endocarditis [NCI p] A1N0
    ("C0700345", "C0006840"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Candidiasis, Vulvovaginal -> Candidiasis [NCI p] A11N0
    ("C0007138", "C0334265"): ("keep", "genuine subtype / clinical form of the target"),  # Carcinoma, Transitional Cell -> Transitional Cell Neoplasm [NCI p] A3N3
    ("C0003811", "C0007222"): ("exclude_wrong", "target is a different disease or too generic"),  # Cardiac Arrhythmia -> Cardiovascular Diseases [NCI g] A1N0
    ("C0302592", "C0007873"): ("keep", "genuine subtype / clinical form of the target"),  # Cervix carcinoma -> Cervical Neoplasm [NCI g] A1N0
    ("C0008149", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Chlamydia Infections -> Bacterial Infection [NCI p] A2N0
    ("C0339959", "C0004626"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Chlamydial pneumonia -> Bacterial Pneumonia [MSH p] A13N0
    ("C0206698", "C0001418"): ("exclude_wrong", "target is a different disease or too generic"),  # Cholangiocarcinoma -> Adenocarcinoma [MSH p] A1N0
    ("C0008513", "C0042164"): ("exclude_wrong", "target is a different disease or too generic"),  # Chorioretinitis -> Uveitis [NCI p] A1N0
    ("C0008526", "C0042164"): ("exclude_wrong", "target is a different disease or too generic"),  # Choroiditis -> Uveitis [NCI g] A2N0
    ("C0409818", "C2316212"): ("keep", "genuine subtype / clinical form of the target"),  # Chronic Infantile Neurological, Cutaneous, -> Cryopyrin-Associated Periodic Syndrome [MSH p] A1N0
    ("C1561643", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Chronic Kidney Diseases -> Kidney Diseases [NCI g] A1N9
    ("C0346421", "C0023418"): ("exclude_wrong", "target is a different disease or too generic"),  # Chronic eosinophilic leukemia -> Leukemia [NCI g] A1N1
    ("C0085695", "C0017152"): ("keep", "genuine subtype / clinical form of the target"),  # Chronic gastritis -> Gastritis [NCI p] A1N0
    ("C0600040", "C0282488"): ("keep", "genuine subtype / clinical form of the target"),  # Chronic interstitial cystitis -> Cystitis, Interstitial [MSH p] A1N0
    ("C0221765", "C0036341"): ("keep", "genuine subtype / clinical form of the target"),  # Chronic schizophrenia -> Schizophrenia [NCI p] A3N1
    ("C2363973", "C0020542"): ("keep", "genuine subtype / clinical form of the target"),  # Chronic thromboembolic pulmonary hypertens -> Hypertension, Pulmonary [NCI p] A2N0
    ("C1333064", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # Classical Hodgkin's Lymphoma -> Lymphoma [NCI g] A2N0
    ("C0343386", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Clostridium difficile infection -> Bacterial Infection [NCI p] A1N3
    ("C0342257", "C0011849"): ("exclude_wrong", "target is a different disease or too generic"),  # Complications of Diabetes Mellitus -> Diabetes Mellitus [MSH p] A1N0
    ("C0010308", "C0020676"): ("keep", "genuine subtype / clinical form of the target"),  # Congenital Hypothyroidism -> Hypothyroidism [NCI p] A2N0
    ("C0018802", "C0018801"): ("keep", "genuine subtype / clinical form of the target"),  # Congestive heart failure -> Heart Failure [NCI p] A6N8
    ("C0011616", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Contact Dermatitis -> Dermatitis [NCI g] A33N0
    ("C0011616", "C0013595"): ("keep", "genuine subtype / clinical form of the target"),  # Contact Dermatitis -> Eczema [NCI p] A33N0
    ("C0010200", "C0035204"): ("exclude_wrong", "target is a different disease or too generic"),  # Coughing -> Respiration Disorders [MSH p] A34N7
    ("C0010418", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Cryptosporidiosis -> Protozoal Infection [NCI p] A1N0
    ("C0948245", "C0242966"): ("exclude_wrong", "target is a different disease or too generic"),  # Cytokine Release Syndrome -> Systemic Inflammatory Response Syndrome [MSH p] A1N0
    ("C0029408", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Degenerative polyarthritis -> Arthritis [NCI p] A36N7
    ("C0005587", "C0005586"): ("keep", "genuine subtype / clinical form of the target"),  # Depression, Bipolar -> Bipolar Disorder [NCI p] A5N0
    ("C0221074", "C0011581"): ("keep", "genuine subtype / clinical form of the target"),  # Depression, Postpartum -> Depressive Disorder [MSH p] A1N2
    ("C0011608", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Dermatitis Herpetiformis -> Dermatitis [NCI p] A12N0
    ("C0011636", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Dermatophytosis -> Dermatitis [NCI p] A1N0
    ("C0011881", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Diabetic Nephropathy -> Kidney Diseases [MSH p] A3N0
    ("C0011974", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Diaper Rash -> Dermatitis [NCI p] A11N0
    ("C0855112", "C0079744"): ("keep", "genuine subtype / clinical form of the target"),  # Diffuse large B-cell lymphoma refractory -> Diffuse Large B-Cell Lymphoma [NCI p] A2N1
    ("C0012546", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Diphtheria -> Bacterial Infection [NCI p] A3N0
    ("C0013404", "C0035204"): ("exclude_wrong", "target is a different disease or too generic"),  # Dyspnea -> Respiration Disorders [MSH p] A1N5
    ("C0013595", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Eczema -> Dermatitis [NCI p] A2N0
    ("C0014013", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Empyema, Pleural -> Bacterial Infection [NCI g] A2N0
    ("C0041472", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Endemic Flea-Borne Typhus -> Bacterial Infection [NCI g] A2N0
    ("C0476089", "C0014170"): ("keep", "genuine subtype / clinical form of the target"),  # Endometrial Carcinoma -> Endometrial Neoplasm [NCI g] A3N1
    ("C0275982", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Enteric campylobacteriosis -> Bacterial Infection [NCI p] A3N0
    ("C0014488", "C1568272"): ("keep", "genuine subtype / clinical form of the target"),  # Epicondylitis -> Tendinopathy [MSH p] A7N0
    ("C0014547", "C0014544"): ("keep", "genuine subtype / clinical form of the target"),  # Epilepsies, Partial -> Epilepsy [NCI p] A5N0
    ("C0085580", "C0020538"): ("keep", "genuine subtype / clinical form of the target"),  # Essential Hypertension -> Hypertension [NCI p] A4N0
    ("C0011606", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Exfoliative dermatitis -> Dermatitis [MSH p] A10N0
    ("C1827301", "C0041296"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Extensively Drug-Resistant Tuberculosis -> Tuberculosis [NCI p] A1N0
    ("C0002390", "C0206062"): ("keep", "genuine subtype / clinical form of the target"),  # Extrinsic allergic alveolitis -> Interstitial Lung Disease [MSH p] A1N0
    ("C0221032", "C0023787"): ("keep", "genuine subtype / clinical form of the target"),  # Familial generalized lipodystrophy -> Lipodystrophy [NCI p] A1N0
    ("C1318543", "C1568272"): ("exclude_wrong", "target is a different disease or too generic"),  # Fibrous histiocytoma of tendon sheath -> Tendinopathy [MSH p] A1N0
    ("C1318543", "C0017525"): ("exclude_wrong", "target is a different disease or too generic"),  # Fibrous histiocytoma of tendon sheath -> Giant Cell Neoplasm [NCI p] A1N0
    ("C0206682", "C0001418"): ("exclude_wrong", "target is a different disease or too generic"),  # Follicular thyroid carcinoma -> Adenocarcinoma [MSH p] A1N0
    ("C0278701", "C0699791"): ("keep", "genuine subtype / clinical form of the target"),  # Gastric Adenocarcinoma -> Gastric Carcinoma [NCI p] A1N0
    ("C1961835", "C0017205"): ("keep", "genuine subtype / clinical form of the target"),  # Gaucher Disease, Type 1 -> Gaucher Disease [MSH p] A3N0
    ("C0343055", "C0033860"): ("keep", "genuine subtype / clinical form of the target"),  # Generalized pustular psoriasis -> Psoriasis [NCI p] A1N0
    ("C0017536", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Giardiasis -> Protozoal Infection [MSH p] A6N0
    ("C0017609", "C0017601"): ("exclude_wrong", "target is a different disease or too generic"),  # Glaucoma, Neovascular -> Glaucoma [MSH p] A1N1
    ("C0017612", "C0017601"): ("keep", "genuine subtype / clinical form of the target"),  # Glaucoma, Open-Angle -> Glaucoma [NCI p] A28N2
    ("C1621958", "C0004114"): ("exclude_wrong", "GBM rolled onto broader astrocytic tumour class"),  # Glioblastoma Multiforme -> Astrocytic Tumor [NCI g] A0N24
    ("C0018081", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Gonorrhea -> Bacterial Infection [NCI p] A23N0
    ("C0018099", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Gout -> Arthritis [NCI p] A9N0
    ("C0018190", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Granuloma Inguinale -> Bacterial Infection [NCI p] A10N0
    ("C3495801", "C0206062"): ("exclude_wrong", "target is a different disease or too generic"),  # Granulomatosis with polyangiitis -> Interstitial Lung Disease [MSH p] A2N0
    ("C1960398", "C0858252"): ("keep", "NCIt level 2"),  # HER2-positive carcinoma of breast -> Breast Adenocarcinoma [NCI g] A5N0
    ("C0023443", "C0023418"): ("exclude_wrong", "target is a different disease or too generic"),  # Hairy Cell Leukemia -> Leukemia [NCI g] A5N0
    ("C0018621", "C2607914"): ("keep", "genuine subtype / clinical form of the target"),  # Hay fever -> Allergic Rhinitis [NCI p] A35N0
    ("C0018671", "C0027651"): ("exclude_wrong", "target is a different disease or too generic"),  # Head and Neck Neoplasms -> Neoplasm [NCI g] A1N0
    ("C0238425", "C0002895"): ("keep", "genuine subtype / clinical form of the target"),  # Hemoglobin SS disease with crisis -> Anemia, Sickle Cell [NCI p] A1N5
    ("C0019147", "C0019151"): ("keep", "genuine subtype / clinical form of the target"),  # Hepatic Coma -> Hepatic Encephalopathy [MSH p] A3N0
    ("C0019156", "C0023895"): ("exclude_wrong", "target is a different disease or too generic"),  # Hepatic Veno-Occlusive Disease -> Liver Diseases [NCI g] A1N0
    ("C0524910", "C0019196"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Hepatitis C, Chronic -> Hepatitis C [NCI p] A23N5
    ("C0019202", "C0023895"): ("exclude_wrong", "target is a different disease or too generic"),  # Hepatolenticular Degeneration -> Liver Diseases [MSH p] A3N0
    ("C0019212", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Hepatorenal Syndrome -> Kidney Diseases [MSH p] A0N2
    ("C0019348", "C0042769"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Herpes simplex infection -> Viral Infection [NCI g] A1N1
    ("C0019829", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # Hodgkin Disease -> Lymphoma [NCI p] A22N3
    ("C0342881", "C0020445"): ("keep", "genuine subtype / clinical form of the target"),  # Homozygous Familial Hypercholesterolemia -> Hyperlipoproteinemia Type II [MSH p] A7N0
    ("C1328504", "C0600139"): ("keep", "genuine subtype / clinical form of the target"),  # Hormone refractory prostate cancer -> Prostate Carcinoma [NCI g] A0N2
    ("C0949083", "C0032285"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Hospital acquired pneumonia -> Pneumonia [NCI p] A8N5
    ("C0020473", "C0242339"): ("keep", "genuine subtype / clinical form of the target"),  # Hyperlipidemia -> Dyslipidemia [MSH p] A14N9
    ("C0178416", "C0002874"): ("keep", "genuine subtype / clinical form of the target"),  # Hypoplastic anemia -> Anemia, Aplastic [MSH p] A1N0
    ("C0021099", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Impetigo -> Bacterial Infection [NCI p] A3N0
    ("C0000810", "C0000786"): ("keep", "subtype; drugs act on uterus"),  # Incomplete spontaneous abortion -> Abortion, Spontaneous [MSH p] A4N0
    ("C0751122", "C0014544"): ("exclude_wrong", "target is a different disease or too generic"),  # Infantile Severe Myoclonic Epilepsy -> Epilepsy [NCI p] A2N0
    ("C0850666", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Infection caused by Helicobacter pylori -> Bacterial Infection [NCI p] A3N0
    ("C0345905", "C0206698"): ("keep", "genuine subtype / clinical form of the target"),  # Intrahepatic Cholangiocarcinoma -> Cholangiocarcinoma [NCI p] A1N2
    ("C0276653", "C0004030"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Invasive Pulmonary Aspergillosis -> Aspergillosis [NCI p] A1N1
    ("C0022073", "C0042164"): ("keep", "NCIt level 2"),  # Iridocyclitis -> Uveitis [NCI g] A12N0
    ("C0162316", "C0002871"): ("exclude_wrong", "target is a different disease or too generic"),  # Iron deficiency anemia -> Anemia [NCI g] A15N0
    ("C3495559", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Juvenile arthritis -> Arthritis [NCI p] A21N3
    ("C0022575", "C0013238"): ("keep", "genuine subtype / clinical form of the target"),  # Keratoconjunctivitis Sicca -> Dry Eye Syndrome [MSH p] A1N0
    ("C0022972", "C0026896"): ("exclude_wrong", "target is a different disease or too generic"),  # Lambert-Eaton Myasthenic Syndrome -> Myasthenia Gravis [MSH p] A2N0
    ("C0280324", "C1168401"): ("keep", "genuine subtype / clinical form of the target"),  # Laryngeal Squamous Cell Carcinoma -> Head and Neck Squamous Cell Carcinoma [NCI g] A1N0
    ("C0022340", "C0027877"): ("keep", "genuine subtype / clinical form of the target"),  # Late-Infantile Neuronal Ceroid Lipfuscinos -> Neuronal Ceroid Lipofuscinosis [MSH p] A1N0
    ("C0238111", "C0014544"): ("exclude_wrong", "target is a different disease or too generic"),  # Lennox-Gastaut syndrome -> Epilepsy [NCI p] A11N0
    ("C0023381", "C0019621"): ("keep", "genuine subtype / clinical form of the target"),  # Letterer-Siwe Disease -> Histiocytosis, Langerhans-Cell [NCI g] A1N0
    ("C0023461", "C0023467"): ("exclude_wrong", "target is a different disease or too generic"),  # Leukemia, Mast-Cell -> Acute Myeloid Leukemia [MSH p] A1N1
    ("C0023474", "C0023418"): ("exclude_wrong", "target is a different disease or too generic"),  # Leukemia, Myeloid, Chronic-Phase -> Leukemia [NCI g] A3N0
    ("C0023480", "C0023418"): ("exclude_wrong", "target is a different disease or too generic"),  # Leukemia, Myelomonocytic, Chronic -> Leukemia [NCI g] A1N13
    ("C0752347", "C0497327"): ("exclude_wrong", "target is a different disease or too generic"),  # Lewy Body Disease -> Dementia [NCI p] A1N0
    ("C0023646", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Lichen Planus -> Dermatitis [NCI p] A1N2
    ("C0348868", "C1541923"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Listerial endocarditis -> Infective Endocarditis [NCI g] A1N0
    ("C0023885", "C0023895"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Liver Abscess -> Liver Diseases [MSH p] A4N0
    ("C0085605", "C0023895"): ("exclude_wrong", "target is a different disease or too generic"),  # Liver Failure -> Liver Diseases [NCI g] A1N0
    ("C0024143", "C0024141"): ("keep", "genuine subtype / clinical form of the target"),  # Lupus Nephritis -> Lupus Erythematosus, Systemic [MSH p] A1N13
    ("C0024286", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Lymphogranuloma Venereum -> Bacterial Infection [NCI g] A9N0
    ("C0024305", "C0024299"): ("exclude_wrong", "target is a different disease or too generic"),  # Lymphoma, Non-Hodgkin -> Lymphoma [NCI p] A1N0
    ("C0079773", "C0079772"): ("keep", "single target dataset (GSE45712)"),  # Lymphoma, T-Cell, Cutaneous -> Lymphoma, T-Cell [NCI g] A34N0
    ("C0555198", "C0017638"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant Glioma -> Glioma [NCI p] A1N0
    ("C0278996", "C0018671"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant Head and Neck Neoplasm -> Head and Neck Neoplasm [NCI p] A5N0
    ("C0024588", "C0020538"): ("keep", "NCIt level 2"),  # Malignant essential hypertension -> Hypertension [NCI g] A2N0
    ("C0153633", "C0006118"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant neoplasm of brain -> Brain Neoplasm [NCI p] A1N0
    ("C0220636", "C0018671"): ("exclude_wrong", "target is a different disease or too generic"),  # Malignant neoplasm of salivary gland -> Head and Neck Neoplasm [NCI g] A1N0
    ("C0024623", "C0038356"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant neoplasm of stomach -> Gastric Neoplasm [NCI p] A10N0
    ("C0007115", "C0040136"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant neoplasm of thyroid -> Thyroid Gland Neoplasm [NCI p] A3N0
    ("C0007847", "C0007873"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant tumor of cervix -> Cervical Neoplasm [NCI p] A2N0
    ("C0007102", "C0009375"): ("keep", "genuine subtype / clinical form of the target"),  # Malignant tumor of colon -> Colon Neoplasm [NCI p] A5N0
    ("C0238462", "C0549473"): ("exclude_wrong", "medullary (C-cell) carcinoma differs from follicular-derived thyroid carcinoma"),  # Medullary carcinoma of thyroid -> Thyroid Gland Carcinoma [NCI p] A3N2
    ("C0085436", "C0025289"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Meningitis, Cryptococcal -> Meningitis [NCI g] A1N0
    ("C0085438", "C0025289"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Meningitis, Fungal -> Meningitis [MSH p] A1N0
    ("C0025293", "C0085437"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Meningitis, Listeria -> Bacterial Meningitis [MSH p] A1N0
    ("C0025295", "C0085437"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Meningitis, Pneumococcal -> Bacterial Meningitis [MSH p] A4N0
    ("C0025295", "C0042769"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Meningitis, Pneumococcal -> Viral Infection [NCI g] A4N0
    ("C0346974", "C0009375"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to colon -> Colon Neoplasm [NCI g] A2N0
    ("C0346974", "C0007102"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to colon -> Malignant Colon Neoplasm [NCI p] A2N0
    ("C0153685", "C0022665"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to kidney -> Kidney Neoplasm [NCI g] A2N0
    ("C0346976", "C0030297"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to pancreas -> Pancreatic Neoplasm [NCI g] A2N0
    ("C0686068", "C0038356"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to stomach -> Gastric Neoplasm [NCI g] A2N0
    ("C0686068", "C0024623"): ("exclude_wrong", "target is a different disease or too generic"),  # Metastatic malignant neoplasm to stomach -> Malignant Gastric Neoplasm [NCI p] A2N0
    ("C0278987", "C0007131"): ("keep", "NCIt level 2: stage subtype of NSCLC"),  # Metastatic non-small cell lung cancer -> Carcinoma, Non-Small-Cell Lung [NCI g] A5N3
    ("C4288754", "C2145472"): ("keep", "genuine subtype / clinical form of the target"),  # Metastatic urothelial carcinoma -> Urothelial Carcinoma [NCI p] A2N7
    ("C2347126", "C0042384"): ("exclude_wrong", "target is a different disease or too generic"),  # Microscopic Polyarteritis -> Vasculitis [NCI g] A2N3
    ("C0276180", "C0042769"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Monkeypox disease -> Viral Infection [NCI p] A1N0
    ("C1527383", "C0036420"): ("keep", "genuine subtype / clinical form of the target"),  # Morphea -> Localized Scleroderma [MSH p] A2N1
    ("C0006848", "C0006840"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Mucocutaneous candidiasis -> Candidiasis [NCI p] A2N0
    ("C0027662", "C0014132"): ("exclude_wrong", "target is a different disease or too generic"),  # Multiple Endocrine Neoplasia -> Endocrine Gland Neoplasms [NCI p] A3N0
    ("C0393664", "C0751967"): ("keep", "genuine subtype / clinical form of the target"),  # Multiple Sclerosis, Acute Relapsing -> Multiple Sclerosis, Relapsing-Remitting [MSH p] A2N3
    ("C0751967", "C0026769"): ("keep", "genuine subtype / clinical form of the target"),  # Multiple Sclerosis, Relapsing-Remitting -> Multiple Sclerosis [NCI p] A11N3
    ("C0026838", "C0026848"): ("exclude_wrong", "target is a different disease or too generic"),  # Muscle Spasticity -> Muscular Diseases [MSH p] A14N1
    ("C0231528", "C0026848"): ("exclude_wrong", "target is a different disease or too generic"),  # Myalgia -> Muscular Diseases [MSH p] A1N0
    ("C0027145", "C0020676"): ("keep", "genuine subtype / clinical form of the target"),  # Myxedema -> Hypothyroidism [NCI p] A2N0
    ("C0751362", "C0027404"): ("keep", "genuine subtype / clinical form of the target"),  # Narcolepsy-Cataplexy Syndrome -> Narcolepsy [MSH p] A3N0
    ("C0027726", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Nephrotic Syndrome -> Kidney Diseases [NCI g] A11N2
    ("C0206754", "C0014132"): ("exclude_wrong", "target is a different disease or too generic"),  # Neuroendocrine Tumors -> Endocrine Gland Neoplasms [NCI p] A1N7
    ("C0027873", "C0029134"): ("exclude_wrong", "target is a different disease or too generic"),  # Neuromyelitis Optica -> Optic Neuritis [MSH p] A3N4
    ("C0268242", "C0028064"): ("keep", "genuine subtype / clinical form of the target"),  # Niemann-Pick Disease, Type A -> Niemann-Pick Disease [NCI p] A1N0
    ("C0268243", "C0028064"): ("keep", "genuine subtype / clinical form of the target"),  # Niemann-Pick Disease, Type B -> Niemann-Pick Disease [NCI p] A1N0
    ("C0220756", "C0028064"): ("exclude_wrong", "NPC is a distinct disease (NPC1/2) from Niemann-Pick A/B"),  # Niemann-Pick Disease, Type C -> Niemann-Pick Disease [NCI p] A1N3
    ("C0028242", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Nocardia Infections -> Bacterial Infection [NCI p] A2N0
    ("C0028840", "C0020538"): ("exclude_wrong", "ocular hypertension is not systemic hypertension"),  # Ocular Hypertension -> Hypertension [NCI p] A16N4
    ("C0406486", "C0035854"): ("keep", "genuine subtype / clinical form of the target"),  # Ocular Rosacea -> Rosacea [MSH p] A3N0
    ("C0917796", "C0029125"): ("keep", "genuine subtype / clinical form of the target"),  # Optic Atrophy, Hereditary, Leber -> Hereditary Optic Atrophy [NCI p] A1N0
    ("C0585442", "C0029463"): ("keep", "genuine subtype / clinical form of the target"),  # Osteosarcoma of bone -> Osteosarcoma [NCI p] A3N0
    ("C0030528", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Paratyphoid Fever -> Bacterial Infection [NCI p] A2N0
    ("C0235480", "C0004238"): ("keep", "genuine subtype / clinical form of the target"),  # Paroxysmal atrial fibrillation -> Atrial Fibrillation [NCI p] A3N0
    ("C0221154", "C0020538"): ("exclude_wrong", "target is a different disease or too generic"),  # Paroxysmal hypertension -> Hypertension [NCI p] A1N0
    ("C0013274", "C0007222"): ("exclude_wrong", "target is a different disease or too generic"),  # Patent ductus arteriosus -> Cardiovascular Diseases [NCI g] A0N1
    ("C0030809", "C0030807"): ("keep", "genuine subtype / clinical form of the target"),  # Pemphigus Vulgaris -> Pemphigus [NCI p] A2N2
    ("C0014869", "C0014868"): ("keep", "genuine subtype / clinical form of the target"),  # Peptic Esophagitis -> Esophagitis [MSH p] A5N0
    ("C0031190", "C0020542"): ("exclude_wrong", "target is a different disease or too generic"),  # Persistent Fetal Circulation Syndrome -> Hypertension, Pulmonary [MSH p] A3N2
    ("C0279543", "C0027022"): ("exclude_wrong", "target is a different disease or too generic"),  # Philadelphia chromosome positive chronic m -> Myeloproliferative Disorders [MSH p] A4N0
    ("C0031946", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Pinta -> Bacterial Infection [NCI g] A2N0
    ("C0032019", "C0006118"): ("exclude_wrong", "target is a different disease or too generic"),  # Pituitary Neoplasms -> Brain Neoplasm [NCI g] A1N0
    ("C0032064", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Plague -> Bacterial Infection [NCI p] A2N0
    ("C0264515", "C0032285"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Pneumonia, Necrotizing -> Pneumonia [MSH p] A1N0
    ("C0031036", "C0042384"): ("exclude_wrong", "target is a different disease or too generic"),  # Polyarteritis Nodosa -> Vasculitis [NCI g] A1N0
    ("C0032463", "C0032461"): ("exclude_wrong", "target signature is Chuvash polycythemia (congenital VHL), not clonal PV"),  # Polycythemia Vera -> Polycythemia [NCI p] A5N9
    ("C0032533", "C0026848"): ("exclude_wrong", "target is a different disease or too generic"),  # Polymyalgia Rheumatica -> Muscular Diseases [MSH p] A1N2
    ("C0162533", "C0023895"): ("exclude_wrong", "target is a different disease or too generic"),  # Porphyrias, Hepatic -> Liver Diseases [MSH p] A1N0
    ("C3805233", "C0026987"): ("keep", "genuine subtype / clinical form of the target"),  # Post-essential thrombocythemia myelofibros -> Myelofibrosis [NCI g] A1N0
    ("C3805232", "C0032461"): ("exclude_wrong", "target is a different disease or too generic"),  # Post-polycythemia vera myelofibrosis -> Polycythemia [NCI g] A1N0
    ("C3805232", "C0032463"): ("exclude_wrong", "target is a different disease or too generic"),  # Post-polycythemia vera myelofibrosis -> Polycythemia Vera [NCI p] A1N0
    ("C0032897", "C0028754"): ("exclude_wrong", "target is a different disease or too generic"),  # Prader-Willi Syndrome -> Obese [MSH p] A2N5
    ("C1292758", "C0079772"): ("exclude_wrong", "target is a different disease or too generic"),  # Precursor T-cell lymphoblastic lymphoma -> Lymphoma, T-Cell [NCI g] A1N0
    ("C0520676", "C0011581"): ("exclude_wrong", "target is a different disease or too generic"),  # Premenstrual Dysphoric Disorder -> Depressive Disorder [MSH p] A5N1
    ("C0001815", "C0027022"): ("exclude_wrong", "target is a different disease or too generic"),  # Primary Myelofibrosis -> Myeloproliferative Disorders [MSH p] A1N15
    ("C0750979", "C0006118"): ("keep", "genuine subtype / clinical form of the target"),  # Primary malignant neoplasm of brain -> Brain Neoplasm [MSH p] A1N0
    ("C0029291", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Psittacosis -> Bacterial Infection [NCI g] A9N0
    ("C0392164", "C0010674"): ("keep", "genuine subtype / clinical form of the target"),  # Pulmonary Cystic Fibrosis -> Cystic Fibrosis [MSH p] A7N0
    ("C0034063", "C0024115"): ("exclude_wrong", "target is a different disease or too generic"),  # Pulmonary Edema -> Lung Diseases [MSH p] A1N1
    ("C0034065", "C0024115"): ("exclude_wrong", "target is a different disease or too generic"),  # Pulmonary Embolism -> Lung Diseases [MSH p] A1N0
    ("C0034067", "C0024117"): ("keep", "genuine subtype / clinical form of the target"),  # Pulmonary Emphysema -> Chronic Obstructive Pulmonary Disease [NCI p] A20N0
    ("C0034069", "C0206062"): ("keep", "genuine subtype / clinical form of the target"),  # Pulmonary Fibrosis -> Interstitial Lung Disease [MSH p] A1N0
    ("C0276688", "C0032285"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Pulmonary cryptococcosis -> Pneumonia [NCI p] A1N0
    ("C0349649", "C0751674"): ("keep", "genuine subtype / clinical form of the target"),  # Pulmonary lymphangioleiomyomatosis -> Lymphangioleiomyomatosis [NCI p] A1N0
    ("C0030246", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Pustulosis of Palms and Soles -> Dermatitis [NCI p] A1N0
    ("C0030246", "C0033860"): ("keep", "genuine subtype / clinical form of the target"),  # Pustulosis of Palms and Soles -> Psoriasis [MSH p] A1N0
    ("C0034186", "C0022658"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Pyelonephritis -> Kidney Diseases [NCI g] A7N3
    ("C0854802", "C0023434"): ("keep", "genuine subtype / clinical form of the target"),  # Recurrent Chronic Lymphoid Leukemia -> Chronic Lymphocytic Leukemia [NCI p] A1N4
    ("C4528668", "C0023467"): ("keep", "genuine subtype / clinical form of the target"),  # Refractory acute myeloid leukemia -> Acute Myeloid Leukemia [NCI p] A1N16
    ("C2981142", "C3463824"): ("keep", "genuine subtype / clinical form of the target"),  # Refractory anemia, without ringed siderobl -> Myelodysplastic Syndrome [NCI p] A1N0
    ("C0238217", "C0022658"): ("exclude_wrong", "target is a different disease or too generic"),  # Renal transplant rejection -> Kidney Diseases [NCI p] A5N1
    ("C0035436", "C0003864"): ("exclude_wrong", "target is a different disease or too generic"),  # Rheumatic Fever -> Arthritis [MSH p] A1N0
    ("C0035436", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Rheumatic Fever -> Bacterial Infection [NCI p] A1N0
    ("C0035457", "C2607914"): ("keep", "genuine subtype / clinical form of the target"),  # Rhinitis, Allergic, Perennial -> Allergic Rhinitis [NCI p] A2N0
    ("C0035460", "C0035455"): ("exclude_wrong", "target is a different disease or too generic"),  # Rhinitis, Vasomotor -> Rhinitis [NCI p] A41N0
    ("C0035793", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Rocky Mountain Spotted Fever -> Bacterial Infection [NCI g] A9N0
    ("C0152486", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Salmonella sepsis -> Bacterial Infection [NCI p] A2N0
    ("C0036285", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Scarlet Fever -> Bacterial Infection [NCI g] A3N0
    ("C0085159", "C0011581"): ("keep", "genuine subtype / clinical form of the target"),  # Seasonal Affective Disorder -> Depressive Disorder [MSH p] A1N0
    ("C0036508", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Seborrheic dermatitis -> Dermatitis [NCI g] A44N0
    ("C0036508", "C0013595"): ("exclude_wrong", "target is a different disease or too generic"),  # Seborrheic dermatitis -> Eczema [NCI p] A44N0
    ("C0268407", "C0002726"): ("keep", "genuine subtype / clinical form of the target"),  # Senile cardiac amyloidosis -> Amyloidosis [NCI p] A2N0
    ("C0011268", "C0497327"): ("keep", "genuine subtype / clinical form of the target"),  # Senile dementia -> Dementia [NCI p] A4N0
    ("C0456103", "C0243026"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Sepsis of the newborn -> Sepsis [MSH p] A7N0
    ("C0341950", "C0032914"): ("keep", "genuine subtype / clinical form of the target"),  # Severe pre-eclampsia -> Pre-Eclampsia [NCI p] A1N0
    ("C0013371", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Shigella Infections -> Bacterial Infection [NCI p] A7N0
    ("C0036980", "C0027051"): ("exclude_wrong", "target is a different disease or too generic"),  # Shock, Cardiogenic -> Myocardial Infarction [MSH p] A1N1
    ("C1527336", "C0003873"): ("exclude_wrong", "target is a different disease or too generic"),  # Sjogren's Syndrome -> Arthritis, Rheumatoid [MSH p] A1N1
    ("C0162627", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Skin Diseases, Bacterial -> Bacterial Infection [MSH p] A11N0
    ("C0152030", "C0021368"): ("exclude_wrong", "target is a different disease or too generic"),  # Skin irritation -> Inflammation [NCI g] A21N0
    ("C0855095", "C0079731"): ("exclude_wrong", "target is a different disease or too generic"),  # Small Lymphocytic Lymphoma -> B-Cell Non-Hodgkin Lymphoma [NCI g] A4N0
    ("C0149925", "C0684249"): ("exclude_wrong", "SCLC is a distinct lung cancer; lung carcinoma signature a poor proxy"),  # Small cell carcinoma of lung -> Lung Carcinoma [NCI g] A5N0
    ("C0037354", "C0042769"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Smallpox -> Viral Infection [NCI p] A1N0
    ("C0242429", "C0031350"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Sore Throat -> Pharyngitis [MSH p] A5N0
    ("C0026847", "C0037928"): ("exclude_wrong", "target is a different disease or too generic"),  # Spinal Muscular Atrophy -> Spinal Cord Diseases [MSH p] A3N0
    ("C0149782", "C0684249"): ("keep", "genuine subtype / clinical form of the target"),  # Squamous cell carcinoma of lung -> Lung Carcinoma [NCI p] A1N0
    ("C1319317", "C0747548"): ("keep", "genuine subtype / clinical form of the target"),  # Squamous cell carcinoma of pharynx -> Pharyngeal Carcinoma [NCI p] A1N0
    ("C1319317", "C1168401"): ("keep", "genuine subtype / clinical form of the target"),  # Squamous cell carcinoma of pharynx -> Head and Neck Squamous Cell Carcinoma [NCI g] A1N0
    ("C4509816", "C0007131"): ("keep", "genuine subtype / clinical form of the target"),  # Squamous non-small cell lung cancer -> Carcinoma, Non-Small-Cell Lung [NCI p] A1N4
    ("C0038157", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Staphylococcal enterocolitis -> Bacterial Infection [NCI g] A1N0
    ("C1318973", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Staphylococcus aureus infection -> Bacterial Infection [NCI g] A49N0
    ("C0271093", "C0024437"): ("exclude_wrong", "target is a different disease or too generic"),  # Stargardt's disease -> Macular Degeneration [NCI p] A1N0
    ("C0038220", "C0036572"): ("exclude_wrong", "target is a different disease or too generic"),  # Status Epilepticus -> Seizure [MSH p] A6N6
    ("C0038395", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Streptococcal Infections -> Bacterial Infection [NCI p] A4N0
    ("C0862889", "C0007117"): ("keep", "genuine subtype / clinical form of the target"),  # Superficial basal cell carcinoma -> Basal Cell Carcinoma [NCI g] A2N0
    ("C0039128", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Syphilis -> Bacterial Infection [NCI g] A12N0
    ("C0036421", "C0011644"): ("keep", "genuine subtype / clinical form of the target"),  # Systemic Scleroderma -> Scleroderma [NCI p] A3N9
    ("C0039520", "C1568272"): ("keep", "genuine subtype / clinical form of the target"),  # Tenosynovitis -> Tendinopathy [MSH p] A10N1
    ("C0039614", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Tetanus -> Bacterial Infection [NCI p] A2N0
    ("C0040021", "C0042384"): ("exclude_wrong", "target is a different disease or too generic"),  # Thromboangiitis Obliterans -> Vasculitis [MSH p] A1N0
    ("C0494475", "C0036572"): ("exclude_wrong", "target is a different disease or too generic"),  # Tonic - clonic seizures -> Seizure [NCI g] A5N0
    ("C0014549", "C0014544"): ("exclude_wrong", "target is a different disease or too generic"),  # Tonic-Clonic Epilepsy -> Epilepsy [NCI g] A8N0
    ("C0040425", "C0031350"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Tonsillitis -> Pharyngitis [MSH p] A9N1
    ("C0040558", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Toxoplasmosis -> Protozoal Infection [NCI p] A3N0
    ("C0040560", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Toxoplasmosis, Congenital -> Protozoal Infection [NCI g] A1N0
    ("C0040921", "C0033740"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Trichomonas Infections -> Protozoal Infection [NCI p] A6N0
    ("C3539878", "C1458155"): ("exclude_wrong", "target is a different disease or too generic"),  # Triple Negative Breast Neoplasms -> Breast Neoplasm [MSH p] A2N4
    ("C0041318", "C0041296"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Tuberculosis, Meningeal -> Tuberculosis [NCI p] A1N0
    ("C0041351", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Tularemia -> Bacterial Infection [NCI p] A4N0
    ("C0041466", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Typhoid Fever -> Bacterial Infection [NCI p] A5N0
    ("C0041473", "C0004623"): ("exclude_infection", "infectious indication: the drug acts on the pathogen"),  # Typhus, Epidemic Louse-Borne -> Bacterial Infection [NCI g] A2N0
    ("C0042109", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Urticaria -> Dermatitis [NCI p] A28N2
    ("C0042373", "C0007222"): ("exclude_wrong", "target is a different disease or too generic"),  # Vascular Diseases -> Cardiovascular Diseases [NCI p] A3N0
    ("C0042571", "C0042594"): ("exclude_wrong", "target is a different disease or too generic"),  # Vertigo -> Vestibular Diseases [MSH p] A4N1
    ("C0042900", "C0011603"): ("exclude_wrong", "target is a different disease or too generic"),  # Vitiligo -> Dermatitis [NCI p] A2N2
    ("C0043515", "C0017185"): ("exclude_wrong", "target is a different disease or too generic"),  # Zollinger-Ellison syndrome -> Digestive System Neoplasm [MSH p] A10N0
    ("C0005283", "C0039730"): ("keep", "genuine subtype / clinical form of the target"),  # beta Thalassemia -> Thalassemia [NCI p] A3N0
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
    (cui -> {node: offset}, node -> label, child node -> set of parent
    nodes, version) for the given CUIs, nodes being "SAB:code" in
    UMLS_SOURCES. offset is 0 when the CUI holds the node's preferred term
    (UMLS_PREFERRED_TTY), 1 otherwise. Node labels come from preferred atoms
    and from the parents endpoint. Parents are fetched up to max_levels steps
    above the CUIs' own codes.
    """
    import requests
    from local import UMLS_API_KEY

    cache = {"version": UMLS_VERSION, "atoms_tty": {}, "parents": {}}
    if UMLS_CACHE_FILE.exists():
        with open(UMLS_CACHE_FILE, encoding="utf-8") as fh:
            cache = json.load(fh)
    cache.pop("atoms", None)                 # older cache layout, without term types
    cache.setdefault("atoms_tty", {})
    sabs = ",".join(UMLS_SOURCES)
    session = requests.Session()
    n_calls = 0

    todo = sorted(c for c in cuis if c not in cache["atoms_tty"])
    for i, cui in enumerate(todo, 1):
        atoms = _uts_get_atoms(session, cui, sabs, UMLS_API_KEY)
        cache["atoms_tty"][cui] = sorted({(a["rootSource"], a["code"].rsplit("/", 1)[1],
                                           a.get("termType", ""), a["name"])
                                          for a in atoms if a.get("rootSource") in UMLS_SOURCES})
        n_calls += 1
        if i % 200 == 0:
            print(f"UTS atoms: {i}/{len(todo)}")
            _save_cache(cache)

    cui2n, label = {}, {}
    for cui in cuis:
        for sab, code, tty, name in cache["atoms_tty"].get(cui, []):
            n = f"{sab}:{code}"
            off = 0 if tty in UMLS_PREFERRED_TTY.get(sab, ()) else 1
            nodes = cui2n.setdefault(cui, {})
            nodes[n] = min(off, nodes.get(n, 1))
            if off == 0:
                label[n] = name

    parents, frontier = {}, set().union(*(set(v) for v in cui2n.values()))
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
                if name:
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


def _nodes(cui2m, cui):
    """{node: offset} of a CUI; MONDO nodes (sets) all have offset 0."""
    v = cui2m.get(cui, {})
    return v if isinstance(v, dict) else {n: 0 for n in v}


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
        for m, off in _nodes(cui2m, cui).items():
            if off == 0:                           # preferred term of the node only
                sig_by_m.setdefault(m, set()).add(cui)

    rows = []
    for cui in sorted(set(pairs["ind_id"]) - sig):
        hits = {}                                   # target CUI -> (level, node)
        for src in order:
            for m, off in _nodes(cui2m, cui).items():
                if _source(m) != src:
                    continue
                reach = {m: 0, **ancestors(m, parents, MAX_LEVELS - off)}
                for a, d in ((a, d + off) for a, d in reach.items()):
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
        if level == 2 and _source(hits[target][1]) in EXCLUDE_LEVEL2_SOURCES:
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
            "target_cui_preferred_labels": "|".join(sorted(
                label.get(n, "") for n, o in _nodes(cui2m, target).items() if o == 0)),
            "tied_target_cuis": "|".join(tied),
            "target_max_psamples": max_ps[target],
            "child_nodes": "|".join(sorted(f"{n}{'' if o == 0 else '(non-preferred)'}"
                                           for n, o in _nodes(cui2m, cui).items())),
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
