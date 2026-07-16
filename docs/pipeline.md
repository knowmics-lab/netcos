# NetCoS — Pipeline & Architecture Reference

The precise, stage-by-stage specification of what the NetCoS pipeline does. This is the "what does the pipeline do at all times" reference: it complements `README.md` (theory + flowchart) and `paper/` (the scientific write-up).

NetCoS (Network-Enhanced Connectivity Score) ranks reference drugs by similarity to a disease's transcriptomic signature, projecting both into a protein–protein / protein–microRNA / protein–metabolite interaction network and propagating signals across it via MITHrIL before computing a Kolmogorov–Smirnov-based connectivity score.

## Two pipelines that share building blocks but invert one step

NetCoS is run through two parallel pipelines. They share the same building blocks — preprocessing, optional MITHrIL propagation, KS-based connectivity score (CS) — but they **invert the drug-collapse step** (where the many perturbation instances of a drug become one signature).

### TSR / LMM pipeline — real-disease data (ALS, IPF)

```
disease & drug raw transcriptomic data
  → preprocessing
  → Linear Mixed Model (lme4 / Julia)
  → meta-analysis (drug side, 6h + 24h time points)   ← DRUG COLLAPSE HAPPENS HERE, BEFORE MITHrIL/CS
  → MITHrIL signal propagation
  → KS-based connectivity score
  → drug ranking
```

Each drug becomes **one** signature (via LMM + meta-analysis) *before* propagation and scoring.

### Validation pipeline — BinChen2017 / ChEMBL

```
disease & drug raw fold-change / LINCS data
  → (optional MITHrIL propagation)
  → KS-based connectivity score computed PER PERTURBATION INSTANCE   ← drug NOT collapsed yet
  → DRUG COLLAPSE HAPPENS HERE, AFTER CS (median / best / sRGES)
  → validation against ChEMBL IC50 (precision/recall, Spearman correlation)
```

Drug collapse happens *after* CS, via `src/validations/chembl/final_step_correlations_IC50_CS.py`.

This pipeline also replicates BinChen2017's own RGES/sRGES numbers as a ground-truth check on the CS algorithm itself, before trusting any ChEMBL correlation built on top of it:

- `src/validations/binchen2017/` — RGES/sRGES replication (ground-truth check on the CS algorithm).
- `src/validations/chembl/` — IC50 correlation + hyperparameter grid search (built on top of the above).

## The CS run registry (`logs/cs_runs.tsv`)

`logs/cs_runs.tsv` is the run registry for both pipelines — every `cs_batch` run is logged there with its hyperparameters and output filename. Resolve "the correct CS file for a given config" via `src/loader.py::load_cs_run()` / `pick_canonical_cs_run()`, **not** by hardcoding paths or globbing `connectivity_score/output/`.

## Config system

There is no CLI flag to switch experiments. Configuration is controlled by **which `conf.py` is active in `src/`**. The `configs/<name>/` folders are **frozen backups of past runs**, not generic templates: each is a snapshot of `conf.py` (and sometimes the specific `src/*.py` step scripts) exactly as they existed when that run produced its results — effectively a versioned record of one pipeline result, kept reproducible even after `src/` moves on. They are named descriptively with dates, e.g. `configs/05 06 2026 Bin Chen 2017 replication + IC50 validation/`.

To reuse one, copy its `conf.py` into `src/` and rename (see `configs/README.txt`) — don't assume every file inside is meant to be copied back; some are historical record only. Treat `src/conf.py` as a swappable artifact, not a stable file. `src/local.py` holds machine-specific paths and is gitignored.

## Repository layout

```
src/                        Python pipeline (numbered step scripts + utilities)
  conf.py                   Active configuration (swapped from configs/)
  local.py                  Machine-local paths (gitignored)
  loader.py                 Single source of truth for resolving/loading CS run files
                             (load_cs_run, pick_canonical_cs_run, cs_out_dir_for)
  sRGES.py                  BinChen2017 sRGES aggregation (faithful port of BinChen-Lab/RGES)
  step{1..8}_*.py           Numbered pipeline steps (disease + drug branches)
  cs_batch.py               Batch entry point for connectivity score
  9_calculate_connectivity_score_wrapper.py   Final scoring entry point
  validations/chembl/       BinChen2017/ChEMBL validation: IC50 correlation, precision/recall, grid search
  validations/binchen2017/  BinChen2017 RGES/sRGES replication (upstream ground-truth check)
modules/                    R algorithm modules (connectivity_score, disease_signature, drug_signature) — abstract-class architecture
tsr/                        R workflows per disease (als_1, als_2, als_NYGC, ipf, colon_tumor)
new_tsr/R/                  Refactored R modules (70+ files, modern OOP)
configs/                    Frozen snapshots/backups of past pipeline runs (see "Config system")
notebooks/                  NetCos_Pipeline.ipynb, ChEMBL validation notebooks
data_exploration/           CS/ranking + EDA notebooks (IPF, ALS_NYGC, ChEMBL)
paper/                      NetCoS methods paper draft (NetCoS_methods_paper.md) + bib tooling
data/                       Reference + test datasets (BinChen2017 for validation, LINCS placeholders) — gitignored
logs/                       cs_runs.tsv = registry of every cs_batch run; grid-search TSVs under chembl_bootstrap_grid_search/
```

## Running the pipeline

- **Full pipeline:** orchestrated by the numbered `step1..7` scripts in `src/`, typically invoked from `cs_batch.py`.
- **ChEMBL validation:** `src/validations/chembl/` — start at `Preprocess_binChen2017_TCGA_LINCS_data.py`.
- **R workflows:** enter the relevant `tsr/<disease>/` folder; main scripts assume that as the working directory.

> Heavy R/Julia jobs can take hours (a single LMM run especially) — plan accordingly.
