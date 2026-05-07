# NetCoS — Project Context for Claude Code

## What this project is
NetCoS (Network-Enhanced Connectivity Score) is a drug-repurposing pipeline that ranks reference drugs by similarity to a disease's transcriptomic signature, projecting both into a protein-protein / protein-microRNA / protein-metabolite interaction network and propagating signals across it via MITHrIL before computing the connectivity score (Kolmogorov-Smirnov–based).

Public mirror: https://github.com/knowmics-lab/netcos
Author: Lorenzo Signorini (knowmics-lab @ University of Catania)

## Pipeline stages
Disease & drug raw transcriptomic data → preprocessing → Linear Mixed Model (lme4 / Julia) → meta-analysis (drug side, 6h + 24h time points) → MITHrIL signal propagation → KS-based connectivity score → drug ranking.

## Validation pipeline stages
Disease & drug raw fold change data -> MITHrIL propagation -> drug ranking -> Validation against known dataset (Precision/Recall)

## Tech stack
- **Python** — pipeline orchestration, connectivity score, ChEMBL validation, plotting. Deps in `requirements.txt`: `joblib`, `matplotlib`, `numpy`, `pandas`, `pyreadr`, `requests`, `scipy`. Package metadata in `setup.py` (name `drug_repurposing`, v0.1.0).
- **R** — LMM fitting, LINCS data handling, signature meta-analysis. Lives under `tsr/` (legacy workflows) and `new_tsr/R/` (refactored OOP version).
- **Julia** — LMM optimization (`tsr/julia/`, `modules/drug_signature/julia/`).
- **MITHrIL** — external Java tool, run as a subprocess (`step2_disease_run_MITHrIL.py`, `step6_drug_run_MITHrIL.py`).

## Repo layout
```
src/                                Python pipeline (numbered step scripts + utilities)
  conf.py                           Active configuration (swapped from configs/)
  local.py                          Machine-local paths (gitignored)
  step{1..7}_*.py                   Numbered pipeline steps (disease + drug branches)
  cs_batch.py                       Batch entry point for connectivity score
  9_calculate_connectivity_score_wrapper.py   Final scoring entry point
  validations/chembl/               ChEMBL IC50 validation, hyperparameter search
modules/                            R algorithm modules (connectivity_score, disease_signature, drug_signature) — abstract-class architecture
tsr/                                R workflows per disease (als_1, als_2, als_NYGC, ipf, colon_tumor)
new_tsr/R/                          Refactored R modules (70+ files, modern OOP)
configs/                            Config templates — copy a file into src/ and rename to conf.py to switch experiments
notebooks/                          NetCos_Pipeline.ipynb, ChEMBL validation notebooks
data exploration/                   EDA notebooks (ALS, IPF, ChEMBL)
data/                               Reference + test datasets (BinChen2017 for validation, LINCS placeholders) — gitignored
instructions/                       Local-only docs, logs, meeting notes — gitignored, mix of English and Italian
```

## Config system — important
There is no CLI flag to switch experiments. Configuration is controlled by which `conf.py` is active in `src/`. Templates live in `configs/<experiment>/conf.py`; to switch, copy that file into `src/` and rename to `conf.py` (see `configs/README.txt`). Treat `src/conf.py` as a swappable artifact, not a stable file. `src/local.py` holds machine-specific paths and is gitignored — never commit it.

## Local-only files (gitignored)
- `data/` — large reference datasets
- `instructions/` — personal logs (`log.txt`, meeting notes), pipeline mechanics in Italian (`come funge mithril.txt`), validation logs, version control notes. **Read these for context** when working on a specific stage; they're often the most up-to-date source of truth.
- `src/local.py`
- `tsr/log_data/`

## Coding conventions
- Python: standard scientific stack (numpy/pandas/scipy idioms). Step scripts are linear and side-effect heavy — they read from `conf.py`, write to disk, and don't return values. When refactoring, preserve the numbered-step file structure unless asked to restructure.
- R: `new_tsr/R/` uses an abstract-class pattern; prefer extending it over modifying `modules/` for new R work.
- Logging: use `src/logger.py` (Python) and `modules/TSRLogger.R` (R) — don't add `print` debugging.
- Don't auto-format files wholesale; the codebase mixes styles intentionally.

## Running the pipeline
- Full pipeline: orchestrated by the numbered `step1..7` scripts in `src/`, typically invoked from `cs_batch.py`.
- ChEMBL validation: `src/validations/chembl/` — start at `Preprocess_binChen2017_TCGA_LINCS_data.py`.
- R workflows: enter the relevant `tsr/<disease>/` folder; main scripts assume that as the working directory.

## What I'd like Claude Code to help with
- Read instructions/log.txt to know what the last task is
- Refactoring the numbered step scripts into a cleaner orchestration layer without losing the step boundaries.
- Writing unit tests around the connectivity score computation in `src/connectivity_score.py`.
- Improving documentation for the pipeline (the README covers theory; `src/` lacks docstrings).
- Reviewing PRs against `knowmics-lab/netcos`.
- Debugging MITHrIL subprocess calls and LINCS preprocessing issues.

## Don'ts
- Don't commit anything under `data/`, `instructions/`, `tsr/log_data/`, or `src/local.py`.
- Don't paraphrase or rename pipeline stage names — they map to terminology in published papers.
- Don't run heavy R/Julia jobs without asking — single LMM runs can take hours.
- Don't edit any file without explicit permission from the user
