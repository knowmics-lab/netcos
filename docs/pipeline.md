# NetCoS — Pipeline & Architecture Reference

The precise, stage-by-stage specification of what the NetCoS pipeline does. It complements `README.md` (theory + flowchart) and `paper/` (the scientific write-up). Where this document and the code disagree, the code is authoritative.

NetCoS (Network-Enhanced Connectivity Score) ranks reference drugs by similarity to a disease's transcriptomic signature, optionally projecting both into a protein–protein / protein–microRNA / protein–metabolite interaction network and propagating signals across it via MITHrIL before computing a Kolmogorov–Smirnov-based connectivity score (CS).

## Two pipelines that share building blocks but invert one step

NetCoS is run through two parallel pipelines. They share the same building blocks — preprocessing, optional MITHrIL propagation, KS-based CS — but they **invert the drug-collapse step** (where the many perturbation instances of a drug become one signature). The `conf.py` flag `drug_collapsed_before_cs` selects which one `cs_batch.py` is running.

### TSR / LMM pipeline — real-disease data (`drug_collapsed_before_cs = True`)

```
disease & drug raw transcriptomic data
  → preprocessing
  → Linear Mixed Model (lme4 / Julia)
  → meta-analysis (drug side, 6h / 24h / pooled 6h_24h)   ← DRUG COLLAPSE HAPPENS HERE, BEFORE MITHrIL/CS
  → (optional) MITHrIL signal propagation
  → KS-based connectivity score
  → drug ranking
```

Each drug becomes **one** signature per timepoint (via LMM + meta-analysis) *before* propagation and scoring. The CS input ids are already `pert_id`s (`<drug>_<timestep>`), so `cs_batch` writes them straight to the `pert_id` column.

### Validation pipeline — BinChen2017 / ChEMBL (`drug_collapsed_before_cs = False`)

```
disease & drug raw fold-change / LINCS data
  → (optional MITHrIL propagation)
  → KS-based connectivity score computed PER PERTURBATION INSTANCE   ← drug NOT collapsed yet
  → DRUG COLLAPSE HAPPENS HERE, AFTER CS (median / best / sRGES)
  → validation against ChEMBL IC50 (precision/recall, Spearman correlation)
```

The CS output carries one row per LINCS perturbation (`LINCS_id`); `pert_id` is mapped in afterwards (`add_pert_id_to_cs`), and drug collapse happens via `src/validations/chembl/final_step_correlations_IC50_CS.py` (method set by `CS_DRUG_COLLAPSE_METHOD`).

This pipeline also replicates BinChen2017's own RGES/sRGES numbers as a ground-truth check on the CS algorithm itself, before trusting any ChEMBL correlation built on top of it:

- `src/validations/binchen2017/` — RGES/sRGES replication (ground-truth check on the CS algorithm); `src/sRGES.py` is the sRGES aggregation it tests.
- `src/validations/chembl/` — IC50 correlation + hyperparameter grid search (built on top of the above).

## Signature representations (what the CS is computed on)

Three `conf.py` switches decide which representation of the signatures enters the CS. `conf.validate_hyperparameters()` rejects the unsupported combinations.

| `cs_mith` | `cs_on_LM` | `CS_ON_PATHWAYS` | Representation | Value column |
|---|---|---|---|---|
| 0 | 0 | False | DEG — gene-level differential expression, all genes | `DE_log2_FC` |
| 0 | 1 | False | DEG_LM — DEG, L1000 landmark genes only | `DE_log2_FC` |
| 1 | 0 | False | mith — MITHrIL perturbation, all nodes | `Perturbation` |
| 1 | 1 | False | mith_LM — MITHrIL perturbation, landmark genes only | `Perturbation` |
| 1 | 0 | True | mith_pw — MITHrIL perturbation aggregated to pathways | `Corrected Accumulator` |

Not allowed: `cs_mith=0` with `CS_ON_PATHWAYS=True` (pathways exist only after MITHrIL), and `cs_on_LM=1` with `CS_ON_PATHWAYS=True`.

## Stage-by-stage: the numbered scripts in `src/`

Disease side:

- `step1_disease_gene_name_to_id_map.py` — maps the disease signature from gene symbols to Entrez ids, removes duplicates, writes the MITHrIL input (`.mi`). Used by the `tsr/` workflows (ALS, IPF); case studies write the `.mi` themselves (see "Running a new disease").
- `step2_disease_run_MITHrIL.py` — runs MITHrIL (Java, subprocess) on the disease input of `conf.disease_run_name`. Without arguments it prints the command; `python step2_disease_run_MITHrIL.py run` executes it (only where `MITH_APP` exists) and appends a row to `logs/mithril_disease_runs.tsv`. The p-values use MITHrIL's default `-number-of-iterations` (2001).
- `step3_mith_out_to_cs_in_disease.py` — gene-level MITHrIL perturbations → disease CS input.

Drug side:

- `step4_map_metanalysis_wrapper_drug.py` — calls `map_metanalysis_gene_wise_to_drug_wise.R` in batches (gene-wise → drug-wise meta-analysis tables).
- `step5_mith2_to_mith3_drug_input.py` — builds the MITHrIL 3 drug input matrix.
- `step6_drug_run_MITHrIL.py` — runs MITHrIL on the drug signatures.
- `step7_drug_mith_out_to_cs_in.py` — gene-level MITHrIL output → drug CS input, one file per signature.
- `step8_DEG_drug_metanalysis_to_cs_in.py` — DEG (no MITHrIL) drug CS input from the per-drug meta-analysis, as `<drug>_<timestep>_signature_gene_id.pkl`.

Both sides:

- `step2b_7b_drug_and_disease_mith_out_to_cs_in.py` — pathway-level CS input from MITHrIL's pathway summary (`*.output.txt`, signed `Corrected Accumulator`). It has a disease runner and a drug runner.

Scoring:

- `cs_batch.py::run_cs_batch_for_conf()` — the single-configuration entry point. Any argument left `None` is read from `conf.py`. It parallelises over drugs (`joblib`), writes one timestamped CS file to `connectivity_score/output/<disease_run_name>/`, and appends one row to `logs/cs_runs.tsv`.
- `run_cs_batch_grid.py` — sweeps the pre-CS switches (above) for one disease on the collapsed-drug pipeline and calls `run_cs_batch_for_conf` per combination, skipping combinations already in the registry.
- `9_calculate_connectivity_score_wrapper.py` — legacy, superseded by `cs_batch.py`. It launches `connectivity_score.py` as a script, which no longer has an entry point.
- `validations/chembl/call_cs_batch_for_different_confs.py` — the ChEMBL-grid wrapper for the validation pipeline. It is specific to the BinChen2017 diseases; do not use it for other diseases.

## Connectivity score methods and sign convention

`CS_METHOD` picks the scoring function from `connectivity_score.calc_connectivity_score_with`:

- `bin_chen` — RGES as defined in Chen et al. 2017 (`bin_chen_connectivity`), p-value from 1 000 Monte Carlo permutations.
- `bin_chen_disease_sorted` — variant that sorts the disease genes (`bin_chen_connectivity_sorted`); `CS_TH = -1.5` is calibrated for this method only.
- `lamb` — declared but not implemented.

Genes (or pathways) are ranked by magnitude by default (`rank_on='magnitude'`).

**Sign convention:** the most negative CS means the drug signature most strongly *reverses* the disease signature. For reversal studies (diseases), candidates are therefore the **most negative** CS. For studies where the input signature is a perturbation to be *mimicked* rather than reversed (e.g. `mir29_antagomir`), candidates are the **most positive** CS — the score is the same, only the reading changes.

`bin_chen` / `bin_chen_disease_sorted` runs logged before the `rank_genes()` sign fix (11/05/2026) use the inverted sign and must not be used or compared.

## Where the disease signature is read from (`DISEASE_SIGNATURE_SOURCE`)

For DEG-level runs, `loader.load_disease_signature()` reads the disease signature from one of three places, chosen by `conf.DISEASE_SIGNATURE_SOURCE`:

| Source | File | Used by |
|---|---|---|
| `"tsr"` | `tsr/output/disease_signature/<DISEASE>/<DISEASE>_signature_gene_id.csv` | R TSR workflows (ALS, IPF); `mir29_antagomir` |
| `"cs_in"` | `connectivity_score/input/disease_signature_2025/<run>/<run>_signature_gene_id.csv` | signatures produced outside `tsr/`: `stan_als_tdp43`, `asthma_GSE152004` |
| `"other"` | `connectivity_score/input/disease_signature_2025/<run>/<run>_signature_gene_id.txt` | BinChen2017 disease signatures |

All three are standardised to `gene_id, DE_log2_FC, adj.p.value`. For `"other"` the column mapping (`id`, `logFC`, `adj_p_value`) is provisional in the code. The inline comment in `conf.py` still lists only `"tsr"` and `"other"`; `"cs_in"` is accepted by the loader. MITHrIL-level runs always read from `connectivity_score/input/disease_signature_2025/<run>/`.

## The CS run registry (`logs/cs_runs.tsv`)

`logs/cs_runs.tsv` is the run registry for both pipelines — every `run_cs_batch_for_conf` call appends a row with its hyperparameters, gene counts at each filter, input directories and output filename. The `cs_run_id` is the output filename stem: `<dd_mm_YYYY_HH_MM>_{mith|DEG}[_LM][_pw]_<CS_METHOD>_connectivity_score`.

Resolve "the correct CS file for a given config" via `src/loader.py::load_cs_run()` / `pick_canonical_cs_run()` (and `cs_out_dir_for()` for the directory), **not** by hardcoding paths or globbing `connectivity_score/output/`.

## Config system

There is no CLI flag to switch experiments. Configuration is controlled by **which `conf.py` is active in `src/`**. Every file path is defined once in `conf.py` and imported (`BASE_DIR`, `IMG_DIR`, `OOUT_DIR`, `LOGS_DIR`, …); scripts do not build paths from literals.

The `configs/<name>/` folders are **frozen backups of past runs**, not generic templates: each is a snapshot of `conf.py` (and sometimes the specific `src/*.py` step scripts) exactly as they existed when that run produced its results. They are named descriptively with dates, e.g. `configs/05 06 2026 Bin Chen 2017 replication + IC50 validation/`.

To reuse one, copy its `conf.py` into `src/` and rename (see `configs/README.txt`) — don't assume every file inside is meant to be copied back; some are historical record only. Treat `src/conf.py` as a swappable artifact, not a stable file. `src/local.py` holds machine-specific paths and is gitignored.

## Analyses built on the pipeline

Each analysis has one slug, mirrored across `src/`, `imgs/` and `other_outputs/`.

- `src/validations/known_drug_benchmark/` — where known ALS/IPF drugs rank under each signature representation (DEG, DEG_LM, mith, mith_LM, mith_pw), with CS method, ranking and drug set held fixed. Descriptive; tests no hypothesis.
- `src/validations/drug_disease_pairs_validation/` — builds a labelled drug–disease pair benchmark (approved vs failed indications) from metadata, without running the CS:
  - `01_` DiSignAtlas / CREEDS disease-signature inventory;
  - `02_` repoDB structure;
  - `03_` repoDB × disease signatures, joined on UMLS CUI;
  - `04_` repoDB drugs → LINCS L1000;
  - `05_` hierarchy roll-up for CUIs without their own signature, through the UMLS source vocabularies (NCIt, MeSH; UTS REST API).
- `src/case_studies/` — single-dataset applications that produce new disease-side inputs:
  - `stan_als_tdp43/` — TDP-43 iPSC motor-neuron bulk RNA-seq: `01_` edgeR QLF + QC from counts, `02_` contrast → DEG CS input (`"cs_in"`) and MITHrIL `.mi`, one disease run per contrast.
  - `mir29_antagomir/` — mouse miR-29 antagomiR log2FC signature: `01_` mouse → human 1:1 orthologs (MGI homology report) → DEG CS input and MITHrIL `.mi`. No p-values in the input (`adj.p.value` = NaN). Mimicking study: read the most positive CS.
  - `asthma_GSE152004/` — asthmatic vs healthy-control nasal airway epithelium RNA-seq (GSE152004, 440 vs 254), a disease of the labelled pair benchmark: `00_` inspection of the GEO download and sample-group map, `01_` disease signature with the `modules/disease_signature/` R classes on NCBI GeneIDs → DEG CS input (`"cs_in"`) and MITHrIL `.mi`.

## Repository layout

```
src/                        Python pipeline — code only, no outputs
  conf.py                   Active configuration (swapped from configs/)
  local.py                  Machine-local paths (gitignored)
  loader.py                 Single source of truth for resolving/loading CS run files
                             (load_cs_run, pick_canonical_cs_run, cs_out_dir_for)
                             and disease signatures (load_disease_signature)
  connectivity_score.py     KS / RGES scoring functions
  sRGES.py                  BinChen2017 sRGES aggregation (faithful port of BinChen-Lab/RGES)
  step{1..8}_*.py, step2b_7b_*.py   Numbered pipeline steps (disease + drug branches)
  cs_batch.py               CS entry point (run_cs_batch_for_conf)
  run_cs_batch_grid.py      Pre-CS hyperparameter sweep for one disease
  validations/chembl/       BinChen2017/ChEMBL validation: IC50 correlation, precision/recall, grid search
  validations/binchen2017/  BinChen2017 RGES/sRGES replication (upstream ground-truth check)
  validations/known_drug_benchmark/       Known-drug positioning across representations
  validations/drug_disease_pairs_validation/   Labelled pair benchmark construction (01_–05_)
  case_studies/             Single-dataset applications (stan_als_tdp43, mir29_antagomir, asthma_GSE152004)
  run_disease_pipeline.sh   Disease side + CS grid for the disease set in conf.py
modules/                    R algorithm modules (connectivity_score, disease_signature, drug_signature) — abstract-class architecture
tsr/                        R workflows per disease (als_1, als_2, als_NYGC, ipf, colon_tumor)
new_tsr/R/                  Refactored R modules (modern OOP)
MITHrIL/                    MITHrIL inputs and outputs (input/, output/ per disease and drug run)
connectivity_score/         CS inputs (input/disease_signature_2025/, input/drug_signature_2025/) and outputs (output/<disease_run>/)
configs/                    Frozen snapshots/backups of past pipeline runs (see "Config system")
notebooks/                  NetCos_Pipeline.ipynb, ChEMBL validation notebooks
data_exploration/           CS/ranking + EDA notebooks (IPF, ALS_NYGC, asthma_GSE152004, ChEMBL)
imgs/                       Every image the project produces, one folder per analysis slug
other_outputs/              Every non-image result (tables, summaries), one folder per analysis slug
paper/                      NetCoS methods paper draft + bib tooling
data/                       Reference + test datasets (BinChen2017, LINCS, repoDB, DiSignAtlas, CREEDS, UMLS cache, dictionaries) — gitignored
logs/                       cs_runs.tsv = registry of every CS run; per-step run logs; grid-search TSVs
```

## Running a new disease (case study)

A new disease from a public RNA-seq dataset is run as a case study, with the drug side already
computed (`connectivity_score/input/drug_signature_2025/all`, drugs collapsed before CS).

1. **Locations.** The slug is the folder name `src/case_studies/<slug>/`. `conf.case_study_paths(slug)`
   derives the standard locations: data in `data/<slug>/`, outputs in `other_outputs/case_studies/<slug>/`
   and `imgs/case_studies/<slug>/`, DEG CS input in `connectivity_score/input/disease_signature_2025/<slug>/`.
   `src/case_studies/<slug>/case_conf.py` adds only the dataset's file names; no case-study variables
   go into `conf.py`.
2. **Download** the counts (for GEO RNA-seq, the NCBI-recomputed `<GSE>_raw_counts_GRCh38.p13_NCBI.tsv.gz`,
   rows = NCBI GeneID) and the series matrix into `data/<slug>/`.
3. **`00_`** inspects the download (dimensions, integer counts, covariates, sample alignment) and writes
   the sample-group map; **`01_`** computes the disease signature with the R classes of
   `modules/disease_signature/` (`LogarithmGeneFilter`, upper-quartile normalisation, voom / limma,
   `topTable(coef = 2)`, as for IPF) keeping NCBI GeneIDs, and writes the DEG CS input and the MITHrIL `.mi`.
   This replaces `tsr/<disease>/` + `step1`.
4. **`conf.py`:** `DISEASE = "<slug>"`, `DISEASE_SIGNATURE_SOURCE = "cs_in"`.
5. **MITHrIL and CS:** `bash src/run_disease_pipeline.sh` runs the case-study scripts, `step2 run`, `step3`,
   `step2b_7b` and `run_cs_batch_grid.py --disease <disease_run_name>`; a stage name as argument
   restarts from that stage, `DRY_RUN=1` prints the stages only.
6. **Exploration:** `data_exploration/data_exploration_<slug>_CS_and_rankings.ipynb`, on the model of the
   IPF / ALS notebooks (section 0: the three disease inputs before CS; then CS distributions, rankings,
   agreement between representations).

## Running the pipeline

- **Full pipeline:** the numbered steps in `src/` in order (disease: step1 → step2 → step3 / step2b_7b; drug: step4 → step5 → step6 → step7 / step2b_7b, or step8 for DEG), then `cs_batch.run_cs_batch_for_conf()` (or `python cs_batch.py` with the active `conf.py`).
- **ChEMBL validation:** `src/validations/chembl/` — start at `Preprocess_binChen2017_TCGA_LINCS_data.py`.
- **R workflows:** enter the relevant `tsr/<disease>/` folder; main scripts assume that as the working directory.

> Heavy R/Julia jobs can take hours (a single LMM run especially), and MITHrIL on the full drug set needs the server — plan accordingly.
