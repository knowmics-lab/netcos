# NetCoS: Network-Enhanced Connectivity Score

The network-enhanced connectivity score compares the gene signature of a specific biological process—such as transcriptomics data from patients with a particular disease—to the gene signatures of reference drugs. It leverages interactomics by projecting this data into a network space encompassing protein-protein, protein-metabolite, and protein-microRNA interactions, and propagates signals across the network. The algorithm then calculates a similarity score between each drug and the process of interest, ultimately ranking the drugs based on their similarity. A visual representation of this pipeline is shown in Figure 1.

<img width="482" alt="image" src="https://github.com/user-attachments/assets/5c572c57-f9ce-46d1-8471-3dc5fff1974f" />

**Figure 1**

### Flowchart of the Network-Enhanced Connectivity Score Algorithm

#### (a) Input and Data Pre-processing:
The algorithm requires two inputs:

- **Disease Signature** – Currently implemented for transcriptomic data [1], this includes paired transcriptomic data from healthy and diseased tissues across patients.
- **Drug Signature** – A collection of paired test-control L1000 experiments for small-molecule perturbagens across cell lines [2].

Raw data undergoes preprocessing (PP), including filtering and removing duplicate gene IDs, ensuring that only unique gene identifiers are retained for further analysis.

#### (b) Linear Mixed Model (LMM):
An LMM is applied separately to drug and disease signatures [3][4]. This model extends traditional linear models by accounting for experimental variability through random effects, specifically RNA plate and cell type. The output is a log fold change representing the difference between the condition and the control data.

#### (c) Meta-Analysis (MA) of Drug Data:
A meta-analysis [5] is performed to account for multiple perturbation time points (6 and 24 hours) and produce a consolidated measurement for each drug, per time point and pooled across both.

#### (d) Signal Propagation via MITHrIL:
The resulting fold change data from the LMM is processed using the MITHrIL pipeline [6], which propagates the signal across an extensive interaction network. This step generates a perturbation score for each network node and an accumulated perturbation for each pathway; the connectivity score can be computed on either.

#### (e) Connectivity Score Calculation:
Using the perturbation scores, a connectivity score is computed between the disease and each drug based on the Kolmogorov-Smirnov statistic [4][7]. The drugs are then ranked according to these scores. Drugs with the most negative connectivity scores are most likely to counteract the disease effect.

Step (d) is optional: the connectivity score can also be computed directly on the differential-expression signatures, which is the baseline NetCoS is compared against. See [`docs/pipeline.md`](docs/pipeline.md) for the full set of signature representations.

## References

[1] M. Prudencio, J. Humphrey, S. Pickles, A.-L. Brown, S. E. Hill, J. M. Kachergus, J. Shi, M. G. Heckman, M. R. Spiegel, C. Cook et al., "Truncated stathmin-2 is a marker of TDP-43 pathology in frontotemporal dementia," *The Journal of Clinical Investigation*, vol. 130, 2020.

[2] L. C. for Transcriptomics (Broad Institute), *L1000 Dataset -small molecule perturbagens- LINCS Trans-Center Project*, 2014.

[3] D. Bates, M. Mächler, B. Bolker, and S. Walker, "Fitting Linear Mixed-Effects Models Using lme4," *Journal of Statistical Software*, vol. 67, pp. 1–48, 2015.

[4] K. K. M. Koudijs, S. Böhringer, and H.-J. Guchelaar, "Validation of transcriptome signature reversion for drug repurposing in oncology," *Briefings in Bioinformatics*, vol. 24, p. bbac490, 2023.

[5] M. Harrer, P. Cuijpers, T. Furukawa, and D. Ebert, *Doing Meta-Analysis with R: A Hands-On Guide*, Chapman and Hall/CRC, 2021.

[6] S. Alaimo, R. Giugno, M. Acunzo, D. Veneziano, A. Ferro, and A. Pulvirenti, "Post-transcriptional knowledge in pathway analysis increases the accuracy of phenotypes classification," *Oncotarget*, vol. 7, p. 54572, 2016.

[7] B. Chen, L. Ma, H. Paik, M. Sirota, W. Wei, M.-S. Chua, S. So, and A. J. Butte, "Reversal of cancer gene expression correlates with drug efficacy and reveals therapeutic targets," *Nature Communications*, vol. 8, p. 16022, 2017.

[8] A.-L. Brown, O. G. Wilkins, M. J. Keuss, S. E. Hill, M. Zanovello, W. C. Lee, A. Bampton, F. C. Y. Lee, L. Masino, Y. A. Qi et al., "TDP-43 loss and ALS-risk SNPs drive mis-splicing and depletion of UNC13A," *Nature*, vol. 603, pp. 131–137, 2022.

## Tech stack

- **Python** — pipeline orchestration, connectivity score, validations, plotting. Dependencies in `requirements.txt` (`joblib`, `matplotlib`, `numpy`, `pandas`, `pyreadr`, `requests`, `scipy`); package metadata in `setup.py` (`drug_repurposing`, v0.1.0).
- **R** — LMM fitting, LINCS data handling, signature meta-analysis, edgeR for count-based case studies. Under `tsr/` (legacy workflows) and `new_tsr/R/` (refactored OOP version).
- **Julia** — LMM optimization (`tsr/julia/`, `modules/drug_signature/julia/`).
- **MITHrIL** — external Java tool, run as a subprocess (`step2_disease_run_MITHrIL.py`, `step6_drug_run_MITHrIL.py`).

## Quickstart

Experiments are switched by which `conf.py` is active in `src/` (there is no CLI flag); the CS entry point is `cs_batch.run_cs_batch_for_conf()`, and every run is logged in `logs/cs_runs.tsv`. Heavy R/Julia jobs can take hours.

## Documentation

- **[`docs/pipeline.md`](docs/pipeline.md)** — precise stage-by-stage pipeline & architecture reference: the two pipelines (TSR/LMM and BinChen2017/ChEMBL) and the drug-collapse inversion, the signature representations (DEG / MITHrIL / pathways), the numbered steps, the CS methods and sign convention, disease-signature sources, the CS run registry, the config system, the analyses built on the pipeline, the repository layout and how to run it.
