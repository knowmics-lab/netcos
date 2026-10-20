# Bozza deposito Zenodo — NetCoS v0.1.0

**Pubblicato il 12/10/2026 — DOI 10.5281/zenodo.23318730 (https://doi.org/10.5281/zenodo.23318730).** I file in `upload/` sono quelli caricati.

Documento di lavoro (non va caricato su Zenodo). Contiene: cosa caricare, i campi del modulo pronti da incollare, cosa resta da decidere prima di pubblicare, cosa c'è e cosa non c'è nell'archivio.

## 1. Come funziona (in breve)

- **Prima della pubblicazione** la bozza su Zenodo è modificabile in tutto. Dal modulo si può prenotare il DOI ("Get a DOI now!"), ma un DOI prenotato **non risolve** finché il record non è pubblicato: per la proposal ERC conviene pubblicare prima della sottomissione, non citare un DOI di bozza.
- **Dopo la pubblicazione** i metadati (autori, titolo, descrizione, affiliazioni) si possono modificare quando vuoi, e il DOI resta lo stesso. I **file** invece no: per cambiarli si carica una *nuova versione*, che ha un suo DOI; un DOI "concept" punta sempre all'ultima versione. Un record si può cancellare solo nei primi 30 giorni.
- Quindi: gli autori si possono sistemare anche dopo; il contenuto dello zip e del PDF va fissato prima di premere *Publish*.

## 2. File da caricare

Cartella `upload/`:

| File | Contenuto |
|---|---|
| `netcos-v0.1.0.zip` | snapshot del codice (319 file, 424 KB): `src/`, `modules/`, `new_tsr/`, `tsr/` (solo codice), `configs/` (snapshot), `docs/pipeline.md`, `README.md`, `LICENSE`, `CITATION.cff`, `requirements.txt`, `setup.py`, file di ambiente conda, e `results/` (3 tabelle che rendono verificabili i run citati nel report) |
| `NetCoS_technical_report_v0.1.pdf` | report tecnico, 7 pagine, 3 figure (sorgente in `report/`, figure in `report/imgs/`) |

Il report si rigenera con `pandoc NetCoS_technical_report_v0.1.md -o NetCoS_technical_report_v0.1.pdf --pdf-engine=xelatex` dalla cartella `report/`.

## 3. Campi del modulo Zenodo

- **Resource type:** Software. (Alternativa: due record separati, il codice come *Software* e il report come *Publication → Report*, collegati con *Related works*. Con un record solo è più semplice; con due il report è citabile come testo a sé.)
- **Title:** `NetCoS: Network-enhanced Connectivity Score for transcriptome-based drug repurposing`
- **Version:** `0.1.0`
- **Publication date:** data di pubblicazione
- **Creators** (in quest'ordine): Signorini, Lorenzo Federico — University of Catania; Catalano, C. — University of Catania; Alaimo, S. — University of Catania; Shasha, D. — New York University; Pulvirenti, A. — University of Catania. Nel modulo conviene inserirli con nome completo e ORCID (il campo cerca l'ORCID e compila il nome); nel PDF e in `CITATION.cff` sono con le iniziali.
- **Description** (da incollare):

  > NetCoS ranks small molecules by how strongly their transcriptional signature reverses a query signature, such as a disease-versus-control contrast. It extends the gene-level connectivity score of transcriptome signature reversion (TSR) with an optional network-propagation step: query and drug signatures are propagated over a directed molecular interaction network with MITHrIL before scoring, and the score can be computed on genes, network nodes or pathways. Drug signatures are derived from LINCS L1000 with a linear mixed model and time-point meta-analysis.
  >
  > This record contains the source code (Python, R, Julia) at version 0.1.0 and a technical report describing the method and the results obtained so far on the public data of Chen et al. (Nat Commun 2017): a replication of the published RGES scores and an exploratory comparison of gene-level and network-propagated scores against drug potency (IC50) in three cancer cell lines. Input data (LINCS L1000, GEO, TCGA, ChEMBL) are not redistributed.

- **Licenses:** MIT (codice). Aggiungere CC BY 4.0 per il report (il nuovo modulo Zenodo accetta più licenze; il PDF lo dichiara già).
- **Keywords:** drug repurposing; connectivity score; transcriptome signature reversion; LINCS L1000; network propagation; MITHrIL; RGES
- **Related works:** `https://github.com/knowmics-lab/netcos` — *is supplemented by* (o *is derived from*), tipo Software.
- **References:** le 8 voci della bibliografia del report.
- **Funding:** TODO (assegno di ricerca / progetto che lo finanzia, se c'è un grant ID).
- **Language:** English.

Dopo la pubblicazione: mettere il DOI in `CITATION.cff` (campo `identifiers`, già predisposto) e nel `README.md` del repo GitHub.

## 4. Da decidere / sistemare prima di pubblicare

**Bloccanti**

1. **ORCID.** Autori decisi (Signorini, Catalano, Alaimo, Shasha, Pulvirenti) e già nel PDF e in `CITATION.cff`; mancano gli ORCID, che si aggiungono anche dopo la pubblicazione nel modulo Zenodo (in `CITATION.cff` invece entrano solo con una nuova versione dello zip).
2. **Data e versione.** `date-released` in `CITATION.cff` è `2026-10-XX`.
3. **Allineare il repo allo zip.** Nello zip ho corretto due commenti che violano la regola delle docstring (nomi/date/riferimenti al processo di revisione); il codice eseguito è identico. Va applicata la stessa modifica nel repo, così lo zip coincide con un tag git (`v0.1.0`):
   - `src/validations/drug_disease_pairs_validation/03_repodb_signature_intersection.py`, riga 38: `Two views, as asked by Pulvirenti (22/09): the pooled` → `Two views: the pooled`
   - `src/validations/drug_disease_pairs_validation/04_repodb_lincs_drug_mapping.py`, righe 122–123: `Every entry was adjudicated by hand (Critic,` / `# 24/09/2026); anything not listed` → `Every entry was adjudicated by hand;` / `# anything not listed`
   - `05_umls_hierarchy_rollup.py` è incluso nello zip ma non è ancora committato.

**Da valutare (non bloccanti)**

4. **Case study delle collaborazioni (consiglio l'opzione b)** (`src/case_studies/`, `configs/28_09_2026_*`): esclusi dallo zip perché riguardano dati non pubblicati di collaboratori e un record Zenodo è permanente. Però `src/conf.py` ha come query attiva `DISEASE = "mir29_antagomir"` e definisce `MIR29_DIR = DATA_DIR / 'Pizzoruso_mir29'`, e `docs/pipeline.md` descrive i due case study. Anche la docstring di `src/conf.py` (righe 8–11) descrive il run miR-29 (RNA-seq di cervello di topo, ortologhi MGI), quindi la configurazione di default dell'archivio è quella del case study non pubblicato, che nessun risultato del report usa. Opzioni: (a) lasciare così (sono già pubblici su GitHub); (b) prima del tag, rimettere in `conf.py` una query neutra (es. IPF) e togliere i due paragrafi da `docs/pipeline.md`; (c) includere i case study dopo averne parlato con i collaboratori.
5. **Percorsi assoluti** nel codice: `MITH_APP="/home/signorini/mithril3/..."` in `src/conf.py` e negli snapshot di `configs/`; `setwd(...)` in `modules/LINCS_export/LINCS_export.R` e `src/map_metanalysis_gene_wise_to_drug_wise.R`; `prefix: C:\Users\los4\...` in fondo a `drug_repurposing.yaml`. Non sono segreti, ma rendono lo snapshot meno portabile.
6. **README.** Il Quickstart non spiega come creare `src/local.py` (escluso perché contiene i percorsi locali) né dove scaricare i dati. Basterebbero 5 righe.
7. **Commenti `#%% gpt 03/11/2025`** in `src/step_data_preprocessing_LINCS_database_filtering.py` (e copia in `configs/`): tracce di uno strumento AI con data; meglio toglierli.
8. **Commenti con nomi** (`# using signorini's mapping`, `# using alaimo's mapping`, `# Catalano dataset`) in `src/conf.py` e negli snapshot: sono attribuzioni di metodo, li lascerei; li segnalo per completezza.

**Fuori dal deposito ma da sapere**

9. **Il repo pubblico su GitHub traccia ancora `sources/`** (22 file: PDF di articoli sotto copyright, due tesi, supplementari), nonostante `/sources/` sia in `.gitignore`: erano stati committati prima. Per questo **non** conviene usare l'integrazione GitHub→Zenodo (archivierebbe anche quei PDF): meglio il caricamento manuale dello zip. Per toglierli da GitHub serve `git rm -r --cached sources/` (e, per toglierli anche dalla storia, una riscrittura della storia). Stesso discorso, più lieve, per `validations/chembl/chembl_input/` (dati ChEMBL, licenza CC BY-SA) e `unict 2024-25.lnk`.
10. Aggiungere `/zenodo_draft/` a `.gitignore` (contiene una copia del codice).

## 5. Manifest dello zip

**Incluso (319 file):** tutti i file di codice tracciati in git sotto `src/` (incluse `validations/binchen2017`, `chembl`, `drug_disease_pairs_validation` 01_–05_, `known_drug_benchmark`), `modules/`, `new_tsr/R/`, `tsr/` (script R e Julia), `configs/` (snapshot BinChen 2017 e main pipeline 18/11/2025), `LINCS_drug_signatures_parallel_computation.R`, `docs/pipeline.md`, `README.md`, `LICENSE`, `requirements.txt`, `setup.py`, `server_env.yaml`, `drug_repurposing.yaml`, più `CITATION.cff` (nuovo) e `results/` (nuovo: `cs_runs_cited.tsv` = righe del registro per i 30 run citati, senza host e percorsi; `ic50_correlation_grid.tsv` = le 240 correlazioni IC50 della §3.2; `rges_replication_summary.tsv` = i tre summary della §3.1).

**Escluso, e perché:**

| Escluso | Motivo |
|---|---|
| `sources/` | PDF e documenti di terzi (copyright) |
| `data/`, `tsr/output/`, `validations/chembl/chembl_input/`, `validations/octad tests/`, `modules/**/*.top.table.tsv` | dati (il deposito è solo codice) |
| `src/case_studies/`, `configs/28_09_2026_*` | collaborazioni con dati non pubblicati (vedi §4.4) |
| `src/local.py`, `configs/**/local.py` | percorsi della macchina locale |
| `notebooks/`, `data_exploration/` | notebook esplorativi con output incorporati |
| `src/older_versions.zip`, `src/validations/chembl/old version/`, `Untitled.ipynb`, `tsr/ipf_old/`, `modules/drug_signature.zip` | versioni superate o vuote |
| `.Rhistory`, `*.bmp`, `*.lnk`, `.gitattributes`, `.gitignore`, `__pycache__/`, `*.egg-info/` | file di sistema/IDE |
| `CLAUDE.md`, `instructions/`, `ai-team-skills/`, `admin/`, `Claude outputs/` | non fanno parte della pipeline |
| `paper/`, `presentations/`, `imgs/`, `other_outputs/`, `logs/` (completi) | documenti e output; il report breve e `results/` li sostituiscono |

## 6. Revisione

Un round di Critic (verdetto *trust-with-caveats*): tutti i numeri ricalcolati dai file sorgente coincidono. Recepite le osservazioni: analisi di sensibilità con le IC50 di tutte le linee cellulari (in HepG2 l'ordine gene-level/MITHrIL non è stabile), asimmetria 24 vs 8 configurazioni, confronto con Chen et al. reso omogeneo (quadrati in Fig. 3), metodo degli IC citato, `results/` aggiunto all'archivio, Fig. 1 ripulita e didascalia corretta, dettagli su `srges`, coppie in conflitto e soglia ≥10 casi e ≥10 controlli. Non fatto: bootstrap appaiato della differenza tra correlazioni (lo script `hyperparameter_grid_search_bootstrap.py` esiste; sarebbe il passo successivo per una versione 0.2).

## 7. Provenienza dei numeri del report

| Affermazione | Fonte verificata |
|---|---|
| Tab. 2, Fig. 2 (RGES: r, ρ, n; firme che deviano > 0.01: 1/100, 10/1 163, 2/244) | `other_outputs/binchen2017_validation/*__vs__lincs_score_1__*__{summary,joined}.tsv` dei run `13_05_2026_11_47_DEG_LM_bin_chen`, `13_05_2026_18_28_DEG_LM_bin_chen`, `14_05_2026_07_03_DEG_LM_bin_chen` (presenti in `logs/cs_runs.tsv`) |
| §3.2: 0.56 / 0.33 / 0.41 (DEG_LM, `srges`, IC50 mediana) | `logs/BinChen2017_chembl_validation_IC50_correlation_runs.tsv`, run `13_05_2026_13_46`, `14_05_2026_02_01`, `14_05_2026_09_57` (DEG_LM) |
| Tab. 3, Fig. 3 (massimo e mediana su 8 configurazioni × 5 rappresentazioni × 3 linee) | stesso registro, righe con run CS dal 13/05/2026 e compound 24/100/58; ricalcolato per questo report |
| Valori di riferimento 0.58 / 0.33 / 0.41 e 24 / 100 / 58 composti | Chen et al. 2017, testo e Fig. 3 (controllato su `sources/BinChen2017.pdf`) |
| Sensibilità alle IC50 di tutte le linee (35 / 103 / 68 composti; HepG2 0.59 vs 0.60, mediane 0.51 vs 0.42/0.42/0.13; MCF7 0.57 vs 0.43; HT-29 0.45 vs 0.40) | stesso registro IC50, righe con 35/103/68 composti |
| IC 95% | calcolati qui (Fisher z, SE = 1.03/√(n−3)); non tengono conto della scelta del massimo |
| §4: 1 578 / 1 141 coppie, 178 malattie; 779 / 739, 157 malattie, 620 farmaci, 65 con entrambe | `other_outputs/validations/drug_disease_pairs_validation/repodb_signature_intersection_summary.tsv`, `repodb_signature_lincs_summary.tsv` (riga `min_psamples = 10`) |
| Metodo (§1–2) | `paper/NetCoS_report_draft.md` §2.1–2.6 e codice |

Nessun risultato da run `bin_chen` precedenti all'11/05/2026.

**Differenza dal report di lavoro.** `paper/NetCoS_report_draft.md` §3.4.4 dice che la griglia include la regola `median` per le istanze; nel registro IC50 le configurazioni sono solo `best` e `srges` (8 per rappresentazione). Il report breve descrive la griglia del registro; nel report di lavoro va corretto il testo. Tabella 3.11 del report di lavoro indica per HepG2/mith_pw `best`/`median`: nel registro `best`/`best` e `best`/`median` danno lo stesso 0.68, quindi è coerente.
