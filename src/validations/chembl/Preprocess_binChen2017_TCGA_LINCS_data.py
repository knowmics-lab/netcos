# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 16:09:56 2026
Preprocessing of LINCS and TCGA data from
Bin Chen 2017 for 
chembl valdiation using NetCos

@author: los4
"""
from pathlib import Path

# this works with conf file related to BinChen2017 chembl validation
from conf import cell_line, DISEASE, MITH_IN_DISEASE, MITH_IN_DRUG, map_name_to_id,\
diseases_of, DATA_DIR, landmark_disease, landmark_drug, LM_flag_disease,\
    LM_flag_drug, BC_DATA, LINCS_BC_DATA, disease_run_name, cell_line_run_name, LOGS_DIR,\
        CS_IN_DISEASE, CS_IN_DRUG
from logger import append_run_metadata
import pandas as pd
import pyreadr
import pickle


from datetime import datetime
import socket

def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def write_cs_disease_input(df, out_dir, disease_run_name):
    ensure_dir(out_dir)
    out = df.rename(columns={"id": "gene_id", "log2FoldChange": "DE_log2_FC","padj":"adj.p.value"}).copy()
    out = out[["gene_id", "DE_log2_FC", "adj.p.value"]].dropna(subset=["gene_id", "DE_log2_FC"])
    out["gene_id"] = out["gene_id"].astype(str)

    out_file = out_dir / f"{disease_run_name}_signature_gene_id.txt"
    out.to_csv(out_file, sep="\t", index=False)
    return out_file, out.shape[0], out["gene_id"].nunique()


def write_cs_drug_inputs(lincs_fc, drug_md, out_dir, overwrite=False):
    """
    Write direct-FoldChange LINCS CS-input files, one pickle per drug/signature.

    Output:
        <drug_id>_signature_gene_id.pkl

    Columns follow step7 CS-input convention:
        signature_id, gene_id, gene, Perturbation, p.value, adj.p.value
    """
    ensure_dir(out_dir)

    lincs_fc = lincs_fc.copy()
    lincs_fc.index = lincs_fc.index.astype(str)

    written = []
    drug_md = drug_md.copy()
    drug_md["id"] = drug_md["id"].astype(str)

    for idx, row in drug_md.iterrows():
        drug_id = row["id"]

        if drug_id not in lincs_fc.columns:
            continue

        out_pkl = out_dir / f"{drug_id}_signature_gene_id.pkl"
        if out_pkl.exists() and not overwrite:
            print(f"skipping existing {out_pkl.name}")
            continue

        df = pd.DataFrame({
            "signature_id": drug_id,
            "gene_id": lincs_fc.index,
            "DE_log2_FC": pd.to_numeric(lincs_fc[drug_id], errors="coerce"),
            "adj.p.value": 1.0,
        }).dropna(subset=["gene_id", "DE_log2_FC"])

        df["gene_id"] = df["gene_id"].astype(str)

        with open(out_pkl, "wb") as fh:
            pickle.dump(df.reset_index(drop=True), fh)

        written.append((drug_id, out_pkl, df.shape[0], df["gene_id"].nunique()))

    return written


#%%
print(' Load LINCS signatures')
lincs_rdata = LINCS_BC_DATA / "lincs_signatures_cmpd_landmark.Rdata"
result = pyreadr.read_r(lincs_rdata) # also works for Rds
print(result.keys()) # 'lincs_signatures'

# LINCS signatures between treated and untreated samples, 
# for 978 landmark genes, 
# for 6511 is_gold drugs (at 6h and 24h timepoints)
# rows: landmark gene_id cols: LINCS perturbagen ids
LINCS_FC_bc = result['lincs_signatures']
del result
print(LINCS_FC_bc.shape)
n_LINCS_perturbagens = LINCS_FC_bc.shape[1]
n_LINCS_genes = LINCS_FC_bc.shape[0]

#%%

    
for current_cell_line, current_disease in diseases_of.items():
    
    print(current_cell_line, current_disease, 'landmark disease?', landmark_disease,  'landmark drug?', landmark_drug)
    print('load disease data')
    tcga = pd.read_excel(BC_DATA / 'SD2.xlsx', sheet_name = current_disease+'_sig.csv')
    n_tcga_total = len(tcga)
    print(tcga.shape, 'all genes')
    
    # Filter disease signature
    # by landmark genes
    
    if landmark_disease:
        filtered_tcga = tcga[tcga.landmark == landmark_disease].copy()
    else:
        filtered_tcga = tcga.copy()
    n_lm_tcga_fc = len(filtered_tcga)
    print(n_lm_tcga_fc, "genes with LINCS landmark =", landmark_disease)
    
    filtered_tcga = filtered_tcga[["id", "log2FoldChange","padj"]].copy()
    filtered_tcga["id"] = filtered_tcga["id"].apply(lambda x: x.split("|")[1])
    
    n_tcga_selected_unique_gene_ids = filtered_tcga["id"].nunique()
    print(n_tcga_selected_unique_gene_ids, 'unique gene ids')

    
    print('convert disease data to MITHrIL input')

    current_disease_run_name = current_disease + LM_flag_disease

    mith_disease_in_filename = current_disease_run_name+'_signature_gene_id.mi'
    disease_mith_input_path = MITH_IN_DISEASE/mith_disease_in_filename
    lincs_filtered_tcga = filtered_tcga[["id", "log2FoldChange"]]
    lincs_filtered_tcga.to_csv(disease_mith_input_path, sep = '\t', index=False, header=False)
    
    
    print("convert disease data to direct CS input")

    current_cs_disease_dir = CS_IN_DISEASE.parent / current_disease_run_name
    disease_cs_input_path, n_cs_disease_rows, n_cs_disease_unique_gene_ids = write_cs_disease_input(
        filtered_tcga, current_cs_disease_dir,   current_disease_run_name)
    
# # These three cancer types have
# # 5 (BRCA), 2 (LIHC) and 12 (COAD) cancer cell lines with
# # the same cell lineage (Fig. 1c and Supplementary Data 1) in LINCS data
# LINCS_cell_lines_of_disease_df = pd.read_excel(DATA_DIR/'BinChen2017'/'SD1.xlsx', sheet_name = current_disease+'_cell_lines.csv')
# LINCS_cell_lines_of_disease=LINCS_cell_lines_of_disease_df.LINCS.dropna()
# print(len(LINCS_cell_lines_of_disease),LINCS_cell_lines_of_disease)
    
    # DRUGS preprocessing
    print('load ',current_cell_line,' Binchen2017 FC data')
    # gene_id gene_symbol 
    landmark_genes = pd.read_csv(LINCS_BC_DATA / 'lincs_landmark.csv') 
    n_landmark_gene_rows = len(landmark_genes)
    n_landmark_gene_ids = landmark_genes["gene_id"].nunique() if "gene_id" in landmark_genes.columns else None
    print(landmark_genes.shape)
    # id for drug at given per time. 66511 drugs (is_gold = 1)
    # relevant columns: id  pert_iname pert_id pert_time cell_id 
    print('load perturbagen metadata')
    drug_md_all = pd.read_csv(LINCS_BC_DATA / 'lincs_sig_info_new.csv') 
    n_drug_md_all = len(drug_md_all)
    n_unique_compounds_all = drug_md_all["pert_id"].nunique() if "pert_id" in drug_md_all.columns else None
    print(drug_md_all.shape)
    print(len(drug_md_all.pert_id.unique()), 'unique compounds')
    
    # filter by cell line
    drug_md = drug_md_all[drug_md_all['cell_id']==current_cell_line].copy()
    n_drug_md_cell_line = len(drug_md)
    n_unique_compounds_cell_line = drug_md["pert_id"].nunique() if "pert_id" in drug_md.columns else None
    n_unique_signature_ids_cell_line = drug_md["id"].astype(str).nunique()
    print(drug_md.shape)

    print('filter LINCS signatures by perturbagens appearing in selected cell line')
    LINCS_FC_bc_filtered = LINCS_FC_bc[ drug_md.id.astype(str)]
    print(LINCS_FC_bc_filtered.shape)
    n_lincs_genes_filtered = LINCS_FC_bc_filtered.shape[0]
    n_lincs_cols_filtered = LINCS_FC_bc_filtered.shape[1]

    # 2418 perturbagens for HEPG2

    print('convert LINCS data to MITHrIL input. Remember to manually remove first tab from header!')

    current_cell_line_run_name = current_cell_line+LM_flag_drug
    mith_in_lincs_filename = 'LINCS_'+current_cell_line_run_name+'.mi'
    drug_mith_input_path = MITH_IN_DRUG / mith_in_lincs_filename
    LINCS_FC_bc_filtered.to_csv(drug_mith_input_path, header = True , sep = '\t', index = True)
    
    print("convert LINCS data to direct CS input")

    current_cs_drug_dir = CS_IN_DRUG.parent / current_cell_line_run_name
    drug_cs_written = write_cs_drug_inputs(LINCS_FC_bc_filtered.iloc[:,0:3], drug_md, current_cs_drug_dir)
    
    n_drug_cs_files = len(drug_cs_written)
    n_drug_cs_rows_total = sum(x[2] for x in drug_cs_written)
    
    
    # save run info
    timestamp=datetime.now().strftime("%d_%m_%Y_%H_%M") 
    preprocessing_run_id =timestamp + "__" + current_cell_line_run_name
    
    metadata_row = {
        "preprocessing_run_id": preprocessing_run_id,
        "timestamp": timestamp,
        "hostname": socket.gethostname(),
        "filter_disease_signature_by_landmark_genes_pre_mith": int(bool(landmark_disease)),
        "filter_drug_signature_by_landmark_genes_pre_mith": int(bool(landmark_drug)),
        "cell_line": current_cell_line,
        "disease": current_disease,
        "cell_line_run_name": current_cell_line_run_name,
        "disease_run_name": current_disease_run_name,
        "tcga_source_file": str(BC_DATA / "SD2.xlsx"),
        "tcga_sheet": current_disease+'_sig.csv',
        "n_tcga_genes_total": n_tcga_total,
        "n_tcga_genes_landmark": n_lm_tcga_fc,
        "n_tcga_selected_unique_gene_ids": n_tcga_selected_unique_gene_ids,
        "lincs_signatures_source_file": str(lincs_rdata),
        "lincs_metadata_source_file": str(LINCS_BC_DATA / "lincs_sig_info_new.csv"),
        "lincs_landmark_source_file": str(LINCS_BC_DATA / "lincs_landmark.csv"),
        "n_LINCS_genes": n_LINCS_genes,
        "n_LINCS_perturbagens": n_LINCS_perturbagens,
        "n_landmark_gene_rows": n_landmark_gene_rows,
        "n_landmark_gene_ids": n_landmark_gene_ids,
        "n_LINCS_drug_metadata_all_rows": n_drug_md_all,
        "n_unique_compounds_all": n_unique_compounds_all,
        "n_drug_metadata_cell_line_rows": n_drug_md_cell_line,
        "n_unique_compounds_cell_line": n_unique_compounds_cell_line,
        "n_unique_signature_ids_cell_line": n_unique_signature_ids_cell_line,
        "n_lincs_filtered_genes": n_lincs_genes_filtered,
        "n_lincs_filtered_cols": n_lincs_cols_filtered,
        "disease_mith_input_file": str(disease_mith_input_path),
        "drug_mith_input_file": str(drug_mith_input_path),
        "disease_cs_input_file": str(disease_cs_input_path),
        "drug_cs_input_dir": str(current_cs_drug_dir),
        "notes": "Pre-MITHrIL preprocessing from BinChen2017 TCGA/LINCS data; one row per disease/cell-line pair.",
    }
    
    append_run_metadata(  LOGS_DIR / "BinChen2017_chembl_validation_preprocessing_runs.tsv", metadata_row)
