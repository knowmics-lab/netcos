# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 16:17:13 2025

@author: L-F-S
"""

from conf import cell_lines_chembl, DATA_DIR, alias_2geneid, LINCS_file, inst_info_file
from tsr_python.LINCS_preprocessing import map_assay_id_and_matrix_idx, main_filter
import pandas as pd
import os
import numpy as np
import h5py


#%% gpt 03/11/2025
from conf import LINCS_DIR, BING_GENES, GENE_INFO_FILE
from typing import List, Optional, Iterable, Dict, Tuple
import re 
###############################################################################
# METADATA LOADING & FILTERING
###############################################################################

def load_inst_info() -> pd.DataFrame:
    """Load LINCS inst_info metadata (samples/experiments)."""
    # inst_info is tab or comma separated depending on distribution; infer sep
    return pd.read_csv(inst_info_file, sep=None, engine='python')

def load_gene_info() -> pd.DataFrame:
    return pd.read_csv(GENE_INFO_FILE, sep="\t")

def load_bing_genes(optional_csv: Optional[str] = None) -> Optional[pd.Series]:
    fn = optional_csv if optional_csv else BING_GENES
    if not os.path.exists(fn):
        return None
    df = pd.read_csv(fn)
    # Accept either a 'gene' column or the first column
    col = "gene" if "gene" in df.columns else df.columns[0]
    return df[col].astype(str)

def _norm_unit(u: str) -> str:
    if not isinstance(u, str):
        return ""
    u = u.strip().lower().replace("µ", "u")  # µm -> um
    # normalize common synonyms to "um"
    return {"um": "um", "μm": "um", "umol": "um", "umolar": "um", "um/l": "um", "um/lit": "um",
            "um/litre": "um", "um/liter": "um", "um/ltr": "um", "umolar": "um",
            "umol/l": "um", "um/litre": "um", "um/liter": "um",
            "umolars": "um", "umolarity": "um",
            "umol·l-1": "um", "uM": "um"}.get(u, u)

def _parse_idose(s: str) -> tuple[float, str] | tuple[None, None]:
    # accept forms like "10 um", "0.04 uM", "1e-1 um"
    if not isinstance(s, str):
        return (None, None)
    m = re.match(r"^\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([a-zA-Zµu]+)\s*$", s.strip())
    if not m:
        return (None, None)
    val = float(m.group(1))
    unit = _norm_unit(m.group(2))
    return (val, unit)

def _dose_mask(df: pd.DataFrame, dose_value: str | float, dose_unit: str | None) -> pd.Series:
    """Return a boolean mask selecting rows with the requested dose."""
    if dose_value is None:
        return pd.Series(True, index=df.index)  # no filtering on dose
    # normalize target
    try:
        tgt = float(dose_value)
    except Exception:
        tgt = None
    unit_norm = _norm_unit(dose_unit) if dose_unit is not None else None

    mask_parts = []

    # 1) numeric columns (pert_dose + optional unit)
    if "pert_dose" in df.columns:
        dose_num = pd.to_numeric(df["pert_dose"], errors="coerce")
        m_num = np.isfinite(dose_num) & np.isfinite(tgt) & np.isclose(dose_num, tgt, rtol=1e-6, atol=1e-12)
        if "pert_dose_unit" in df.columns and unit_norm:
            m_unit = df["pert_dose_unit"].astype(str).map(_norm_unit) == unit_norm
            m_num = m_num & m_unit
        mask_parts.append(m_num.fillna(False))

    # 2) combined string column (pert_idose: e.g., "10 um")
    if "pert_idose" in df.columns:
        parsed = df["pert_idose"].astype(str).map(_parse_idose)
        idose_val = parsed.map(lambda x: x[0])
        idose_unit = parsed.map(lambda x: x[1])
        m_idose = idose_val.notna()
        if np.isfinite(tgt):
            m_idose = m_idose & np.isclose(idose_val.astype(float), tgt, rtol=1e-6, atol=1e-12)
        if unit_norm:
            m_idose = m_idose & (idose_unit == unit_norm)
        mask_parts.append(m_idose.fillna(False))

    if not mask_parts:
        # no dose columns at all → keep everything (and let summary guide you)
        return pd.Series(True, index=df.index)

    mask = mask_parts[0]
    for m in mask_parts[1:]:
        mask = mask | m
    return mask.fillna(False)

def _print_filter_summary(df, cell_id, times, dose_value, dose_unit, full_before_filter=None):
    print("\n[LINCS filter summary]")
    print(f" cell={cell_id}  times={tuple(map(str, times))}  dose={dose_value} {dose_unit or ''}".rstrip())
    # show pert_type composition at the final kept df
    pt_counts = df["pert_type"].value_counts(dropna=False).to_dict()
    print(f" pert_type counts (kept): {pt_counts}")
    for t in sorted({str(x) for x in df['pert_time'].astype(str)}):
        d_t = df[df['pert_time'].astype(str) == t]
        n_dmso = (d_t['pert_iname'] == 'DMSO').sum()
        d_drug = d_t[d_t['pert_type'] == 'trt_cp']
        print(f"  time={t}h:  drug rows={len(d_drug):5d}  unique_drugs={d_drug['pert_iname'].nunique():4d}  DMSO rows={n_dmso:5d}")
    # help: what doses exist among drugs BEFORE applying the dose mask (for the chosen cell/times)
    if full_before_filter is not None:
        base_drugs = full_before_filter[full_before_filter["pert_type"] == "trt_cp"].copy()
        if not base_drugs.empty:
            # Prefer pert_idose histogram if available; else show pert_dose + unit
            if "pert_idose" in base_drugs.columns:
                top = base_drugs["pert_idose"].astype(str).str.strip().value_counts().head(12)
                print("  [pre-filter] top drug pert_idose values:", ", ".join([f"{k} ({v})" for k, v in top.items()]))
            if "pert_dose" in base_drugs.columns:
                topd = base_drugs["pert_dose"].astype(str).str.strip().value_counts().head(8)
                
                mode_series = (
                base_drugs.get("pert_dose_unit", pd.Series([], dtype=object))
                .astype(str).str.lower()
                .mode()
                )
                unit_mode = mode_series.iat[0] if len(mode_series) > 0 else "NA"
                
                print(f"  [pre-filter] dose unit (mode among drugs): {unit_mode}")
                print("  [pre-filter] top numeric doses:", ", ".join([f"{k} ({v})" for k, v in topd.items()]))
    print("")

def filter_inst(
    inst: pd.DataFrame,
    cell_id: str | None = None,   # None or "*" => keep all cells
    times=("6","24"),
    dose_value="10",
    dose_unit="um",
) -> pd.DataFrame:
    """
    Robust filtering:
      • optionally filter by cell_id; if None/"*", keep all cells
      • always filter by times
      • keep drugs (trt_cp) + DMSO controls
      • apply dose filter to drugs using pert_dose/pert_dose_unit and/or pert_idose
      • print per-time × per-cell summary
    """
    df = inst.copy()

    # normalize "all cells"
    want_all_cells = (cell_id is None) or (str(cell_id).strip() == "*")

    # times filter (always)
    df = df[df.get("pert_time", "").astype(str).isin([str(t) for t in times])]

    # optional cell filter
    if not want_all_cells:
        df = df[df.get("cell_id", "").astype(str).str.upper() == str(cell_id).upper()]

    # keep a copy (for dose diagnostics before masking)
    base_subset = df.copy()

    # identify drugs vs DMSO
    pert_type = df.get("pert_type", pd.Series(index=df.index, dtype=object)).astype(str)
    pert_iname = df.get("pert_iname", pd.Series(index=df.index, dtype=object)).astype(str)

    is_drug = pert_type.eq("trt_cp")
    is_dmso = pert_type.eq("ctl_vehicle") & pert_iname.eq("DMSO")

    # dose mask (ONLY for drugs)
    m_dose = _dose_mask(df, dose_value=dose_value, dose_unit=dose_unit)
    is_drug = is_drug & m_dose

    # keep dose-matched drugs + DMSO (same cell/time already enforced)
    kept = df[is_drug | is_dmso].copy()

    # sanity: require at least one DMSO
    if not (kept.get("pert_iname", "").astype(str) == "DMSO").any():
        raise RuntimeError(
            "No DMSO controls after filtering. "
            "Try dose_value=None to inspect doses and units, then choose accordingly."
        )

    # ---- summary printing ----
    print("\n[LINCS filter summary]")
    cells_shown = "(all cells)" if want_all_cells else str(cell_id).upper()
    print(f" cells={cells_shown}  times={tuple(map(str, times))}  dose={dose_value} {dose_unit or ''}".rstrip())

    # overall pert_type composition
    if "pert_type" in kept.columns:
        print(f" pert_type counts (kept): {kept['pert_type'].value_counts(dropna=False).to_dict()}")

    # per time × per cell breakdown
    times_present = sorted({str(x) for x in kept.get("pert_time", pd.Series([], dtype=object)).astype(str)})
    cells_present = kept.get("cell_id", pd.Series([], dtype=object)).astype(str)
    for t in times_present:
        d_t = kept[kept.get("pert_time", "").astype(str) == t]
        for c in sorted(d_t.get("cell_id", pd.Series([], dtype=object)).astype(str).unique()):
            d_tc = d_t[cells_present == c]
            n_dmso = (d_tc.get("pert_iname", "").astype(str) == "DMSO").sum()
            d_drug = d_tc[d_tc.get("pert_type", "").astype(str) == "trt_cp"]
            print(f"  time={t}h  cell={c}:  drug rows={len(d_drug):5d}  unique_drugs={d_drug.get('pert_iname','').nunique():4d}  DMSO rows={n_dmso:5d}")

    # dose hints BEFORE masking (for chosen times/cells)
    base_drugs = base_subset[base_subset.get("pert_type", "").astype(str) == "trt_cp"].copy()
    if not base_drugs.empty:
        if "pert_idose" in base_drugs.columns:
            top = base_drugs["pert_idose"].astype(str).str.strip().value_counts().head(12)
            if not top.empty:
                print("  [pre-filter] top drug pert_idose values:", ", ".join([f"{k} ({v})" for k, v in top.items()]))
        if "pert_dose" in base_drugs.columns:
            topd = base_drugs["pert_dose"].astype(str).str.strip().value_counts().head(8)
            
            mode_series = (
            base_drugs.get("pert_dose_unit", pd.Series([], dtype=object))
            .astype(str).str.lower()
            .mode()
            )
            unit_mode = mode_series.iat[0] if len(mode_series) > 0 else "NA"

            
            if not topd.empty:
                print(f"  [pre-filter] dose unit (mode among drugs): {unit_mode}")
                print("  [pre-filter] top numeric doses:", ", ".join([f"{k} ({v})" for k, v in topd.items()]))

    print("")
    return kept.reset_index(drop=True)

#%% main gpt 03/11/2025
cell_id= 'HEPG2' # or None
times: Iterable[str] = ("6", "24")
dose_value: Optional[str] = "10"
dose_unit: Optional[str] = "um"
restrict_to_bing = True # temp, onoly ~1k landmark genes. easier computations
n_jobs: int = 1,

inst = load_inst_info()


inst = filter_inst(inst, cell_id=cell_id, times=times, dose_value=dose_value, dose_unit=dose_unit)


if restrict_to_bing:
    bing = load_bing_genes()
    if bing is not None:
        # gene_basename in split files might be gene symbol; if your split files use IDs, align here.
        keep = set(s.lower() for s in bing.astype(str))
    else:
        keep = None
else:
    keep = None
# genes = list_split_rds_genes()
# if keep:
#     genes = [g for g in genes if g.lower() in keep]
# if gene_whitelist:
#     wl = set(x.lower() for x in gene_whitelist)
#     genes = [g for g in genes if g.lower() in wl]

#%%


#%% load sample ids 
with h5py.File(LINCS_file, "r") as f:
    # Explore the internal structure
    print(list(f.keys()))  # Typically includes '0' or 'DATA'
    print(f.name)
    
    matrix = f['/0/DATA/0/matrix']
    
    print(list(f['0/META/COL'].keys()))
    # Read sample IDs (row IDs in the matrix)
    # A sample is a single expression profile 
    # for one cell line, exposed to one drug (or control)
    # at a given dose, for one time point.
    sample_ids = f['0/META/COL/id'][:].astype(str)
    
    gene_ids = f['0/META/ROW/id'][:].astype(str)
    print(matrix.shape)  # (1319138, 12328)
        
    # Create mapping between assay id and index in matrix
    print('Loading metadata')
    inst_df = pd.read_csv(inst_info_file, sep='\t', low_memory=False)
    id_to_index, inst_df = map_assay_id_and_matrix_idx(inst_df, sample_ids)
    
## FILTERING TESTS:
#%% filtering
# Catalano's filtering, for 3 cell lines
from conf import cell_lines_chembl
cell_lines_chembl = [cell_line.upper() for cell_line in cell_lines_chembl]
pert_times=[6, 24]
doses=[10]
min_cell_lines=0 

filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines_chembl, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 100249 matching samples.
# Found 5787 pert_ids (=perturbagen id = drug id)
len(filtered_inst_df.pert_iname.unique())
# 5301 pert_inames
#%% no filtering (should return all data). It does

cell_lines = []
pert_times=[]
doses=[]
min_cell_lines=0
 
filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 1319138 matching samples.
# Found 43526 pert_ids
len(filtered_inst_df.pert_iname.unique())
#%% 
cell_lines = ["HEPG2"]
pert_times=[6,24]
doses=[10]
min_cell_lines=0
 
filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 1393 matching samples.
# Found 122 pert_ids
#%%

cell_lines = []
pert_times=[6,24]
doses=[10]
min_cell_lines=5 

filtered_inst_df, sorted_indices, sorted_idx = \
    main_filter(inst_df,id_to_index, min_cell_lines=min_cell_lines,\
                cell_lines=cell_lines, pert_times=pert_times,\
                doses=doses, only_int=True)
# Found 370258 matching samples.
# Found 5155 pert_ids
len(filtered_inst_df.pert_iname.unique())
# but 4844 pert_inames..

# and if I change min_cell_lines filtering to work on pert_iname instead of pert_id 
# I get a bit more even, which makes sense because there
# are less pert_inames...
#%%  
metadata_cols = ['rna_plate', 'rna_well', 'pert_id', 'pert_iname', 'pert_type',
       'pert_dose', 'pert_dose_unit', 'pert_time', 'pert_time_unit', 'cell_id',
       'assay_index']
  
pert_types=['ctl_vehicle', 'trt_cp', 'trt_lig', 'ctl_vector', 'trt_sh',
       'ctl_untrt', 'trt_oe', 'trt_oe.mut']

inst_df['pert_id'][inst_df['pert_type']=='ctl_vehicle'].unique()
# array(['DMSO', 'PBS', 'H2O'])
inst_df['pert_id'][inst_df['pert_type']=='ctl_vector'].unique()
# quindi c'e anche PBS...
# ma