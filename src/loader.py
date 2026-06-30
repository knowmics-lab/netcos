#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 10:24:39 2025

@ author: L-F-S
"""
import os
from pathlib import Path

import pandas as pd
import pickle
from conf import CS_DIR, TSR_OUT_DRUG, TSR_OUT_DISEASE, CS_IN_DISEASE, CS_IN_DRUG,\
    CHEMBL_INPUT_DATA_DIR, BC_DISEASE_DATA, DISEASE_SIGNATURE_SOURCE, cs_log_filename
import conf as _conf  # for conf-default fallbacks in load_cs_run
from preprocessing_utils import get_drugs_list, get_signature_ids_list_from_cs_input


# ---------------------------------------------------------------------------
# NetCoS CS-run discovery / picking
# ---------------------------------------------------------------------------
#
# Shared by:
#   - src/validations/binchen2017/compare_rges_per_disease.py
#       (per-signature RGES replication: NetCoS CS vs BinChen RGES)
#   - src/validations/binchen2017/test_sRGES_replication.py
#       (per-drug sRGES replication: NetCoS CS -> sRGES vs BinChen sRGES)

def cs_out_dir_for(disease, landmark_disease=False, cs_dir=None):
    """Return the per-disease NetCoS CS output directory.

    Parameters
    ----------
    disease : str
        Disease symbol (e.g. 'LIHC', 'BRCA', 'COAD').
    landmark_disease : bool
        If True, points to '<disease>_LM' which holds the landmark-filtered
        disease-signature runs.
    cs_dir : Path | None
        Override the conf-supplied CS root. Default conf.CS_DIR.
    """
    base = Path(cs_dir) if cs_dir is not None else CS_DIR
    suffix = "_LM" if landmark_disease else ""
    return base / "output" / (disease + suffix)


def _parse_cs_run_id_ts(cs_run_id):
    """Parse the DD_MM_YYYY_HH_MM prefix of a cs_run_id into a Timestamp.

    cs_run_ids are always of the form ``DD_MM_YYYY_HH_MM_<rest>`` (see
    src/logger.py + src/conf.py::make_cs_filename), which is a more reliable
    date source than the ``timestamp`` column in cs_runs.tsv (that column
    historically mixes ``dd/mm/yyyy HH:MM`` with ISO-8601).
    """
    parts = str(cs_run_id).split("_")
    if len(parts) < 5:
        return pd.NaT
    try:
        d, m, y, h, mn = parts[:5]
        return pd.Timestamp(year=int(y), month=int(m), day=int(d),
                            hour=int(h), minute=int(mn))
    except (ValueError, TypeError):
        return pd.NaT


def pick_canonical_cs_run(
    disease,
    landmark_disease_pref=False,
    *,
    cs_log_path=None,
    cs_dir=None,
    cs_method="bin_chen",
    mith=0,
    cs_on_LM=1,
    cs_on_pathways=0,
    require_file_on_disk=True,
    verbose=True,
):
    """Pick the most recent NetCoS CS run matching the canonical Bin Chen
    RGES-replication config.

    Defaults match the canonical Bin Chen replication setup:
        CS_METHOD = 'bin_chen' (NOT 'bin_chen_disease_sorted'),
        mith = 0, cs_on_LM = 1, CS_ON_PATHWAYS = 0.

    Robustness:
      1) Sort key is the DD_MM_YYYY_HH_MM prefix parsed from the cs_run_id
         itself (via _parse_cs_run_id_ts). Falls back to the cs_runs.tsv
         ``timestamp`` column if the prefix doesn't parse.
      2) If ``require_file_on_disk`` (default), drop matching rows whose
         .tsv has been deleted from the local CS output dir — this prevents
         picking an orphan row left behind after manual cleanup. Skipped
         rows are listed when verbose=True.

    Parameters
    ----------
    disease : str
        Disease symbol (e.g. 'LIHC').
    landmark_disease_pref : bool
        Prefer the landmark-disease '<disease>_LM' rows.
    cs_log_path : Path | str | None
        Path to logs/cs_runs.tsv. Defaults to conf.cs_log_filename.
    cs_dir : Path | None
        CS output root. Defaults to conf.CS_DIR. Only used for the
        require_file_on_disk filesystem check.
    cs_method, mith, cs_on_LM, cs_on_pathways : filter values
        Pre-CS hyperparameters that identify the canonical Bin Chen run.
    require_file_on_disk : bool
        If True (default), skip rows whose .tsv is missing from disk.
    verbose : bool
        Print informative messages about which rows were skipped / picked.

    Returns
    -------
    str | None
        The filename (cs_run_id + '.tsv') of the most recent matching run,
        or None if no match is found.
    """
    log_path = Path(cs_log_path) if cs_log_path is not None else Path(cs_log_filename)
    if not log_path.exists():
        if verbose:
            print(f"[pick_canonical_cs_run] cs_runs.tsv not found at {log_path}")
        return None

    try:
        df = pd.read_csv(log_path, sep="\t")
    except Exception as e:
        if verbose:
            print(f"[pick_canonical_cs_run] could not read {log_path}: {e}")
        return None

    suffix = "_LM" if landmark_disease_pref else ""
    target_disease_run = disease + suffix
    mask = (
        (df["disease_run_id"] == target_disease_run)
        & (df["mith"].astype(int) == int(mith))
        & (df["cs_on_LM"].astype(int) == int(cs_on_LM))
        & (df["CS_ON_PATHWAYS"].astype(int) == int(cs_on_pathways))
        & (df["CS_METHOD"] == cs_method)
    )
    hits = df.loc[mask].copy()
    if hits.empty:
        if verbose:
            print(f"[pick_canonical_cs_run] no rows in {log_path.name} match "
                  f"disease_run_id={target_disease_run} mith={mith} "
                  f"cs_on_LM={cs_on_LM} CS_ON_PATHWAYS={cs_on_pathways} "
                  f"CS_METHOD={cs_method!r}")
        return None

    out_dir = cs_out_dir_for(disease, landmark_disease_pref, cs_dir=cs_dir)
    candidates, skipped = [], []
    for _, row in hits.iterrows():
        run_id = str(row["cs_run_id"])
        fname = run_id + ".tsv"
        if require_file_on_disk and not (out_dir / fname).exists():
            skipped.append(fname)
            continue
        ts = _parse_cs_run_id_ts(run_id)
        if pd.isna(ts):
            ts = pd.to_datetime(row.get("timestamp"),
                                 errors="coerce", dayfirst=True)
        candidates.append((ts, fname))

    if verbose and skipped:
        print(f"[pick_canonical_cs_run] skipped {len(skipped)} matching "
              f"cs_runs.tsv row(s) whose .tsv is missing from "
              f"{out_dir}:")
        for f in skipped:
            print(f"    {f}")

    if not candidates:
        if verbose:
            print(f"[pick_canonical_cs_run] no canonical CS .tsv exists on disk "
                  f"for disease={disease!r} landmark_disease_pref={landmark_disease_pref}. "
                  f"Either re-run cs_batch.py with the canonical config or pass "
                  f"the .tsv path explicitly.")
        return None

    candidates.sort(key=lambda x: (pd.isna(x[0]), x[0]))
    chosen_ts, chosen_fname = candidates[-1]
    if verbose:
        print(f"[pick_canonical_cs_run] picked: {chosen_fname}  (ts={chosen_ts})")
    return chosen_fname


def load_cs_run(
    *,
    # ----- filters: None means "pull current value from conf.py" -----
    disease=None,
    landmark_disease=None,
    drug_run_id=None,
    mith=None,
    cs_on_LM=None,
    CS_ON_PATHWAYS=None,
    CS_METHOD=None,
    # ----- explicit overrides (bypass the filter-based picker) -----
    selected_cs_run_id=None,
    cs_file_path=None,
    # ----- knobs -----
    cs_log_path=None,
    cs_dir=None,
    score_col=None,
    require_file_on_disk=True,
    verbose=True,
):
    """One-stop loader for a NetCoS CS .tsv produced by cs_batch.py.

    Call with **no arguments** to load the CS run that matches the
    hyperparameters currently set in src/conf.py — i.e. the run you'd produce
    by running cs_batch.py right now with the active conf::

        df, path, run_id = load_cs_run()                  # use conf.py defaults

    Common partial overrides::

        df, path, run_id = load_cs_run(CS_METHOD='bin_chen_disease_sorted')
        df, path, run_id = load_cs_run(disease='BRCA', landmark_disease=False)
        df, path, run_id = load_cs_run(selected_cs_run_id='13_05_2026_..._connectivity_score')
        df, path, run_id = load_cs_run(cs_file_path='/abs/path/to/run.tsv')

    Filter semantics:
        Each filter kwarg (disease, landmark_disease, drug_run_id, mith,
        cs_on_LM, CS_ON_PATHWAYS, CS_METHOD) accepts an explicit value, or
        ``None`` to inherit the corresponding conf.py default. Filters are
        used to query cs_runs.tsv via pick_canonical_cs_run, which:
          - sorts matches by the DD_MM_YYYY_HH_MM prefix of cs_run_id
            (more reliable than the timestamp column);
          - skips rows whose .tsv has been deleted from the local CS output
            dir (set ``require_file_on_disk=False`` to disable).
        The most recent surviving row wins.

    Override semantics (in priority order)::
        cs_file_path  >  selected_cs_run_id  >  filter-based picker

    Returns:
        (df, resolved_path, cs_run_id)
        - df            : pandas.DataFrame (full CS .tsv contents)
        - resolved_path : pathlib.Path of the file actually loaded
        - cs_run_id     : str (cs_run_id without the '.tsv' suffix)

    If ``score_col`` is supplied, the function validates it exists and
    coerces it to numeric (rows with non-numeric values are dropped).
    """
    # 1) Direct file path bypass.
    if cs_file_path is not None:
        path = Path(cs_file_path)
        if not path.exists():
            raise FileNotFoundError(f"cs_file_path does not exist: {path}")
        df = pd.read_csv(path, sep="\t")
        return _validate_cs_df(df, score_col, path), path, path.stem

    # 2) Resolve conf.py defaults for any filter left as None.
    if disease is None:
        disease = _conf.DISEASE
    if landmark_disease is None:
        landmark_disease = bool(getattr(_conf, "landmark_disease", False))
    if drug_run_id is None:
        drug_run_id = getattr(_conf, "cell_line_run_name", None)
    if mith is None:
        mith = int(getattr(_conf, "cs_mith", 0))
    if cs_on_LM is None:
        cs_on_LM = int(getattr(_conf, "cs_on_LM", 1))
    if CS_ON_PATHWAYS is None:
        CS_ON_PATHWAYS = int(bool(getattr(_conf, "CS_ON_PATHWAYS", False)))
    if CS_METHOD is None:
        CS_METHOD = str(getattr(_conf, "CS_METHOD", "bin_chen"))
    if selected_cs_run_id is None:
        selected_cs_run_id = getattr(_conf, "selected_cs_run_id", None)

    out_dir = cs_out_dir_for(disease, landmark_disease, cs_dir=cs_dir)

    # 3) Explicit cs_run_id bypass (no filter validation; just open the file).
    if selected_cs_run_id is not None:
        run_id = str(selected_cs_run_id)
        fname = run_id if run_id.endswith(".tsv") else run_id + ".tsv"
        path = out_dir / fname
        if not path.exists():
            raise FileNotFoundError(
                f"selected_cs_run_id resolves to a path that does not exist: {path}"
            )
        if verbose:
            print(f"[load_cs_run] selected_cs_run_id override -> {path}")
        df = pd.read_csv(path, sep="\t")
        bare_run_id = run_id[:-4] if run_id.endswith(".tsv") else run_id
        return _validate_cs_df(df, score_col, path), path, bare_run_id

    # 4) Filter-based pick from cs_runs.tsv.
    fname = pick_canonical_cs_run(
        disease,
        landmark_disease_pref=landmark_disease,
        cs_log_path=cs_log_path,
        cs_dir=cs_dir,
        cs_method=CS_METHOD,
        mith=mith,
        cs_on_LM=cs_on_LM,
        cs_on_pathways=CS_ON_PATHWAYS,
        require_file_on_disk=require_file_on_disk,
        verbose=verbose,
    )
    if fname is None:
        raise FileNotFoundError(
            "load_cs_run: no CS run found matching:\n"
            f"  disease={disease}  landmark_disease={landmark_disease}\n"
            f"  drug_run_id={drug_run_id}\n"
            f"  mith={mith}  cs_on_LM={cs_on_LM}  CS_ON_PATHWAYS={CS_ON_PATHWAYS}\n"
            f"  CS_METHOD={CS_METHOD!r}\n"
            f"Either change those parameters, pass selected_cs_run_id explicitly, "
            f"or pass cs_file_path explicitly."
        )
    path = out_dir / fname
    df = pd.read_csv(path, sep="\t")
    bare_run_id = fname[:-4] if fname.endswith(".tsv") else fname
    return _validate_cs_df(df, score_col, path), path, bare_run_id


def _validate_cs_df(df, score_col, source_path):
    """Helper for load_cs_run: optionally validate & numeric-coerce score_col."""
    if score_col is None:
        return df
    if score_col not in df.columns:
        raise ValueError(
            f"score_col {score_col!r} not in {source_path} columns "
            f"(available: {list(df.columns)})"
        )
    df = df.copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col]).copy()
    if df.empty:
        raise ValueError(
            f"After numeric coercion of {score_col!r} in {source_path}, "
            f"no rows remain."
        )
    return df


def standardize_disease_signature(df, source):
    """
    Convert source-specific disease DEG files to the schema expected by CS:
        gene_id, DE_log2_FC, adj.p.value
    """
    df = df.copy()

    if source == "tsr":
        rename_map = {
            "gene_id": "gene_id",
            "DE_log2_FC": "DE_log2_FC",
            "adj.p.value": "adj.p.value",
        }

    elif source == "binchen":
        rename_map = {
            "id": "gene_id",              # change to actual column name
            "logFC": "DE_log2_FC",          # change to actual DEG column name
            "adj_p_value": "adj.p.value",  # change to actual p-adjust column name
        }

    else:
        raise ValueError(f"Unknown disease signature source: {source}")

    df = df.rename(columns=rename_map)

    required = ["gene_id", "DE_log2_FC", "adj.p.value"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing standardized disease columns after loading {source}: {missing}")

    out = df[required].copy()
    out["gene_id"] = out["gene_id"].astype(str)
    out["DE_log2_FC"] = pd.to_numeric(out["DE_log2_FC"], errors="coerce")
    out["adj.p.value"] = pd.to_numeric(out["adj.p.value"], errors="coerce")
    return out.dropna(subset=["gene_id", "DE_log2_FC"])

def load_disease_signature(disease, cs_input_dir=None, mith=False, pathway=False, source=None):
    """
    Load disease signature and return a standardized dataframe.

    DEG output schema:
        gene_id, DE_log2_FC, adj.p.value

    MITHrIL output schema is left unchanged.

    cs_input_dir: optional Path. If None, falls back to module-level
        CS_IN_DISEASE (i.e. the path computed from conf.py defaults).
        The multi-conf wrapper (call_cs_batch_for_different_confs.py)
        passes its own per-combo CS_IN_DISEASE so non-default
        landmark_disease reads from the right subdir
        (LIHC_LM/ instead of LIHC/, etc.). Mirrors the cs_input_dir arg
        of load_single_signature_cs_input.
    """
    source = DISEASE_SIGNATURE_SOURCE if source is None else source
    signature = "_signature.csv" if not pathway else "_pathways.tsv"
    if cs_input_dir is None:
        cs_input_dir = CS_IN_DISEASE

    if mith:
        filename = cs_input_dir / f"{disease}_mith3{signature}"
        return pd.read_csv(filename, sep="\t")

    else:
        if source == "tsr":
            filename = TSR_OUT_DISEASE / f"{disease}_signature_gene_id.csv"
            df = pd.read_csv(filename, sep=";", decimal=",", dtype={"gene_id": "str"})
            return standardize_disease_signature(df, source="tsr")

        if source == "other":
            filename = cs_input_dir / f"{disease}_signature_gene_id.txt"
            df = pd.read_csv(filename, sep='\t')
            return standardize_disease_signature(df, source="binchen")

    raise ValueError("source must be one of: 'tsr', 'binchen'")
    
def load_disease_signature_old(DISEASE, mith=False, pathway=False):
    '''
    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    signature = "_signature.csv" if pathway == False else "_pathways.tsv"
    
    if not mith:
        filename=os.path.join(TSR_OUT_DISEASE,DISEASE+'_signature_gene_id.csv')
        return pd.read_csv(filename, sep=';',decimal=',', dtype={'gene_id':'str'})
   
    else:
        filename=os.path.join(CS_IN_DISEASE,DISEASE+'_mith3'+signature)
        return pd.read_csv(filename, sep='\t')

def load_single_drug_signature(drug, mith=False, pkl=True):
    '''
    Loads DEG drug wise drug signature (with unique gene id rows)

    inputs:
        DISEASE: str disease symbol (directory name)

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.

    returns a pandas.Dataframe where rows are genes.
    '''
    if not mith:
        
        filename=os.path.join(TSR_OUT_DRUG,'LINCS','metanalysis_drug_wise_filtered',drug+'_metanalysis')
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t',  dtype={'gene_id':'str'})
        else:
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            return data
                
    else:
        
        filename=os.path.join(CS_IN_DRUG,drug)
        
        if not pkl:
            return pd.read_csv(filename+'.csv', sep='\t')
        else:
            with open(filename+'.pkl', 'rb') as f:
                data=pickle.load(f)
            return data
           

def load_drug_signatures_old(mith=False, pkl=True):
    '''
    inputs:

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.
    returns a pandas.Dataframe which is a concatenation by row
    of all drug signature data. rows are genes signature for given drug 
    (genes will appear multiple time, once for every drug)
    '''
    
    print('loading all drug signatures...')
    
    drugs_list=get_drugs_list(mith)
    
    if not mith:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list): 
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    else:
        DEG_drug_list=[]
        for n, drug_file in enumerate(drugs_list):
            drug=drug_file.split('_')[0]
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    print(n+1, 'drugs loaded')
    return pd.concat(DEG_drug_list)

def load_all_drug_signatures(mith=False, pkl=True):
    '''
    inputs:

        mith: bool
            def. False: loading DEG FC signature. if True, loading MIthRIL 
            perturbation signature.
    returns a pandas.Dataframe which is a concatenation by row
    of all drug signature data. rows are genes signature for given drug 
    (genes will appear multiple time, once for every drug)
    '''
    
    print('loading all drug signatures...')
    
    drugs_list=get_signature_ids_list_from_cs_input(mith=mith)
    
    if not mith:
        DEG_drug_list=[]
        for n, drug in enumerate(drugs_list): 
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    else:
        DEG_drug_list=[]
        for n, drug in enumerate(drugs_list):
            DEG_drug_list.append(load_single_drug_signature(drug, mith, pkl))
    
    print(n+1, 'drugs loaded')
    return pd.concat(DEG_drug_list)

def load_single_signature_cs_input(signature_id, cs_input_dir, mith= True, pathway=False):
    if not pathway:
        if mith:
            filename = os.path.join(cs_input_dir, f"{signature_id}.pkl")
        else:
            filename = os.path.join(cs_input_dir, f"{signature_id}_signature_gene_id.pkl")
        with open(filename, "rb") as f:
            return pickle.load(f)
    else:
        filename = os.path.join(cs_input_dir, f"{signature_id}_pathways.tsv")
        return pd.read_csv(filename, sep='\t')
        
    
def load_landmark_gene_ids(path):
    landmark_genes = pd.read_csv(path / 'lincs_landmark.csv') 
    return list(landmark_genes.gene_id.sort_values())

###############################################################################
#
# Connectivity Score and Validation data
#
###############################################################################
def add_pert_id_to_cs( lincs_metadata_path,cs_df, cs_id_col='LINCS_id', metadata_id_col='id', metadata_pert_col='pert_id'):
    """
    Adds a 'pert_id' column to a CS dataframe by mapping LINCS_id via LINCS metadata.
    
    Parameters
    ----------
    cs_df : pd.DataFrame
        Connectivity score dataframe containing LINCS_id column
    lincs_metadata_path : str
        Path to lincs_sig_info_new.csv
    cs_id_col : str
        Column in cs_df (default 'LINCS_id')
    metadata_id_col : str
        Column in metadata corresponding to LINCS_id (default 'id')
    metadata_pert_col : str
        Column in metadata for pert_id (default 'pert_id')
    """

    if not metadata_pert_col in cs_df.columns:
        # load only what we need
        meta = pd.read_csv(lincs_metadata_path, usecols=[metadata_id_col, metadata_pert_col], dtype='str')
    
        # build mapping dict
        id_to_pert = dict(zip(meta[metadata_id_col], meta[metadata_pert_col]))
    
        # map
        cs_df[metadata_pert_col] = cs_df[cs_id_col].map(id_to_pert)
    
        # optional sanity check
        n_missing = cs_df['pert_id'].isna().sum()
        if n_missing > 0:
            print(f"Warning: {n_missing} LINCS_id values could not be mapped to pert_id")

    return cs_df

    
def load_drug_rankings(path, filename=None, pert_time='all',mith='mith', lincs_metadata_path =None):
    
    path= Path(path)
    
    if not filename:
        filename=mith+'_connectivity_score.tsv'
    df= pd.read_csv(path/filename, sep='\t', header=0, dtype='str')

    if lincs_metadata_path:
        df = add_pert_id_to_cs(lincs_metadata_path, df) 
    
    if not pert_time == 'all':
        return df[df.perturbation_time==pert_time]
        
    return df

def translate_cl(cell_line):
    if cell_line == 'HT29':
        return 'HT-29'
    return cell_line

# ---------------------------------------------------------------------------
# BinChen 2017 sRGES replication loaders
# Used by src/validations/binchen2017/test_sRGES_replication.py.
# All functions default to the disease / paths active in conf.py.
# ---------------------------------------------------------------------------

def binchen_disease_dir(disease=None):
    """Return the BinChen per-disease data directory (BC_DISEASE_DATA / disease)."""
    disease = disease if disease is not None else _conf.DISEASE
    return BC_DISEASE_DATA / disease


def add_lincs_metadata_to_cs(
    cs_df,
    lincs_metadata_path=None,
    *,
    cs_id_col="LINCS_id",
    metadata_id_col="id",
    extra_cols=("pert_iname", "cell_id", "pert_dose", "pert_time", "is_gold"),
):
    """Enrich a NetCoS CS dataframe with the LINCS metadata columns required
    by sRGES. Joins on LINCS_id <-> id in lincs_sig_info{_new}.csv."""
    if lincs_metadata_path is None:
        lincs_metadata_path = _conf.LINCS_METADATA_PATH
    usecols = [metadata_id_col, *extra_cols]
    meta = (pd.read_csv(lincs_metadata_path, usecols=usecols, dtype="str")
              .drop_duplicates(subset=[metadata_id_col]))
    out = cs_df.merge(meta, left_on=cs_id_col, right_on=metadata_id_col, how="left")
    if metadata_id_col in out.columns and metadata_id_col != cs_id_col:
        out = out.drop(columns=[metadata_id_col])
    if "pert_dose" in out.columns:
        out["pert_dose"] = pd.to_numeric(out["pert_dose"], errors="coerce")
    if "pert_time" in out.columns:
        out["pert_time"] = pd.to_numeric(out["pert_time"], errors="coerce")
    n_missing = out["pert_iname"].isna().sum() if "pert_iname" in out.columns else 0
    if n_missing:
        print(f"[warn] add_lincs_metadata_to_cs: {n_missing} LINCS_id rows had "
              f"no metadata match (dropped from sRGES input).")
        out = out.dropna(subset=["pert_iname", "cell_id", "pert_dose", "pert_time"])
    return out


def load_netcos_cs(cs_file_path, lincs_metadata_path=None):
    """Load a NetCoS CS .tsv and enrich it with per-signature LINCS metadata."""
    if lincs_metadata_path is None:
        lincs_metadata_path = _conf.LINCS_METADATA_PATH
    df = pd.read_csv(cs_file_path, sep="\t", dtype={"LINCS_id": str})
    if "LINCS_id" not in df.columns or "connectivity_score" not in df.columns:
        raise ValueError(
            f"NetCoS CS file missing required columns "
            f"LINCS_id/connectivity_score: {cs_file_path}"
        )
    return add_lincs_metadata_to_cs(df, lincs_metadata_path)


def load_netcos_cs_combined(cs_file_paths, lincs_metadata_path=None):
    """Concatenate multiple per-cell-line NetCoS CS .tsv files and enrich.

    BinChen's full sRGES expects one row per (drug-signature, cell), so when
    the NetCoS pipeline emits one CS file per cell line they must be stacked
    before sRGES aggregation. Per-(sig, cell) duplicates are intentionally
    preserved.
    """
    if lincs_metadata_path is None:
        lincs_metadata_path = _conf.LINCS_METADATA_PATH
    parts = []
    for p in cs_file_paths:
        d = pd.read_csv(p, sep="\t", dtype={"LINCS_id": str})
        d["__cs_source_file__"] = Path(p).name
        parts.append(d)
    if not parts:
        raise ValueError("No CS files provided")
    df = pd.concat(parts, ignore_index=True)
    return add_lincs_metadata_to_cs(df, lincs_metadata_path)


def load_binchen_all_lincs_score(disease=None):
    """Load Bin Chen's per-instance RGES file (their per-signature output).

    Note: the column is named `cmap_score` in this file, but semantically it
    is the per-instance RGES (matches `RGES` in lincs_score_1.csv exactly
    for rows that overlap). It is renamed to `RGES` on load.
    """
    disease = disease if disease is not None else _conf.DISEASE
    path = binchen_disease_dir(disease) / "all_lincs_score.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path).rename(columns={"cmap_score": "RGES"})
    df["pert_dose"] = pd.to_numeric(df["pert_dose"], errors="coerce")
    df["pert_time"] = pd.to_numeric(df["pert_time"], errors="coerce")
    return df


def load_reference_full_sRGES(disease=None):
    """Load Bin Chen's per-drug sRGES (SD6-equivalent local file)."""
    disease = disease if disease is not None else _conf.DISEASE
    path = binchen_disease_dir(disease) / "lincs_cancer_sRGES.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def load_reference_fig3_sRGES(disease=None):
    """Load Bin Chen's per-drug sRGES restricted to drugs with IC50
    (SD5 / Fig. 3)."""
    disease = disease if disease is not None else _conf.DISEASE
    path = binchen_disease_dir(disease) / "rges_ic50_normalized.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def cor_weights_for_disease(disease=None):
    """Build the lincs_cell_id -> cor map used by the sRGES_all_cmpds.R cor
    weighting. Returns None if either input file is missing on disk."""
    # Lazy import: keep sRGES.py optional for callers that never need cor.
    from sRGES import load_cell_line_weights_binchen
    disease = disease if disease is not None else _conf.DISEASE
    cl_path = binchen_disease_dir(disease) / f"cell_line_{disease}_tacle.csv"
    ccle_path = BC_DISEASE_DATA / "raw" / "cell_line_lincs_ccle.csv"
    if not cl_path.exists() or not ccle_path.exists():
        return None
    return load_cell_line_weights_binchen(str(cl_path), str(ccle_path))


# ---------------------------------------------------------------------------
# Bin Chen 2017 ChEMBL IC50 loader (used by validations/chembl/...)
# ---------------------------------------------------------------------------

def load_IC50(ic50_file, cancer_type, cell_line, IC50_ONLY=True, IC_50_binchen_SD5=True):#, median_IC50=False):
    
    if IC_50_binchen_SD5==True:
        log_dict = {}
        # 'standard_inchi' is loaded so downstream merges (e.g. the
        # CS-vs-IC50 join in final_step_correlations_IC50_CS.py) can use
        # InChI as primary join key, matching Bin Chen 2017's convention.
        # Falls back gracefully if SD8 doesn't carry it.
        sd8_cols = pd.read_excel(ic50_file, sheet_name=cancer_type, nrows=0).columns
        wanted_cols = ['pert_iname', 'pert_id', 'standard_value', 'standard_units',
                       'standard_type', 'cell_line', 'activity',
                       'standard_value_median', 'standard_inchi']
        use_cols = [c for c in wanted_cols if c in sd8_cols]
        df = pd.read_excel(ic50_file,
                           sheet_name=cancer_type, header=0,
                           usecols=use_cols)
        #filter for cell line
        cell_line = translate_cl(cell_line)
        if cell_line is not None:
            log_dict['IC50_rows'] = len(df)
            df = df[df['cell_line'].astype(str).str.upper()==cell_line.upper()]
            log_dict['IC50_rows_filtered_by_cell'] = len(df)
        
        # IC50 filtering
        if IC50_ONLY:
            df=df[df.standard_type=='IC50']
            log_dict['IC50_rows_filtered_by_IC50'] = len(df)
    
        
        return df.sort_values(by='standard_value_median'), log_dict
    
    # load downloaded IC50 data 
    # problema: mappa chemblid to pert id e poi fallo
    df=pd.read_csv(CHEMBL_INPUT_DATA_DIR/(cell_line+"_activity.tsv"),sep='t')
    return df