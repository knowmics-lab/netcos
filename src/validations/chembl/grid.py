# -*- coding: utf-8 -*-
"""
Shared hyperparameter grid for the ChEMBL validation protocol.

Two grids live here:

1. Pre-MITHrIL and pre-CS grid
   The same combos that call_cs_batch_for_different_confs.py sweeps when
   computing the connectivity score. Each upstream combo defines a unique
   row in logs/cs_runs.tsv and a unique CS output file.

2. Validation hyperparameter grid
   The combos that hyperparameter_grid_search_bootstrap.py sweeps on top of
   each available CS file. These are the post-CS knobs from conf.py
   (drug-collapse methods, IC50 filtering, cell-line filtering, thresholds).

Both wrappers import from this module so the grid stays in one place.

How to add a new validation hyperparameter
------------------------------------------
Example: turning ic50_threshold (currently fixed at 10.0) into a hyperparameter.

  1. Add a list here, e.g. ``ic50_thresholds = [10.0, 1.0, 100.0]``.
  2. In hyperparameter_grid_search_bootstrap.py:
     a. Import ``ic50_thresholds`` from this module.
     b. Add it to the ``itertools.product(...)`` call that builds
        ``downstream_combos``.
     c. Unpack the new element inside ``build_eval_configs`` and pass it to
        ``EvalConfig(ic50_threshold=...)``. The field is already defined on
        EvalConfig, so no schema change is needed for thresholds.
  3. (Optional) If the new param has cross-param constraints (e.g. only
     valid for certain CS_METHODs), add the check to ``validate_hyperparameters``
     in conf.py and have ``build_eval_configs`` skip combos that fail it.

For non-threshold params, also add the field to EvalConfig and to the
``combo_summary`` dict so it appears in the summary TSV.

@author: L-F-S
"""

# ===================== Pre-MITHrIL and pre-CS =====================
# Mirrors the conf.py hyperparameters. Single source of truth for the CS sweep.
cell_lines        = ['HEPG2', 'MCF7', 'HT29']
landmark_diseases = [False, True]
landmark_drugs    = [True]                                    # BinChen drug data is already LM only
CS_METHODs        = [ 'bin_chen_disease_sorted']#'bin_chen']#
cs_miths          = [0, 1]
cs_on_LMs         = [0, 1]
CS_ON_PATHWAYSs   = [False, True]


# ==================== Validation hyperparameters ===================
# Validation hyperparameters evaluated by hyperparameter_grid_search_bootstrap.py on top
# of each CS file. Names mirror the conf.py validation parameters.
cs_drug_collapse_methods   = ['srges','best']#['best', 'median', 'srges']
ic50_drug_collapse_methods = ['best', 'median', 'srges']
ic50_onlys                 = [True]#, False]

# Connectivity-score threshold for the "predicted positive" classification.
# NOTE: -1.5 is calibrated for CS_METHOD='bin_chen_disease_sorted'. For
# CS_METHOD='bin_chen' (RGES-style scores), a different threshold is required;
# conf.validate_hyperparameters refuses (bin_chen, -1.5) so those grid
# combos get filtered out. Run RGES_study.py to find a reasonable bin_chen
# threshold and add it here.
cs_thresholds              = [-1.5]

# 'per_cell_line' is a sentinel meaning "filter IC50 to the cell line of the
# current upstream combo". None means "use all cell lines in the IC50 file".
# The validation wrapper resolves 'per_cell_line' to the upstream cell_line
# at iteration time.
cell_line_IC50_options     = ['per_cell_line', None]
