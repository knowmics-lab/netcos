"""
Generic case-study analysis steps, run by case_study_pipeline.py.

Each module takes the study specification (the module
case_studies/<slug>/case_conf.py, passed as `cc`) and reads every
dataset-specific name from it:

- weights.py                 step 'weights': CCLE / LINCS cell-line similarity
                             to the disease samples, sRGES cell-line weights
- per_drug_scores.py         step 'srges': per-drug sRGES / wsRGES / best from
                             the per-instance connectivity scores
- repodb_auroc.py            step 'validation', when cc.REPODB is declared
- cs_checks.py               step 'validation', when cc.CS_CHECKS is declared
- known_drugs.py             step 'validation', when cc.KNOWN_DRUGS is declared
- common.py                  house figure style and configuration colours
"""
