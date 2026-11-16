#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis pipeline of a case study: cell-line weights, per-drug sRGES scores
and validation, the same code for every study.

The study is conf.DISEASE (= slug = folder name in src/case_studies/); its
specification is src/case_studies/<slug>/case_conf.py, which declares the
dataset-specific names. The steps run are conf.CASE_STUDY_STEPS, in this
order:

- 'weights'    case_study.weights: CCLE / LINCS cell-line similarity to the
               disease samples (case_conf.SAMPLES, NCBI_TPM) and the sRGES
               cell-line weights.
- 'srges'      case_study.per_drug_scores: per-drug sRGES_weighted,
               sRGES_unweighted and best from the per-instance CS runs.
- 'validation' every validation the study declares in its case_conf:
               REPODB       -> case_study.repodb_auroc
               CS_CHECKS    -> case_study.cs_checks
               KNOWN_DRUGS  -> case_study.known_drugs
               A field that is not declared means that validation is not run.
               CS_CHECKS runs before KNOWN_DRUGS, which also reports the
               known drugs in the 1-based rankings of cs_checks.

Before any step, validate_case_conf checks that the case_conf declares every
field the chosen steps need and stops with a readable message otherwise;
conf.CS_METHOD must be 'bin_chen'. Every figure uses the house style
(lfs_style); case_conf.FIGURES holds only figure text (FIGURES_TEXT_KEYS).
Inputs and outputs of each
step are listed in the module docstrings; tables go to
other_outputs/case_studies/<slug>/, figures to imgs/case_studies/<slug>/.

Usage
-----
    cd src && PYTHONPATH=. python case_study_pipeline.py
"""

import importlib.util
import sys

import conf
from case_study.common import apply_style

STEPS = ('weights', 'srges', 'validation')


SAMPLES_KEYS = {  # format: keys it needs besides 'format' and 'groups'
    'samples_table': ['group_column'],
    'series_matrix': ['file', 'characteristic', 'control_value'],
    'csv': ['file', 'sample_column', 'group_column', 'columns'],
}
STRATA_KEYS = ['name', 'from_column', 'prefixes', 'classes', 'label', 'short_label']
FIGURES_KEYS = {  # step or report: FIGURES text keys it needs
    'weights': ['dataset', 'disease', 'case_axis', 'control_axis', 'panel_c_title'],
    'known_drugs': ['known_drugs_label'],
}
FIGURES_TEXT_KEYS = ['dataset', 'disease', 'group_labels', 'case_axis', 'control_axis', 'panel_c_title',
                     'known_drugs_label']  # every key FIGURES may hold
REPODB_KEYS = ['cui', 'negative_classes', 'negative_sets']
CS_CHECKS_KEYS = ['cell_line']               # optional: n_check, n_jobs
CS_CHECKS_OPTIONAL = ['n_check', 'n_jobs']


def validate_case_conf(cc, steps):
    """List every missing or malformed field the chosen steps need; stop if there is any."""
    name = f'case_conf of {cc.P.slug}'
    err = []
    fig = getattr(cc, 'FIGURES', None)
    if fig is None:
        err.append('FIGURES is missing (figure text of every step)')
        fig = {}
    extra = [k for k in fig if k not in FIGURES_TEXT_KEYS]
    if extra:
        err.append(f'FIGURES holds only figure text {FIGURES_TEXT_KEYS}; unknown keys {extra} '
                   f'(figure style is the house style for every study)')

    def need(field, d, keys, why):
        miss = [k for k in keys if k not in d]
        if miss:
            err.append(f'{field} lacks {miss} ({why})')

    if 'weights' in steps:
        S = getattr(cc, 'SAMPLES', None)
        if S is None:
            err.append("SAMPLES is missing (step 'weights')")
        else:
            fmt = S.get('format')
            if fmt not in SAMPLES_KEYS:
                err.append(f"SAMPLES['format'] = {fmt!r}: expected one of {list(SAMPLES_KEYS)}")
            else:
                need('SAMPLES', S, SAMPLES_KEYS[fmt], f'format {fmt!r}')
            g = list(S.get('groups', {}))
            if len(g) != 2 or g[0] == 'control' or 'control' not in g:
                err.append(f"SAMPLES['groups'] = {g}: expected two groups, the case group first and 'control' second")
            if 'strata' in S:
                need("SAMPLES['strata']", S['strata'], STRATA_KEYS, 'sample classes')
        if not hasattr(cc, 'NCBI_TPM'):
            err.append("NCBI_TPM is missing (step 'weights')")
        need('FIGURES', fig, FIGURES_KEYS['weights'], "text of the step 'weights' panel")
    if 'validation' in steps:
        if hasattr(cc, 'REPODB'):
            need('REPODB', cc.REPODB, REPODB_KEYS, 'repoDB AUROC')
        if hasattr(cc, 'CS_CHECKS'):
            need('CS_CHECKS', cc.CS_CHECKS, CS_CHECKS_KEYS, 'CS checks')
            extra = [k for k in cc.CS_CHECKS if k not in CS_CHECKS_KEYS + CS_CHECKS_OPTIONAL]
            if extra:
                err.append(f'CS_CHECKS: unknown keys {extra} (expected {CS_CHECKS_KEYS + CS_CHECKS_OPTIONAL})')
        if hasattr(cc, 'KNOWN_DRUGS_REPORT'):
            err.append('KNOWN_DRUGS_REPORT is no longer used: every study with KNOWN_DRUGS gets the same report')
        if hasattr(cc, 'KNOWN_DRUGS'):
            if not isinstance(cc.KNOWN_DRUGS, dict) or not cc.KNOWN_DRUGS:
                err.append('KNOWN_DRUGS must be a non-empty dict {drug name: class}')
            need('FIGURES', fig, FIGURES_KEYS['known_drugs'], 'text of the known-drug panel')
    if err:
        sys.exit(f'{name}: fix these fields before running steps {steps}:\n  - ' + '\n  - '.join(err))


def load_case_conf(slug):
    path = conf.CASE_STUDIES_SRC / slug / 'case_conf.py'
    if not path.exists():
        sys.exit(f'conf.DISEASE = {slug!r}: no case-study specification {path}; '
                 f'set conf.DISEASE to the slug of a folder of src/case_studies/ that has a case_conf.py')
    spec = importlib.util.spec_from_file_location(f'case_conf_{slug}', path)
    cc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cc)
    return cc


def run_validation(cc):
    ran = []
    if hasattr(cc, 'REPODB'):
        from case_study import repodb_auroc
        repodb_auroc.run(cc)
        ran.append('repoDB AUROC')
    rank1based = None
    if hasattr(cc, 'CS_CHECKS'):
        from case_study import cs_checks
        rank1based = cs_checks.run(cc)
        ran.append('CS checks')
    if hasattr(cc, 'KNOWN_DRUGS'):
        from case_study import known_drugs
        known_drugs.run(cc, rank1based)
        ran.append('known drugs')
    print(f'validation: {", ".join(ran) if ran else "nothing declared in case_conf"}')


if __name__ == '__main__':
    assert conf.CS_METHOD == 'bin_chen', f"conf.CS_METHOD is {conf.CS_METHOD!r}; only 'bin_chen' is used"
    slug = conf.DISEASE
    cc = load_case_conf(slug)
    assert cc.P.slug == slug
    steps = list(conf.CASE_STUDY_STEPS)
    unknown = [s for s in steps if s not in STEPS]
    if unknown:
        sys.exit(f'conf.CASE_STUDY_STEPS: unknown steps {unknown}; expected some of {list(STEPS)}')
    validate_case_conf(cc, steps)
    apply_style()
    print(f'case study {slug}; steps: {", ".join(s for s in STEPS if s in steps)}')

    if 'weights' in steps:
        print('\n=== step weights')
        from case_study import weights
        weights.run(cc)
    if 'srges' in steps:
        print('\n=== step srges')
        from case_study import per_drug_scores
        per_drug_scores.run(cc)
    if 'validation' in steps:
        print('\n=== step validation')
        run_validation(cc)
