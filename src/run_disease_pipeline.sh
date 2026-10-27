#!/usr/bin/env bash
# =====================================================================
# Disease side + connectivity score for the disease set in src/conf.py
#
#   case-study scripts  src/case_studies/<DISEASE>/0*_*.py  (if the folder exists:
#                       download inspection, disease signature -> DEG CS input + MITHrIL .mi)
#   step2               MITHrIL on the disease signature
#   step3               gene-level MITHrIL perturbations -> disease CS input
#   step2b_7b           pathway-level MITHrIL output -> disease CS input
#                       (its drug runner also converts drug pathway files not yet converted)
#   cs_grid             run_cs_batch_grid.py --disease <disease_run_name>
#
# The disease, its run name and every path come from src/conf.py (DISEASE,
# landmark_disease, DISEASE_SIGNATURE_SOURCE); set them there first.
# A disease without src/case_studies/<DISEASE>/ (e.g. ipf, als_NYGC through tsr/ + step1)
# must already have its MITHrIL input .mi and DEG CS input; the case_study stage is skipped.
# The CS configurations are the GRID block of run_cs_batch_grid.py.
#
# Usage (from anywhere; runs in src/):
#   bash src/run_disease_pipeline.sh                 # all stages
#   bash src/run_disease_pipeline.sh step2           # from step2 on
#   bash src/run_disease_pipeline.sh cs_grid         # only the CS grid
#   DRY_RUN=1 bash src/run_disease_pipeline.sh       # print the stages, run nothing
# Stage names: case_study step2 step3 step2b_7b cs_grid
# =====================================================================

set -euo pipefail
cd "$(dirname "$0")"

STAGES=(case_study step2 step3 step2b_7b cs_grid)
FROM=${1:-case_study}
DRY_RUN=${DRY_RUN:-0}

if [[ ! " ${STAGES[*]} " =~ " ${FROM} " ]]; then
    echo "unknown stage '${FROM}'; stages: ${STAGES[*]}" >&2
    exit 1
fi

DISEASE=$(python -c "import conf; print(conf.DISEASE)")
RUN_NAME=$(python -c "import conf; print(conf.disease_run_name)")
SOURCE=$(python -c "import conf; print(conf.DISEASE_SIGNATURE_SOURCE)")
echo "conf.py: DISEASE=${DISEASE}  disease_run_name=${RUN_NAME}  DISEASE_SIGNATURE_SOURCE=${SOURCE}"

run() {   # run <stage> <command...>
    local stage=$1; shift
    echo; echo "===== ${stage}: $* ====="
    if [[ "${DRY_RUN}" == "1" ]]; then return 0; fi
    "$@"
    echo "===== ${stage} OK ====="
}

started=0
for stage in "${STAGES[@]}"; do
    [[ "${stage}" == "${FROM}" ]] && started=1
    [[ ${started} -eq 0 ]] && continue
    case ${stage} in
        case_study)
            if [[ -d "case_studies/${DISEASE}" ]]; then
                for script in case_studies/"${DISEASE}"/0[0-9]_*.py; do
                    run case_study python -u "${script}"
                done
            else
                echo; echo "===== case_study: no src/case_studies/${DISEASE}/, skipped ====="
            fi ;;
        step2)     run step2     python -u step2_disease_run_MITHrIL.py run ;;
        step3)     run step3     python -u step3_mith_out_to_cs_in_disease.py ;;
        step2b_7b) run step2b_7b python -u step2b_7b_drug_and_disease_mith_out_to_cs_in.py ;;
        cs_grid)   run cs_grid   python -u run_cs_batch_grid.py --disease "${RUN_NAME}" ;;
    esac
done
echo; echo "Done: ${RUN_NAME}. CS runs in logs/cs_runs.tsv, outputs in connectivity_score/output/${RUN_NAME}/"
