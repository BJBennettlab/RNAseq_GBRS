#!/bin/bash
# Check every mate rule of the run (scripts/validate_gbrs.py; docs/validation.md): finished samples, sex from
# expression, duplicate genotypes and, when DNA_ALLELEPROBS or STAR_COUNTS is set, identity against the mice's DNA
# genotypes and against a STAR count matrix. Reports in WORKDIR/validation/<rule>/; the job fails when a check fails.
#
#   ./run.sh validate     (./run.sh all queues it after combine)
#SBATCH --job-name=gbrs_validate
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurmout/gbrs_validate_%j.out
set -euo pipefail
umask 007
HERE="${PIPELINE_ROOT:?submit through ./run.sh validate}"
source "$HERE/config.sh"

# numpy comes from the GBRS container, or from python3 on PATH without one.
binds=$(printf '%s\n' "$WORKDIR" "$REF_DIR" "$HERE" "$(dirname "$SAMPLE_LIST")" \
    ${DNA_ALLELEPROBS:+"$(dirname "$DNA_ALLELEPROBS")"} ${STAR_COUNTS:+"$(dirname "$STAR_COUNTS")"} | sort -u | paste -sd, -)
py() {
    if [[ -n "$GBRS_SIF" ]]; then
        apptainer exec --cleanenv --bind "$binds" "$GBRS_SIF" python3 "$@"
    else
        python3 "$@"
    fi
}
if [[ -n "$GBRS_SIF" ]]; then load_software "$MODULE_APPTAINER"; else load_software "$MODULE_PYTHON"; fi

status=0
for rule in $MATE_RULES; do
    echo "validate $rule $(date -Is)"
    py "$HERE/scripts/validate_gbrs.py" --sheet "$SAMPLE_SHEET" --samples "$GBRS_OUT" --rule "$rule" \
        --grid "$REF_DIR/ref.genome_grid.GRCm39.tsv" --out "$WORKDIR/validation/$rule" --sample-list "$SAMPLE_LIST" \
        --same-mouse "$SAME_MOUSE_SCORE" ${DNA_ALLELEPROBS:+--dna "$DNA_ALLELEPROBS"} \
        ${STAR_COUNTS:+--star "$STAR_COUNTS"} || status=1
done
echo "validate done $(date -Is), status $status"
exit $status
