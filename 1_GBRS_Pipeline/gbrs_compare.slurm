#!/bin/bash
# Compare this run's results with another result for the same libraries, sample by sample
# (scripts/compare_gbrs_runs.py; docs/validation.md, Comparing with earlier results).
#
#   ./run.sh compare --rule=either --other='/old/run/{sample}/{sample}.diploid.genes.expected_read_counts' --name=old
#SBATCH --job-name=gbrs_compare
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurmout/gbrs_compare_%j.out
set -euo pipefail
umask 007
HERE="${PIPELINE_ROOT:?submit through ./run.sh compare}"
source "$HERE/config.sh"
: "${COMPARE_RULE:?}" "${COMPARE_OTHER:?}" "${COMPARE_NAME:?}"
binds=$(printf '%s\n' "$WORKDIR" "$HERE" "$(dirname "${COMPARE_OTHER%%\{sample\}*}x")" | sort -u | paste -sd, -)
args=(--sheet "$SAMPLE_SHEET" --a "$GBRS_OUT/{sample}/$COMPARE_RULE/{sample}.diploid.genes.expected_read_counts"
      --b "$COMPARE_OTHER" --out "$WORKDIR/comparison/${COMPARE_RULE}_vs_${COMPARE_NAME}.tsv")
if [[ -n "$GBRS_SIF" ]]; then
    load_software "$MODULE_APPTAINER"
    apptainer exec --cleanenv --bind "$binds" "$GBRS_SIF" python3 "$HERE/scripts/compare_gbrs_runs.py" "${args[@]}"
else
    load_software "$MODULE_PYTHON"
    python3 "$HERE/scripts/compare_gbrs_runs.py" "${args[@]}"
fi
