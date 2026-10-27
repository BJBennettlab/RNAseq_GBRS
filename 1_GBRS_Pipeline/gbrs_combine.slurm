#!/bin/bash
# Join every finished sample into gene-by-sample tables, one folder per mate rule:
# WORKDIR/combined/<rule>/ (counts, diplotypes, founder counts, qc.tsv, SHA256SUMS.txt). Refuses unfinished or
# changed samples and reads every written cell back (scripts/combine_gbrs.py).
#
#   ./run.sh combine      (./run.sh all queues it after the samples)
#SBATCH --job-name=gbrs_combine
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurmout/gbrs_combine_%j.out
set -euo pipefail
umask 007
HERE="${PIPELINE_ROOT:?submit through ./run.sh combine}"
source "$HERE/config.sh"
load_software "$MODULE_PYTHON"
for rule in $MATE_RULES; do
    echo "combine $rule $(date -Is)"
    python3 "$HERE/scripts/combine_gbrs.py" "$SAMPLE_SHEET" "$GBRS_OUT" "$WORKDIR/combined/$rule" "$rule"
done
echo "combine done $(date -Is)"
