#!/bin/bash
# GBRS workflow entry point. Settings: run_settings.sh (copy from examples/).
#
#   ./run.sh check        settings, software and reference check
#   ./run.sh containers   pull the pinned GBRS container to GBRS_SIF
#   ./run.sh references   download and verify the GBRS reference into REF_DIR
#   ./run.sh sheet        check SAMPLE_LIST against READS_DIR, write samples.tsv
#                         (add --allow-unlisted to run a subset, e.g. a pilot)
#   ./run.sh submit       one Slurm array task per sample (--array=LIST for chosen rows,
#                         --dependency=afterok:JOB to wait for another job)
#   ./run.sh status       finished, running and failed samples
#   ./run.sh combine      join finished samples into gene-by-sample tables, one folder per mate rule
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_ROOT="$HERE"
source "$HERE/config.sh"
cmd=${1:-}; shift || true

need_python() { load_software "$MODULE_PYTHON"; }
rows() { [[ -s "$SAMPLE_SHEET" ]] || { echo "Run ./run.sh sheet first" >&2; exit 1; }; echo $(( $(wc -l < "$SAMPLE_SHEET") - 1 )); }

case "$cmd" in
check)
    fail=0
    [[ "$WORKDIR" != *CHANGEME* ]] || { echo "Set WORKDIR"; fail=1; }
    [[ "$SLURM_ACCOUNT" != CHANGEME ]] || { echo "Set SLURM_ACCOUNT"; fail=1; }
    [[ -d "$READS_DIR" ]] || { echo "READS_DIR not found: $READS_DIR"; fail=1; }
    [[ -f "$SAMPLE_LIST" ]] || { echo "SAMPLE_LIST not found: $SAMPLE_LIST"; fail=1; }
    [[ -f "$REF_DIR/reference_manifest.json" ]] || { echo "No verified reference in $REF_DIR (./run.sh references)"; fail=1; }
    if [[ -n "$GBRS_SIF" ]]; then
        [[ -s "$GBRS_SIF" ]] || { echo "Container missing: $GBRS_SIF (./run.sh containers)"; fail=1; }
    else
        for t in gbrs emase bowtie samtools; do command -v $t >/dev/null || { echo "$t not on PATH"; fail=1; }; done
    fi
    case "$WORKDIR" in "$READS_DIR"*) echo "WORKDIR must not be inside READS_DIR"; fail=1;; esac
    (( fail == 0 )) && echo "Settings look complete" || exit 1
    ;;
containers)
    load_software "$MODULE_APPTAINER"
    [[ -e "$GBRS_SIF" ]] && { echo "$GBRS_SIF exists"; exit 0; }
    mkdir -p "$(dirname "$GBRS_SIF")"
    apptainer pull --disable-cache "$GBRS_SIF" "$GBRS_IMAGE"
    sha256sum "$GBRS_SIF" | tee "$GBRS_SIF.sha256"
    ;;
references)
    need_python; python3 "$HERE/scripts/fetch_references.py" "$REF_DIR" ;;
sheet)
    need_python; mkdir -p "$WORKDIR"
    python3 "$HERE/scripts/sample_sheet.py" "$SAMPLE_LIST" "$READS_DIR" "$SAMPLE_SHEET" "$@" ;;
submit)
    n=$(rows); array="1-$n%$THROTTLE"; extra=()
    for a in "$@"; do
        case "$a" in
            --array=*) array="${a#--array=}" ;;
            --dependency=*|--begin=*|--qos=*|--nice=*) extra+=("$a") ;;  # passed to sbatch
            *) echo "Unknown submit option: $a" >&2; exit 1 ;;
        esac
    done
    if [[ -s "$WORKDIR/jobs.tsv" && "$array" == "1-$n%$THROTTLE" ]]; then
        echo "A full run is already recorded in $WORKDIR/jobs.tsv; use --array=LIST to rerun chosen rows" >&2; exit 1
    fi
    mkdir -p "$WORKDIR/slurmout" "$GBRS_OUT"
    cd "$WORKDIR"
    job=$(sbatch --parsable --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" --cpus-per-task="$CPUS" \
        --mem="$MEMORY" --time="$TIME_LIMIT" --array="$array" --export=ALL,PIPELINE_ROOT="$HERE" \
        "${extra[@]}" "$HERE/1_GBRS_Pipeline/gbrs_sample.slurm")
    printf '%s\t%s\t%s\n' "$job" "$array" "$(date -Is)" >> "$WORKDIR/jobs.tsv"
    echo "Submitted $job for rows $array"
    ;;
status)
    n=$(rows); done_n=$(find "$GBRS_OUT" -mindepth 2 -maxdepth 2 -name done.json 2>/dev/null | wc -l)
    partial=$(find "$GBRS_OUT" -mindepth 1 -maxdepth 1 -name '*.partial' 2>/dev/null | wc -l)
    echo "samples $n, finished $done_n, unfinished folders $partial"
    [[ -s "$WORKDIR/jobs.tsv" ]] && cut -f1 "$WORKDIR/jobs.tsv" | paste -sd, | xargs -r -I{} sacct -j {} -X -n -o State | sort | uniq -c
    ;;
combine)
    need_python
    for rule in $MATE_RULES; do
        python3 "$HERE/scripts/combine_gbrs.py" "$SAMPLE_SHEET" "$GBRS_OUT" "$WORKDIR/combined/$rule" "$rule"
    done ;;
*)
    sed -n '2,12p' "$0"; exit 1 ;;
esac
