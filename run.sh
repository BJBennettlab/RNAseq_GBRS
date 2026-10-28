#!/bin/bash
# GBRS workflow entry point. Settings: run_settings.sh (copy from examples/), or a file named by PIPELINE_CONFIG.
# Everything that reads data or computes runs as a Slurm job; this script checks settings, lists files and submits.
#
#   ./run.sh sheet      match SAMPLE_LIST to its reads by exact name and write WORKDIR/samples.tsv
#                       (--allow-unlisted to run a subset of READS_DIR, such as a pilot)
#   ./run.sh check      settings, inputs, container and reference
#   ./run.sh all        setup when needed, every sample, combine and validate, each waiting for the last
#   ./run.sh status     samples finished and unfinished, and the state of every recorded job
#
# Single steps:
#   ./run.sh setup      pull the pinned container and download and verify the reference
#   ./run.sh submit     the per-sample array (--array=LIST for chosen rows of samples.tsv)
#   ./run.sh combine    gene-by-sample tables in WORKDIR/combined/<rule>/
#   ./run.sh validate   checks in WORKDIR/validation/<rule>/ (docs/validation.md)
#   ./run.sh compare --rule=RULE --other='PATH/{sample}/...expected_read_counts' --name=NAME
#                       per-sample comparison with another result for the same libraries
# submit, combine, validate and compare pass --dependency=, --begin=, --qos= and --nice= to sbatch.
# Every submitted job is recorded in WORKDIR/jobs.tsv.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_ROOT="$HERE"
source "$HERE/config.sh"
cmd=${1:-}; shift || true
# Recorded in every sample's done.json: the git commit (with -dirty for uncommitted changes), or the REVISION
# file written when the pipeline was copied without its .git folder.
PIPELINE_REVISION=$(git -C "$HERE" describe --always --dirty --abbrev=12 2>/dev/null || cat "$HERE/REVISION" 2>/dev/null || echo unknown)
export PIPELINE_REVISION

usage() { sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; }
rows() { [[ -s "$SAMPLE_SHEET" ]] || { echo "Run ./run.sh sheet first" >&2; exit 1; }; echo $(( $(wc -l < "$SAMPLE_SHEET") - 1 )); }
list_has() { head -1 "$SAMPLE_LIST" | tr '\t' '\n' | tr -d '\r' | grep -qx "$1"; }
reference_ready() { [[ -f "$REF_DIR/reference_manifest.json" && -e "$REF_DIR/ref.fa.fai" ]]; }
container_ready() { [[ -z "$GBRS_SIF" || -s "$GBRS_SIF" ]]; }

# sb STAGE SCRIPT DETAIL [sbatch options...]: submit from WORKDIR, record the job, print its ID.
sb() {
    local stage=$1 script=$2 detail=$3 job; shift 3
    # A file keeps the log folder from being removed as empty before a queued job starts and writes to it.
    mkdir -p "$WORKDIR/slurmout" && touch "$WORKDIR/slurmout/.keep"
    job=$(cd "$WORKDIR" && sbatch --parsable --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" \
        --export=ALL,PIPELINE_ROOT="$HERE" "$@" "$HERE/$script")
    job=${job%%;*}
    printf '%s\t%s\t%s\t%s\t%s\n' "$stage" "$job" "$detail" "$PIPELINE_REVISION" "$(date +%Y-%m-%dT%H:%M:%S%z)" \
        >> "$WORKDIR/jobs.tsv"
    echo "$job"
}

# Options for submit, combine and validate.
ARRAY=""; PASS=(); RULE=""; OTHER=""; NAME=""
parse_options() {
    for a in "$@"; do
        case "$a" in
            --array=*) [[ "$cmd" == submit ]] || { echo "--array is for submit" >&2; exit 1; }; ARRAY="${a#--array=}" ;;
            --rule=*|--other=*|--name=*)
                [[ "$cmd" == compare ]] || { echo "${a%%=*} is for compare" >&2; exit 1; }
                case "$a" in --rule=*) RULE="${a#*=}" ;; --other=*) OTHER="${a#*=}" ;; *) NAME="${a#*=}" ;; esac ;;
            --dependency=*|--begin=*|--qos=*|--nice=*) PASS+=("$a") ;;
            --allow-unlisted) ;;
            *) echo "Unknown option: $a" >&2; exit 1 ;;
        esac
    done
}

do_check() {  # do_check [before-setup]
    local fail=0
    [[ "$WORKDIR" != *CHANGEME* ]] || { echo "Set WORKDIR"; fail=1; }
    [[ "$SLURM_ACCOUNT" != CHANGEME ]] || { echo "Set SLURM_ACCOUNT"; fail=1; }
    if [[ ! -f "$SAMPLE_LIST" ]]; then
        echo "SAMPLE_LIST not found: $SAMPLE_LIST"; fail=1
    elif ! { list_has r1 && list_has r2; }; then
        [[ -d "$READS_DIR" ]] || { echo "READS_DIR not found: $READS_DIR (or give r1 and r2 columns)"; fail=1; }
        case "$WORKDIR/" in "$READS_DIR"/*) echo "WORKDIR must not be inside READS_DIR"; fail=1;; esac
    fi
    for f in "$DNA_ALLELEPROBS" "$STAR_COUNTS"; do
        [[ -z "$f" || -s "$f" ]] || { echo "Not found: $f"; fail=1; }
    done
    for rule in $MATE_RULES; do
        [[ "$rule" == either || "$rule" == both ]] || { echo "MATE_RULES may contain either and both, not $rule"; fail=1; }
    done
    if [[ "${1:-}" != before-setup ]]; then
        reference_ready || { echo "No verified reference in $REF_DIR (./run.sh setup)"; fail=1; }
        container_ready || { echo "Container missing: $GBRS_SIF (./run.sh setup)"; fail=1; }
    fi
    if [[ -z "$GBRS_SIF" ]]; then
        for t in gbrs emase bowtie samtools; do command -v $t >/dev/null || { echo "$t not on PATH"; fail=1; }; done
    fi
    [[ ! -s "$SAMPLE_SHEET" ]] || echo "samples.tsv: $(rows) samples"
    (( fail == 0 )) && echo "Settings look complete" || return 1
}

submit_samples() {  # submit_samples ARRAY [sbatch options...]
    local array=$1; shift
    mkdir -p "$GBRS_OUT"
    sb samples 1_GBRS_Pipeline/gbrs_sample.slurm "$array" --cpus-per-task="$CPUS" --mem="$MEMORY" \
        --time="$TIME_LIMIT" --array="$array" "$@"
}

case "$cmd" in
check)
    do_check ;;
sheet)
    load_software "$MODULE_PYTHON"; mkdir -p "$WORKDIR"
    python3 "$HERE/scripts/sample_sheet.py" "$SAMPLE_LIST" "$READS_DIR" "$SAMPLE_SHEET" "$@" ;;
setup)
    j=$(sb setup 0_Setup/gbrs_setup.slurm "$GBRS_IMAGE")
    echo "Submitted setup job $j; log in $WORKDIR/slurmout/" ;;
submit)
    parse_options "$@"
    n=$(rows); full="1-$n%$THROTTLE"; array=${ARRAY:-$full}
    if [[ -z "$ARRAY" && -s "$WORKDIR/jobs.tsv" ]] && awk -F'\t' '$1 == "samples" || $1 ~ /^[0-9]+$/ {f = 1} END {exit !f}' "$WORKDIR/jobs.tsv"; then
        echo "A sample array is already recorded in $WORKDIR/jobs.tsv; use --array=LIST to rerun chosen rows" >&2; exit 1
    fi
    j=$(submit_samples "$array" "${PASS[@]+"${PASS[@]}"}")
    echo "Submitted $j for rows $array" ;;
combine)
    parse_options "$@"
    j=$(sb combine 1_GBRS_Pipeline/gbrs_combine.slurm "$MATE_RULES" "${PASS[@]+"${PASS[@]}"}")
    echo "Submitted combine job $j" ;;
validate)
    parse_options "$@"
    j=$(sb validate 1_GBRS_Pipeline/gbrs_validate.slurm "$MATE_RULES" "${PASS[@]+"${PASS[@]}"}")
    echo "Submitted validate job $j" ;;
compare)
    parse_options "$@"
    [[ " $MATE_RULES " == *" $RULE "* ]] || { echo "--rule must be one of: $MATE_RULES" >&2; exit 1; }
    [[ "$OTHER" == *"{sample}"* ]] || { echo "--other needs a path template with {sample}" >&2; exit 1; }
    [[ "$NAME" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "--name needs letters, digits, dot, dash or underscore" >&2; exit 1; }
    j=$(COMPARE_RULE=$RULE COMPARE_OTHER=$OTHER COMPARE_NAME=$NAME sb compare 1_GBRS_Pipeline/gbrs_compare.slurm \
        "$RULE vs $NAME: $OTHER" "${PASS[@]+"${PASS[@]}"}")
    echo "Submitted compare job $j; result in $WORKDIR/comparison/${RULE}_vs_${NAME}.tsv" ;;
all)
    mkdir -p "$WORKDIR"
    exec 8>"$WORKDIR/.submit.lock"
    flock -n 8 || { echo "Another submission is running for this WORKDIR" >&2; exit 1; }
    [[ ! -s "$WORKDIR/jobs.tsv" ]] || { echo "$WORKDIR/jobs.tsv exists: inspect those jobs, then use submit, combine or validate" >&2; exit 1; }
    do_check before-setup
    [[ -s "$SAMPLE_SHEET" ]] || { echo "Run ./run.sh sheet first" >&2; exit 1; }
    n=$(rows); after=()
    if ! reference_ready || ! container_ready; then
        j=$(sb setup 0_Setup/gbrs_setup.slurm "$GBRS_IMAGE"); after=("--dependency=afterok:$j")
        echo "setup     $j"
    fi
    j=$(submit_samples "1-$n%$THROTTLE" "${after[@]+"${after[@]}"}"); echo "samples   $j ($n)"
    j=$(sb combine 1_GBRS_Pipeline/gbrs_combine.slurm "$MATE_RULES" --dependency=afterok:$j); echo "combine   $j"
    j=$(sb validate 1_GBRS_Pipeline/gbrs_validate.slurm "$MATE_RULES" --dependency=afterok:$j); echo "validate  $j"
    echo "Jobs recorded in $WORKDIR/jobs.tsv; ./run.sh status shows progress" ;;
status)
    count() { find "$GBRS_OUT" -mindepth "$1" -maxdepth "$1" -name "$2" 2>/dev/null | wc -l | tr -d ' '; }
    n=$(rows); done_n=$(count 2 done.json); partial=$(count 1 '*.partial'); failed=$(count 1 '*.failed_*')
    echo "samples $n, finished $done_n, running or stopped $partial, earlier failed attempts $failed"
    if [[ -s "$WORKDIR/jobs.tsv" ]]; then
        # Current format: stage, job, detail, revision, time. Older files start with the job ID.
        awk -F'\t' '{print ($1 ~ /^[0-9]+$/) ? $1 : $2}' "$WORKDIR/jobs.tsv" | while read -r job; do
            printf '%s\t' "$job"; sacct -j "$job" -X -n -o State 2>/dev/null | sort | uniq -c | tr -s ' ' | paste -sd, -
        done
    fi
    for rule in $MATE_RULES; do
        s="$WORKDIR/validation/$rule/summary.json"
        [[ -s "$s" ]] && echo "validation $rule: $(grep -o '"passed": [a-z]*' "$s")"
    done ;;
*)
    usage; exit 1 ;;
esac
