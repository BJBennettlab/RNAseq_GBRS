#!/bin/bash
# GBRS for one sample: one Slurm array task per row of WORKDIR/samples.tsv.
#
#   bowtie 1 (R1 and R2 separately, all best-stratum alignments, up to 3 mismatches) -> bam2emase,
#   then for each rule in MATE_RULES (either: a pair counts where either mate aligned, the 2025 rule;
#   both: only where both mates aligned, the GBRS authors' step) -> combine mates -> compress
#   -> quantify (multiway) -> reconstruct (sex- and generation-specific transition
#   probabilities) -> quantify with the reconstructed genotypes (diploid)
#   -> interpolate -> export -> plot
#
# BAMs and uncompressed EMASE files stay in node-local scratch and are removed.
# The compressed alignment file and every GBRS output are kept in
# WORKDIR/samples/<sample>/<rule>/. done.json records inputs, tool versions, alignment
# counts and output checksums; a sample with done.json is never overwritten.
#
#   ./run.sh submit        (or ./run.sh submit --array=1,5 for chosen rows)
#SBATCH --job-name=gbrs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=slurmout/gbrs_%A_%a.out

set -euo pipefail
umask 007

find_pipeline_root() {
    local c
    for c in "${PIPELINE_ROOT:-}" "${SLURM_SUBMIT_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." 2>/dev/null && pwd)"; do
        if [[ -n "$c" && -f "$c/config.sh" ]]; then echo "$c"; return 0; fi
    done
    echo "ERROR: cannot locate config.sh; run through ./run.sh or export PIPELINE_ROOT" >&2
    return 1
}
HERE="$(find_pipeline_root)"
source "$HERE/config.sh"
# The container brings bowtie, GBRS and EMASE; SAMtools and Python come from the site or PATH.
if [[ -n "$GBRS_SIF" ]]; then
    load_software "$MODULE_APPTAINER" "$MODULE_SAMTOOLS" "$MODULE_PYTHON"
else
    load_software "$MODULE_BOWTIE" "$MODULE_SAMTOOLS" "$MODULE_PYTHON"
fi

task=${SLURM_ARRAY_TASK_ID:?run as a Slurm array task}
row=$(awk -F'\t' -v n="$((task + 1))" 'NR == n' "$SAMPLE_SHEET")
[[ -n "$row" ]] || { echo "ERROR: task $task has no row in $SAMPLE_SHEET" >&2; exit 1; }
IFS=$'\t' read -r sample r1 r2 r1_bytes r2_bytes sex generation <<< "$row"
echo "task $task sample $sample sex $sex generation $generation host $HOSTNAME start $(date -Is)"

out="$GBRS_OUT/$sample"
[[ -e "$out/done.json" ]] && { echo "ERROR: $out/done.json exists; this sample is finished" >&2; exit 1; }
[[ "$(stat -c %s "$r1")" == "$r1_bytes" && "$(stat -c %s "$r2")" == "$r2_bytes" ]] \
    || { echo "ERROR: $sample reads changed since the sample sheet was written" >&2; exit 1; }

tprob="${TPROB_TEMPLATE//\{sex\}/$sex}"; tprob="${tprob//\{generation\}/$generation}"
[[ "$tprob" == /* ]] || tprob="$REF_DIR/$tprob"
[[ -s "$tprob" ]] || { echo "ERROR: transition probabilities not found: $tprob" >&2; exit 1; }

scratch="${TMPDIR:-/tmp}/gbrs_${sample}_${SLURM_JOB_ID:-$$}"
# Keep an earlier failed attempt for inspection instead of mixing its files in.
[[ -e "$out.partial" ]] && mv "$out.partial" "$out.failed_$(date +%Y%m%dT%H%M%S)"
mkdir -p "$scratch" "$out.partial"
trap 'rm -rf "$scratch"' EXIT
cpus=${SLURM_CPUS_PER_TASK:-$CPUS}
binds=$(printf '%s\n' "$WORKDIR" "$REF_DIR" "$HERE" "$(dirname "$r1")" "$(dirname "$r2")" "$scratch" "$(dirname "$tprob")" | sort -u | paste -sd, -)
tool() {
    if [[ -n "$GBRS_SIF" ]]; then
        apptainer exec --cleanenv --bind "$binds" "$GBRS_SIF" "$@"
    else
        "$@"
    fi
}
p="$out.partial/$sample"

for mate in R1 R2; do
    reads=$r1; [[ $mate == R2 ]] && reads=$r2
    echo "align $mate $(date -Is)"
    tool bowtie -p "$cpus" -q -a --best --strata --sam -v 3 -x "$BOWTIE_INDEX" "$reads" 2> "$p.bowtie_$mate.log" \
        | samtools view -b -o "$scratch/$mate.bam" -
    echo "bam2emase $mate $(date -Is)"
    [[ -s "$scratch/$mate.bam" ]] || { echo "ERROR: alignment of $mate wrote nothing" >&2; exit 1; }
    tool gbrs bam2emase -i "$scratch/$mate.bam" -m "$REF_DIR/emase.fullTranscripts.info" \
        -h "$HAPLOTYPES" -o "$scratch/$mate.emase.h5"
    [[ -s "$scratch/$mate.emase.h5" ]] || { echo "ERROR: bam2emase of $mate wrote nothing" >&2; exit 1; }
    rm -f "$scratch/$mate.bam"
done

for rule in $MATE_RULES; do
    d="$out.partial/$rule"; q="$d/$sample"; mkdir -p "$d"
    echo "combine mates ($rule) $(date -Is)"
    case "$rule" in
        both)   tool emase get-common-alignments -i "$scratch/R1.emase.h5" -i "$scratch/R2.emase.h5" \
                    -o "$scratch/merged.emase.h5" ;;
        either) tool python3 "$HERE/scripts/merge_mates.py" either "$scratch/R1.emase.h5" "$scratch/R2.emase.h5" \
                    "$scratch/merged.emase.h5" ;;
        *)      echo "ERROR: MATE_RULES may contain either and both, not $rule" >&2; exit 1 ;;
    esac
    [[ -s "$scratch/merged.emase.h5" ]] || { echo "ERROR: combining mates ($rule) wrote nothing" >&2; exit 1; }
    tool gbrs compress -i "$scratch/merged.emase.h5" -o "$q.compressed.emase.h5"
    [[ -s "$q.compressed.emase.h5" ]] || { echo "ERROR: compress ($rule) wrote nothing" >&2; exit 1; }
    rm -f "$scratch/merged.emase.h5"
    echo "rule $rule: quantify, reconstruct, interpolate, export, plot"
    "$HERE/1_GBRS_Pipeline/gbrs_downstream.sh" "$q.compressed.emase.h5" "$q" "$sample" "$sex" "$generation"
done
rm -f "$scratch"/*.emase.h5

versions="$(tool bash -c 'pip list 2>/dev/null | grep -i -E "^(gbrs|emase) " | tr -s " " | paste -sd";"; bowtie --version | head -1' | paste -sd"|")|$(samtools --version | head -1)"
python3 - "$out.partial" "$sample" "$r1" "$r2" "$sex" "$generation" "$tprob" "${GBRS_IMAGE:-}" "$versions" "$MATE_RULES" <<'PY'
import hashlib, json, re, sys
from pathlib import Path
d, sample, r1, r2, sex, gen, tprob, image, versions, rules = sys.argv[1:]
d = Path(d)
def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()
stats = {}
for mate in ('R1', 'R2'):
    log = (d / f'{sample}.bowtie_{mate}.log').read_text()
    m = re.search(r'reads processed: (\d+)', log); a = re.search(r'at least one reported alignment: (\d+)', log)
    if not m or not a: raise SystemExit(f'bowtie log for {mate} lacks its summary')
    stats[mate] = {'reads': int(m.group(1)), 'aligned': int(a.group(1))}
if stats['R1']['reads'] != stats['R2']['reads']: raise SystemExit('R1 and R2 read counts differ')
needed = ['multiway.genes.tpm', 'genotypes.tsv', 'genoprobs.npz', 'diploid.genes.expected_read_counts',
          'interpolated.genoprobs.npz', 'gbrs.interpolated.genoprobs.tsv', 'compressed.emase.h5']
for rule in rules.split():
    for n in needed:
        p = d / rule / f'{sample}.{n}'
        if not p.is_file() or p.stat().st_size == 0: raise SystemExit(f'missing output {rule}/{p.name}')
record = {'sample': sample, 'reads': {'r1': r1, 'r2': r2}, 'sex': sex, 'generation': int(gen),
          'transition_probabilities': tprob, 'container': image, 'versions': versions, 'bowtie': stats,
          'mate_rules': rules.split(),
          'outputs': {str(p.relative_to(d)): sha(p) for p in sorted(d.rglob('*')) if p.is_file()}}
(d / 'done.json').write_text(json.dumps(record, indent=1) + '\n')
PY
mv "$out.partial" "$out"
echo "done $sample $(date -Is)"
