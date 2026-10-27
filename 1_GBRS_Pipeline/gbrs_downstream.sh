#!/bin/bash
# GBRS steps after the two mates are combined and compressed: quantify (multiway), reconstruct, quantify with the
# reconstructed genotypes (diploid), interpolate, export and plot. Called by gbrs_sample.slurm for each mate rule;
# can also finish a sample from an existing compressed alignment file.
#
# Some gbrs commands print an error and still exit 0, so every step's output files are checked before the next step.
# GBRS_DATA points at the reference directory: reconstruct reads chromosome lengths from GBRS_DATA/ref.fa.fai and
# plot reads the founder colours from it.
#
#   PIPELINE_ROOT=/path/to/RNAseq_GBRS gbrs_downstream.sh COMPRESSED_EMASE OUTBASE SAMPLE SEX GENERATION
set -euo pipefail
HERE="${PIPELINE_ROOT:?set PIPELINE_ROOT to the repository directory}"
source "$HERE/config.sh"
[[ $# -eq 5 ]] || { sed -n '2,10p' "$0"; exit 1; }
c=$1; q=$2; sample=$3; sex=$4; generation=$5

tprob="${TPROB_TEMPLATE//\{sex\}/$sex}"; tprob="${tprob//\{generation\}/$generation}"
[[ "$tprob" == /* ]] || tprob="$REF_DIR/$tprob"
need() { for f in "$@"; do [[ -s "$f" ]] || { echo "ERROR: $(basename "$f") was not written" >&2; exit 1; }; done; }
need "$c" "$tprob" "$REF_DIR/ref.fa.fai"
mkdir -p "$(dirname "$q")"
binds=$(printf '%s\n' "$REF_DIR" "$HERE" "$(dirname "$c")" "$(dirname "$q")" "$(dirname "$tprob")" | sort -u | paste -sd, -)
tool() {
    if [[ -n "$GBRS_SIF" ]]; then
        apptainer exec --cleanenv --env GBRS_DATA="$REF_DIR" --bind "$binds" "$GBRS_SIF" "$@"
    else
        GBRS_DATA="$REF_DIR" "$@"
    fi
}

echo "quantify $(date -Is)"
tool gbrs quantify -i "$c" -g "$REF_DIR/emase.gene2transcripts.tsv" -L "$REF_DIR/emase.pooled.fullTranscripts.info" \
    -M "$EMASE_MODEL" --report-alignment-counts -o "$q"
need "$q.multiway.genes.tpm" "$q.multiway.genes.expected_read_counts"
echo "reconstruct $(date -Is)"
tool gbrs reconstruct -e "$q.multiway.genes.tpm" -t "$tprob" -x "$REF_DIR/gbrs_emissions_all_tissues.avecs.npz" \
    -g "$REF_DIR/ref.gene_pos.ordered_ensBuild_105.npz" -o "$q"
need "$q.genotypes.tsv" "$q.genoprobs.npz"
echo "quantify with genotypes $(date -Is)"
tool gbrs quantify -i "$c" -G "$q.genotypes.tsv" -g "$REF_DIR/emase.gene2transcripts.tsv" \
    -L "$REF_DIR/emase.pooled.fullTranscripts.info" -M "$EMASE_MODEL" --report-alignment-counts -o "$q"
need "$q.diploid.genes.expected_read_counts" "$q.diploid.genes.tpm"
echo "interpolate, export, plot $(date -Is)"
tool gbrs interpolate -i "$q.genoprobs.npz" -g "$REF_DIR/ref.genome_grid.GRCm39.tsv" \
    -p "$REF_DIR/ref.gene_pos.ordered_ensBuild_105.npz" -o "$q.interpolated.genoprobs.npz"
need "$q.interpolated.genoprobs.npz"
tool gbrs export -i "$q.interpolated.genoprobs.npz" -s "$HAPLOTYPES" -g "$REF_DIR/ref.genome_grid.GRCm39.tsv" \
    -o "$q.gbrs.interpolated.genoprobs.tsv"
need "$q.gbrs.interpolated.genoprobs.tsv"
tool gbrs plot -i "$q.interpolated.genoprobs.npz" -o "$q.plotted.genome.pdf" -n "$sample"
need "$q.plotted.genome.pdf"
echo "downstream done $(date -Is)"
