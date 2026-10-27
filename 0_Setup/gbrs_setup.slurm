#!/bin/bash
# Setup job: pull the digest-pinned GBRS container to GBRS_SIF and download the GBRS reference into REF_DIR
# (Zenodo record 15091310: published MD5 checked, SHA-256 recorded for every file). Parts already in place are
# checked, not downloaded again. Apptainer does not run on Ceres login nodes, so this is a Slurm job.
#
#   ./run.sh setup        (./run.sh all submits it when either part is missing)
#SBATCH --job-name=gbrs_setup
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=slurmout/gbrs_setup_%j.out
set -euo pipefail
umask 007
HERE="${PIPELINE_ROOT:?submit through ./run.sh setup}"
source "$HERE/config.sh"

if [[ -n "$GBRS_SIF" ]]; then
    load_software "$MODULE_APPTAINER"
    if [[ ! -s "$GBRS_SIF" ]]; then
        echo "pull $GBRS_IMAGE $(date -Is)"
        mkdir -p "$(dirname "$GBRS_SIF")"
        export APPTAINER_TMPDIR="${TMPDIR:-/tmp}"
        partial="${GBRS_SIF%.sif}.partial.sif"
        rm -f "$partial"
        apptainer pull --disable-cache "$partial" "$GBRS_IMAGE"
        mv "$partial" "$GBRS_SIF"
        printf '%s\n' "$GBRS_IMAGE" > "$GBRS_SIF.image"
        (cd "$(dirname "$GBRS_SIF")" && sha256sum "$(basename "$GBRS_SIF")" > "$(basename "$GBRS_SIF").sha256")
    fi
    [[ ! -s "$GBRS_SIF.image" || "$(cat "$GBRS_SIF.image")" == "$GBRS_IMAGE" ]] \
        || { echo "ERROR: $GBRS_SIF was pulled from $(cat "$GBRS_SIF.image"), not GBRS_IMAGE" >&2; exit 1; }
    [[ ! -s "$GBRS_SIF.sha256" ]] || (cd "$(dirname "$GBRS_SIF")" && sha256sum -c "$(basename "$GBRS_SIF").sha256")
    apptainer exec --cleanenv "$GBRS_SIF" bash -c \
        'pip list 2>/dev/null | grep -i -E "^(gbrs|emase|numpy|tables) "; bowtie --version | head -1'
fi
load_software "$MODULE_PYTHON"
python3 "$HERE/scripts/fetch_references.py" "$REF_DIR"
echo "setup done $(date -Is)"
