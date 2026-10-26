# RNAseq_GBRS

Genome Reconstruction by RNA-Seq (GBRS) pipeline for allele-specific expression (ASE) quantification in Collaborative Cross (CC), CC-RIX, and Diversity Outbred (DO) mice. GBRS leverages the known founder haplotype structure of CC mice to reconstruct local genotypes from RNA-seq reads and quantify founder strain-specific expression.

Key contributions of this repository are the novel implementations of GBRS for CC mice (`3_CC_Implementation/`) and CCRIX mice (`4_CCRIX_Implementation/`), including generation of CC-specific transition probability matrices and validation of DO.G0-based genotype reconstruction for CCRIX samples.

## Overview

GBRS uses RNA-seq reads to simultaneously perform:
1. **Genotype reconstruction** — infer founder haplotype probabilities at each genomic locus from read alignment patterns
2. **Allele-specific expression** — quantify expression separately for each of the 8 CC founder alleles

This approach is based on [Choi et al. 2021 (GBRS)](https://github.com/churchill-lab/gbrs) and the [EMASE](https://github.com/churchill-lab/emase) framework.

## Repository Structure

```
RNAseq_GBRS/
├── 1_GBRS_Pipeline/          # Slurm job: alignment → quantification → genotype reconstruction, one task per sample
├── scripts/                  # reference download, sample sheet check, result joining
├── config.sh, run.sh         # settings and the entry point (copy examples/ceres.sh to run_settings.sh)
├── 2_GBRS_Analysis/          # R Markdown: ASE visualization and interpretation
├── 3_CC_Implementation/      # CC-specific transition probability generation and reference data
└── 4_CCRIX_Implementation/   # CCRIX validation workflow and results (uses DO.G0 transprob)
```

## Workflow

```
Preprocessed FASTQs (from RNAseq_STAR)
    │
    ▼
1_GBRS_Pipeline (one Slurm array task per sample)
    ├── Bowtie 1 alignment of R1 and R2 to the pooled-founder transcriptome
    ├── BAM → EMASE format
    ├── Alignments common to both mates (emase get-common-alignments)
    ├── EMASE quantification (multiway)
    ├── Genotype reconstruction (GBRS reconstruct)
    ├── Allele-specific expression quantification with the reconstructed genotypes
    ├── Genotype interpolation, export and genome plot
    └── run.sh combine: gene-by-sample counts, diplotypes and founder counts
    │
    ▼
2_GBRS_Analysis             ← uses transition probability files from 3_CC_Implementation/
    └── ASE visualization and interpretation
```

## Dependencies

| Stage | Tools |
|-------|-------|
| HPC pipeline | GBRS and EMASE 1.0 with Bowtie 1 (pinned container), SAMtools, Apptainer, Slurm |
| R analysis | R, GBRS R package, ape, ggplot2 |

## Running the pipeline

The input is the preprocessed reads from [RNAseq_STAR](https://github.com/BJBennettlab/RNAseq_STAR)
(`01-HTS_Preproc/<sample>/<sample>_R1.fastq.gz` and `_R2.fastq.gz`) and a tab-separated sample list with the
columns `sample`, `sex` (M or F) and `generation`. On Ceres:

```bash
git clone https://github.com/BJBennettlab/RNAseq_GBRS.git
cd RNAseq_GBRS
cp examples/ceres.sh run_settings.sh   # set WORKDIR, READS_DIR, SAMPLE_LIST, SLURM_ACCOUNT
./run.sh containers                    # GBRS 1.0 container, pinned by digest
./run.sh references                    # Zenodo record 15091310, MD5-checked
./run.sh sheet                         # every sample's reads found by exact name
./run.sh check
./run.sh submit
./run.sh status
./run.sh combine                       # after every sample has finished
```

On another Slurm cluster, set `SLURM_PARTITION`, the module names (or `SOFTWARE_MODE=path`) and keep the container.
Samples are always matched by their full name: a sample named `MOUSE_1` never picks up `MOUSE_10`'s reads. Each
sample's folder holds a `done.json` with its inputs, sex and generation, the transition-probability file, tool
versions, alignment counts and output checksums. `run.sh combine` refuses to join an unfinished or changed sample.

| Output in `WORKDIR` | Contents |
| --- | --- |
| `samples/<sample>/` | multiway and diploid quantification, `genotypes.tsv`, `genoprobs.npz`, interpolated and exported genotype probabilities, genome plot, compressed alignments, bowtie logs, `done.json` |
| `combined/gbrs_counts.tsv` | gene by sample expected read counts (diploid quantification, column `total`) |
| `combined/gbrs_diplotypes.tsv` | gene by sample reconstructed diplotype |
| `combined/gbrs_founder_counts.tsv` | gene by sample and founder expected read counts |
| `combined/qc.tsv`, `SHA256SUMS.txt` | read pairs, aligned reads per mate, sex, generation, versions; checksums |

## Reference Data

Reference files required for EMASE and GBRS (pooled-founder reference transcriptome, GBRS_DATA, etc.) are archived on Zenodo:

> Lloyd, Michael. (2025). *Reference Data for EMASE and GBRS*. Zenodo. https://doi.org/10.5281/zenodo.15091310

Dawson Budke is a contributor to this dataset.

## Related Repository

Standard STAR alignment and DE analysis pipeline: [RNAseq_STAR](https://github.com/BJBennettlab/RNAseq_STAR)

## Author

Dawson Budke — dpbudke@ucdavis.edu
