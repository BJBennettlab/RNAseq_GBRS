# RNAseq_GBRS

Genome Reconstruction by RNA-Seq (GBRS) pipeline for allele-specific expression (ASE) quantification in Collaborative Cross (CC), CC-RIX, and Diversity Outbred (DO) mice. GBRS leverages the known founder haplotype structure of CC mice to reconstruct local genotypes from RNA-seq reads and quantify founder strain-specific expression.

Key contributions of this repository are the novel implementations of GBRS for CC mice (`3_CC_Implementation/`) and CCRIX mice (`4_CCRIX_Implementation/`), including generation of CC-specific transition probability matrices and validation of DO.G0-based genotype reconstruction for CCRIX samples.

## Overview

GBRS uses RNA-seq reads to simultaneously perform:
1. **Genotype reconstruction** — infer founder haplotype probabilities at each genomic locus from read alignment patterns
2. **Allele-specific expression** — quantify expression separately for each of the 8 CC founder alleles

This approach is based on [Choi et al. 2021 (GBRS)](https://github.com/churchill-lab/gbrs) and the [EMASE](https://github.com/churchill-lab/emase) framework.
[docs/method.md](docs/method.md) explains what GBRS adds over a single-reference count such as
[RNAseq_STAR](https://github.com/BJBennettlab/RNAseq_STAR), which mice and settings it needs, and the two ways of
combining the reads of a pair.

## Repository Structure

```
RNAseq_GBRS/
├── run.sh, config.sh         # entry point and settings (copy examples/ceres.sh to run_settings.sh)
├── 0_Setup/                  # Slurm job: pinned container and verified reference
├── 1_GBRS_Pipeline/          # Slurm jobs: one task per sample, then combine and validate
├── scripts/                  # sample sheet, reference download, mate combining, joining, validation, DNA export
├── docs/                     # method, validation, run time
├── 2_GBRS_Analysis/          # R Markdown: ASE visualization and interpretation
├── 3_CC_Implementation/      # CC-specific transition probability generation and reference data
└── 4_CCRIX_Implementation/   # CCRIX validation workflow and results (uses DO.G0 transprob)
```

## Workflow

```
Paired-end reads (usually RNAseq_STAR's HTStream output) + sample list (sample, sex, generation)
    │
    ▼
./run.sh all     (every step is a Slurm job; each waits for the one before)
    ├── setup      pinned GBRS container, reference from Zenodo with checksums
    ├── samples    one array task per sample:
    │                Bowtie 1 alignment of R1 and R2 to the pooled-founder transcriptome → EMASE format
    │                → combine the two reads by each mate rule ("either" and/or "both")
    │                → quantify → reconstruct the genome → quantify with the reconstructed genotypes
    │                → interpolate, export and plot the genotype probabilities
    ├── combine    gene-by-sample counts, diplotypes and founder counts for each mate rule
    └── validate   finished samples, sex, duplicate genotypes, DNA and STAR identity
    │
    ▼
2_GBRS_Analysis   ← uses transition probability files from 3_CC_Implementation/
    └── ASE visualization and interpretation
```

## Running the pipeline

On Ceres:

```bash
git clone https://github.com/BJBennettlab/RNAseq_GBRS.git
cd RNAseq_GBRS
cp examples/ceres.sh run_settings.sh   # set WORKDIR, READS_DIR, SAMPLE_LIST, SLURM_ACCOUNT
./run.sh sheet                         # finds every sample's reads by exact name, writes WORKDIR/samples.tsv
./run.sh check
./run.sh all                           # setup if needed, every sample, combine, validate
./run.sh status
```

The sample list is tab-separated with the columns `sample`, `sex` (M or F) and `generation`. By default each
sample's reads are `READS_DIR/<sample>/<sample>_R1.fastq.gz` and `_R2.fastq.gz`, the layout of RNAseq_STAR's
HTStream stage; `r1` and `r2` columns name the read files directly instead. Optional columns for the checks:
`mouse_id` (the mouse's ID in the DNA genotype file), `genotype_group` (samples expected to share a genotype) and
`star_library`. For CC and CC-RIX mice, see the table in [docs/method.md](docs/method.md).

Samples are always matched by their full name: a sample named `MOUSE_1` never picks up `MOUSE_10`'s reads, and a
read file listed for two samples stops the run. Each sample's folder holds a `done.json` with its inputs, sex and
generation, transition-probability file, tool versions, pipeline revision, alignment counts and output checksums.
A finished sample is never overwritten; rerun failed rows with `./run.sh submit --array=LIST`.

On another Slurm cluster, start from `examples/other_slurm.sh`. Run time, memory and storage are in
[docs/runtime.md](docs/runtime.md).

| Output in `WORKDIR` | Contents |
| --- | --- |
| `samples/<sample>/<rule>/` | multiway and diploid quantification, `genotypes.tsv`, `genoprobs.npz`, interpolated and exported genotype probabilities, genome plot, compressed alignments |
| `samples/<sample>/` | bowtie logs and `done.json` |
| `combined/<rule>/gbrs_counts.tsv` | gene by sample expected read counts (diploid quantification, column `total`) |
| `combined/<rule>/gbrs_diplotypes.tsv` | gene by sample reconstructed diplotype |
| `combined/<rule>/gbrs_founder_counts.tsv` | gene by sample and founder expected read counts |
| `combined/<rule>/qc.tsv`, `SHA256SUMS.txt` | read pairs, aligned reads per mate, sex, generation, versions; checksums |
| `validation/<rule>/` | `summary.json` and one table per check ([docs/validation.md](docs/validation.md)) |
| `jobs.tsv`, `slurmout/` | every submitted job with its pipeline revision; job logs |

## Dependencies

| Stage | Tools |
|-------|-------|
| HPC pipeline | Slurm, Apptainer, SAMtools, Python 3; GBRS, EMASE and Bowtie 1 from the pinned container |
| DNA genotypes for the identity check (optional) | R with the mice's R/qtl2 allele probabilities (`scripts/export_dna_alleleprobs.R`) |
| R analysis | R, GBRS R package, ape, ggplot2 |

## Reference Data

Reference files required for EMASE and GBRS (pooled-founder reference transcriptome, GBRS_DATA, etc.) are archived on Zenodo:

> Lloyd, Michael. (2025). *Reference Data for EMASE and GBRS*. Zenodo. https://doi.org/10.5281/zenodo.15091310

Dawson Budke is a contributor to this dataset.

## Related Repository

Standard STAR alignment and DE analysis pipeline: [RNAseq_STAR](https://github.com/BJBennettlab/RNAseq_STAR)

## Author

Dawson Budke — dpbudke@ucdavis.edu
