# 1_GBRS_Pipeline

| Job | Started by | Does |
| --- | --- | --- |
| `gbrs_sample.slurm` | `./run.sh submit` (or `all`) | GBRS for one sample per array task |
| `gbrs_downstream.sh` | `gbrs_sample.slurm` | the steps after the reads are combined, with an output check after each |
| `gbrs_combine.slurm` | `./run.sh combine` (or `all`) | gene-by-sample tables for each mate rule (`scripts/combine_gbrs.py`) |
| `gbrs_validate.slurm` | `./run.sh validate` (or `all`) | the checks in `docs/validation.md` (`scripts/validate_gbrs.py`) |

Run them from the repository root after `./run.sh sheet` has written `WORKDIR/samples.tsv` (see the main README).

## Steps in each task

| Step | Command | Output kept in `WORKDIR/samples/<sample>/` (per rule in `<rule>/`) |
| --- | --- | --- |
| Align R1 and R2 separately to the pooled-founder transcriptome | `bowtie -q -a --best --strata --sam -v 3` | `<sample>.bowtie_R1.log`, `_R2.log` |
| Convert each BAM to an EMASE alignment file | `gbrs bam2emase -h A,B,C,D,E,F,G,H` | (temporary) |
| Combine the two mates, once per rule in `MATE_RULES` | `either`: `scripts/merge_mates.py`; `both`: `emase get-common-alignments` | (temporary) |
| Compress | `gbrs compress` | `<rule>/<sample>.compressed.emase.h5` |
| Quantify (multiway) | `gbrs quantify -M 4 --report-alignment-counts` | `<sample>.multiway.*` |
| Reconstruct the genome | `gbrs reconstruct` with `tranprob.DO.G<generation>.<sex>.npz` | `<sample>.genotypes.tsv`, `.genoprobs.npz`, `.genotypes.npz` |
| Quantify with the reconstructed genotypes | `gbrs quantify -G` | `<sample>.diploid.*` (allele-specific expression) |
| Interpolate, export and plot | `gbrs interpolate`, `export`, `plot` | `<sample>.interpolated.genoprobs.npz`, `.gbrs.interpolated.genoprobs.tsv`, `.plotted.genome.pdf` |

The steps from quantification onward are in `gbrs_downstream.sh`, which checks every step's output files (some
gbrs commands report an error and still exit 0) and can finish a sample from an existing compressed alignment file.
BAMs and uncompressed EMASE files are written to node-local scratch and removed. A task writes into
`<sample>.partial/` and renames it to `<sample>/` only after every output exists and `done.json` is written; a
failed attempt is kept as `<sample>.failed_<time>/`. A finished sample is never overwritten.

## Mate rules

A read pair is aligned as two separate reads. `either` lets the pair count toward a transcript and founder where
either mate aligned (the rule behind the lab's 2025 GBRS results); `both` requires both mates to align there and
drops pairs with one unaligned mate (the GBRS authors' paired-end step). Sample identity is the same under both;
counts and allele assignments differ somewhat, so keep one rule within a study and say which one was used.

## Settings that change results

`MATE_RULES`, `EMASE_MODEL` (4), `HAPLOTYPES`, `TPROB_TEMPLATE` with each sample's sex and generation, the reference directory
(recorded with checksums in `reference_manifest.json`) and the container image (pinned by digest). `done.json`
records all of them for each sample.

## Resources

On Ceres, 16 CPUs and 64 GB per sample; 1.5 to 2.6 hours for a 20-million-pair library, mostly alignment
([docs/runtime.md](../docs/runtime.md)). `THROTTLE` limits how many samples run at once.
