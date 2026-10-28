# What GBRS does and when to use it

## GBRS and STAR

[RNAseq_STAR](https://github.com/BJBennettlab/RNAseq_STAR) aligns reads to one reference genome (C57BL/6J,
GRCm39) and counts reads per gene. It works for any mouse, and with another reference for any organism.

GBRS (Genome Reconstruction by RNA-Seq, with EMASE) aligns reads to the transcripts of all eight founder strains of
the Collaborative Cross (CC) and Diversity Outbred (DO) populations: A/J (A), C57BL/6J (B), 129S1/SvImJ (C),
NOD/ShiLtJ (D), NZO/HlLtJ (E), CAST/EiJ (F), PWK/PhJ (G) and WSB/EiJ (H). It shares each read among the founders and
transcripts it fits, and from the pattern of expression along the genome it reconstructs which two founders each
mouse carries at every position.

For DO, CC and CC-RIX mice this gives three things a single-reference count does not:

1. Counts that depend less on the B6 reference. Reads from founders whose sequence differs from B6 (CAST, PWK and WSB
   differ most) align to their own founder's transcripts instead of mismatching against B6. With one reference,
   such reads align less often, which lowers the count for mice carrying those founders at a gene and can look like
   a genetic effect on expression.
2. Allele-specific expression: expected read counts for each founder allele of each gene
   (`diploid.genes.expected_read_counts`, columns A to H).
3. A genotype from RNA: founder probabilities along the genome (`gbrs.interpolated.genoprobs.tsv`) and a diplotype
   for each gene. These let each library be checked against its mouse's DNA genotype and against every other
   library ([validation.md](validation.md)).

Limits: GBRS works only for mice descended from the eight founders. It counts annotated transcripts (Ensembl 105 in
the current reference), so it does not find new splice junctions. It aligns with Bowtie 1, which allows up to three
mismatches and no gaps.

For DO, CC and CC-RIX studies, run both from the same preprocessed reads. STAR counts are the lab's common format
across studies and carry the standard alignment QC; GBRS counts serve eQTL and allele-specific analyses in these
populations and the identity checks. Which counts a paper uses is a study decision. For other mice, such as
knockouts on a B6 background, use STAR only.

## Inputs and settings

**Reads.** Paired-end FASTQ, usually RNAseq_STAR's HTStream output (PhiX removed, PCR duplicates removed,
adapters, poly-A tails, Ns and low-quality ends trimmed, reads under 50 bp dropped). JAX's own pipeline starts from
untrimmed reads; either works, but keep one choice within a study. A sample list with `r1` and `r2` columns can
name any read files.

**Sex** selects the sex-specific transition probabilities, which differ for the X chromosome. Take it from the
study's records. The validation step checks it against Xist and Y-chromosome gene expression.

**Generation** selects the transition probabilities, which describe how often the founder at a position changes
between neighbouring positions; that rises with each generation of outbreeding.

| Mice | `TPROB_TEMPLATE` | `generation` in the sample list |
| --- | --- | --- |
| DO | `transition_probabilities/tranprob.DO.G{generation}.{sex}.npz` (default) | the mouse's DO generation (files exist for 0 to 100) |
| CC | `transition_probabilities/tranprob.CC.G20.{sex}.npz` | 20 (not used) |
| CC-RIX | default | 0 (`4_CCRIX_Implementation/README.md`) |

When repeating an earlier study, keep its choices and record them. The lab's DO liver rerun of October 2026 used generation
43 for every mouse, as the 2025 run did, although the recorded generations are 42 to 44.

**Mate rule.** Bowtie aligns the two reads of a pair separately; `MATE_RULES` says how the two are combined.

| Rule | A read pair counts toward a transcript and founder | Used by |
| --- | --- | --- |
| `either` | where either read aligned | the lab's 2025 GBRS results |
| `both` | only where both reads aligned (`emase get-common-alignments`) | the GBRS authors and JAX's pipeline |

One run can produce both from the same alignment. On one DO liver library, `either` counted 18.7 million reads and `both`
16.9 million; 91% of genes had the same diplotype and founder counts correlated at 0.975. Sample identity is the
same under both. Use one rule within a study and report which.

**Reference.** Zenodo record 15091310 (Lloyd 2025): pooled founder transcriptome and Bowtie index, emission and
transition probabilities, gene positions and a genome grid of 74,165 points on GRCm39. `./run.sh setup` downloads it,
checks the published MD5 of each file and records SHA-256 for all of them in `reference_manifest.json`.

## Steps

See `1_GBRS_Pipeline/README.md` for the commands of each step. In order: align each read of the pair, convert to
EMASE alignment files, combine the two reads by the mate rule, compress, quantify with all eight founders allowed
(multiway), reconstruct the genome, quantify again with each gene limited to its two reconstructed founders
(diploid), interpolate the genotype probabilities onto the genome grid, export and plot.

## Versions

| GBRS | Container | Notes |
| --- | --- | --- |
| 1.1.0 | `quay.io/jaxcompsci/gbrs_py3:v1.1.0-338c782` (default, pinned by digest) | March 2026: reworked file handling and memory use, bug fixes; used by JAX's pipeline. On four Bennett lab DO libraries (October 2026) its diploid counts, genotype probabilities and genotypes were byte-identical to 1.0.0's, with the same run time and memory. |
| 1.0.0 | `quay.io/jaxcompsci/gbrs_py3:v1.0.1` | The version behind the lab's 2025 results and the lab's DO liver rerun of October 2026. Rerun from the same reads, it reproduces that library's December 2025 counts exactly. A study already processed with it can keep it. |

JAX's Nextflow pipeline ([cs-nf-pipelines](https://github.com/TheJacksonLaboratory/cs-nf-pipelines), workflow
`gbrs`) runs the same steps with the same reference and the `both` rule. This repository keeps the lab's Slurm setup,
the `either` rule, exact-name sample matching, the CC work and the validation step.
