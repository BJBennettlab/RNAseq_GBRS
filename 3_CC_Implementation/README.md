# 3_CC_Implementation

Novel implementation of GBRS for Collaborative Cross (CC) and CCRIX mice. The standard GBRS tool was originally developed for Diversity Outbred (DO) mice; this directory contains the work required to extend it to CC and reciprocal intercross (CCRIX) populations, including generation of CC-specific transition probability matrices.

## Contents

### Scripts

| File | Description |
|------|-------------|
| `CC.trans.prob.Rmd` | Generates CC founder transition probability matrices used by GBRS for genotype reconstruction. Adapting this to CCRIX requires generating CCRIX-specific transprob files (in progress). |
| `cc_parser.py` | Helper script for parsing CC-format haplotype and concordance files. |

### Transition Probability Matrices

Used by genome reconstruction for CC mice. The same files are in the GBRS Zenodo reference (Lloyd 2025), so a CC run sets
`TPROB_TEMPLATE="transition_probabilities/tranprob.CC.G20.{sex}.npz"` in its settings (the sample list's `generation`
column is then not used; give 20).

| File | Description |
|------|-------------|
| `tranprob.CC.G20.F.h5` | CC transition probability matrix — females, generation 20 (HDF5) |
| `tranprob.CC.G20.M.h5` | CC transition probability matrix — males, generation 20 (HDF5) |
| `tranprob.CC.G20.F.npz` | Same as above (NumPy compressed format) |
| `tranprob.CC.G20.M.npz` | Same as above (NumPy compressed format) |

### Reference Data (`Data/CC_indexes/`)

Haplotype concordance and strain index files used during CC implementation and validation.

| Directory | Contents |
|-----------|----------|
| `CC001/` | `concordance_results.csv`, `strain_index.csv` |
| `CC019/` | `concordance_results.csv`, `strain_index.csv` |
| `CC040/` | `concordance_results.csv`, `strain_index.csv` |
| (root) | `filtered_haplotypes_CC001/019/040.csv` |

## Usage

Open `CC.trans.prob.Rmd` in RStudio and knit to regenerate transition probability matrices:

```r
rmarkdown::render("CC.trans.prob.Rmd")
```

To use newly generated files, point `TPROB_TEMPLATE` at them (absolute path, with `{sex}` for M or F).

## Notes

- CCRIX implementation is complete — see `4_CCRIX_Implementation/` for the full CCRIX validation workflow and results. Custom CCRIX transition probability files were not required; DO.G0 transprob files from the GBRS Zenodo reference (Lloyd 2025) are used instead.
