# Write the mice's DNA founder probabilities for GBRS's identity check (scripts/validate_gbrs.py --dna,
# setting DNA_ALLELEPROBS).
#
# Inputs, both RDS files from an R/qtl2 analysis of the mice's DNA genotypes (for example GigaMUGA arrays):
#   ALLELEPROB.rds   founder allele probabilities, the output of qtl2::genoprob_to_alleleprob(): one array per
#                    chromosome, mice x 8 founders (A to H) x markers
#   CROSS_OR_PMAP    the cross2 object (its pmap is used) or the physical map itself, positions in Mb
# Output: a gzipped TSV with the columns mouse_id, chr, pos_mb, A to H for every Nth autosomal marker (default
# every 10th; a few thousand markers are enough to tell mice apart). mouse_id is the row name in qtl2; give each
# sample's mouse_id in the GBRS sample list when it differs from the sample name.
#
# Usage: Rscript scripts/export_dna_alleleprobs.R ALLELEPROB.rds CROSS_OR_PMAP.rds OUT.tsv.gz [EVERY]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript export_dna_alleleprobs.R ALLELEPROB.rds CROSS_OR_PMAP.rds OUT.tsv.gz [EVERY]")
out_path <- args[3]
every <- if (length(args) > 3) as.integer(args[4]) else 10L
if (file.exists(out_path)) stop(out_path, " exists")
apr <- readRDS(args[1])
map <- readRDS(args[2])
if (!is.null(map$pmap)) map <- map$pmap

founders <- LETTERS[1:8]
out <- gzfile(out_path, "w")
writeLines(paste(c("mouse_id", "chr", "pos_mb", founders), collapse = "\t"), out)
n_markers <- 0L
for (chr in intersect(names(apr), as.character(1:19))) {
    a <- apr[[chr]]
    if (!identical(dimnames(a)[[2]], founders)) stop("chromosome ", chr, ": founders are not A to H")
    markers <- dimnames(a)[[3]]
    keep <- seq(1L, length(markers), by = every)
    pos <- unname(map[[chr]][markers[keep]])
    if (anyNA(pos)) stop("chromosome ", chr, ": markers missing from the physical map")
    for (k in seq_along(keep)) {
        p <- matrix(a[, , keep[k]], nrow = dim(a)[1])
        values <- apply(p, 1, function(x) paste(sprintf("%.4f", x), collapse = "\t"))
        writeLines(sprintf("%s\t%s\t%.6f\t%s", rownames(a), chr, pos[k], values), out)
    }
    n_markers <- n_markers + length(keep)
}
close(out)
cat(sprintf("wrote %d mice x %d markers to %s\n", dim(apr[[1]])[1], n_markers, out_path))
