#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript triangular_run.R <vcf_file> <popmap_file> <p1> <p2> <diff>", call. = FALSE)
}

vcf_file   <- args[1]
popmap_file <- args[2]
p1 <- args[3]
p2 <- args[4]
diff <- as.numeric(args[5])

# Set up libraries
library(triangulaR)
library(vcfR)

# Read input files
vcfR <- read.vcfR(vcf_file)
popmap <- read.table(popmap_file)
colnames(popmap) <- c("id", "pop")

# Run allele frequency difference
vcfR.diff <- alleleFreqDiff(vcfR = vcfR, pm = popmap, p1 = p1, p2 = p2, difference = diff)
save(vcfR.diff, file = paste0("vcfR.diff_", diff, "_", basename(vcf_file), ".RData"))

# Run hybrid index
hi.het <- hybridIndex(vcfR = vcfR.diff, pm = popmap, p1 = p1, p2 = p2)
#save(hi.het, file = paste0("hi_het_", diff, "_", basename(vcf_file), ".RData"))
write.csv(hi.het, file = paste0("hi_het_", diff, "_", basename(vcf_file), ".csv"))
