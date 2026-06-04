library(bgchm)
setwd("./8_Cline/TraitRegions/BGChm/manyloci")


P0 <- read.table("NA.auto.Cline.TraitRegions_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("NA.auto.Cline.TraitRegions_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("NA.auto.Cline.TraitRegions_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

# Handle missing data
ploidy_P1<-matrix(2, nrow = nrow(P1), ncol = ncol(P1))
ploidy_P1[P1 == -1 | P1 == -2] <- 0 

ploidy_P0<-matrix(2, nrow = nrow(P0), ncol = ncol(P0))
ploidy_P0[P0 == -1 | P0 == -2] <- 0 

ploidy_hyb<-matrix(2, nrow = nrow(hyb), ncol = ncol(hyb))
ploidy_hyb[hyb == -1 | hyb == -2] <- 0 

ploidy <- list(ploidy_hyb = ploidy_hyb, ploidy_P0 = ploidy_P0, ploidy_P1 = ploidy_P1)

P0[P0 == -1 | P0 == -2] <- NA
P1[P1 == -1 | P1 == -2] <- NA
hyb[hyb == -1 | hyb == -2] <- NA

L<-dim(hyb)[2]
rset<-sample(1:L,1000,replace=FALSE)

hyb_s <- hyb[, rset]
P0_s  <- P0[, rset]
P1_s  <- P1[, rset]

pldat_s <- list(
  ploidy_hyb = ploidy$ploidy_hyb[, rset],
  ploidy_P0  = ploidy$ploidy_P0[, rset],
  ploidy_P1  = ploidy$ploidy_P1[, rset]
)

## ---- Estimate hybrid index ----
h_out <- est_hi(
  Gx = hyb[, rset],
  G0 = P0[, rset],
  G1 = P1[, rset],
  ploidy = "mixed",
  model = "genotype",
  pldat = pldat_s
)

write.table(
  h_out$hi,
  file = paste0("NA_h_est.txt"),
  row.names = FALSE,
  quote = FALSE
)

## ---- Estimate genotype classes ----
gc_out <- est_genocl(
  Gx = hyb[, rset],
  G0 = P0[, rset],
  G1 = P1[, rset],
  H  = h_out$hi[, 1],
  ploidy = "mixed",
  model = "genotype",
  hier = TRUE,
  n_iters = 4000,
  pldat = pldat_s
)
save(gc_out, file = "NA_bgchm1.RData")

cat("Mean SDc:", mean(gc_out$SDc, na.rm = TRUE), "\n")
cat("Mean SDv:", mean(gc_out$SDv, na.rm = TRUE), "\n")
