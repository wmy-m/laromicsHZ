library(bgchm)

P0 <- read.table("EU.auto.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("EU.auto.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("EU.auto.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("EU.auto.Cline_input.RData")

### Clear environment-------------------------------------------------------

rm(list = ls())


P0 <- read.table("EU.ChrW.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("EU.ChrW.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("EU.ChrW.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("EU.ChrW.Cline_input.RData")

### Clear environment-------------------------------------------------------

rm(list = ls())


P0 <- read.table("EU.ChrZ.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("EU.ChrZ.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("EU.ChrZ.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("EU.ChrZ.Cline_input.RData")

### Clear environment-------------------------------------------------------

rm(list = ls())


P0 <- read.table("NA.auto.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("NA.auto.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("NA.auto.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("NA.auto.Cline_input.RData")

### Clear environment-------------------------------------------------------

rm(list = ls())


P0 <- read.table("NA.ChrW.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("NA.ChrW.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("NA.ChrW.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("NA.ChrW.Cline_input.RData")
### Clear environment-------------------------------------------------------

rm(list = ls())


P0 <- read.table("NA.ChrZ.Cline_P0.txt", header = TRUE, sep = "")
row.names(P0) <- P0[, 1]  # Set the first column as row names
P0 <- P0[, -1]            # Remove the first column (sample names)

P1 <- read.delim("NA.ChrZ.Cline_P1.txt", header = TRUE, sep = "")
row.names(P1) <- P1[, 1]  # Set the first column as row names
P1 <- P1[, -1]

hyb <- read.delim("NA.ChrZ.Cline_hybrids.txt", header = TRUE, sep = "")
row.names(hyb) <- hyb[, 1]  # Set the first column as row names
hyb <- hyb[, -1]

hyb <- as.matrix(hyb)
P0 <- as.matrix(P0)
P1 <- as.matrix(P1)

genotype <- rbind(P0, P1, hyb)

save.image("NA.ChrZ.Cline_input.RData")
