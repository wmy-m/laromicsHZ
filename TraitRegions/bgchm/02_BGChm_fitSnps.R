library(bgchm)
setwd("./8_Cline/TraitRegions/BGChm/manyloci")

#The same code is used for EU

#load dataset
load("NA_bgchm1.RData")

# read in prior estimates of hybrid index
h<-read.table("NA_h_est.txt",header=TRUE)

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

# enter point estimates of cline standard deviations inferred from the subset of loci
sdc<-1.179983 
sdv<-0.4350847 

#EU point estimates
#sdc<-1.255208 
#sdv<-0.7166359
# batch size
bsize<-1000

# batch number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

# bounds for this batch
# if the last batch is incomplete, you could set some maximum
# upper bound (ub) here
lb<-(k-1)*bsize + 1
ub<-lb+bsize-1
#ub <-7953 #for the last set (total snp number)

# number of loci and hybrids for this analysis
L<-length(lb:ub)
N<-nrow(h)

# subset genotype objects
sGhyb<-hyb[,lb:ub]
sGP0<-P0[,lb:ub]
sGP1<-P1[,lb:ub]
sPloidy <- list(
  ploidy_hyb = ploidy$ploidy_hyb[, lb:ub],
  ploidy_P0  = ploidy$ploidy_P0[, lb:ub],
  ploidy_P1  = ploidy$ploidy_P1[, lb:ub]
)

# in this example, I am just saving the center and gradient,
# point estimates and CIs, for each locus
# one snp at a time, L rows for loci, 3 columns for point est and CIs
onev<-matrix(NA,nrow=L,ncol=3)
onec<-matrix(NA,nrow=L,ncol=3)

# loop over loci
for (i in 1:L) {

  pldat_i <- list(
  ploidy_hyb = sPloidy$ploidy_hyb[, i],
  ploidy_P0  = sPloidy$ploidy_P0[, i],
  ploidy_P1  = sPloidy$ploidy_P1[, i]
)

  out <- est_genocl(
    Gx = sGhyb[, i],
    G0 = sGP0[, i],
    G1 = sGP1[, i],
    H  = h[, 1],
    model = "genotype",
    ploidy = "mixed",
    pldat = pldat_i,
    hier = FALSE,
    SDc = sdc,
    SDv = sdv
  )

  onev[i, ] <- out$gradient
  onec[i, ] <- out$center
}

# save the R object with batch ID included
out<-paste("NA.auto.clinesOut",k,".rda",sep="")
save(list=ls(),file=out)
