# To check the current library paths:
.libPaths()

# To add a new directory to the library paths:
.libPaths(c("/scicore/home/marque0000/wu0006/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))

library(bgchm)
setwd("/scicore/home/marque0000/GROUP/wu0006/8_Cline")

#load dataset
load("NA.ChrZ.Cline_input.RData")

# read in prior estimates of hybrid index
h<-read.table("NA.ChrZ.Cline_input_h_est.txt",header=TRUE)

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
sdc<-0.03409006 
sdv<-0.3134507 

# batch size
bsize<-5000

# batch number
myargs<-commandArgs(trailingOnly=TRUE)
k<-as.numeric(myargs[1])

# bounds for this batch
# if the last batch is incomplete, you could set some maximum
# upper bound (ub) here
lb<-(k-1)*bsize + 1
#ub<-lb+bsize-1
ub <- 38343 #for the last set (total snp number)

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
out<-paste("NA.ChrZ.clinesOut",k,".rda",sep="")
save(list=ls(),file=out)
