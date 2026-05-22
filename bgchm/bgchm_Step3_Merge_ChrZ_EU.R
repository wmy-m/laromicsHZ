library(bgchm)

## get a list of all of the rda files
cf<-list.files(pattern="ChrZ.clinesOut")

## load the first rda file and create our storage objects
load(cf[[1]])
Lt<-dim(hyb)[2] ## number of loci
## object for gradient = v and center = c
## one row per locus, 3 columns for point estimate and CIs
v_est<-matrix(NA,nrow=Lt,ncol=3)
c_est<-matrix(NA,nrow=Lt,ncol=3)
name_est<-matrix(NA,nrow=Lt,ncol=1)

## Add genomic position
locus_names <- colnames(genotype)

## using lb and ub stored in rda file
v_est[lb:ub,]<-onev
c_est[lb:ub,]<-onec
name_est[lb:ub] <- locus_names[lb:ub]

## loop over all of the rda files, adding the estimates
for(cfi in 2:length(cf)){
  load(cf[cfi])
  v_est[lb:ub,]<-onev
  c_est[lb:ub,]<-onec
  name_est[lb:ub] <- locus_names[lb:ub]
}
## save object
save(list=ls(),file="EU.combinedClines.ChrZ.rda")

## now impose s2z constraints, make plots, etc.
