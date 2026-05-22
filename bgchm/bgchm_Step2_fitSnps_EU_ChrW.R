library(bgchm)

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

#Handling missing data, making ploidy file and setting missing data to NA
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

p_out<-est_p(G0=P0,G1=P1,model="genotype",ploidy="mixed",pldat = ploidy,HMC=FALSE)
h_out<-est_hi(Gx=hyb,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="mixed",pldat = ploidy)

plot(sort(h_out$hi[,1]),ylim=c(0,1),pch=19,xlab="Individual (sorted by HI)",ylab="Hybrid index (HI)")
segments(1:100,h_out$hi[order(h_out$hi[,1]),3],1:100,h_out$hi[order(h_out$hi[,1]),4])

## fit a hierarchical genomic cline model for all 51 loci using the estimated
## hybrid indexes and parental allele frequencies (point estimates)
## use 4000 iterations and 2000 warmup to make sure we get a nice effective sample size
gc_out<-est_genocl(Gx=hyb,p0=p_out$p0[,1],p1=p_out$p1[,1],H=h_out$hi[,1],model="genotype",ploidy="mixed",pldat = ploidy,hier=TRUE,n_iters=4000)

## how variable is introgression among loci? Lets look at the cline SDs
## these are related to the degree of coupling among loci overall
gc_out$SDc
gc_out$SDv

#> gc_out$SDc
#50%        5%       95% 
#0.4711462 0.1099172 1.3617573 
#> gc_out$SDv
#50%        5%       95% 
#0.6322142 0.2578985 1.2318464 

## examine a plot of the joint posterior distribution for the SDs
pp_plot(objs=gc_out,param1="sdv",param2="sdc",probs=c(0.5,0.75,0.95),colors="black",addPoints=TRUE,palpha=0.1,pdf=FALSE,pch=19)

## impose sum-to-zero constraint on log/logit scale
## not totally necessary, but this is mostly a good idea
sz_out<-sum2zero(hmc=gc_out$gencline_hmc,transform=TRUE,ci=0.90)

## plot genomic clines for the 51 loci, first without the sum-to-zero constraint
## then with it... these differ more for some data sets than others
gencline_plot(center=gc_out$center[,1],v=gc_out$gradient,pdf=FALSE)
#gencline_plot(center=sz_out$center[,1],v=sz_out$gradient,pdf=FALSE)
which(gc_out$gradient[,2] > 1) ## index for loci with credibly steep clines
sum(gc_out$gradient[,2] > 1) ## number of loci with credibly steep clines


## summarize loci with credible deviations from genome-average gradients, here the focus is
## specifically on steep clines indicative of loci introgressing less than the average
which(sz_out$gradient[,2] > 1) ## index for loci with credibly steep clines
sum(sz_out$gradient[,2] > 1) ## number of loci with credibly steep clines

## last, lets look at interspecific ancestry for the same data set, this can
## be especially informative about the types of hybrids present
q_out<-est_Q(Gx=hyb,p0=p_out$p0[,1],p1=p_out$p1[,1],model="genotype",ploidy="mixed",pldat = ploidy)

## plot the results
tri_plot(hi=q_out$hi[,1],Q10=q_out$Q10[,1],pdf=FALSE,pch=19)

save.image("EU.ChrW.combinedClines.rda")
