### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------

library(hzar)
require(coda)

### Load data and set parameters  -----------------------------------------

data <- read.csv("./EU.HI.csv", header = TRUE)
meta <- read.table("./EU.meta.txt",header=TRUE)
colnames(meta) <- c("locality", "distance","nSample")
meta$distance <- NULL
data <- merge(meta, data, by = "locality", all.y = TRUE)

#task_id=1
#n_chunks=1
Nrun=3
Nchain=3
chainLength <- 100000

# Register parallel processing if available
if(require(doMC)){
  cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  if (is.na(cores)) cores <- 4
  registerDoMC(cores = cores)
} else {
  registerDoSEQ()
}

mainSeed <- list(A = c(596, 528, 124, 978, 544, 99),
                 B = c(528, 124, 978, 544, 99, 596),
                 C = c(124, 978, 544, 99, 596, 528))

# Initialize an empty list to hold the analysis results for all loci
gull <- list()

# Get column names that end with '_a'
locus_columns <- setdiff(colnames(data), c("locality", "distance", "id","nSample"))
locus_subset <- locus_columns

if (is.null(locus_subset)) {
  cat("No loci assigned to this task\n")
  quit(save = "no")
}

# Define all models to fit
model_cfg <- data.frame(
  scaling = c("fixed", "free", "fixed", "free", "fixed", "free"),
  tails   = c("none", "none", "left", "left", "right", "right"),
  modelID = paste0("model", 1:6),
  stringsAsFactors = FALSE
)

Nmods=dim(model_cfg)[1]

# Loop over each locus/trait

#locus_id="EU.auto.Cline.EUorb3"

recheck_modelfit=NULL
Converge_modelfit=NULL

for (locus_id in locus_subset) {
  
  
  Sys.time() 
  gull[[locus_id]] <- list()
  
  # Initialize storage
  gull[[locus_id]]$obs <- list()
  gull[[locus_id]]$models <- list()
  gull[[locus_id]]$fitRs <- list()
  gull[[locus_id]]$runs <- list()
  gull[[locus_id]]$analysis <- list()
  
  # Prepare observed molecular data
  gull[[locus_id]]$obs <- hzar.doMolecularData1DPops(
    data$distance,
    data[[locus_id]],
    data$nSample
  )
  
  
  # --- Load all models dynamically ---
  for (i in seq_len(nrow(model_cfg))) {
    cfg <- model_cfg[i, ]
    gull[[locus_id]]$models[[cfg$modelID]] <- 
      hzar.makeCline1DFreq(gull[[locus_id]]$obs, cfg$scaling, cfg$tails)
  }
  
  # --- Add box constraints to all models ---
  gull[[locus_id]]$models <- lapply(
    gull[[locus_id]]$models,
    hzar.model.addBoxReq,
    -30, 2000
  )
  
  # --- Initialize MCMC fit requests for all models ---
  gull[[locus_id]]$fitRs$init <- lapply(
    gull[[locus_id]]$models,
    hzar.first.fitRequest.old.ML,
    obsData = gull[[locus_id]]$obs,
    verbose = FALSE
  )
  
  # --- Set MCMC parameters dynamically ---
  for (i in seq_len(nrow(model_cfg))) {
    id <- model_cfg$modelID[i]
    gull[[locus_id]]$fitRs$init[[id]]$mcmcParam$chainLength <- chainLength
    gull[[locus_id]]$fitRs$init[[id]]$mcmcParam$burnin <- chainLength %/% 10
    # Assign seeds from mainSeed (optional, can cycle through A, B, C)
    gull[[locus_id]]$fitRs$init[[id]]$mcmcParam$seed[[1]] <- mainSeed[[c("A","B","C")[ (i-1) %% 3 + 1 ]]]
  }
  
  # --- Run initial MCMC fits ---
  gull[[locus_id]]$runs$init <- lapply(
    gull[[locus_id]]$fitRs$init,
    hzar.doFit
  )
  
  # --- Request next run and replicate request 3 times ---
  gull[[locus_id]]$fitRs$chains <- lapply(
    gull[[locus_id]]$runs$init,
    hzar.next.fitRequest
  )
  gull[[locus_id]]$fitRs$chains <- hzar.multiFitRequest(
    gull[[locus_id]]$fitRs$chains,
    each = Nrun,
    baseSeed = NULL
  )
  
  # --- Run this and do 3 chains in parallel ---
  gull[[locus_id]]$runs$chains <- hzar.doChain.multi(
    gull[[locus_id]]$fitRs$chains,
    doPar = TRUE,
    inOrder = FALSE,
    count = Nchain
  )
  
  
  ## -- evaluate convergence -- 
  ## get all dependnet and indepndnet chains;:
  
  runPerMod=Nrun*Nchain
  mcmc_list = lapply(unlist(gull[[locus_id]]$runs$chains, recursive = F), function(x) x$mcmcRaw)
  mcmc_list_perMod = split(mcmc_list, rep(1:Nmods, each = runPerMod))
  mcmc_list_perMod = lapply(mcmc_list_perMod,mcmc.list )
  
  
  #Heidel = lapply(mcmc_list_perMod, heidel.diag)
  #Heidel_out = unlist(lapply(unlist(Heidel, recursive = F), function(x) all(x[,1] == 1 )))
  
  Gelman = lapply(mcmc_list_perMod, gelman.diag)
  Gelman_out = unlist(lapply(Gelman, function(x) x$mpsrf < 1.2))
  
  ESS = lapply(mcmc_list_perMod, effectiveSize)
  ESS_out = all(unlist(ESS) >=200)
  
  
  #  if( !all(c(all(Heidel_out), all(Gelman_out), ESS_out)) ) {
  if( !all(c(all(Gelman_out), ESS_out)) ) {
    Converge_modelfit = c(Converge_modelfit, locus_id)
    #failed_mod=names(Heidel_out[!Heidel_out])
  }
  
  # --- Convert to dataGroups for analysis ---
  gull[[locus_id]]$analysis$initDGs <- list(nullModel = hzar.dataGroup.null(gull[[locus_id]]$obs))
  for (id in names(gull[[locus_id]]$runs$init)) {
    gull[[locus_id]]$analysis$initDGs[[id]] <- hzar.dataGroup.add(gull[[locus_id]]$runs$init[[id]])
  }
  gull[[locus_id]]$analysis$oDG <- hzar.make.obsDataGroup(gull[[locus_id]]$analysis$initDGs)
  gull[[locus_id]]$analysis$oDG <- hzar.copyModelLabels(
    gull[[locus_id]]$analysis$initDGs,
    gull[[locus_id]]$analysis$oDG
  )
  gull[[locus_id]]$analysis$oDG <- hzar.make.obsDataGroup(
    lapply(gull[[locus_id]]$runs$chains, hzar.dataGroup.add),
    gull[[locus_id]]$analysis$oDG
  )
  
  # --- Compute AICc and select best model ---
  gull[[locus_id]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(gull[[locus_id]]$analysis$oDG)
  bestModel <- rownames(gull[[locus_id]]$analysis$AICcTable)[
    which.min(gull[[locus_id]]$analysis$AICcTable$AICc)
  ]
  print(paste("Best model for", locus_id, ":", bestModel))
  gull[[locus_id]]$analysis$model.selected <- gull[[locus_id]]$analysis$oDG$data.groups[[bestModel]]
  # --- Optional: print maximum likelihood parameters ---
  #print(gull[[locus_id]]$analysis$model.selected$ML.cline$param.free)
  
  # test if the best model failed convergence test
  best_model_index <- match(bestModel, model_cfg$modelID)
  best_chain <- mcmc_list_perMod[[best_model_index]]
  gelman_res <- gelman.diag(best_chain)
  ess_res <- effectiveSize(best_chain)
  gelman_fail <- any(gelman_res$mpsrf >= 1.2)
  ess_fail <- any(ess_res < 200)
  
  if (gelman_fail || ess_fail) {
    recheck_modelfit <- c(recheck_modelfit, locus_id)
  }
  
}

fit_modelSelected = lapply(gull, function(x) x$analysis$model.selected)
#hzar.plot.cline(fit_modelSelected$EU.auto.Cline.EUorb3,pch=NA,col='darkorange',add=T,xlim=c(0,2000),lwd=1.5)

save(fit_modelSelected, file = paste0("EU.HZAR_v2_traitHI.RData"))
write.table(recheck_modelfit, file = paste0("EU.HZAR_v2_traitHI_convergence_failed.bestmodel.txt"))
#write.table(Converge_modelfit, file = paste0("EU.HZAR_v2_traitHI_convergence_failed.txt"))

#[1] "Best model for EU.auto.Cline.EUorb1 : model2"
#[1] "Best model for EU.auto.Cline.EUorb1ADF0.5 : model2"
#[1] "Best model for EU.auto.Cline.EUorb2 : model2"
#[1] "Best model for EU.auto.Cline.EUorb2ADF0.5 : model2"
#[1] "Best model for EU.auto.Cline.EUorb3 : model2"
#[1] "Best model for EU.auto.Cline.EUorb3ADF0.5 : model6"
#[1] "Best model for EU.auto.Cline.NAorb1 : model2"
#[1] "Best model for EU.auto.Cline.NAorb1ADF0.5 : model2"
#[1] "Best model for EU.auto.Cline.NAorb2 : model2"
#[1] "Best model for EU.auto.Cline.NAorb2ADF0.5 : model2"
#[1] "Best model for EU.auto.Cline.NAorb3 : model2"
#[1] "Best model for EU.auto.Cline.NAorb3ADF0.5 : model6"
#[1] "Best model for EU.auto.Cline.NAtip1 : model2"
#[1] "Best model for EU.auto.Cline.NAtip1ADF0.5 : model2"
