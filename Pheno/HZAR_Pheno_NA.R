library(hzar)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

# Individual-level trait file
traits <- read.table("./NA.Phenotype.txt", header = TRUE, sep = "\t") 
locations <- read.table("./NA.geodistance.txt",header=TRUE, sep = "\t")

run_hzar <- function(trait_name,
                               traits_df,
                               locations_df,
                               chainLength = 1e5,
                               mainSeed = list(
                                 A = c(978, 544, 99, 596, 528, 124)
                               )) {
  
  library(dplyr)
  library(hzar)
  library(ggplot2)
  
  gull_trait <- list()
  gull_trait[[trait_name]] <- list(
    obs = list(),
    models = list(),
    fitRs = list(),
    runs = list(),
    analysis = list()
  )
  
  # Prepare observation data
  gull_trait[[trait_name]]$obs <- hzar.doNormalData1DRaw(
    hzar.mapSiteDist(locations_df$LocalityID, locations_df$distance),
    traits_df$Location,
    traits_df[[trait_name]]  # trait column
  )
  
  # Model loader helper
  gull_load_model <- function(scaling, tails, id = paste(scaling, tails, sep = ".")) {
    gull_trait[[trait_name]]$models[[id]] <<- 
      hzar.makeCline1DNormal(gull_trait[[trait_name]]$obs, tails)
    
    if (grepl("fixed", scaling, ignore.case = TRUE)) {
      hzar.meta.fix(gull_trait[[trait_name]]$models[[id]])$muL <<- TRUE
      hzar.meta.fix(gull_trait[[trait_name]]$models[[id]])$muR <<- TRUE
      hzar.meta.fix(gull_trait[[trait_name]]$models[[id]])$varL <<- TRUE
      hzar.meta.fix(gull_trait[[trait_name]]$models[[id]])$varR <<- TRUE
    }
    
    frame <- gull_trait[[trait_name]]$obs$frame
    hzar.meta.init(gull_trait[[trait_name]]$models[[id]])$muL <<- frame[1, "mu"]
    hzar.meta.init(gull_trait[[trait_name]]$models[[id]])$varL <<- frame[1, "var"]
    hzar.meta.init(gull_trait[[trait_name]]$models[[id]])$muR <<- frame[nrow(frame), "mu"]
    hzar.meta.init(gull_trait[[trait_name]]$models[[id]])$varR <<- frame[nrow(frame), "var"]
  }
  
  # All model configurations
  modelConfigs <- expand.grid(
    scaling = c("fixed", "free"),
    tails   = c("none", "left", "right", "both"),
    stringsAsFactors = FALSE
  )
  
  # Build models
  for (i in seq_len(nrow(modelConfigs))) {
    cfg <- modelConfigs[i, ]
    modelID <- paste(cfg$scaling, cfg$tails, sep = ".")
    gull_load_model(cfg$scaling, cfg$tails, modelID)
    gull_trait[[trait_name]]$models[[modelID]] <- 
      hzar.model.addBoxReq(gull_trait[[trait_name]]$models[[modelID]], -30, 6310000)
  }
  
  # Compile fit requests
  gull_trait[[trait_name]]$fitRs$init <- lapply(
    gull_trait[[trait_name]]$models,
    function(model) hzar.first.fitRequest.gC(model, obsData = gull_trait[[trait_name]]$obs, verbose = FALSE)
  )
  
  # Initial fits
  gull_trait[[trait_name]]$runs$init <- lapply(
    names(gull_trait[[trait_name]]$fitRs$init),
    function(mid) hzar.doFit(gull_trait[[trait_name]]$fitRs$init[[mid]])
  )
  names(gull_trait[[trait_name]]$runs$init) <- names(gull_trait[[trait_name]]$fitRs$init)
  
  # Prepare MCMC chains
  gull_trait[[trait_name]]$fitRs$chains <- lapply(
    gull_trait[[trait_name]]$runs$init,
    hzar.next.fitRequest
  )
  
  gull_trait[[trait_name]]$fitRs$chains <- hzar.multiFitRequest(
    gull_trait[[trait_name]]$fitRs$chains,
    each = 3,
    baseSeed = mainSeed$A
  )
  
  gull_trait[[trait_name]]$runs$chains <- hzar.doChain.multi(
    gull_trait[[trait_name]]$fitRs$chains,
    doPar = TRUE,
    inOrder = FALSE,
    count = 3
  )
  
  # Convert to hzar.dataGroup and analyze
  gull_trait[[trait_name]]$analysis$oDG <- hzar.make.obsDataGroup(
    lapply(gull_trait[[trait_name]]$runs$chains, hzar.dataGroup.add)
  )
  
  gull_trait[[trait_name]]$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(
    gull_trait[[trait_name]]$analysis$oDG
  )
  
  # Identify best model
  bestIndex <- which.min(gull_trait[[trait_name]]$analysis$AICcTable$AICc)
  gull_trait[[trait_name]]$analysis$model.selected <-
    gull_trait[[trait_name]]$analysis$oDG$data.groups[[bestIndex]]
  
  bestModelName <- paste(modelConfigs$scaling[bestIndex],
                         modelConfigs$tails[bestIndex],
                         sep=".")
    
  return(list(
    trait_name = trait_name,
    gull_trait = gull_trait[[trait_name]],
    bestModelName = bestModelName
  ))
}


traits_to_run <- c("Tarsus_R", "BillShapeRatio")  

trait <- traits_to_run[task_id]
cat("Running trait:", trait, "\n")
result <- run_hzar(trait, traits, locations)

saveRDS(result, file = paste0(trait, "_result.RDS"))

