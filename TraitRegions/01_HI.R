### Clear environment-------------------------------------------------------

rm(list = ls())

### Working directory and dependencies--------------------------------------
library(vcfR)
library(triangulaR)

### Extracting hybrid index for each trait ---------------------------------
# Read input files
popmap <- read.table("./EUpopmap.cline.txt")
#popmap <- read.table("./NApopmap.cline.txt")

colnames(popmap) <- c("id", "pop")

vcf_files <- list.files("./", pattern = "^EU.*\\.vcf\\.gz$", full.names = TRUE)
#vcf_files <- list.files("./", pattern = "^NA.*\\.vcf\\.gz$", full.names = TRUE)

# Get HI for each trait
for (vcf_path in vcf_files) {
  
  # Read VCF
  vcfR <- read.vcfR(vcf_path)
  
  # Run hybrid index
  hi.het <- hybridIndex(vcfR = vcfR, pm = popmap, p1 = "GLAUC", p2 = "OCCI")
  #p1 = "GLAUC", p2 = "OCCI"
  #p1 = "ARGE", p2 = "CACH"
  # Create output file name based on VCF name
  trait_name <- tools::file_path_sans_ext(basename(vcf_path))
  trait_name <- sub("\\.vcf$", "", trait_name)  # remove .vcf if needed
  
  output_file <- paste0("data.", trait_name, ".csv")
  
  # Write results
  write.csv(hi.het, file = output_file, row.names = FALSE)
  
}

### Formatting input file --------------------------------------------------
locality <- read.table("Samplelist.txt")
colnames(locality) <- c("id", "locality")

dist <- read.table("EU.geodistance.txt", header=TRUE)  #In Pheno folder
colnames(dist) <- c("locality","Raw", "distance")
dist$Raw=NULL

#dist <- read.table("NA.geodistance.txt", header=TRUE)  #In Pheno folder
#colnames(dist) <- c("locality", "distance")


# Folder with all your HI CSVs
csv_files <- list.files(pattern = "^data\\.EU.*\\.csv$", full.names = TRUE)

# Initialize empty list to store data frames
hi_list <- list()

for (file in csv_files) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df <- df[, c("id", "hybrid.index")]  
  
  # Name the HI column based on the trait
  trait_name <- sub("^data\\.(NA.*)\\.csv$", "\\1", basename(file))
  colnames(df)[2] <- trait_name
  
  hi_list[[trait_name]] <- df
}

# Merge all data frames
hi_merged <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), hi_list)
hi_merged2 <- merge(locality, hi_merged, by = "id", all.y = TRUE)
hi_merged3 <- merge(dist, hi_merged2, by = "locality", all.y = TRUE)
#hi_merged3$raw <- NULL

write.csv(hi_merged3, file = "NA.HI.AFD0.5.csv", row.names = FALSE)

