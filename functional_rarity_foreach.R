

###---------------------------------------#
###--- Load data and run FR analysis ------
###---------------------------------------#


### Get the variable
# 
# if(Sys.info()[["sysname"]]=="Linux") {
#   args <- commandArgs()
#   print(args)
#   i <- args[6]
#   cat("Variable", i, "\n")
#   # spp <- gsub("_", " ", args[6]) # get original species name without underscore
# }


### Load libraries
if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("maptools")) { install.packages("maptools", dependencies = T) ; library(maptools)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("plyr")) { install.packages("plyr", dependencies = T) ; library(plyr)}
if (!require("devtools")) { install.packages("devtools", dependencies = T) ; library(devtools)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = T) ; library(doParallel)}
if (!require("foreach")) { install.packages("foreach", dependencies = T) ; library(foreach)}
if (!require("FD")) { install.packages("FD", dependencies = T) ; library(FD)}
if (!require("cluster")) { install.packages("cluster", dependencies = T) ; library(cluster)}
if (!require("vegan")) { install.packages("vegan", dependencies = T) ; library(vegan)}
if (!require("ape")) { install.packages("ape", dependencies = T) ; library(ape)}


path <- "/mnt/domisch/domisch/data/na_fish"
source(paste0(path, "/scripts/Xtree.R"))
source(paste0(path, "/scripts/fw_fish_FD_dendro.R"))

### Source functions for PE (Rosauer et al)
source(paste0(path, "/scripts/phylomatchr.R"))
source(paste0(path, "/scripts/FastXtreePhylo.R"))


### Load trait data
load("/mnt/domisch/domisch/data/na_fish/traits_eco/fish_traits_and_tree_with_weights.RData")



# ======================================================#
INDIR <- "/mnt/domisch/domisch/data/na_fish/stacked_outputs/IN"
# ======================================================#

### Choose if grids or volume used in the calculation

### GRIDS
# source(paste0(path, "/scripts/phyloendemism_sami.R"))
# OUTDIR <- "/mnt/domisch/domisch/data/na_fish_all/traits_eco/OUT_FR_grids"

### VOLUME
source(paste0(path, "/scripts/phyloendemism_volume.R"))
OUTDIR <- "/mnt/domisch/domisch/data/na_fish/traits_eco/OUT_FR_volume"
dir.create(OUTDIR)


cl <- makePSOCKcluster(10, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


foreach(i=seq(1:4014), .errorhandling = "stop", .packages = c("ape", "vegan", "FD", "cluster"), .verbose=F) %dopar% {

  cat("File", i, "\n")
  
  ### Check if file exists, else skip
  if(file.exists(paste0(INDIR, "/subfile_", i, ".RData"))) {  # if input file exists 
    
    if(!file.exists(paste0(OUTDIR, "/FR_per_grid_", i, ".RData"))) { # ..and if output file does not exist
      
      load(paste0(INDIR, "/subfile_", i, ".RData"))
      
      ### Run FD analysis for each grid

      ### Calculates the full tree, and calculates abundance (probability) -weighted FD for each grid
      trait_weights <- c(rep(1, 7)) # eco 
      ### Remove species not present in the tree anymore (=out of domain, extinct)
      toss_these <- names(tmp_sub_agg)[!(names(tmp_sub_agg) %in% traits_final$sequence_original)]
      tmp_sub_agg <- tmp_sub_agg[,!(names(tmp_sub_agg) %in% toss_these)]
      my_rownames <- row.names(tmp_sub_agg)
      tmp_sub_agg <- as.data.frame(sapply(tmp_sub_agg, as.numeric), row.names = my_rownames)
    
      ### Phylogentic endemism (Rosauer et al.)
      FR_per_grid_out <- phyloendemism(tmp_sub_agg, mytree, weighted = TRUE)
      
      # ### Save as .txt file to merge in bash later more efficiently
      i_padd <- sprintf("%04d", as.numeric(i) ) # add leading zeroes
      FR_per_grid_out <- cbind(FR_per_grid_out, data.frame(grid_id_long=row.names(FR_per_grid_out))) # already written in "comm" column
      write.table(FR_per_grid_out, paste0(OUTDIR, "/FR_per_grid_", i_padd, ".txt"), row.names = F, col.names = F, quote = F)
      
      rm(tmp_sub_agg, FR_per_grid_out); gc()
    } else {
      cat("Output file exists", i, "\n")
    } 
  } else {
    cat("Skipping file", i, "\n")
  }
  
} # close foreach loop
stopCluster(cl)
