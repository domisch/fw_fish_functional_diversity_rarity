

###==========================#
###---- Trait analyses ----
###==========================#


usePackage <- function(p){
  if (!is.element(p, installed.packages()[,1])) install.packages(p, dep = TRUE) 
  library(p, character.only = TRUE)
}


usePackage("raster"); usePackage("maptools"); usePackage("foreach"); usePackage("doParallel"); usePackage("data.table"); usePackage("plyr"); usePackage("devtools")
usePackage("FD"); usePackage("cluster"); usePackage("vegan"); 
# install_github("ibartomeus/fundiv")
usePackage("fundiv")

path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()




### Define functions
# source(paste0(path, "/traits/Xtree.R"))

options(stringsAsFactors = FALSE)
"%ni%" <- Negate("%in%")



### Customized FD_dendro function, avoids plotting the tree but saves it as an object). This one uses lapply and foreach and runs faster by magnitudes
source("fw_fish_FD_dendro.R")


INDIR <- paste0(path, "/stacked_outputs/IN")
dir.create(INDIR)
# OUTDIR <- paste0(path, "/traits/OUT")
# dir.create(OUTDIR)


###-----------------------------------------------------#
###--- Prepare data at grid cell level for traits ------
###-----------------------------------------------------#


### Get cell numbers of all stream grid cells

### Load trait data
load(paste0(path, "/traits_eco/fish_traits_and_tree_with_weights.RData"))
# loads "traits_final" and "mytree"
# traits are weighted according to each category

### Load the massive table
tmp <- fread("stacked_outputs/sprior_807_probs_per_grid.txt")
setnames(tmp, c("probs", "range_size", "grid_id_long", "species"))

setkey(tmp,grid_id_long)

### Need to sort the table by grid_id_long
tmp <- arrange(tmp, grid_id_long) # plyr

length(unique(tmp$grid_id_long)) # 4122473 grid cells
min <- min(tmp$grid_id_long) # 4547
max <- max(tmp$grid_id_long) # 60214729

### Divide into chunks, has to be split by grid_id (not by length of file)

### Create vector of splitting numbers
mysplit <- c(seq.int(from=min, to=max, by=15000))[-1]
# mysplit <- c(seq.int(from=min, to=max, by=2000))[-1]
offset <- 1



for(i in seq_along(mysplit)) { # original, from 1...n

cat("File", i, "\n")

### FIRST RUN, uses the "min" information
if (i==1) {
  tmp_sub <- subset(tmp, grid_id_long >= min & grid_id_long <= mysplit[i])
} else {
  tmp_sub <- subset(tmp, grid_id_long >= mysplit[i-offset] & grid_id_long <= mysplit[i])
}

# tmp_sub <- subset(tmp, grid_id_long >= mysplit[i-offset] & grid_id_long <= mysplit[i])

### If has any rows, could be in the ocean, writes "skipped" info out
if (nrow(tmp_sub) != 0) {

### Get data into wide format
tmp_sub$probs <- as.numeric(tmp_sub$probs)

## New version using reshape
tmp_names <- c("grid_id_long", "species", "probs")
tmp_sub <-  tmp_sub[,..tmp_names] 
tmp_sub <- tmp_sub[!duplicated(tmp_sub),]
### get into wide format
tmp_sub_agg <- reshape(tmp_sub, idvar = "grid_id_long", timevar = "species", direction = "wide")
names(tmp_sub_agg) <- gsub("probs.", "", names(tmp_sub_agg))


### Get the grid_id_long
tmp_sub_agg_id <- as.data.frame(with( tmp_sub, tapply(grid_id_long, list(grid_id_long), unique)))
names(tmp_sub_agg_id) <- "V1"

### Columns (=species) need to be in the same order as te rows in "trait_final"
species_order <- row.names(traits_final)

### Not all species modelled in each grid - add empty columns
### Which species do not occur in the subset?
missing <- species_order[species_order %ni% names(tmp_sub_agg) ]

### Create dummy data (missing species)
tmp_dt <- data.table(matrix(ncol = length(missing), nrow = nrow(tmp_sub_agg)))
names(tmp_dt) <- missing

### Add the dummy data
tmp_sub_agg <- cbind(tmp_sub_agg, tmp_dt)

### Set column order (species names) as in traits
setcolorder(tmp_sub_agg, species_order) # data.table

### Replace all NA with zero (probabilities)
tmp_sub_agg[is.na(tmp_sub_agg)] <- 0

### Fuction needs rownames
tmp_sub_agg <- as.data.frame(tmp_sub_agg)
row.names(tmp_sub_agg) <- tmp_sub_agg_id$V1
tmp_sub_agg$grid_id_long <-  NULL 

save(tmp_sub_agg, file=paste0(INDIR, "/subfile_", i, ".RData"))
rm(tmp_sub_agg, tmp_sub_agg_id, tmp_dt); gc()
} else {
  cat("Skipping file", i, "\n")
  write.table(NULL, paste0(INDIR, "/subfile_", i, "_skipped.txt"))
    }
}

stopCluster(cl)





###---------------------------------------#
###--- Load data and run FD analysis ------
###---------------------------------------#

>>>>>> ---> see the FD_batch_per_grid and FR_batch_per_grid <<<<<<

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


# install_github("ibartomeus/fundiv")
require(fundiv)
# zip -r  IN.zip IN
path <- "/mnt/domisch/data/fw_fish"
source(paste0(path, "/traits/Xtree.R"))
source(paste0(path, "/traits/my_FD_dendro.R"))

### Load trait data
load(paste0(path, "/traits/fish_traits_and_tree_with_weights.RData")) 



library(foreach); library(doParallel); library(fundiv)
n_cores=4

cl <- makePSOCKcluster(n_cores, outfile="") # outfile=""
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


INDIR <- paste0(path, "/traits/IN")
OUTDIR <- paste0(path, "/traits/OUT")
TABLEDIR <- paste0(path, "/traits/OUT_tables")

mysplit <- c(seq.int(from=min, to=max, by=10000))[-1]
length(mysplit)

for(i in 523:length(mysplit)) {  # resume at other position in case the computation stopped 
  cat("File", i, "\n")
  
  ### Check if file exists, else skip (Ocean grid cells)
  if(file.exists(paste0(INDIR, "/subfile_", i, ".RData"))) {  # if input file exists 
    
    if(!file.exists(paste0(OUTDIR, "/FDw_per_grid_", i, ".RData"))) { # ..and if output file does not exist
      
      load(paste0(INDIR, "/subfile_", i, ".RData"))
      
      ### Run FD analysis for each grid
      mylist <- setNames(split(tmp_sub_agg, seq(nrow(tmp_sub_agg))), rownames(tmp_sub_agg))
      
            ### Caluclates the full tree, and calculates abundance (probability) -weighted FD for each grid
      trait_weights <- c(rep(1, 8), rep(1/11, 11), rep(1/25, 25))
      
      FD_per_grid_out <-  fw_fish_FD_dendro(S = traits_final[-1], A = tmp_sub_agg, w=trait_weights,
                                       Distance.method= "gower",  Cluster.method = "average", ord = "podani",
                                       Weigthedby = "abundance")
      
      save(FD_per_grid_out, file=paste0(OUTDIR, "/FDw_per_grid_", i, ".RData"))
      
      rm(tmp_sub_agg, FD_per_grid_out); gc()
    } else {
      cat("Output file exists", i, "\n")
    } 
  } else {
    cat("Skipping file", i, "\n")
  }
} # close for loop
stopCluster(cl)
#




###---- Reload and save as txt for bash 'cat' ----

path <- "/data/domisch/data/fw_fish"
if (!require("doParallel")) { install.packages("doParallel", dependencies = T) ; library(doParallel)}
if (!require("foreach")) { install.packages("foreach", dependencies = T) ; library(foreach)}
n_cores=30

mysplit <- c(seq.int(from=min, to=max, by=10000))[-1]


OUTDIR <- paste0(path, "/traits_eco/OUT")
TABLEDIR <- paste0(path, "/traits_eco/OUT_tables")
dir.create(TABLEDIR)


cl <- makePSOCKcluster(n_cores, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


### functional rarity
foreach(i=seq_along(mysplit), .verbose=T, .errorhandling=c('stop')) %dopar% {
  # for(i in seq_along(mysplit)) {
  cat("File", i, "\n")
  if(file.exists(paste0(OUTDIR, "/FRw_per_grid_", i, ".RData"))) { # ..and if output file does not exist
    load(paste0(OUTDIR, "/FRw_per_grid_", i, ".RData"))
    ### Save as .txt file to merge in bash later more efficiently
    i_padd <- sprintf("%04d", i) # add leading zeroes
    
    ### Add grid id
    names(FR_per_grid_out) <- "functional_endemism"
    grid_id <- as.data.frame(row.names(FR_per_grid_out)); names(grid_id) <- "grid_id"
    FR_per_grid_out <- cbind(grid_id, FR_per_grid_out)
      
    write.table(FR_per_grid_out, paste0(TABLEDIR, "/FRw_per_grid_functional_rarity_", i_padd, ".txt"), row.names = F, col.names = F, quote = F)
    
    rm(FR_per_grid_out) # gc()
  }
}


stopCluster(cl)

### Write colnames once
FR_per_grid_out <- names(FR_per_grid_out)
write.table(FR_per_grid_out, paste0(TABLEDIR, "/colnames.txt"), row.names = F, col.names = T, quote = F)

# comm: vector with the name of the community
# n_sp: vector listing the number of species for each community
# n_tr: vector listing the number of traits used
# FDpg: vector listing FDpg (petchey and gaston) for each community
# FDw: vector listing FD weighthed by species relative abundances/biomass in each community
# FDwcomm: vector listing FD weighthed by species abundances/biomass across all communities
# qual.FD: vector repeating the quality of the dendogram representation. clustering performance is assessed by the correlation with the cophenetic distance







###======================================#
###--- Append FD txt tables in bash  ----
###======================================#

DIR=/mnt/domisch/data/fw_fish
cd $DIR


## eco
for FILE in $DIR/traits_eco/OUT/FDw_per_grid_*; do
echo FILE $FILE
cat $FILE >> $DIR/traits_eco/eco_FDw_per_grid_all_cat.txt
done

cat $DIR/traits_eco/eco_FDw_per_grid_all_cat.txt | wc -l



###=====================================================#
###--- Append functional rarity txt tables in bash  ----
###=====================================================#


DIR=/mnt/domisch/data/fw_fish
cd $DIR


### eco
# rangesize
for FILE in $DIR/traits_eco/OUT_FR_grids//FR_per_grid_*; do
echo FILE $FILE
cat $FILE >> $DIR/traits_eco/eco_funrar_per_grid_all_cat_rangesize.txt
done

cat $DIR/traits_eco/eco_funrar_per_grid_all_cat_rangesize.txt | wc -l


# volume
for FILE in $DIR/traits_eco/OUT_FR_volume/FR_per_grid_*; do
echo FILE $FILE
cat $FILE >> $DIR/traits_eco/eco_funrar_per_grid_all_cat_volume.txt
done

cat $DIR/traits_eco/eco_funrar_per_grid_all_cat_volume.txt | wc -l








###============================#
###---- Get FD into raster ----
###============================#


if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("maptools")) { install.packages("maptools", dependencies = T) ; library(maptools)}
if (!require("foreach")) { install.packages("foreach", dependencies = T) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = T) ; library(doParallel)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("plyr")) { install.packages("plyr", dependencies = T) ; library(plyr)}
# if (!require("rgdal")) { install.packages("rgdal", dependencies = T) ; library(rgdal)}



library(raster); library(maptools); library(foreach); library(doParallel); library(data.table); library(plyr)
# path <- "/project/fas/sbsc/sd566/fw_fish/"
path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))
getwd()


source("raster.as.data.table_coord.R")


### Load as data.table
tmp <- fread(paste0(path, "/traits_eco/eco_FDw_per_grid_all_cat.txt"))
setnames(tmp, c("grid_id_long", "n_sp", "n_tr", "FDpg", "FDw", "FDwcomm", "qual.FD"))
setkey(tmp, grid_id_long)
tmp <- tmp[!duplicated(tmp),]

### Load raster template ---
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key

### Assign the predicted values to the raster 
FD_all <- merge(domain_cells, tmp, by="grid_id_long", all.x=T)

### Sort data to match the spatial configuration
FD_all <- arrange(FD_all, seq_id) # plyr
FD_all <- FD_all[!duplicated(FD_all$seq_id),] 

### Create raster and save. Make integer to reduce disk space
number_species_r <- setValues(domain_raster, FD_all$n_sp)
number_traits_r <- setValues(domain_raster, FD_all$n_tr)
FDpg_r <- setValues(domain_raster, FD_all$FDpg)
FDw_r <- setValues(domain_raster, FD_all$FDw)
FDwcomm_r <- setValues(domain_raster, FD_all$FDwcomm)

### Set nodata value
NAvalue(number_species_r) <- 0
NAvalue(number_traits_r) <- 0
NAvalue(FDpg_r) <- 0
NAvalue(FDw_r) <- 0
NAvalue(FDwcomm_r) <- 0


### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))

number_species_r <- mask(number_species_r, mask_mexico, inverse=T, updatevalue=NA)
number_traits_r <- mask(number_traits_r, mask_mexico, inverse=T, updatevalue=NA)
FDpg_r <- mask(FDpg_r, mask_mexico, inverse=T, updatevalue=NA)
FDw_r <- mask(FDw_r, mask_mexico, inverse=T, updatevalue=NA)
FDwcomm_r <- mask(FDwcomm_r, mask_mexico, inverse=T, updatevalue=NA)



number_species_r <- trim(number_species_r, padding=0, values=NA)
number_traits_r <- trim(number_traits_r, padding=0, values=NA)
FDpg_r <- trim(FDpg_r, padding=0, values=NA)
FDw_r <- trim(FDw_r, padding=0, values=NA)
FDwcomm_r <- trim(FDwcomm_r, padding=0, values=NA)


writeRaster(number_species_r, paste0(path, "/stacked_outputs/eco_number_species_mex_masked.tif"), overwrite=T)
writeRaster(number_traits_r, paste0(path, "/stacked_outputs/eco_number_traits_mex_masked.tif"), overwrite=T)
writeRaster(FDpg_r, paste0(path, "/stacked_outputs/eco_FDpg_mex_masked.tif"), overwrite=T)
writeRaster(FDw_r, paste0(path, "/stacked_outputs/eco_FDw_mex_masked.tif"), overwrite=T)
writeRaster(FDwcomm_r, paste0(path, "/stacked_outputs/eco_FDwcomm_mex_masked.tif"), overwrite=T)


# again with lakes averaged
# writeRaster(number_species_r, paste0(path, "/stacked_outputs/lakes_avg/eco_number_species_mex_masked_lakes_avg.tif"), overwrite=T)
# writeRaster(number_traits_r, paste0(path, "/stacked_outputs/lakes_avg/eco_number_traits_mex_masked_lakes_avg.tif"), overwrite=T)
# writeRaster(FDpg_r, paste0(path, "/stacked_outputs/lakes_avg/eco_FDpg_mex_masked_lakes_avg.tif"), overwrite=T)
# writeRaster(FDw_r, paste0(path, "/stacked_outputs/lakes_avg/eco_FDw_mex_masked_lakes_avg.tif"), overwrite=T)
# writeRaster(FDwcomm_r, paste0(path, "/stacked_outputs/lakes_avg/eco_FDwcomm_mex_masked_lakes_avg.tif"), overwrite=T)






###---- Functonal rarity (Rosauer et al.) to raster ----


### Load as data.table

### rangesize
tmp <- fread(paste0(path, "/traits_eco/eco_funrar_per_grid_all_cat_rangesize.txt"))

### volume
tmp <- fread(paste0(path, "/traits_eco/eco_funrar_per_grid_all_cat_volume.txt"))



setnames(tmp, c("funrar_rosauer", "grid_id_long"))
setkey(tmp, grid_id_long)

tmp <- tmp[!duplicated(tmp),]


### Load raster template ---
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key

### Assign the predicted values to the raster 
FR_all <- merge(domain_cells, tmp, by="grid_id_long", all.x=T)


### Sort data to match the spatial configuration
FR_all <- arrange(FR_all, seq_id) # plyr
FR_all <- FR_all[!duplicated(FR_all$seq_id),]

### Create raster and save. Make integer to reduce disk space
FR_all_r <- setValues(domain_raster, FR_all$funrar_rosauer)

### Set nodata value
NAvalue(FR_all_r) <- 0


### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))

FR_all_r <- mask(FR_all_r, mask_mexico, inverse=T, updatevalue=NA)
# FR_all_r <- trim(FR_all_r, padding=0, values=NA)


### rangesize:
writeRaster(FR_all_r, paste0(path, "/stacked_outputs/eco_funrar_rosauer_rangesize_mex_masked.tif"), overwrite=T)



### volume:
writeRaster(FR_all_r, paste0(path, "/stacked_outputs/eco_funrar_rosauer_volume_mex_masked.tif"), overwrite=T)

