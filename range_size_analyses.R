
###-----------------------------------------------#
###--- Evaluate range size / restrictedness ------
###-----------------------------------------------#


library(raster)
library(maptools)
library(foreach)
library(doParallel)
library(parallel)

setwd("/mnt/domisch/data/fw_fish")
path=getwd()

cl <- makePSOCKcluster(15, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers()
rasterOptions(tmpdir="/data/domisch/R_temp_delete")
rasterOptions(chunksize=1000,maxmemory=1000)


### With rangemap
S <- read.table(paste0(path, "/models_done_sprior_769.txt"), h=F)
fun_rescale <- function(x) { scales::rescale(x, to=c(0, 1), from=range(0, 7)) } 

rs_models <- foreach(i=S[[1]], .combine=rbind, .inorder=T ,.errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  cat("#-------- Checking species", i, "\n")
  rasterOptions(tmpdir="/data/domisch/data/fw_fish/R_temp_delete")
  r <- raster(paste0(path, "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
  r <- calc(r, fun_rescale)
  rs <- cellStats(r, sum, na.rm=T)
  rs
}


rs_models_df <- cbind(S, rs_models)


names(rs_models_df) <- c("species", "rs_model")
write.csv(rs_models_df, "rangesize_models_769.csv", row.names = F)

stopCluster(cl)



### Rangemaps
S <- read.table(paste0(path, "/sequence_trait_lookup_827_Jan19.txt"), h=F)

rs_models <- foreach(i=S[[1]], .combine=rbind, .inorder=T ,.errorhandling="pass", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  cat("#-------- Checking species", i, "\n")
  rasterOptions(tmpdir="/home/fas/sbsc/sd566/scratch/fw_fish/R_temp_delete")
  r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif"))
  rs <- cellStats(r, sum, na.rm=T)
  rs
  }


rs_rangemaps_df <- cbind(S, rs_models)
names(rs_rangemaps_df) <- c("species", "rs")
write.csv(rs_rangemaps_df, "rangesize_rangemaps_827_Jan19.csv", row.names = F)




### Rangemaps polygons
S <- read.table(paste0(path, "/sequence_trait_lookup_827_Jan19.txt"), h=F)
rs_polygons <- foreach(i=S[[1]], .combine=rbind, .inorder=T ,.errorhandling="pass", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  cat("#-------- Checking species", i, "\n")
  rasterOptions(tmpdir="/home/fas/sbsc/sd566/scratch/fw_fish/R_temp_delete")
  s <- shapefile(paste0(getwd(), "/species_folders/", i, "/", i, "_rangemap.shp"))
  r_tmp <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif"))
  r <- rasterize(s, r_tmp, field=1)
  rs <- cellStats(r, sum, na.rm=T)
  rs
}


rs_rangemaps_raw_blobs <- cbind(S, rs_polygons)
names(rs_rangemaps_raw_blobs) <- c("species", "rs_raw_blobs")
write.csv(rs_rangemaps_raw_blobs, "rangesize_rangemaps_raw_blobs_827_Jan19.csv", row.names = F)


### Merge all
rs_sprior_df <- read.csv("rangesize_models_772.csv", h=T)

names(rs_sprior_df) <- c("species", "rs_sprior")
rs_rangemaps_df <- read.csv("rangesize_rangemaps_827_Jan19.csv", h=T)
names(rs_rangemaps_df) <- c("species", "rs_rangemap")


rs_polygons_df <- read.csv("rangesize_rangemaps_raw_blobs_846_Sep18.csv", h=T)
 rs_all <- merge(rs_sprior_df, rs_rangemaps_df, by="species", all.y=T)



### For 6 species no rangemap overlap with stream network (Canada)
rs_all$rs_rangemap <- ifelse(rs_all$rs_rangemap ==0, NA, rs_all$rs_rangemap)
write.csv(rs_all, "rangesize_all.csv", row.names = F)

### For those species that lack modeled data, fill in the rangemap data (and keep track of these)
rs_all$filled_data <- ifelse(is.na(rs_all$rs_sprior),1, 0)

### Remove those 6 species that have no rangemap data
rs_all <- subset(rs_all, !is.na(rs_rangemap))
sum(rs_all$filled_data) # 55 species have rangemap data filled in (and Canada-species were removed)

### Copy the rangemap data to the sprior and maxent columns
dim(subset(rs_all, is.na(rs_sprior))) # 29

rs_all$rs_sprior <- ifelse(rs_all$filled_data==1, rs_all$rs_rangemap, rs_all$rs_sprior)
subset(rs_all, filled_data == 1) 

#### Add the raw polygon rangesize information
rs_all_with_poly <- merge(rs_all, rs_polygons_df, by="species")
write.csv(rs_all_with_poly, "rangesize_all_772_827_filled_missing_spp.csv", row.names = F)


### Create frequency distribution of rangesize
png(paste0("rangesize_frequency841_2018.png"), res=200, width=2000)
par(mfrow=c(1,2))
### With rangemap
# png(paste0("rangesize_frequency_sprior867.png"))
hist(rs_all$rs_sprior, breaks=100, xlab="Rangesize in ~1km pixels", main="Rangesize: models with offset")
dev.off()

# ### Maxent
# # png(paste0("rangesize_frequency_maxent867.png"))
# hist(rs_all$rs_maxent, breaks=100, xlab="Rangesize in ~1km pixels", main="Rangesize: models without offset")
# # dev.off()

### Rangemaps
# png(paste0("rangesize_frequency_rangemaps867.png"))
hist(rs_all$rs_sprior, breaks=100, xlab="Rangesize in ~1km pixels", main="Rangesize: rangemaps")
dev.off()






###====================================================#
###---- Get top 10% range-restriced species ----
###====================================================#

rs_all <- read.csv("rangesize_all_772_826_filled_missing_spp.csv", h=T)

# volume
rs_all <- read.csv("/mnt/domisch/data/fw_fish/traits_eco/master_table_spp_fam_comm_ED_RS_volume_residuals_Aug21.csv", h=T) # Aug21: 807 species
S <- read.table(paste0(path, "/sequence_trait_lookup_807_within_domain_for_stacking_Aug21.txt"), h=F)

rs_all <- rs_all[rs_all$sequence_original %in% S$V1, ]
str(rs_all) # 807 species


### Get 10% most restricted
thres_sprior_10 <- quantile(rs_all$rs_sprior, probs=0.1, na.rm=T)

### Remove the NA n the specific colum
sub_sprior_10 <- subset(rs_all, rs_all$rs_sprior <= thres_sprior_10)
sub_sprior_10 <- sub_sprior_10[!is.na(sub_sprior_10$rs_sprior),]

# volume
thres_sprior_10 <- quantile(rs_all$volume.x, probs=0.1, na.rm=T)
sub_sprior_10 <- subset(rs_all, rs_all$volume.x <= thres_sprior_10)
sub_sprior_10 <- sub_sprior_10[!is.na(sub_sprior_10$volume.x),]


### Export the tables
write.csv(sub_sprior_10, "stacked_outputs/rangesize_sprior_10%.csv", row.names = F)
write.csv(sub_sprior_10, "stacked_outputs/volume_sprior_10%.csv", row.names = F)


### Stack the small-ranged species for each dataset
library(raster); library(maptools)
e <- extent(-145, -52, 5, 60)
fun_rescale <- function(x) { scales::rescale(x, to=c(0, 1), from=range(0, 7)) } 


library(raster); library(maptools); library(foreach); library(doParallel)
setwd("/mnt/domisch/domisch/data/fw_fish")
path=getwd()
cl <- makePSOCKcluster(20, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers()


### With rangemap
rs_sprior10 <- foreach(i=sub_sprior_10[["sequence_original"]], .errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  rasterOptions(tmpdir="/data/domisch/R_temp_delete")
  tmp <- subset(sub_sprior_10, sequence_original==i)
  if (tmp$filled_data == 1) {
    # if (tmp$use_rangemap == 1) { 
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif")) 
    r <- extend(r, e, value=NA)
    r
   } else {  
  cat("#-------- Loading modeled results of", paste0(i), "\n")
  r <- raster(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
  r <- calc(r, fun_rescale) # rescale only if modelled map is loaded
  r <- extend(r, e, value=NA)
  r
    }
  }
gc()
stopCluster(cl)
rs_sprior10 <- stack(unlist(rs_sprior10))



# volume
rs_sprior10 <- foreach(i=sub_sprior_10[["sequence_original"]], .errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  rasterOptions(tmpdir="/data/domisch/R_temp_delete")
  tmp <- subset(sub_sprior_10, sequence_original==i)
  if (tmp$filled_data == 1) {
    # if (tmp$use_rangemap == 1) { 
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif")) 
    r <- extend(r, e, value=NA)
    r
  } else {  
    cat("#-------- Loading modeled results of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
    r <- calc(r, fun_rescale) # rescale only if modelled map is loaded
    r <- extend(r, e, value=NA)
    r
  }
}
gc()
stopCluster(cl)
rs_sprior10 <- stack(unlist(rs_sprior10))




beginCluster(20)

### Get richness for the 10% subset of species with smallest ranges
rs_sprior10_sum <- calc(rs_sprior10, sum, na.rm=T)
writeRaster(rs_sprior10_sum, "stacked_outputs/volume_sprior10_sum.tif", overwrite=T)
endCluster()


### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))
rs_sprior10_sum_mex <- mask(rs_sprior10_sum, mask_mexico, inverse=T, updatevalue=0)
rs_sprior10_sum_mex <- trim(rs_sprior10_sum_mex, padding=0, values=NA)
writeRaster(rs_sprior10_sum_mex, paste0(getwd(),"/stacked_outputs/vol_sprior10_sum_mex_masked.tif"), overwrite=T)







###===================================#
###--- Get all rasters as a table  ---
###===================================#



library(raster)
library(maptools)
library(foreach)
library(doParallel)
setwd("/mnt/domisch/data/fw_fish")
path=getwd()

library(parallel); library(doParallel)
cl <- makePSOCKcluster(5, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers()

rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))

### Stream layer
e <- extent(-145, -52, 5, 60)
grid_id_short <- raster(paste0(path, "/global_layers/additional_layers/global_grid_ID.tif"))
grid_id_short <- crop(grid_id_short, e)
names(grid_id_short) <- "grid_id_short"
writeRaster(grid_id_short, paste0(path, "/global_layers/additional_layers/NA_grid_ID_streams_short.tif"), overwrite=T)

### Continuous land layer = long ID to match the grid cells unsing setValues()
grid_id_long <- raster(paste0(path, "/global_layers/additional_layers/grid_id_global_all_cells.tif"))
grid_id_long <- crop(grid_id_long, e)
grid_id_long <- mask(grid_id_long, grid_id_short)
# grid_id_long <- mask(grid_id_long, lakes_5000skm, inverse=T, updatevalue=NA))
names(grid_id_long) <- "grid_id_long"
writeRaster(grid_id_long, paste0(path, "/global_layers/additional_layers/NA_grid_ID_streams_long.tif"), overwrite=T)



### Create grid-id layer using the identical procedure as the species rasters
### Load a random species map as a template
i="Acantharchus_pomotis"
grid_id_long <- raster(paste0(path, "/species_folders/", i, "/rangemap_r_corr_final.tif"))
grid_id_long <- extend(grid_id_long, e, value=NA)
values(grid_id_long) <- 1:ncell(grid_id_long) # fill values
str_order <- raster(paste0(path, "/layers_NA/str_order_lakes0.tif"))
values(str_order) <- ifelse(values(str_order) >=0, 1, NA )
grid_id_long <- mask(grid_id_long, str_order) # mask stream network
names(grid_id_long) <- "grid_id_long"
writeRaster(grid_id_long, paste0(path, "/global_layers/additional_layers/NA_grid_ID_streams_long.tif"), overwrite=T)


library(raster); library(data.table)
source("D:/codes/github/na_fish/fw_fish/raster.as.data.table.R")
setwd("D:/projects/na_fish_manuscript/traits_from_Oct17/upload_GEE")
path=getwd()
rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=getwd())

e <- extent(-145, -52, 5, 60)
grid_id_long <- raster(paste0(path, "/stream_order_lakes0.tif"))
grid_id_long <- crop(grid_id_long, e)
values(grid_id_long) <- 1:ncell(grid_id_long)

str_order <- raster(paste0(path, "/stream_order_lakes0.tif"))
str_order <- crop(str_order, e)
grid_id_long <- mask(grid_id_long, str_order) # maskvalue=F
names(grid_id_long) <- "grid_id_long"
writeRaster(grid_id_long, paste0(path, "/NA_grid_ID_streams_long_nb.asc"), overwrite=T) # OK
tail(values(grid_id_long), 100)


save(grid_id_long, file="grid_id_long_with_lakes.RData")
load("grid_id_long_with_lakes.RData")



### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/mask_mexico2.shp"))

grid_id_long <- mask(grid_id_long, mask_mexico, inverse=T, updatevalue=NA)
grid_id_long <- trim(grid_id_long, padding=0, values=NA)
writeRaster(grid_id_long, paste0(path, "/NA_grid_ID_streams_long_mex_masked.tif"), overwrite=T)
test <- as.data.table(grid_id_long, xy=T)
test <- na.omit(test)
test
save(grid_id_long, file="grid_id_long_with_lakes_mex_masked.RData")
writeRaster(grid_id_long, paste0(path, "/NA_grid_ID_streams_long_nb_mex_masked.asc"), overwrite=T)
na.omit(tail(values(grid_id_long), 1000) )



# grid_id_long <- raster(paste0(path, "/global_layers/additional_layers/NA_grid_ID_streams_long.tif"))
grid_id_long <- raster(paste0(path, "/global_layers/additional_layers/NA_grid_ID_streams_long.tif"))
e <- extent(-145, -52, 5, 60) # North America
fun_rescale <- function(x) { scales::rescale(x, to=c(0, 1), from=range(0, 7)) } 




### With range map (spatial prior)
S <- read.table(paste0(path, "/models_done_sprior.txt"), h=F)


foreach(i=S[[1]], .inorder=F ,.errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  cat("#-------- Running species", i, "\n")
  rasterOptions(tmpdir="/home/fas/sbsc/sd566/scratch/fw_fish/R_temp_delete")
  r <- raster(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
  r <- calc(r, fun_rescale)
  rs <- cellStats(r, sum, na.rm=T)
  r_extended <- extend(r, e, value=NA)
  r_extended_df <- as.data.frame(r_extended)
  names(r_extended_df) <- "probs"
  r_extended_df$range_size <- ifelse(r_extended_df$probs >0, rs, 0)
  save(r_extended_df, file=paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/probs_range_size_per_grid.RData"))

  r_extended_df <- subset(r_extended_df, probs>0)
  r_extended_df$grid_id_long <- rownames(r_extended_df)
  r_extended_df$species <- paste(i)
  write.table(r_extended_df, paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/probs_range_size_per_grid.txt"), row.names=F, col.names=F, quote=F)
  rm(r_extended_df)
  }




### Check missing species
missing <- list()
for (i in S[[1]]) {
# remove underscore
# i <- gsub("_", " ", i) # get original species name without underscore..
#   print(i)
  if(!file.exists(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/probs_range_size_per_grid.RData"))) {
    missing[[i]] <- i
    cat(gsub(" ", "_", i), "missing", "\n") # show progress
  }
}

S <- do.call(rbind, missing); dim(S)[1]
S <- droplevels(as.data.frame(S))
write.table(missing, "data_missing_sprior.txt", row.names=F, col.names=F, quote=F)



### Maxent
S <- read.table(paste0(path, "/models_done_maxent.txt"), h=F)

foreach(i=S[[1]], .inorder=F ,.errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  cat("#-------- Running species", i, "\n")
  rasterOptions(tmpdir="/home/fas/sbsc/sd566/scratch/fw_fish/R_temp_delete")
  r <- raster(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_maxent/bin_sum.tif"))
  values(r) <- ifelse(values(r) >=4, values(r), 0)
  r <- calc(r, fun_rescale)
  rs <- cellStats(r, sum, na.rm=T)
  r_extended <- extend(r, e, value=NA)
  r_extended_df <- as.data.frame(r_extended)
  names(r_extended_df) <- "probs"
  r_extended_df$range_size <- ifelse(r_extended_df$probs >0, rs, 0)
   
  r_extended_df<- subset(r_extended_df, probs>0)
  r_extended_df$grid_id_long <- rownames(r_extended_df)
  ### Add the species name in the table
  r_extended_df$species <- paste(i)
  write.table(r_extended_df, paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_maxent/probs_range_size_per_grid.txt"), row.names=F, col.names=F, quote=F)
rm(r_extended_df)
  }


### Check missing species
 missing <- list()
for (i in S[[1]]) {
# remove underscore
# i <- gsub("_", " ", i) # get original species name without underscore..
#   print(i)
  if(!file.exists(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_maxent/probs_range_size_per_grid.RData"))) {
    missing[[i]] <- i
    cat(gsub(" ", "_", i), "missing", "\n") # show progress
  }
}
S <- do.call(rbind, missing); dim(S)[1]
S <- droplevels(as.data.frame(S))
write.table(missing, "data_missing_maxent.txt", row.names=F, col.names=F, quote=F)


gc()




### Rangemaps
S <- read.table(paste0(path, "/species_sequence_873.txt"), h=F)
# S <- read.table(paste0(path, "/data_missing_rangemaps.txt"), h=F)
# S <- as.data.frame(t(S))

foreach(i=S[[1]], .inorder=F ,.errorhandling="pass", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
   cat("#-------- Running species", i, "\n")
   rasterOptions(tmpdir="/home/fas/sbsc/sd566/scratch/fw_fish/R_temp_delete")
  r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif"))
  rs <- cellStats(r, sum, na.rm=T)
  r_extended <- extend(r, e, value=NA)
  r_extended_df <- as.data.frame(r_extended)
  names(r_extended_df) <- "probs"
  r_extended_df$range_size <- ifelse(r_extended_df$probs >0, rs, 0)
  save(r_extended_df, file=paste0(getwd(), "/species_folders/", i, "/rangemap_count_range_size_per_grid.RData"))

  r_extended_df<- subset(r_extended_df, probs>0)
  r_extended_df$grid_id_long <- rownames(r_extended_df)
  ### Add the species name in the table
  r_extended_df$species <- paste(i)
  write.table(r_extended_df, paste0(getwd(), "/species_folders/", i, "/rangemap_count_range_size_per_grid.txt"), row.names=F, col.names=F, quote=F)
rm(r_extended_df)
  }



### Check missing species
 missing <- list()
for (i in S[[1]]) {
# remove underscore
# i <- gsub("_", " ", i) # get original species name without underscore..
#   print(i)
  # if(!file.exists(paste0(getwd(), "/species_folders/", i, "/rangemap_count_range_size_per_grid.RData"))) {
  	if(!file.exists(paste0(getwd(), "/species_folders/", i, "/rangemap_count_range_size_per_grid.txt"))) {
    missing[[i]] <- i
    cat(gsub(" ", "_", i), "missing", "\n") # show progress
  }
}



missing <- as.data.frame(do.call(rbind, missing)); dim(missing)[1]
S <- as.data.frame(do.call(rbind, missing)); dim(S)[1]
# S <- as.data.frame(rownames(S))
# S <- as.data.frame(droplevels(S))
write.table(S, "data_missing_rangemaps.txt", row.names=F, col.names=F, quote=F)









###========================================================#
###--- Append tables in bash due to RAM limitation in R ---#
###========================================================#

DIR=/mnt/domisch/data/fw_fish
cd $DIR

cat $DIR/sequence_trait_lookup_807_within_domain_for_stacking_Aug21.txt  | wc -l 


### With rangemap
while read SPP_VAR ; do
echo Running species $SPP_VAR
tmp=$DIR/species_folders/${SPP_VAR}/spatial_priors/binary_best_model/probs_range_size_per_grid.txt
cat $tmp >> $DIR/stacked_outputs/sprior_769_probs_per_grid.txt
done < $DIR/sequence_trait_lookup_769_for_getting_model_eval.txt


cat $DIR/stacked_outputs/sprior_769_probs_per_grid.txt | wc -l




### Maxent
while read SPP_VAR ; do
echo Running species $SPP_VAR
tmp=$DIR/species_folders/${SPP_VAR}/spatial_priors/binary_maxent/probs_range_size_per_grid.txt
cat $tmp >> $DIR/stacked_outputs/maxent_873_probs_per_grid.txt
done < $DIR/models_done_maxent.txt


### Run the 42 missing rangemaps for sprior and maxent 

fromdos $DIR/add_these_rangemaps_to_sprior_stack_Aug21.txt

while read SPP_VAR ; do
echo Running species $SPP_VAR
tmp=$DIR/species_folders/${SPP_VAR}/rangemap_count_range_size_per_grid.txt
cat $tmp >> $DIR/stacked_outputs/rangemap_38_sprior_count_range_size_per_grid_Aug21.txt
done < $DIR/add_these_rangemaps_to_sprior_stack_Aug21.txt


### Add the 42 missing species to sprior and maxent data
cat $DIR/stacked_outputs/rangemap_38_sprior_count_range_size_per_grid_Aug21.txt >> $DIR/stacked_outputs/sprior_769_probs_per_grid.txt 

cat $DIR/stacked_outputs/sprior_769_probs_per_grid.txt # --> renamed to $DIR/stacked_outputs/sprior_807_probs_per_grid.txt 








###========================================================#
###---- Range-size aggregation and raster layers in R -----
###========================================================#

usePackage <- function(p){
  if (!is.element(p, installed.packages()[,1])) install.packages(p, dep = TRUE) 
  library(p, character.only = TRUE)
}


usePackage("raster"); usePackage("maptools"); usePackage("foreach"); usePackage("doParallel"); usePackage("data.table"); usePackage("plyr")
path <- "/mnt/domisch/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()


### Raster-to-data.table function
source("raster.as.data.table_coord.R")


###---- Range size ----

### Load as data.table
tmp <- fread("stacked_outputs/sprior_807_probs_per_grid.txt")
setnames(tmp, c("probs", "range_size", "grid_id_long", "species"))
setkey(tmp, grid_id_long)


### Evaluate pixel-based range size (average)
med_RS_perGrid <- tmp[, as.double(median(range_size)),by = grid_id_long] # median, fix bug in data.table
setnames(med_RS_perGrid, c("grid_id_long", "med_RS"))

### Evaluate pixel-based range size (weighted average)
avg_RS_weightedperGrid <- tmp[,lapply(.SD,weighted.mean,w=1/range_size), by=grid_id_long, .SDcols=c("range_size")]
setnames(avg_RS_weightedperGrid, c("grid_id_long", "avg_RS_weighted"))


### Add volume for estimating median volume per species
vol <- read.csv("traits_eco/master_table_spp_fam_comm_ED_RS_volume_residuals.csv")
vol <- subset(vol, select=c(sequence_original,  volume.x))
names(vol) <- c("species", "volume")
tmp <- merge(tmp, vol, by="species")

### Get median and averaged volume-based range-size (use identical names to follow code below..)
med_RS_perGrid <- tmp[, as.double(median(volume)),by = grid_id_long]
setnames(med_RS_perGrid, c("grid_id_long", "med_RS"))


avg_RS_weightedperGrid <- tmp[,lapply(.SD,weighted.mean,w=1/volume), by=grid_id_long, .SDcols=c("volume")]
setnames(avg_RS_weightedperGrid, c("grid_id_long", "avg_RS_weighted"))



### Add residual from rangesize ~ volume and map median residuals
# --> map the discrepancies of range size ~ volume 
vol <- read.csv("traits_eco/master_table_spp_fam_comm_ED_RS_volume_residuals.csv")
vol <- subset(vol, select=c(sequence_original, residuals,	volume_DIV_rs_sprior_m2))
names(vol) <- c("species", "residuals", "volume_DIV_rs_sprior_m2")
vol$volume_DIV_rs_sprior_m2_log <- log(vol$volume_DIV_rs_sprior_m2)
tmp <- merge(tmp, vol, by="species")

### Get median and averaged volume-based range-size (use identical names to follow code below..)
med_RS_perGrid <- tmp[, as.double(median(residuals)),by = grid_id_long]
setnames(med_RS_perGrid, c("grid_id_long", "med_RS"))


### Set keys of the fish-tables
setkey(med_RS_perGrid, grid_id_long)
setkey(avg_RS_weightedperGrid, grid_id_long)



### Load raster / domain template
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key


### Assign the predicted values to the raster 
med_RS_NAwide <- merge(domain_cells, med_RS_perGrid, by="grid_id_long", all.x=T)
avg_RS_weighted_NAwide <- merge(domain_cells, avg_RS_weightedperGrid, by="grid_id_long", all.x=T)
### Sort data to match the spatial configuration
med_RS_NAwide <- arrange(med_RS_NAwide, seq_id) 
avg_RS_weighted_NAwide <- arrange(avg_RS_weighted_NAwide, seq_id) 


### Create raster and save. Make integer to reduce disk space
med_RS_NAwide_r <- setValues(domain_raster, med_RS_NAwide$med_RS)
avg_RS_weighted_NAwide_r  <- setValues(domain_raster, avg_RS_weighted_NAwide$avg_RS_weighted)

### Set nodata value
NAvalue(med_RS_NAwide_r) <- 0
NAvalue(avg_RS_weighted_NAwide_r) <- 0



### With range map
writeRaster(med_RS_NAwide_r, paste0(path, "/stacked_outputs/med_RS_sprior807.tif"), overwrite=T)
writeRaster(avg_RS_weighted_NAwide_r, paste0(path, "/stacked_outputs/avg_RS_weighted_sprior807.tif"), overwrite=T)

# # volume:
writeRaster(med_RS_NAwide_r, paste0(path, "/stacked_outputs/med_RS_volume_sprior807.tif"), overwrite=T)
writeRaster(avg_RS_weighted_NAwide_r, paste0(path, "/stacked_outputs/avg_RS_weighted_volume_sprior807.tif"), overwrite=T)

# # volume residuals:
writeRaster(med_RS_NAwide_r, paste0(path, "/stacked_outputs/med_volume_residuals_807.tif"), overwrite=T)



### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))
### Border to the US kept original, rest simplified to speed up

med_RS_NAwide_r_mex <- mask(med_RS_NAwide_r, mask_mexico, inverse=T, updatevalue=NA)
avg_RS_weighted_NAwide_r_mex <- mask(avg_RS_weighted_NAwide_r, mask_mexico, inverse=T, updatevalue=NA)

writeRaster(med_RS_NAwide_r_mex, paste0(path, "/stacked_outputs/med_RS_sprior807_mex_masked.tif"), overwrite=T)
writeRaster(avg_RS_weighted_NAwide_r_mex, paste0(path, "/stacked_outputs/avg_RS_weighted_sprior807_mex_masked.tif"), overwrite=T)

# # volume:
writeRaster(med_RS_NAwide_r_mex, paste0(path, "/stacked_outputs/med_RS_volume_sprior807_mex_masked.tif"), overwrite=T)
writeRaster(avg_RS_weighted_NAwide_r_mex, paste0(path, "/stacked_outputs/avg_RS_weighted_volume_sprior807_mex_masked.tif"), overwrite=T)
# 

# # volume residuals:
writeRaster(med_RS_NAwide_r_mex, paste0(path, "/stacked_outputs/med_volume_residuals_mex_masked_807.tif"), overwrite=T)


# ### Maxent
# writeRaster(med_RS_NAwide_r, paste0(path, "/stacked_outputs/med_RS_maxent873.tif"), overwrite=T)
# writeRaster(avg_RS_weighted_NAwide_r, paste0(path, "/stacked_outputs/avg_RS_weighted_maxent873.tif"), overwrite=T)


# ### Rangemaps
# writeRaster(med_RS_NAwide_r, paste0(path, "/stacked_outputs/med_RS_rangemaps873.tif"), overwrite=T)
# writeRaster(avg_RS_weighted_NAwide_r, paste0(path, "/stacked_outputs/avg_RS_weighted_rangemaps873.tif"), overwrite=T)
# 





###=================================================================#
###----- Range-size rarity  and  Average range size rarity -----
###=================================================================#


library(raster)
library(maptools)
library(foreach)
library(doParallel)
library(data.table)
library(plyr)
path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()



### Raster-to-data.table function
source("raster.as.data.table_coord.R")


### Load raster template
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key


### Load range size
tmp_rs <- fread("stacked_outputs/sprior_807_probs_per_grid.txt")
setnames(tmp_rs, c("probs", "range_size", "grid_id_long", "species"))
setkey(tmp_rs, grid_id_long)




### Add volume for estimating median volume  per species --
vol <- read.csv("traits_eco/master_table_spp_fam_comm_ED_RS_volume_residuals.csv", h=T)
vol <- subset(vol, select=c(sequence_original, volume.x))
names(vol) <- c("species", "volume")
tmp_rs <- merge(tmp_rs, vol, by="species")

### Get inverse range size
tmp_rs$inverse_range_size <- 1/tmp_rs$volume

## Aggregate data: get rangesize rarity
rs_rar <- tmp_rs[, as.double(sum(inverse_range_size)),by = grid_id_long] 
setnames(rs_rar, c("grid_id_long", "range_size_rarity"))
setkey(rs_rar, grid_id_long)


### Get richness
rich <- tmp_rs[, as.double(length(unique(species))),by = grid_id_long] 
setnames(rich, c("grid_id_long", "richness"))
setkey(rich, grid_id_long)

### Merge the number of species to range-size rarity (aggregated per grid) table
rs_rar <- merge(rs_rar, rich, by="grid_id_long", all.x=T)  

### Divide rangesize rarity by richness -> average rangesize rarity
rs_rar$avg_range_size_rarity <- rs_rar$range_size_rarity / rs_rar$richness

### Merge to domain cells (raster)
tmp_merge <- merge(domain_cells, rs_rar, by="grid_id_long", all.x=T)

### Sort by seq_id to match the raster....
tmp_merge <- arrange(tmp_merge, seq_id)



### Get inverse range size
tmp_rs$inverse_range_size <- 1/tmp_rs$range_size


## Aggregate data: get rangesize rarity
rs_rar <- tmp_rs[, as.double(sum(inverse_range_size)),by = grid_id_long] 
setnames(rs_rar, c("grid_id_long", "range_size_rarity"))
setkey(rs_rar, grid_id_long)


### Get richness
rich <- tmp_rs[, as.double(length(unique(species))),by = grid_id_long] 
setnames(rich, c("grid_id_long", "richness"))
setkey(rich, grid_id_long)

### Merge the number of species to range-size rarity (aggregated per grid) table
rs_rar <- merge(rs_rar, rich, by="grid_id_long", all.x=T)  

### Divide rangesize rarity by richness -> average rangesize rarity
rs_rar$avg_range_size_rarity <- rs_rar$range_size_rarity / rs_rar$richness

### Merge to domain cells (raster)
tmp_merge <- merge(domain_cells, rs_rar, by="grid_id_long", all.x=T)

### Sort by seq_id to match the raster
tmp_merge <- arrange(tmp_merge, seq_id)


### Create raster and save
range_size_rarity_r <- setValues(domain_raster, tmp_merge$range_size_rarity)
average_range_size_rarity_r <- setValues(domain_raster, tmp_merge$avg_range_size_rarity)
richness_r <- setValues(domain_raster, tmp_merge$richness)


### Set nodata value
NAvalue(range_size_rarity_r) <- 0
NAvalue(average_range_size_rarity_r) <- 0
NAvalue(richness_r) <- 0

### Write raster
writeRaster(range_size_rarity_r, paste0(path, "/stacked_outputs/range_size_rarity.tif"), overwrite=T)
writeRaster(average_range_size_rarity_r, paste0(path, "/stacked_outputs/average_range_size_rarity.tif"), overwrite=T)
writeRaster(richness_r, paste0(path, "/stacked_outputs/richness.tif"), overwrite=T)

### Write raster (volume)
writeRaster(range_size_rarity_r, paste0(path, "/stacked_outputs/range_size_rarity_volume.tif"), overwrite=T)
writeRaster(average_range_size_rarity_r, paste0(path, "/stacked_outputs/average_range_size_rarity_volume.tif"), overwrite=T)
writeRaster(richness_r, paste0(path, "/stacked_outputs/richness_volume.tif"), overwrite=T)





### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))

range_size_rarity_r <- mask(range_size_rarity_r, mask_mexico, inverse=T, updatevalue=NA)
average_range_size_rarity_r <- mask(average_range_size_rarity_r, mask_mexico, inverse=T, updatevalue=NA)
richness_r <- mask(richness_r, mask_mexico, inverse=T, updatevalue=NA)

range_size_rarity_r <- trim(range_size_rarity_r, padding=0, values=NA)
average_range_size_rarity_r <- trim(average_range_size_rarity_r, padding=0, values=NA)
richness_r <- trim(richness_r, padding=0, values=NA)


writeRaster(range_size_rarity_r, paste0(getwd(), "/stacked_outputs/range_size_rarity_mex_masked.tif"), overwrite=T)
writeRaster(average_range_size_rarity_r, paste0(getwd(), "/stacked_outputs/average_range_size_rarity_mex_masked.tif"), overwrite=T)
writeRaster(richness_r, paste0(getwd(), "/stacked_outputs/richness_mex_masked.tif"), overwrite=T)


### Write raster (volume)
writeRaster(range_size_rarity_r, paste0(path, "/stacked_outputs/range_size_rarity_volume_mex_masked.tif"), overwrite=T)
writeRaster(average_range_size_rarity_r, paste0(path, "/stacked_outputs/average_range_size_rarity_volume_mex_masked.tif"), overwrite=T)
writeRaster(richness_r, paste0(path, "/stacked_outputs/richness_volume_mex_masked.tif"), overwrite=T)




