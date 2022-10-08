

###==========================================================================#
###---- Get the 10% most distinct species and plot their distributions ------
###==========================================================================#


if (!require("classInt")) { install.packages("classInt", dependencies = T) ; library(classInt)}
if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = T) ; library(rgdal)}
if (!require("dismo")) { install.packages("dismo", dependencies = T) ; library(raster)}
if (!require("XML")) { install.packages("XML", dependencies = T) ; library(XML)}
if (!require("maps")) { install.packages("maps", dependencies = T) ; library(maps)}
if (!require("sp")) { install.packages("sp", dependencies = T) ; library(sp)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}
if (!require("rasterVis")) { install.packages("rasterVis", dependencies = T) ; library(rasterVis)}
if (!require("gridExtra")) { install.packages("gridExtra", dependencies = T) ; library(gridExtra)}


path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()



path_OUT <- paste0(path, "/traits_eco")
rs_all <- read.csv(paste0(path_OUT, "/master_table_spp_fam_comm_ED_RS_volume_residuals.csv"), h=T)


### Get 10% most restricted
thres_sprior_10 <- quantile(rs_all$ED, probs=0.9, na.rm=T)

### Remove the NA n the specific colum
sub_sprior_10 <- subset(rs_all, rs_all$ED >= thres_sprior_10)
sub_sprior_10 <- sub_sprior_10[!is.na(sub_sprior_10$ED),]


### Export the tables
write.csv(sub_sprior_10, paste0(path_OUT, "/ED_sprior_10%.csv"), row.names = F)



add_rangemaps <- read.table("add_these_rangemaps_to_sprior_stack.txt", h=F)
names(add_rangemaps) <- "species"
add_rangemaps$use_rangemap <- 1
sub_sprior_10 <- merge(sub_sprior_10, add_rangemaps, by.x="sequence_original", by.y="species", all.x=T)
sub_sprior_10$use_rangemap <- ifelse(is.na(sub_sprior_10$use_rangemap), 0, sub_sprior_10$use_rangemap)



### Stack the small-ranged species for each dataset
e <- extent(-145, -52, 5, 60)
fun_rescale <- function(x) { scales::rescale(x, to=c(0, 1), from=range(0, 7)) } 


library(raster); library(maptools); library(foreach); library(doParallel)
setwd("/mnt/domisch/domisch/data/fw_fish")
path=getwd()
cl <- makePSOCKcluster(2, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers()

rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=paste0("/data/domisch/R_temp_delete/"))


rs_sprior10 <- foreach(i=sub_sprior_10[["sequence_original"]], .errorhandling="stop", .verbose=T, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  rasterOptions(tmpdir="/data/domisch/data/fw_fish/R_temp_delete")
  rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=paste0(path, "/R_temp_delete/"))
  tmp <- subset(sub_sprior_10, sequence_original ==i)
  if (tmp$filled_data == 1) { 
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif")) 
    # r <- mask(r, lakes_5000skm, inverse=T, updatevalue=0)
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


rs_sprior10 <- stack(unlist(rs_sprior10))
stopCluster(cl)
gc()
  
  rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=paste0(path, "/R_temp_delete/"))
  
  beginCluster(2)
  
  ### Get richness for the 10% subset of species with smallest ranges
  rs_sprior10_sum <- calc(rs_sprior10, sum, na.rm=T)
  writeRaster(rs_sprior10_sum, paste0(path_OUT, "/ED_sprior10_sum.tif"), overwrite=T)
  ### Mask Mexico
  mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))
  
  rs_sprior10_sum <- mask(rs_sprior10_sum, mask_mexico, inverse=T, updatevalue=0)
  rs_sprior10_sum <- trim(rs_sprior10_sum, padding=0, values=NA)
  writeRaster(rs_sprior10_sum, paste0(path_OUT, "/ED_sprior10_sum_mex_masked.tif"), overwrite=T)
  
  endCluster()
  


  
  
  


###=====================================================#
###---- Map the richness of top RS and ED combined -----
###=====================================================#

path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()


if (!require("classInt")) { install.packages("classInt", dependencies = T) ; library(classInt)}
if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = T) ; library(rgdal)}
if (!require("dismo")) { install.packages("dismo", dependencies = T) ; library(raster)}
if (!require("XML")) { install.packages("XML", dependencies = T) ; library(XML)}
if (!require("maps")) { install.packages("maps", dependencies = T) ; library(maps)}
if (!require("sp")) { install.packages("sp", dependencies = T) ; library(sp)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}
if (!require("rasterVis")) { install.packages("rasterVis", dependencies = T) ; library(rasterVis)}
if (!require("gridExtra")) { install.packages("gridExtra", dependencies = T) ; library(gridExtra)}


rm(top_RS_ED_sum, species)
path_OUT <- paste0(path, "/traits_eco")
species <- read.csv(paste0(path_OUT, "/eco_top10_ED_VOL.csv"), h=T) # volume


### Stack the species 
e <- extent(-145, -52, 5, 60)
fun_rescale <- function(x) { scales::rescale(x, to=c(0, 1), from=range(0, 7)) } 

library(raster); library(maptools); library(foreach); library(doParallel)
setwd("/mnt/domisch/domisch/data/fw_fish")
path=getwd()
cl <- makePSOCKcluster( nrow(species), outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers()



top_RS_ED <- foreach(i=species[["sequence_original"]], .errorhandling="stop",.verbose=T, .final=stack, .packages=c("raster", "scales")) %dopar% {
  # k <- gsub("_", " ", i)
  rasterOptions(tmpdir="/data/domisch/data/fw_fish/R_temp_delete")
  tmp <- subset(species, sequence_original ==i)
  if (tmp$filled_data == 1) {
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif"))
    # r <- mask(r, lakes_5000skm, inverse=T, updatevalue=0)
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




### Get richness for the subset of species with smallest ranges
beginCluster(4)
top_RS_ED_sum <- calc(top_RS_ED, sum, na.rm=T)

### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))

top_RS_ED_sum <- mask(top_RS_ED_sum, mask_mexico, inverse=T, updatevalue=NA)
writeRaster(top_RS_ED_sum, paste0(path_OUT, "/top_10_VOL_ED_combined_mex_masked.tif"), overwrite=T)
endCluster()



