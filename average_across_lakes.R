

###=============================#
### ---- Average lake data -----
###=============================#

### Combinen Lake Erie and Lake St. Clair
shapefile("lakes_5000skm_stclair_erie_merged.shp")

### Due to 500km buffer in SDMs, lakes may have stripes
library(raster)
path="/mnt/domisch/data/fw_fish"
rasterOptions(tmpdir="/data/domisch/fw_fish/R_temp_delete")

lakes <- shapefile(paste0(path, "/lakes_5000km/lakes_5000skm_stclair_erie_merged.shp"))
  
mex_masked <- raster(paste0(path, "/stacked_outputs/sprior_sum807_mex_masked.tif"))
e <- extent(mex_masked)
lakes_mex_masked <- crop(lakes, e, snap="out")

shapefile(lakes_mex_masked, paste0(path, "/lakes_5000km/lakes_5000skm_stclair_erie_merged_mex_masked.shp"), overwrite=T)


setwd(paste0(path, "/stacked_outputs"))
rast_list <- list.files(path=getwd(), pattern="*_mex_masked.tif", full.names=F)
rasters <- sapply(rast_list, raster)

library(foreach); library(doParallel)
cl <- makePSOCKcluster(length(rasters), outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


### Run through each map
foreach(i=names(rasters), .errorhandling="stop",.verbose=T, .packages=c("raster")) %dopar% {
  # print((i))
  cat(i, "\n")
  r <- raster(paste0(path, "/stacked_outputs/", i))
  # r <- raster(paste0(path, "/stacked_outputs/tmp/", i))
  lake_val <- extract(r, lakes_mex_masked, fun=mean, na.rm=T, sp=T)
  i_tmp <- gsub(".tif", "", i)
  r_up <- rasterize(lake_val, r, fun="first", field=i_tmp, update=T, na.rm=T, updateValue="all", background=NA)
  
  filename=paste0(i_tmp, "_lakes_avg.tif")
  writeRaster(r_up, paste0(path, "/stacked_outputs/lakes_avg/", filename), overwrite=T)
  
}
stopCluster(cl)


