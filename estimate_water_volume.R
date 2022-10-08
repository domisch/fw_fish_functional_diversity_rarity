


###========================================================#
###---- Get average stream flow for each stream order -----
###========================================================#


### Download FLO1K data
wget https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10598674/FLO1K.ts.1960.2015.qma.nc
wget https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10598158/FLO1K.ts.1960.2015.qav.nc
wget https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10609456/FLO1K.ts.1960.2015.qmi.nc

### Average and create a tif
cdo timmean FLO1K.ts.1960.2015.qav.nc   mean_flow_1km_1960_2015.nc
cdo timmean FLO1K.ts.1960.2015.qma.nc   max_flow_1km_1960_2015.nc
cdo timmean FLO1K.ts.1960.2015.qmi.nc   min_flow_1km_1960_2015.nc

gdal_translate -of GTiff  mean_flow_1km_1960_2015.nc   mean_flow_1km_1960_2015.tif  -co COMPRESS=LZW
gdal_translate -of GTiff  max_flow_1km_1960_2015.nc   max_flow_1km_1960_2015.tif  -co COMPRESS=LZW
gdal_translate -of GTiff  min_flow_1km_1960_2015.nc   min_flow_1km_1960_2015.tif  -co COMPRESS=LZW


### cdo timmean: Processed 52254720000 values from 1 variable over 56 timesteps [1106.20s 17GB]
# units: m3/s
# missing_value=-1


library(raster)
library(data.table)
library(foreach)
library(doParallel)

path <- "/data/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))
getwd()




### Extract the multi-year annual flow by each stream order
flo1k_mean <- raster(paste0(path, "/FLO1k/mean_flow_1km_1960_2015.tif"))
flo1k_min <- raster(paste0(path, "/FLO1k/min_flow_1km_1960_2015.tif"))
flo1k_max <- raster(paste0(path, "/FLO1k/max_flow_1km_1960_2015.tif"))


### Load layers
strord <- raster(paste0(path, "/layers_NA/str_order_lakes0.tif"))

flo1k_mean <- crop(flo1k_mean, extent(strord))
flo1k_min <- crop(flo1k_min, extent(strord))
flo1k_max <- crop(flo1k_max, extent(strord))

### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))
### Border to the US kep original, rest simplified to speed up


strord <- mask(strord, mask_mexico, inverse=T, updatevalue=NA)

flo1k_mean <- mask(flo1k_mean, mask_mexico, inverse=T, updatevalue=NA)
flo1k_min <- mask(flo1k_min, mask_mexico, inverse=T, updatevalue=NA)
flo1k_max <- mask(flo1k_max, mask_mexico, inverse=T, updatevalue=NA)


cl <- makePSOCKcluster(9, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

### Get mean values per stream order
flow_all <- foreach(i=seq(1:9), .combine=rbind, .packages="raster")  %dopar% { 
  new_i <- i-1

  mask <- strord
  values(mask) <- ifelse(values(mask)==new_i, 1, NA)
  
  flo1k_mean_tmp <- mask(flo1k_mean, mask)
  avg_mean <- cellStats(flo1k_mean_tmp, "mean")
  
  flo1k_min_tmp <- mask(flo1k_min, mask)
  avg_min <- cellStats(flo1k_min_tmp, "mean")
  
  flo1k_max_tmp <- mask(flo1k_max, mask)
  avg_max <- cellStats(flo1k_max_tmp, "mean")
  
  
  flo1k_mean_tmp <- mask(flo1k_mean, mask)
  sd_mean <- cellStats(flo1k_mean_tmp, "sd")
  
  flo1k_min_tmp <- mask(flo1k_min, mask)
  sd_min <- cellStats(flo1k_min_tmp, "sd")
  
  flo1k_max_tmp <- mask(flo1k_max, mask)
  sd_max <- cellStats(flo1k_max_tmp, "sd")
  
  
  str_ord <- new_i
  res <- cbind(str_ord, avg_mean, sd_mean, avg_min, sd_min, avg_max, sd_max)
  return(res)
  }

stopCluster(cl)


write.csv(flow_all, paste0(path, "/flow1k_mean_min_max_average_sd_per_stream_ord.csv"), row.names = F, quote = F)




###===================================#
###----- Water volume estimates ------
###===================================#

### Calculate for each species an estimate of water volume occupied, i.e. Sum of pixel x typical flow rate of stream order that pixel has
### Load species table
myspp <- read.csv(paste0(path, "/rangesize_all_spp.csv"))



library(foreach); library(doParallel)
cl <- makePSOCKcluster(25, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

rm(volume, volume_tmp, r, tmp_flow)
### sprior with offset
estimated_volume <- foreach(i=myspp[["species"]], .errorhandling="stop",.verbose=T, .combine=rbind, .packages=c("raster")) %dopar% {
  # k <- gsub("_", " ", i)
  rasterOptions(tmpdir="/data/domisch/data/fw_fish/R_temp_delete")
  tmp <- subset(myspp, species ==i)
  if (tmp$filled_data == 1) { 
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/rangemap_r_corr_final.tif")) 
    values(r) <- ifelse(values(r)>0, 1, NA)
    r <- extend(r, strord, value=NA) # else "different extent error.."
    tmp_flow <- crop(flo1k_mean, r)
    tmp_flow <- extend(tmp_flow, strord, value=NA)
    volume_tmp <- mask(tmp_flow, r)
    volume <- cellStats(volume_tmp, "sum")
    volume <- data.frame(species=i, volume=volume)
    volume
   
  } 
  else {  
    cat("#-------- Loading modeled results of", paste0(i), "\n")
    r <- raster(paste0(getwd(), "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
    values(r) <- ifelse(values(r)>0, 1, NA)
    r <- extend(r, strord, value=NA)
    tmp_flow <- crop(flo1k_mean, r)
    tmp_flow <- extend(tmp_flow, strord, value=NA)
    volume_tmp <- mask(flo1k_mean, r)
    volume <- cellStats(volume_tmp, "sum")
    volume <- data.frame(species=i, volume=volume)
    volume
  }
  
}

stopCluster(cl)
head(estimated_volume)
write.csv(estimated_volume, paste0(path, "/estimated_volume_per_species.csv"), quote = F, row.names = F)
