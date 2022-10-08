

###========================================#
###---- Water volume per stream order -----
###========================================#

### Calculate the summed (1km pixel) area of all pixels in a stream order and multiply with the mean flow
flo1k_mean <- raster(paste0(path, "/FLO1k/mean_flow_1km_1960_2015.tif"))


### Load layers
strord <- raster(paste0(path, "/layers_NA/str_order_lakes0.tif"))

flo1k_mean <- crop(flo1k_mean, extent(strord))

### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))

strord <- mask(strord, mask_mexico, inverse=T, updatevalue=NA)
flo1k_mean <- mask(flo1k_mean, mask_mexico, inverse=T, updatevalue=NA)
writeRaster(flo1k_mean, paste0(path, "/layers_NA/flo1k_mean_mex_masked.tif"))


cl <- makePSOCKcluster(9, outfile="")
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

### Get mean values per stream order
flow_all <- foreach(i=seq(1:9), .combine=rbind, .packages="raster")  %dopar% { 
  new_i <- i-1
  
  mask <- strord
  values(mask) <- ifelse(values(mask)==new_i, 1, NA)
  number_pixels <- cellStats(mask, "sum")
  
  
  flo1k_mean_tmp <- mask(flo1k_mean, mask)
  avg_flow_mean <- cellStats(flo1k_mean_tmp, "mean")

  volume_per_strord <- number_pixels*avg_flow_mean
  
  myres <- data.frame(str_ord=new_i,
                      number_pixels=number_pixels,
                      avg_flow_mean=avg_flow_mean,
                      volume_per_strord=volume_per_strord)
  myres
}

### Export table
write.csv(flow_all, paste0(path, "/estimated_volume_in_stream_orders.csv"), quote = F, row.names = F)

