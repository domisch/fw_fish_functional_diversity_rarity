

###========================================#
###--- Create stream network template --
###========================================#

### Create new stream network template without Mexico
r <- raster(paste0(getwd(), "/global_layers/additional_layers/lentic_lotic01.tif"))
e <- extent(-145, -52, 25, 60) # North America
r <- crop(r, e)
writeRaster(r, paste0(getwd(), "/stacked_outputs/NA_str_network_NoData255.tif"), overwrite=T)

### Mask Mexico
mask_mexico <- shapefile(paste0(getwd(), "/shape_mexico/mask_mexico2.shp"))

r <- mask(r, mask_mexico, inverse=T, updatevalue=NA)

writeRaster(r, paste0(getwd(), "/stacked_outputs/NA_str_network_mex_masked_NoData255.tif"), NAFlag=255, overwrite=T)