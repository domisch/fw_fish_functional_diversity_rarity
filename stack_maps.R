

### North America extent
NA_extent=extent(template_NorthAmerica) # from prepare_env_layers.R



stacked_out <- foreach(i=myspp_sequence[["species"]], .errorhandling="stop",.verbose=T, .combine=stack, .packages=c("raster")) %dopar% {
  rasterOptions(tmpdir="/mnt/domisch/data/fw_fish/R_temp_delete")
  spp_i <- subset(myspp, species ==i)
  if (spp_i$filled_data == 1) {
    cat("#-------- Loading rangemap of", paste0(i), "\n")
    r <- raster(paste0(path, "/species_folders/", i, "/rangemap_corrected.tif"))
    r <-  extend(r, NA_extent, value=NA)
    
  } else {
    cat("#-------- Loading modeled results of", paste0(i), "\n")
    r <- raster(paste0(path, "/species_folders/", i, "/spatial_priors/binary_best_model/bin_sum.tif"))
    r <-  extend(r, NA_extent, value=NA)
  }
}
  