


###====================================================================#
###--- Get median elevation, latitude and stream order per species ----
###====================================================================#

# https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
# https://drsimonj.svbtle.com/quick-plot-of-all-variables 

# quit(save = "no")
# R
# Sys.sleep(3600)
# Sys.sleep(5)
path <- "/mnt/domisch/data/fw_fish"
setwd(path)
getwd()

path_FIGURES <- paste0(path, "/figures")
dir.create(path_FIGURES)
path_OUT <- paste0(path, "/traits_eco")



if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = T) ; library(rgdal)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("hexbin")) { install.packages("hexbin", dependencies = T) ; library(hexbin)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))



source(paste0(path, "/scripts/raster.as.data.table_coord.R"))



### Load layers
SR <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/sprior_sum807_mex_masked_lakes_avg_gdal.tif"))
FDpg <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDpg_mex_masked_lakes_avg_gdal.tif"))
FDw <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDw_mex_masked_lakes_avg_gdal.tif"))
FR_rs <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_rangesize_mex_masked_lakes_avg_gdal.tif"))
FR_vol <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_volume_mex_masked_lakes_avg_gdal.tif"))
RS_med <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/med_RS_sprior807_mex_masked_lakes_avg_gdal.tif"))
RS_rarity <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_mex_masked_lakes_avg_gdal.tif"))
median_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_median_ED_per_grid_mex_masked_lakes_avg_gdal.tif"))

### These need to cover the entire range:
dem <- raster(paste0(path, "/layers_NA/dem_site.tif"))
str_order <- raster(paste0(path, "/layers_NA/str_order_lakes0.tif"))
flow <- raster(paste0(path, "/FLO1k/mean_flow_1km_1960_2015.tif"))
flow <- crop(flow, str_order)


### Get elevation per stream and lake grid cell
na_mask <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/hydrographic_network.tif")) # NoData 255, lakes = 0, stream=1

### Crop to mask 
dem_mask <- crop(dem, na_mask)
str_order_mask <- crop(str_order, na_mask)
flow_mask <- crop(flow, na_mask)

### Get same coverage as the mask
dem_mask <- mask(dem_mask, na_mask)
str_order_mask <- mask(str_order_mask, na_mask)
flow_mask <- mask(flow_mask, na_mask)

my_stack <- stack(na_mask, dem_mask, str_order_mask, flow_mask)


### Get as datatable and export
my_stack_dt <- as.data.table(my_stack, na.rm=T)
save(my_stack_dt, file=paste0(path, "/dem_str_order_flow_stream_lake_per_grid_mex_masked.RData")) 
### --> continue in functional_diversity.R



### Make sure that the grid-ID matches to species table
e <- extent(-145, -52, 5, 60)

### Load raster/domain template
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes_mex_masked.RData"))
grid_id_long <- extend(grid_id_long, e)
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
names(domain_cells) <- "grid_id_long"
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, grid_id_long) # set key


domain_raster <- extend(domain_raster, e, value=NA)
SR <- extend(SR, e, value=NA)
FDpg <- extend(FDpg, e, value=NA)
FDw <- extend(FDw, e, value=NA)
dem <- extend(dem, e, value=NA)
str_order <- extend(str_order, e, value=NA)
flow <- extend(flow, e, value=NA)
FR_rs <- extend(FR_rs, e, value=NA)
FR_vol <- extend(FR_vol, e, value=NA)
RS_med <- extend(RS_med, e, value=NA)
RS_rarity <- extend(RS_rarity, e, value=NA)




### Convert rasters to data.table
mystack <- stack(domain_raster, SR, FDpg, FDw, dem, str_order, flow, FR_rs, FR_vol, RS_med, RS_rarity)
mystack_dt <- as.data.table.raster(mystack,  xy=T)
mystack_dt$seq_id <- seq.int(1:nrow(mystack_dt))
save(mystack_dt, file=paste0(path, "/stacked_outputs/eco_datatable_domain_raster_SR_FDpg_FDw_dem_str_order_latitude.RData"))



### Check correlation 
names(mystack_dt)
cor.test(mystack_dt$sprior_sum807_mex_masked_lakes_avg_gdal, mystack_dt$eco_FDw_mex_masked_lakes_avg_gdal, method="spearman", use="complete.obs", exact = F) # 0.65
cor.test(mystack_dt$sprior_sum807_mex_masked_lakes_avg_gdal, mystack_dt$y, method="spearman", use="complete.obs", exact = F) # 0.42
cor.test(mystack_dt$eco_FDw_mex_masked_lakes_avg_gdal, mystack_dt$y, method="spearman", use="complete.obs", exact = F) # 0.11
cor.test(mystack_dt$sprior_sum807_mex_masked_lakes_avg_gdal, mystack_dt$dem_site, method="spearman", use="complete.obs", exact = F) # 0.34
cor.test(mystack_dt$eco_FDw_mex_masked_lakes_avg_gdal, mystack_dt$dem_site, method="spearman", use="complete.obs", exact = F) # 0.35



### Load large species-per-grid table
tmp_rs <- fread("stacked_outputs/sprior_807_probs_per_grid.txt")
setnames(tmp_rs, c("probs", "range_size", "grid_id_long", "species"))
setkey(tmp_rs, grid_id_long)

### Merge elevation, latitude and stream order to large table
tmp <- merge(tmp_rs, mystack_dt, by="grid_id_long", all.y=T) # masked mex. for SR etc, not env. layers (=dem, flow etc)
gc()

### Check if any duplicates
tmp <- tmp[!duplicated(tmp),]
gc()



### Get median values per species
med_latitude_per_species <- tmp[, as.double(median(y, na.rm = TRUE)),by = species]
names(med_latitude_per_species) <- c("species", "median_latitude")
med_elevation_per_species <- tmp[, as.double(median(dem_site, na.rm = TRUE)),by = species]
names(med_elevation_per_species) <- c("species", "median_elevation")
med_flow_per_species <- tmp[, as.double(median(mean_flow_1km_1960_2015, na.rm = TRUE)),by = species]
names(med_flow_per_species) <- c("species", "median_flow")


### Get the mode (majority) stream order
tmp_stream_ord <- subset(tmp, select=c(species, str_order_lakes0))
tmp_stream_ord <- na.omit(tmp_stream_ord)
mode_str_order_per_species <- setkey(tmp_stream_ord[, list(freq = .N), by=list(species, str_order_lakes0)], 
                                    species, freq)[J(unique(species, na.rm = TRUE)), mult="last"]
mode_str_order_per_species <- subset(mode_str_order_per_species, select=-c(freq))



### Combine and export
out <- merge(med_latitude_per_species, med_elevation_per_species, by="species")
out <- merge(out, mode_str_order_per_species, by="species")
write.csv(out, paste0(path, "/stacked_outputs/lakes_avg_gdal/median_latitude_elevation_mode_stream_order_per_spp.csv"), row.names = F, quote = F )
write.csv(med_flow_per_species, paste0(path, "/stacked_outputs/lakes_avg_gdal/median_flow_per_spp.csv"), row.names = F, quote = F )


