


###-------------------------------------------------#
###---------- Prepare environmental data -----------
###-------------------------------------------------#


library(raster)

### Set variables for the cluster:

if(Sys.info()[["sysname"]]=="Linux") path_var="/lustre/scratch/client/fas/sbsc/fw_fish/"

### Set the temporary folder for the raster-package (will be deleted after each species) 
dir.create(paste(path_var, "R_temp_var_KEEP/", sep=""))
rasterOptions(chunksize=1000,maxmemory=1000, progress="text", tmpdir=paste(path_var, "R_temp_var_KEEP/", sep=""))
### Delete the files in the R-tmp folder
unlink(paste(path_var, "R_temp_var_KEEP/*.gr*", sep=""), recursive=T, force = T)


str_order <- raster(paste(path_var, "additional_layers/stream_order_lakes0.tif", sep=""))
lentic_lotic01 <- raster(paste(path_var, "additional_layers/lentic_lotic01.tif", sep=""))
flow_length <- raster(paste(path_var, "flow_accumulation/flow_length.tif", sep=""))
flow_acc <- raster(paste(path_var, "flow_accumulation/flow_acc.tif", sep=""))

### Watershed topography
slope_site <- raster(paste(path_var, "additional_layers/slope_site.tif", sep="")) #  site-level slope
dem_site <- raster(paste(path_var, "additional_layers/dem_site.tif", sep=""))
slope_range <- raster(paste(path_var, "slope/slope_range.tif", sep="")) #  site-level slope
dem_range <- raster(paste(path_var, "elevation/dem_range.tif", sep=""))
bifurcation_ratio <-  raster(paste(path_var, "additional_layers/bifurcation_ratio.tif", sep=""))

### Climate
hydro_avg_01 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_01.tif", sep=""))
hydro_avg_07 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_07.tif", sep=""))
hydro_avg_10 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_10.tif", sep="")) # temp. warmest quarter
hydro_avg_11 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_11.tif", sep="")) # temp. coldest quarter
hydro_avg_12 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_12.tif", sep=""))
hydro_avg_15 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_15.tif", sep=""))
hydro_avg_16 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_16.tif", sep="")) # precip. wettest quarter
hydro_avg_17 <- raster(paste(path_var, "hydroclim_average+sum/hydro_avg_17.tif", sep="")) # precip. driest quarter

### Land use
# Class 1: Evergreen/deciduous needleleaf trees
# Class 2: Evergreen broadleaf trees
# Class 3: Deciduous broadleaf trees
# Class 4: Mixed/other trees
# Class 5: Shrubs
# Class 6: Herbaceous vegetation
# Class 7: Cultivated and managed vegetation
# Class 8: Regularly flooded shrub/herbaceous vegetation
# Class 9: Urban/built-up
# Class 10: Snow/ice
# Class 11: Barren lands/sparse vegetation
# Class 12: Open water
lc_wavg_01 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_01.tif", sep=""))
lc_wavg_02 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_02.tif", sep=""))
lc_wavg_03 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_03.tif", sep=""))
lc_wavg_04 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_04.tif", sep=""))
lc_wavg_05 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_05.tif", sep=""))
lc_wavg_06 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_06.tif", sep=""))
lc_wavg_07 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_07.tif", sep=""))
lc_wavg_08 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_08.tif", sep=""))
lc_wavg_09 <- raster(paste(path_var, "landcover_weighted_average/lc_wavg_09.tif", sep=""))
  

### Geology
# geo_wsum_05 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_05.tif", sep=""))
# geo_wsum_11 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_11.tif", sep=""))
# geo_wsum_14 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_14.tif", sep=""))
geo_wsum_20 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_20.tif", sep=""))
# geo_wsum_28 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_28.tif", sep=""))
geo_wsum_42 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_42.tif", sep=""))
geo_wsum_47 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_47.tif", sep=""))
geo_wsum_48 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_48.tif", sep=""))
geo_wsum_56 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_56.tif", sep=""))
# geo_wsum_70 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_70.tif", sep=""))
# geo_wsum_72 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_72.tif", sep=""))
geo_wsum_76 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_76.tif", sep=""))
geo_wsum_80 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_80.tif", sep=""))
# geo_wsum_88 <- raster(paste(path_var, "geology_weighted_sum/geo_wsum_88.tif", sep=""))



### Soil - average
# soil_min_01.tif     Soil organic carbon (ORCDRC) across sub-catchment
# soil_min_02.tif	    Soil pH in H2O (PHIHOX) across sub-catchment
# soil_min_03.tif	    Sand content mass fraction (SNDPPT) across sub-catchment
# soil_min_04.tif	    Silt content mass fraction (SLTPPT) across sub-catchment
# soil_min_05.tif	    Clay content mass fraction (CLYPPT) across sub-catchment
# soil_min_06.tif	    Coarse fragments (> 2 mm fraction) volumetric (CFRVOL) across sub-catchment
# soil_min_07.tif	    Cation exchange capacity (CEC) across sub-catchment
# soil_min_08.tif	    Bulk density of the fine earth fraction (BLD) across sub-catchment
# soil_min_09.tif	    Depth to bedrock (R horizon) up to maximum 240 cm (BDRICM) across sub-catchment
# soil_min_10.tif	    Predicted probability of occurence (0-100%) of R horizon (BDRLOG) across sub-catchment

soil_avg_01 <- raster(paste(path_var, "soil_average/soil_avg_01.tif", sep=""))
soil_avg_02 <- raster(paste(path_var, "soil_average/soil_avg_02.tif", sep=""))
soil_avg_03 <- raster(paste(path_var, "soil_average/soil_avg_03.tif", sep=""))
soil_avg_04 <- raster(paste(path_var, "soil_average/soil_avg_04.tif", sep=""))
soil_avg_05 <- raster(paste(path_var, "soil_average/soil_avg_05.tif", sep=""))
soil_avg_06 <- raster(paste(path_var, "soil_average/soil_avg_06.tif", sep=""))
soil_avg_07 <- raster(paste(path_var, "soil_average/soil_avg_07.tif", sep=""))
soil_avg_08 <- raster(paste(path_var, "soil_average/soil_avg_08.tif", sep=""))
soil_avg_09 <- raster(paste(path_var, "soil_average/soil_avg_09.tif", sep=""))
soil_avg_10 <- raster(paste(path_var, "soil_average/soil_avg_10.tif", sep=""))




### Load global grid ID
grid_id_global <- raster(paste(path_var, "/additional_layers/grid_id_global_all_cells.tif", sep=""))

### Load 1km (aggregated from 30m) Landsat layer  
proportion_freshwater <- raster(paste(path_var, "/additional_layers/land_frequency_GFC2014_landsat_30m.tif", sep=""))
proportion_freshwater <- crop(proportion_freshwater, str_order, snap="in") # HydroSHEDS extent
names(proportion_freshwater) <- "proportion_freshwater"

# compareRaster(bio1_weighted, dem_range)

### Stack variables
layers_global <- stack(str_order,
                   lentic_lotic01,
                   flow_length,
                   flow_acc,
                   slope_site,
                   slope_range,
                   dem_range,
                   dem_site,
                   bifurcation_ratio,
                   hydro_avg_01, 
                   hydro_avg_07,
                   hydro_avg_10,
                   hydro_avg_11,
                   hydro_avg_12,
                   hydro_avg_15,
                   hydro_avg_16,
                   hydro_avg_17,
                   lc_wavg_01, 
                   lc_wavg_02,
                   lc_wavg_03,
                   lc_wavg_04,
                   lc_wavg_05,
                   lc_wavg_06,
                   lc_wavg_07,
                   lc_wavg_08,
                   lc_wavg_09, 
                   geo_wsum_20,
                   geo_wsum_42,
                   geo_wsum_47, 
                   geo_wsum_48, 
                   geo_wsum_56, 
                   geo_wsum_76, 
                   geo_wsum_80,
                   soil_avg_01,
                   soil_avg_02,
                   soil_avg_03,
                   soil_avg_04,
                   soil_avg_05,
                   soil_avg_06,
                   soil_avg_07,
                   soil_avg_08,
                   soil_avg_09,
                   soil_avg_10,
                   proportion_freshwater,
                   grid_id_global)


### Remove raster layers to save space   
rm(list=ls(pattern="hydro"))
rm(list=ls(pattern="lc"))
rm(list=ls(pattern="geo"))
rm(list=ls(pattern="geo"))
rm(list=ls(pattern="soil"))
rm(list=ls(pattern="slope"))
rm(list=ls(pattern="dem"))
rm(list=ls(pattern="flow"))
rm("str_order", "lentic_lotic01", "bifurcation_ratio")


### load the mask for North America (from HydroSHEDS, North America and Central America merged)
load(paste(path_var, "/temp_mask_na_ca.Rdata", sep=""))

### Align the extent to the rasters
temp_mask_na_ca <- alignExtent(temp_mask_na_ca, layers_global, snap='in')
### Crop the layer to the extent
template_NorthAmerica <- crop(layers_global[["stream_order_lakes0"]], temp_mask_na_ca, snap='in')

### Save data
save(layers_global, grid_id_global, template_NorthAmerica, file=paste(path_var, "layers_global.RData", sep=""))


### Clip Mexico
library(raster)
mex <- shapefile("/data/domisch/adm_shapefile/gadm36_levels/gadm36_0.shp")
mex <- subset(mex, NAME_0=="Mexico")
save(mex, file="/data/domisch/data/fw_fish/mexico_shape.RData")



