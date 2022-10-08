

### Get the argumnets from bash --> species name is object number 6
### R --slave -f  prepare_na_fish.R  --args $spp   --no-save  --no-restore  --silent
args <- commandArgs()
print(args)
spp <- gsub("_", " ", args[6]) # get original species name without underscore..



###----------------------------------------------#
### ----- PREPARE MODEL EXTENT AND DATA ---------#
###----------------------------------------------#

if (!require("maptools")) { install.packages("maptools", dependencies = TRUE) ; library(maptools)}
if (!require("foreign")) { install.packages("foreign", dependencies = TRUE) ; library(foreign)}
if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("sp")) { install.packages("sp", dependencies = TRUE) ; library(sp)}
if (!require("rgeos")) { install.packages("rgeos", dependencies = TRUE) ; library(rgeos)}
if (!require("dismo")) { install.packages("dismo", dependencies = TRUE) ; library(dismo)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = TRUE) ; library(rgdal)}
if (!require("rasterVis")) { install.packages("rasterVis", dependencies = TRUE) ; library(rasterVis)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("PresenceAbsence")) { install.packages("PresenceAbsence", dependencies = TRUE) ; library(PresenceAbsence)}
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)}
if (!require("arm")) { install.packages("arm", dependencies = TRUE) ; library(arm)}
if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(ncdf4) }
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("snow")) { install.packages("snow", dependencies = TRUE) ; library(snow)}
if (!require("doSNOW")) { install.packages("doSNOW", dependencies = TRUE) ; library(doSNOW)}
if (!require("gdistance")) { install.packages("gdistance", dependencies = TRUE) ; library(gdistance)}
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("reshape")) { install.packages("reshape", dependencies = TRUE) ; library(reshape)}




### Set variables for the cluster:
if(Sys.info()[["sysname"]]=="Linux") path="/lustre/scratch/client/fas/sbsc/na_fish_all"
if(Sys.info()[["sysname"]]=="Linux") n_cores=8

### Create folder for species

path_spp <- paste0(path, "/species_folders/", spp)
dir.create(path_spp)

### Create and set the temporary folder for the raster-package (will be deleted after each species) 
dir.create(paste0(path_spp, "/R_temp/"))

### Delete the files in the R-tmp folder
unlink(paste0(path_spp, "/R_temp/*.gr*"), recursive=T, force = T)
### Delete the old predictor-correlation-files
unlink(paste0(path_spp, "/*.csv*"), recursive=F, force = T)


### Set raster-package options depending on the platform
rasterOptions(chunksize=1000,maxmemory=1000, tmpdir=paste0(path_spp, "/R_temp/"))
if(Sys.info()[["sysname"]]=="Linux") rasterOptions(tmpdir=paste0(path_spp, "/R_temp/"))
# progress="text"

### Load layers
load(paste0(path, "/global_layers/layers_global.RData"))




###----------------------------------------------#
### -------- MODEL EXTENT AND RANGEMAP ----------
###----------------------------------------------#


### Load the fish range maps of central and north america
fish_poly <- readShapePoly(paste0(path, "/rangemaps/na_fish.shp"))
proj4string(fish_poly) <- "+proj=longlat +ellps=WGS84"

### Subset range map
rangemap <- subset(fish_poly, latin == spp) 
proj4string(rangemap) <- "+proj=longlat +ellps=WGS84"
writePolyShape(rangemap, paste0(path_spp, "/", spp, "_rangemap.shp"))
cat(showWKT(proj4string(rangemap)),file=paste0(path_spp, "/", spp, "_rangemap.prj"))


## Crop the layers to the extent of the rangemap unit + add some area
ext_vars <- extent(rangemap) # +5 
ext_vars@xmin <- ext_vars@xmin - 5
ext_vars@xmax <- ext_vars@xmax + 5
ext_vars@ymin <- ext_vars@ymin - 5
ext_vars@ymax <- ext_vars@ymax + 5


### Align the extent to the rasters
ext_vars <- alignExtent(ext_vars, template_NorthAmerica, snap='in')
### Create a raster template
template_domain <- crop(template_NorthAmerica, ext_vars, snap="in") # is the stream network, rangemap + 5x5 degree rectangle


rangemap$presence <- 1
rangemap_r <- rasterize(rangemap, template_domain, field="presence", small=T, na.rm=T, background=0) 
e <- extent(-145, -52, 5, 60) # North America
rangemap_r <- extend(rangemap_r, e, value=NA)
writeRaster(rangemap_r, paste0(path_spp, "/rangemap_raw_no_correction.tif"), datatype="INT1U", overwrite=T)


### Delete the files in the R-tmp folder
unlink(paste0(path_spp, "/R_temp/*.gr*"), recursive=T, force = T)

### Keep the files you want to keep in a vector before deleting everything
KEEP = c(
"spp_summary", # species list for loop
"spp",  # current species iteratir
"path",  # main directory
"n_cores") # number of cores for cluster (raster package)


REMOVE = ls()
  
### Delete everything but the objects in "KEEP"
rm(list=(REMOVE[is.na(match(REMOVE, KEEP))]))

### Garbage collection - free memory 
gc()




