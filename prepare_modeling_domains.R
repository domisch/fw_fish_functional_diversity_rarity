

### R --slave -f  prepare_fw_fish.R  --args $spp   --no-save  --no-restore  --silent
args <- commandArgs()
print(args)
spp <- args[6]


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
rasterOptions(tmpdir=paste0(path_spp, "/R_temp/"))
if(Sys.info()[["sysname"]]=="Linux") rasterOptions(tmpdir=paste0(path_spp, "/R_temp/"))

### Load layers
load(paste0(path, "/global_layers/layers_global.RData"))



###----------------------------------------------#
### -------- MODEL EXTENT AND RANGEMAP ----------
###----------------------------------------------#


### Load points (all_fish_NA), all data cleaned, names checked and points snapped to entire NA stream network 

load(paste0(path, "/all_snapping_1934-2016.RData")) 
all_fish_NA_sp <- na_fish_1km_year
all_fish_NA_sp$raw_name_in_shp <- gsub(" ", "_", all_fish_NA_sp$raw_name_in_shp) 
coordinates(all_fish_NA_sp) <- c("longitude", "latitude")
proj4string(all_fish_NA_sp) <- "+proj=longlat +ellps=WGS84"


### Load the fish range maps of central and north america
fish_poly <- readShapePoly(paste0(path, "/rangemaps/fish_poly_pageburr_natureserv.shp"))
proj4string(fish_poly) <- "+proj=longlat +ellps=WGS84"

### Subset range map
rangemap <- subset(fish_poly, latin == spp) 
proj4string(rangemap) <- "+proj=longlat +ellps=WGS84"
writePolyShape(rangemap, paste0(path_spp, "/", spp, "_rangemap.shp"))
cat(showWKT(proj4string(rangemap)),file=paste0(path_spp, "/", spp, "_rangemap.prj"))


### Crop the layers to the extent of the rangemap unit + add 5 degrees
ext_vars <- extent(rangemap) # +5 
ext_vars@xmin <- ext_vars@xmin - 5
ext_vars@xmax <- ext_vars@xmax + 5
ext_vars@ymin <- ext_vars@ymin - 5
ext_vars@ymax <- ext_vars@ymax + 5


### Align the extent to the rasters
ext_vars <- alignExtent(ext_vars, layers_global, snap='in')
### Create a raster template
template_domain <- crop(template_NorthAmerica, ext_vars, snap="in") # is the stream network, rangemap + 5x5 degree rectangle



###-------------------------------------------------#
###-------- Prepare species data for models ---------
###-------------------------------------------------#


### Subset points, get only focal species for mapping:
pts_temp <- subset(all_fish_NA_sp, all_fish_NA_sp$raw_name_in_shp == spp)
writePointsShape(pts_temp, paste0(path_spp, "/", spp, "_points"))
cat(showWKT(proj4string(pts_temp)),file=paste0(path_spp, "/", spp, "_points.prj"))



### Distance from range map to all other cells along the stream network --> connectivity 
rangemap$presence <- 1
rangemap_r <- rasterize(rangemap, template_domain, field="presence", small=T, na.rm=T, background=0)  
### Prepare the distance layer 
rangemap_r <- mask(rangemap_r, template_domain)




### Use global level 8 watershed raster (balance between detail and connectivity)
### Crop to the modeling domain
wshed_id <- raster(paste0(path, "/hybas_lev00_all_8.tif"))
tmp_ex <- alignExtent(template_domain, wshed_id, snap="in") 
wshed_id <- crop(wshed_id, tmp_ex, snap="in"); rm(tmp_ex)

if (extent(wshed_id)== extent(template_domain)) {
  ### mask to get streams
  wshed_id <- mask(wshed_id, template_domain)
} else {
  ### if watershed-id is smaller that domain, grow with NA to make identical extent 
  wshed_id <- extend(wshed_id, extent(template_domain), value=NA)
  wshed_id <- mask(wshed_id, template_domain)
}


### Check how many rangemap cells fall within each watershed
rangemap_in_whsed <- as.data.frame(zonal(rangemap_r, wshed_id, fun='sum', digits=0, na.rm=TRUE))
### Number of stream cells in each watershed
size_whsed <- as.data.frame(zonal(wshed_id, wshed_id, fun='count', digits=0, na.rm=TRUE)) # watersheds have been masked to streams
### Add to main df
rangemap_in_whsed$size_whsed <- size_whsed$count; rm(size_whsed)

### Get new rangemap: if majority of cells in each watersehd is covered by rangemap, keep these
rangemap_in_whsed$rangemap_corr <- ifelse(rangemap_in_whsed$sum >= rangemap_in_whsed$size_whsed/2, 1, 0)

### Replace the "incorrect" rangemap values with zero (and in cae majority is covered, write 1 1)
### First get the watershed_id for each rangemap==1 cell
rangemap_r_corr <- rangemap_r
values(rangemap_r_corr) <- ifelse(values(rangemap_r)==1, values(wshed_id), 0)
# writeRaster(rangemap_r_corr, paste0(path_spp, "/rangemap_r_corr_temp.tif"), overwrite=T)
rangemap_r_corr_tmp <- subs(rangemap_r_corr, rangemap_in_whsed[c("zone","rangemap_corr")], by="zone")
### Add the rest of the stream network
### Fill NA's with zero
rangemap_r_corr_tmp[is.na(rangemap_r_corr_tmp)] <- 0
### Mask stream network again
rangemap_r_corr_final <- mask(rangemap_r_corr_tmp, template_domain)
writeRaster(rangemap_r_corr_final, paste0(path_spp, "/rangemap_r_corr_final.tif"), overwrite=T)

rm(rangemap_r_corr)
rm(rangemap_r_corr_tmp)
gc()
### Calculate distance
dist_rmap_corr <- gridDistance(rangemap_r_corr_final, origin=1, omit=NA)
# writeRaster(dist_rmap_corr, paste0(path_spp, "/dist_rmap_corr.tif"), overwrite=T)


### Correct for points: extend the connectivity to these reaches as well
### Get connectivity for points
pts_temp$presence <- 1 # Convert presences to raster and set to 1
pts_r <- rasterize(pts_temp, template_domain, field="presence", na.rm=T, background=0) 

### Only for additional area: add to distance calculation
pts_r_corr <- pts_r # copy for output
### Replace values only if presence points were within non-connected
values(pts_r_corr) <- ifelse(values(pts_r) == 1 & values(rangemap_r_corr_final) == 0, 1, values(rangemap_r_corr_final))
writeRaster(pts_r_corr, filename=paste0(path_spp, "/rangemap_points_corrected.tif"), overwrite=T)
gc()
dist_rmap_corr <- gridDistance(pts_r_corr, origin=1, omit=NA)
writeRaster(dist_rmap_corr, paste0(path_spp, "/dist_rmap_corr.tif"), overwrite=T) #overwrite 
gc()




### Get a stream-ID dataframe (=id_short, only stream cells)
my_stream <- dist_rmap_corr # was dist_rmap_crop, change the stream cell -extent, only cells are kept that are connected to rangemap 
### To go back (using all cells, use " my_stream" from above..) 
stream_cells <- raster::as.data.frame(my_stream, xy=T) # only stream cells that are connected to range map
stream_cells <- as.data.frame(na.omit(stream_cells)) # separated due to update in raster-package...?
stream_cells <- stream_cells[-3]
### Add the short and long ID
stream_cells$id_short <- seq(1:nrow(stream_cells))
colnames(stream_cells)[3] <- "id_short"
max(stream_cells$id_short)



grid_id <- template_domain
grid_id[] <- seq(1:ncell(dist_rmap_corr))
names(grid_id) <- "grid_id"


### Get grid_ID of entire raster for the streams...
temp <-  raster::as.data.frame(raster::extract(grid_id, stream_cells[c("x", "y")])) # df=T
colnames(temp) <- "id_long"
stream_cells <- cbind(stream_cells, temp); rm(temp)


temp <- as.data.frame(raster::extract(grid_id_global, stream_cells[c("x", "y")])) # df=T
colnames(temp) <- "grid_id_global"
stream_cells <- cbind(stream_cells, temp); rm(temp)




###-------------------------------------------------#
###-------- Prepare species data for models ---------
###-------------------------------------------------#


### Crop points to extent
all_fish_sp <- crop(all_fish_NA_sp, ext_vars, snap = "in")
length(unique(all_fish_sp$raw_name_in_shp)) # numer species sampled within this extent


### Get the presences and temp.absences of other pecies
pts_sp <- subset(all_fish_sp, all_fish_sp$raw_name_in_shp == spp) # all_fish_sp is cropped
pts_abs <- subset(all_fish_sp, all_fish_sp$raw_name_in_shp != spp)

### Convert to raster
pres_stream <- rasterize(pts_sp, my_stream, field="raw_name_in_shp", small=T, fun="count", na.rm=T)
abs_stream <- rasterize(pts_abs, my_stream, field="raw_name_in_shp", small=T, fun="count", na.rm=T) 


### Get the number of years a single site (grid) was visited
### Use the "..."  argument to accept na.rm=T
n_year_stream <- rasterize(all_fish_sp, my_stream, field="year_num", small=T, fun=function(x, ...) {length(unique(na.omit(x)))}, na.rm=T)

### Mask the stream network
pres_stream <- mask(pres_stream, my_stream) 
abs_stream <- mask(abs_stream, my_stream) 
n_year_stream <- mask(n_year_stream, my_stream) 


### Change names
names(pres_stream) <- "pres_stream"
names(abs_stream) <- "abs_stream"
names(n_year_stream) <- "n_year_stream"

### Add the presences and absences to main df
tmp_stack <- stack(pres_stream, abs_stream, n_year_stream)
names(tmp_stack) <- c("pres", "abs", "n_year")

library(foreach)
library(doParallel)

cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
getDoParWorkers() 

tmp <- foreach(rasname = iter(names(tmp_stack)), .combine='cbind', .packages = "raster") %dopar% {
  tmp_extr <- extract(tmp_stack[[rasname]], stream_cells$id_long, df=T)
  tmp_extr[2] # return only the second column with the values
}
stopCluster(cl)

stream_cells <- cbind(stream_cells, tmp); rm(tmp)




###--------------------------------------------#
###----------- Rangemap distance --------------
###--------------------------------------------#

start=0 # 
end=150000 # in meters, 150 km 

### Subset the data
dist_table <- extract(dist_rmap_corr, stream_cells$id_long, df=T)
dist_table<- dist_table[2]
colnames(dist_table) <- "dist_m"
dist_table$id_long <- stream_cells$id_long

within_buff <- subset(dist_table, dist_table$dist_m > 0 & dist_table$dist_m < end)
colnames(within_buff)[1] <- "dist_buff_m"


### Apply the scaling from 0-1
library(scales); detach(package:arm, unload=T)
within_buff$scale_linear <- rescale(within_buff$dist_buff_m, to=c(0, 1))
within_buff$scale_quad <- rescale((within_buff$dist_buff_m)^2, to=c(0, 1))
within_buff$scale_exp <- rescale(exp((within_buff$dist_buff_m / 10000)), to=c(0, 1))
within_buff$scale_logit <- rescale(logit(within_buff$dist_buff_m, min=0, max=end), to=c(0, 1))



### Merge back
stream_cells <- merge(stream_cells, within_buff, by="id_long", all.x=T) ; rm(within_buff)
### Add the raw distance
stream_cells <- merge(stream_cells, dist_table, by="id_long", all.x=T) ; rm(dist_table)

stream_cells$scale_linear <- ifelse(stream_cells$dist_m  == 0, 0, stream_cells$scale_linear)
stream_cells$scale_linear <- ifelse(stream_cells$dist_m  >= end, 1, stream_cells$scale_linear)


stream_cells$scale_quad <- ifelse(stream_cells$dist_m  == 0, 0, stream_cells$scale_quad)
stream_cells$scale_quad <- ifelse(stream_cells$dist_m  >= end, 1, stream_cells$scale_quad)


stream_cells$scale_exp <- ifelse(stream_cells$dist_m  == 0, 0, stream_cells$scale_exp)
stream_cells$scale_exp <- ifelse(stream_cells$dist_m  >= end, 1, stream_cells$scale_exp)

stream_cells$scale_logit <- ifelse(stream_cells$dist_m  == 0, 0, stream_cells$scale_logit)
stream_cells$scale_logit <- ifelse(stream_cells$dist_m  >= end, 1, stream_cells$scale_logit)


### Invert (close is 1, far away 0...1)
stream_cells$scale_linear <- (stream_cells$scale_linear  * (-1)) +1
stream_cells$scale_quad <- (stream_cells$scale_quad  * (-1)) +1
stream_cells$scale_exp <- (stream_cells$scale_exp  * (-1)) +1
stream_cells$scale_logit <- (stream_cells$scale_logit  * (-1)) +1



library(foreach)
library(doParallel)
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
getDoParWorkers() 



env <- foreach(rasname = iter(names(layers_global)), .combine='cbind', .packages = "raster") %dopar% {
  tmp_extr <- extract(layers_global[[rasname]], stream_cells$grid_id_global, df=T)
  tmp_extr[2] # return only the second column with the values
}
stopCluster(cl)

### Lakes will be coded as 10 (highest stream order is 9)
env$stream_order_lakes0 <- ifelse(env$stream_order_lakes0 == 0, 10, env$stream_order_lakes0)


### Standardize by two standard deviations
detach("package:scales", unload=T) # has also the "rescale" function
library(arm)
rescale_fun <- function(x) {rescale(x, "full")}
env_s <- as.data.frame(apply(subset(env, select=-c(lentic_lotic01, grid_id_global_all_cells)), 2, rescale_fun)) # arm -package

### Add the other non-transformed layers (numeric for correlation analysis)
env_s$lentic_lotic01 <- as.factor(env$lentic_lotic01)


### Remove geology columns that have only NA's (=no coverage in the modeling domain)
env_s <- env_s[, colSums(is.na(env_s)) != nrow(env_s)]


### Environmental and species data fit?
if (dim(env_s)[1] != dim(stream_cells)[1]) { cat("Check data")}



### ---- Subset data 70/30 % for validation -----------
`%ni%` <- Negate(`%in%`)
temp  <- subset(stream_cells, stream_cells$presences>=1)
write.table(nrow(temp), paste0(path_spp, "/number_unique_records.txt"), row.names=F, col.names=F, quote=F)


temp2 <- temp[sample.int(nrow(subset(temp, temp$presences>=1)), round(sum(temp$presences>=1) * 0.3, 0), replace=F), ]


stream_cells$validation <- ifelse(stream_cells$id_long %in% temp2$id_long, 1, 0) # get 30% of presences
stream_cells$validation <- ifelse(stream_cells$presences == 0, 1, stream_cells$validation) # get absences

stream_cells$fitting <- ifelse(stream_cells$id_long %ni% temp2$id_long, 1, 0) # get 70% of presences
stream_cells$fitting <- ifelse(stream_cells$presences == 0, 1, stream_cells$fitting) # get absences


### Create the df for the range-wide predictions
data <- cbind(stream_cells, env_s)

### Create the df for the validation data. 30% of the presences and all absences
vdata <- droplevels(subset(data, data$validation == 1))
coordinates(vdata) <- c("x", "y")

### Create the df for the model, i.e. only those cells with observation AND NOT selected for the validation
fdata <- droplevels(subset(data, data$fitting == 1)) 

# mean(fdata$id_long %in% vdata$id_long)
### Make occurrences binary pres/abs
fdata$obs_glm <- ifelse(fdata$presences >1, 1, fdata$presences) 
### Clean data frames
data <- droplevels(data)


save(fdata, data, vdata, neighbours, number_neigh, file=(paste0(path_spp, "/", spp, ".RData")))


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



