


###=======================================================#
### ---- Bivariate plot of rangesize and distictness -----
###=======================================================#


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


path <- "/mnt/domisch//data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0("/data/domisch/R_temp_delete"))
getwd()



### Code modified from here:
# http://rfunctions.blogspot.de/2015/03/bivariate-maps-bivariatemap-function.html

col.matrix<-colmat(nquantiles=10, upperleft="blue", upperright="yellow", bottomleft="green", bottomright="red", xlab="My x label", ylab="My y label")
save(colmat, col.matrix, file=paste0(path, "/col.matrix_bivariate.RData"))
load(paste0(path, "/col.matrix_bivariate.RData"))


### Customized function:
bivariate.map_customized<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  
  quanmean<-getValues(rasterx)
  quanmean[quanmean==0] <- NA # avoid zero inflation, change zero to NA
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  brks <- brks + seq_along(brks) * 0.0000001 # add some jitter 
  # brks<-with(temp, unique(quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))) # keep only unique
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
 
  quanvar<-getValues(rastery)
  quanmean[quanmean==0] <- NA # avoid zero inflation, change zero to NA
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, unique(quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles)))))
  brks <- brks + seq_along(brks) * 0.0000001  # add some jitter 
  # brks<-with(temp, unique(quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))) # keep only unique
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}



### Load raster layers
range_size_VOLUME <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/med_RS_volume_sprior807_mex_masked_lakes_avg_gdal.tif"))
eco_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_median_ED_per_grid_mex_masked_lakes_avg_gdal.tif"))

# volume FD
bivmap_rangesize_VOLUME_eco_ED <- bivariate.map_customized(range_size_VOLUME, eco_ED, colormatrix=col.matrix, nquantiles=10); gc()
writeRaster(bivmap_rangesize_VOLUME_eco_ED, paste0(path, "/stacked_outputs/lakes_avg_gdal/bivmap_rangesize_VOLUME_eco_ED_mex_masked_lakes_avg.tif"), overwrite=T)


