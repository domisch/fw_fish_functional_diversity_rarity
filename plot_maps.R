




###==============================================#
###---- Plot all maps aggregated by factor x ----
###==============================================#


### Packages and options
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
if (!require("ggpolypath")) { install.packages("ggpolypath", dependencies = T) ; library(ggpolypath)}
if (!require("ggmap")) { install.packages("ggmap", dependencies = TRUE) ; library(ggmap)} # for fortifying shapefiles
if (!require("proto")) { install.packages("proto", dependencies = TRUE) ; library(proto)} # required for fixing holes in polygons
if (!require("ggpolypath")) { install.packages("ggpolypath", dependencies = TRUE) ; library(ggpolypath)} 
if (!require("mapproj")) { install.packages("mapproj", dependencies = TRUE) ; library(mapproj)} 
if (!require("BAMMtools")) { install.packages("BAMMtools", dependencies = TRUE) ; library(BAMMtools)} 
if (!require("broom")) { install.packages("broom", dependencies = TRUE) ; library(broom)} 
if (!require("RStoolbox")) { install.packages("RStoolbox", dependencies = TRUE) ; library(RStoolbox)} 

path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))
getwd()


col <- rev(rainbow(100, start = 0, end = 0.7)) 
mycol=col
# myres=10000
myres=4e6





### Load layers
### richness and range size maps
sprior_sum807 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/sprior_sum807_mex_masked_lakes_avg_gdal.tif"))
average_range_size_rarity <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_mex_masked_lakes_avg_gdal.tif"))
range_size_rarity <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/range_size_rarity_mex_masked_lakes_avg_gdal.tif"))
average_range_size_rarity_volume <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_volume_mex_masked_lakes_avg_gdal.tif"))
range_size_rarity_volume <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/range_size_rarity_volume_mex_masked_lakes_avg_gdal.tif"))
med_RS_sprior807_volume <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/med_RS_volume_sprior807_mex_masked_lakes_avg_gdal.tif"))
med_RS_sprior807 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/med_RS_sprior807_mex_masked_lakes_avg_gdal.tif"))
volume_sprior10_sum <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/vol_sprior10_sum_mex_masked_lakes_avg_gdal.tif"))
eco_FDw <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_FDw_mex_masked_lakes_avg_gdal.tif"))
eco_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_median_ED_per_grid_mex_masked_lakes_avg_gdal.tif"))
eco_ED10 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/ED_sprior10_sum_mex_masked_lakes_avg_gdal.tif"))
eco_funrar_VOL <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_volume_mex_masked_lakes_avg_gdal.tif"))
volume_residual <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/med_volume_residuals_mex_masked_lakes_avg_gdal.tif"))
volume_DIV_rs_sprior_m2 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/med_volume_DIV_rs_sprior_m2_mex_masked_lakes_avg_gdal.tif"))
eco_funrar_RS <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_rangesize_mex_masked_lakes_avg_gdal.tif"))
eco_funrar_top10 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/ED_sprior10_sum_mex_masked_lakes_avg_gdal.tif"))
top_10_RS_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/top_10_RS_ED_combined_mex_masked_lakes_avg_gdal.tif"))
top_10_VOL_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/top_10_VOL_ED_combined_mex_masked_lakes_avg_gdal.tif"))
top_10_VOL_ED <- raster(paste0(path, "/traits_eco/top_10_RS_VOL_ED_BINARY_combined_mex_masked.tif")) # August 2022, USE THIS






myfact <- 12
agg_type <- "mean"


### Create directory
map_pdf_DIR <- paste0(path, "/stacked_outputs/pdf")
dir.create(map_pdf_DIR)
pdf_OUT <- paste0(map_pdf_DIR, "/factor_",myfact,"_", agg_type, "_poly")
dir.create(pdf_OUT)


### Aggregate by fact
sprior_sum807_agg <- aggregate(sprior_sum807, fact=myfact, agg_type, na.rm=T)
average_range_size_rarity_agg <- aggregate(average_range_size_rarity, fact=myfact, agg_type, na.rm=T)
average_range_size_rarity_agg <- aggregate(range_size_rarity, fact=myfact, agg_type, na.rm=T)

average_range_size_rarity_volume_agg <- aggregate(average_range_size_rarity_volume, fact=myfact, agg_type, na.rm=T)
range_size_rarity_807_volume_agg <- aggregate(range_size_rarity_volume, fact=myfact, agg_type, na.rm=T)
med_RS_sprior807_volume_agg <- aggregate(med_RS_sprior807_volume, fact=myfact, agg_type, na.rm=T)


med_RS_sprior807_agg <- aggregate(med_RS_sprior807, fact=myfact, agg_type, na.rm=T)
rangemaps_807_corrected_agg <- aggregate(rangemaps_807_corrected, fact=myfact, agg_type, na.rm=T)
rangesize_sprior10_sum_agg <- aggregate(rangesize_sprior10_sum, fact=myfact, agg_type, na.rm=T)
volume_sprior10_sum_agg <- aggregate(volume_sprior10_sum, fact=myfact, agg_type, na.rm=T)


eco_FDw_agg <- aggregate(eco_FDw, fact=myfact, agg_type, na.rm=T)
eco_ED_agg <- aggregate(eco_ED, fact=myfact, agg_type, na.rm=T)
eco_ED10_agg <- aggregate(eco_ED10, fact=myfact, agg_type, na.rm=T)

eco_funrar_VOL_agg <- aggregate(eco_funrar_VOL, fact=myfact, agg_type, na.rm=T)
writeRaster(eco_funrar_VOL_agg, paste0(pdf_OUT, "/eco_funrar_VOL_agg_12_mean.tif"))

eco_funrar_RS_agg <- aggregate(eco_funrar_RS, fact=myfact, agg_type, na.rm=T)

volume_residual_agg <- aggregate(volume_residual, fact=myfact, agg_type, na.rm=T)
volume_DIV_rs_sprior_m2_agg <- aggregate(volume_DIV_rs_sprior_m2, fact=myfact, agg_type, na.rm=T)

eco_funrar_top10_agg <- aggregate(eco_funrar_top10, fact=myfact, agg_type, na.rm=T)

top_10_RS_ED_agg <- aggregate(top_10_RS_ED, fact=myfact, agg_type, na.rm=T)
top_10_VOL_ED_agg <- aggregate(top_10_VOL_ED, fact=myfact, agg_type, na.rm=T)
writeRaster(top_10_VOL_ED_agg, "/tmp/top_10_VOL_ED_agg_12x_max.tif")

### Load the raster template
sprior_sum807 <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/sprior_sum807_mex_masked_lakes_avg_gdal.tif"))
sprior_sum807_agg <- aggregate(sprior_sum807, fact=myfact, agg_type, na.rm=T)



myres=ncell(sprior_sum807_agg)
e <- extent(sprior_sum807_agg)
e@xmin <- -145
e@xmax <- -50
e@ymin <- 25
e@xmax <- 60 # keep slighly bigger than final extent to avoid cut-off lines

### Load shape of North America and project
shape_na <- shapefile(paste0(path, "/layers_NA/NA_dissolved/na_dissolved.shp"))
newproj <-"+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"

shape_na_proj <-  spTransform(shape_na, newproj)
shape_na_fortify <- fortify(shape_na_proj)



### small insets
e_small <- e # copy
e_small@xmin <- -88.5
e_small@xmax <- -86.5
e_small@ymin <- 33.5
e_small@ymax <- 35

### ONLY FOR TOP 10% ED_VOL
# e_small <- e # copy   --> small_2
# e_small@xmin <- -120.5
# e_small@xmax <- -119.1
# e_small@ymin <- 38.6
# e_small@ymax <- 40.4

### ONLY FOR TOP 10% ED_VOL
# e_small <- e # copy  --> small_3
# e_small@xmin <- -98.19
# e_small@xmax <- -97.49
# e_small@ymin <- 29.45
# e_small@ymax <- 30.2


dir.create(paste0(path, "/stacked_outputs/inset_shp"))
e_small_poly <- as(e_small, 'SpatialPolygons')
crs(e_small_poly) <- crs(sprior_sum807)
e_small_poly_proj <- sp::spTransform(e_small_poly, newproj)
e_small_poly_proj_fortify <- fortify(e_small_poly_proj)


mynames<- c("sprior_sum807_agg",
            "average_range_size_rarity_agg",
            "med_RS_sprior807_volume_agg",
            "med_RS_sprior807_agg",
            "volume_sprior10_sum_agg",
            "eco_FDw_agg",
            "eco_ED_agg", 
            "eco_ED10_agg",
            "eco_funrar_VOL_agg", 
            "eco_funrar_RS_agg",
            "eco_funrar_top10_agg",
            "top_10_RS_ED_agg",
            "top_10_VOL_ED_agg")


i="top_10_VOL_ED_agg"





  for (i in mynames) {
  
  cat("\t Plotting", i, "by factor", myfact, "with aggregation type", agg_type, "\n")
      
    ### Rangesize FDw
    tmp <- get(i)

    
    tmp <- crop(tmp, e) # shape cropped above
    
    ##--- LOG --- # uncomment this to get log-maps
    # tmp <- calc(tmp, fun=function(x) log(x))
     
    ### Reproject to equal area
    tmp_proj <- projectRaster(tmp, crs=newproj, method="ngb", na.rm=T) # bilinear may produce negative outputs?
    # tmp_proj <- projectRaster(tmp, crs=newproj, method="bilinear", na.rm=T)
    myres <- ncell(tmp_proj)
    tmp_proj_trimmed <- trim(tmp_proj, padding=0, values=NA)
    
    ### Avoid cut -error due to small values
   if(names(tmp_proj)=="eco_funrar_rosauer_rangesize_mex_masked_lakes_avg_gdal") {tmp_proj <- tmp_proj*100000000} 
    
    min = raster::cellStats(tmp_proj, "min", na.rm=T)
    max = raster::cellStats(tmp_proj, "max", na.rm=T)
    diff <- max - min
    std = raster::cellStats(tmp_proj, sd, na.rm=T)
    
    myval <- values(tmp_proj)
    

    my_int = getJenksBreaks(myval, 11) # see here https://stackoverflow.com/questions/5304057/partition-into-classes-jenks-vs-kmeans
    
    ### Add jitter
    my_int = my_int + seq_along(my_int) * 0.0000001
      
    
    ### "Cut" the raster object based on jenks breaks
    tmp_proj_cut <- cut(tmp_proj,  breaks = my_int, include.lowest = T) # labels=T, dig.lab=15
    pol <- rasterToPolygons(tmp_proj_cut, dissolve=T)
    pol_fortify <- fortify(pol, region="layer")
    pol_fortify$id <- as.numeric(pol_fortify$id)
    pol_fortify$id <- as.factor(pol_fortify$id)
    
    
    ### Keep track of the integer codes on the raster:
    tmp_proj_df <- as.data.table(tmp_proj, xy=T, na.rm=T)
    names(tmp_proj_df) <- c("x", "y", "var")
    tmp_proj_df$valueDiscr <- cut(tmp_proj_df$var, breaks = my_int, include.lowest = T)
    tmp_proj_df$id <- cut(tmp_proj_df$var, breaks = my_int, include.lowest = T, labels = F)
    
    ### merge to main dataframe
    pol_fortify$id <- as.numeric(as.character(pol_fortify$id ))
    my_lookup <-  subset(tmp_proj_df, select=c(valueDiscr, id)) 
    my_lookup <- my_lookup[!duplicated(my_lookup),]
    pol_fortify2 <- merge(pol_fortify, my_lookup, by="id", all.x=T)
    n_breaks <- length(unique(pol_fortify2$valueDiscr))
    
  
    ### Plot
    p <- ggplot(data = pol_fortify2, # the input data
                aes(x = long, y = lat, fill = valueDiscr, group = group)) + # define variables #id
      geom_polypath() + # plot the DAs
      ylab("Latitude")+xlab("Longitude") + # axis lables
      scale_fill_manual(values = colorRamps::matlab.like(n_breaks)) +
      guides(fill = guide_legend(reverse=T)) + 
      theme(panel.grid.major = element_blank(), # remove grey plot backgound
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      theme(axis.title.x = element_text(face="bold", size=1),
            axis.text.x  = element_text(size=10)) +
      theme(axis.title.y = element_text(face="bold", size=1),
            axis.text.y  = element_text(size=10)) +
      theme(axis.line = element_line(colour="black", size=1))+
      theme(axis.ticks.length=unit(.15, "cm")) 
    
    p <- p + coord_equal()
    ### Add the border as polypath (else ends up in weird lines)
    p <- p + geom_polygon(data=shape_na_fortify, aes(x=long, y=lat, group=group), size=0.8, col="black",  fill=NA)
    p <- p + geom_polygon(data=e_small_poly_proj_fortify, aes(x=long, y=lat, group=group), size=1.5, col="black",  fill=NA)
        p <- p +  coord_cartesian(xlim= c( -2350000, 3100000), ylim= c( -2000000, 1950000), expand=T)
    
    ### export pdf
    pdf(paste0(pdf_OUT, "/", i, myfact, "x_", agg_type, "_BINARY.pdf"), width=15, height=10)
    plot(p)
    dev.off()
    
    
    ### export svg
    svg(paste0(pdf_OUT, "/", i, myfact, "x_", agg_type, "_BINARY.svg"), width=15, height=10) # NEED TO SPECIFY WIDTH & HEIGHT when using polygons!
    plot(p)
    dev.off()
    
    rm(p); gc()
    
  }





###--- Plot NORMAL MAPS smaller inset / window ----
e_small <- e # copy
e_small@xmin <- -88.5
e_small@xmax <- -86.5
e_small@ymin <- 33.5
e_small@ymax <- 35



### inset for richness and FD
e_small <- e # copy
e_small@xmin <- -113
e_small@xmax <- -111
e_small@ymin <- 33
e_small@ymax <- 34.5


pdf_OUT_SMALL_1 <- paste0(path, "/stacked_outputs/pdf/small_SR_FD_geom_poly")
dir.create(pdf_OUT_SMALL_1)


mynames<- c(
  "sprior_sum807",
  "eco_FDw",
  "average_range_size_rarity",
  "med_RS_sprior807_volume",
  "med_RS_sprior807",
  "volume_sprior10_sum",
  "eco_ED10",
  "eco_funrar_VOL",
  "eco_funrar_RS",
  "eco_funrar_top10",
  "top_10_RS_ED",
  "top_10_VOL_ED")



dir.create(paste0(path, "/stacked_outputs/inset_shp_SR_FD"))
e_small_poly <- as(e_small, 'SpatialPolygons')
crs(e_small_poly) <- crs(sprior_sum807)
e_small_poly_fortify <- fortify(e_small_poly)
shapefile(e_small_poly, paste0(path, "/stacked_outputs/inset_shp_SR_FD/small_SR_FD.shp"), overwrite=T)


### Get stream network as a background
str_net <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/hydrographic_network.tif")) # lakes=0, stream=1, NoData=255
str_net <- crop(str_net, e_small) # shape cropped above
values(str_net) <- ifelse(values(str_net)==1, 0, values(str_net)) # lakes should have same colour as streams


for (i in mynames) {
  
  cat("\t Plotting", i, "by factor_", myfact, "_SMALL 1 ", agg_type, "\n")

  tmp <- get(i)
  tmp <- crop(tmp, e_small) # shape cropped above
  
  ### add stream network in background
  values(tmp) <- ifelse( is.na(values(tmp)), values(str_net),values(tmp) )
  tmp <- trim(tmp, padding=0, values=NA)
  
  
  ### Reproject to equal area
  tmp_proj <- tmp # do not project --> missing pixels
  myres <- ncell(tmp_proj)
  tmp_proj_trimmed <- trim(tmp_proj, padding=0, values=NA)
  
  
  min = raster::cellStats(tmp_proj, "min", na.rm=T)
  max = raster::cellStats(tmp_proj, "max", na.rm=T)
  diff <- max - min
  std = raster::cellStats(tmp_proj, sd, na.rm=T)

  myval <- values(tmp_proj)
  my_int = getJenksBreaks(myval, 11)
  
    ### Add jitter
  my_int = my_int + seq_along(my_int) * 0.0000001

  
  tmp_proj_df <- as.data.table(tmp_proj, xy=T, na.rm=T)
  names(tmp_proj_df) <- c("x", "y", "var")
  
  tmp_proj_df$valueDiscr <- cut(tmp_proj_df$var, breaks = my_int, include.lowest = T)
  n_breaks <- length(unique(tmp_proj_df$valueDiscr))

  
  p <- ggplot(data = tmp_proj_df, aes(x=x,y=y)) +
    geom_tile(aes(fill = valueDiscr), color=NA) +
    scale_fill_manual(values = colorRamps::matlab.like(n_breaks)) + # , na.value="#000080"
    guides(fill = guide_legend(reverse=T)) +
    theme(panel.grid.major = element_blank(), # remove grey plot backgound
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    theme(axis.title.x = element_text(face="bold", size=1),
          axis.text.x  = element_text(size=10)) +
    theme(axis.title.y = element_text(face="bold", size=1),
          axis.text.y  = element_text(size=10)) +
    theme(axis.line = element_line(colour="black", size=1))+
    theme(axis.ticks.length=unit(.15, "cm"))
  p <- p + coord_equal() 
  
  
  ## Add black frame
  p <- p + geom_path(data=e_small_poly_fortify, aes(x=long, y=lat, group=group), size=0.5, col="black") 
  
  pdf(paste0(pdf_OUT_SMALL_1, "/", i, "_SMALL_WGS_1", ".pdf"), width=25, height=10)
  plot(p)
  dev.off()
  
  svg(paste0(pdf_OUT_SMALL_1, "/", i, "_SMALL_WGS_1", ".svg"), width=25, height=10)
  plot(p)
  dev.off()
  rm(p, p); gc()
  
}












###---- Plot bivariate map aggregated by factor x -----


### Load map
bivmap_rangesize_VOLUME_eco_ED <- raster(paste0(path_rast, "/bivmap_rangesize_VOLUME_eco_ED_mex_masked_lakes_avg.tif"))

### Aggreagte
bivmap_rangesize_VOLUME_eco_ED_agg <- aggregate(bivmap_rangesize_VOLUME_eco_ED, fact=myfact, agg_type, na.rm=T)


mynames <- c("bivmap_rangesize_VOLUME_eco_ED_agg")


for (i in mynames) {
  
  cat("\t Plotting", i, "by factor", myfact, "with aggregation type", agg_type, "\n")

  tmp <- get(i)
  tmp <- crop(tmp, e) # shape cropped above
    
  ### Reproject to equal area
  tmp_proj <- projectRaster(tmp, crs=newproj, method="ngb") 
  myres <- ncell(tmp_proj)

  pol <- rasterToPolygons(tmp_proj, dissolve=T)
  myname <- names(pol)
  pol_fortify <- fortify(pol, region=myname)
  pol_fortify$id <- as.numeric(pol_fortify$id)

      
   p <- ggplot(data = pol_fortify, # the input data
                  aes(x = long, y = lat, fill = id, group = group)) + # define variables #id
    geom_polypath(color = NA, lwd = 0) +
    theme(strip.background = element_blank(), # remove title background
                  strip.text.x = element_blank()) + # remove title
    scale_fill_gradientn( colours = as.vector(col.matrix), na.value="transparent")+
      theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent"))+ # get rid of legend panel bg
    theme(axis.title.x = element_text(face="bold", size=1),
          axis.text.x  = element_text(size=10)) +
    theme(axis.title.y = element_text(face="bold", size=1),
          axis.text.y  = element_text(size=10)) +
    theme(axis.line = element_line(colour="black", size=1))+
    theme(axis.ticks.length=unit(.15, "cm"))
  
  p <- p + coord_equal()
  ### Add the border as polypath (else ends up in weird lines)
  p <- p + geom_polygon(data=shape_na_fortify, aes(x=long, y=lat, group=group), size=0.8, col="black",  fill=NA)
  p <- p +  coord_cartesian(xlim= c( -2350000, 3100000), ylim= c( -2000000, 1950000), expand=T)
  
  
pdf(paste0(pdf_OUT, "/", i, myfact, "x_", agg_type, ".pdf"), width=15, height=10)
plot(p)
dev.off()


svg(paste0(pdf_OUT, "/", i, myfact, "x_", agg_type, ".svg"),  width=15, height=10)
plot(p)
dev.off()

png(paste0(pdf_OUT, "/", i, myfact, "x_", agg_type, ".png"), width=15, height=10, "in", res=1800) #width=1500, height=1000,
plot(p)
dev.off()


rm(p); gc()


}


###--- Plot BIVARIATE smaller inset / window ----
### inset for bivariate
e_small <- e # copy
e_small@xmin <- -108.5
e_small@xmax <- -107
e_small@ymin <- 33.5
e_small@ymax <- 35


pdf_OUT_SMALL_1 <- paste0(path, "/stacked_outputs/pdf/small_1_geom_tile")
dir.create(pdf_OUT_SMALL_1)


dir.create(paste0(path, "/stacked_outputs/inset_shp"))
e_small_poly <- as(e_small, 'SpatialPolygons')
crs(e_small_poly) <- crs(sprior_sum807)
e_small_poly_fortify <- fortify(e_small_poly)


shapefile(e_small_poly, paste0(path, "/stacked_outputs/inset_shp/inset_bivariate.shp"), overwrite=T)


mynames <- c("bivmap_rangesize_VOLUME_eco_ED")
load(paste0(path, "/col.matrix_bivariate.RData")) # col.matrix


for (i in mynames) {
  
  cat("\t Plotting", i, "by factor_", myfact, "_SMALL 1 ", agg_type, "\n")
  
  tmp <- get(i)
  
  tmp <- crop(tmp, e_small) # shape cropped above
  tmp <- trim(tmp, padding=0, values=NA)
  
  ### Reproject to equal area
  tmp_proj <- tmp # use this for the manuscript
  myres <- ncell(tmp_proj)
    
 tmp_proj_df <- as.data.table(tmp_proj, xy=T, na.rm=T)
 names(tmp_proj_df) <- c("x", "y", "var")
    
 
    p <- ggplot(data = tmp_proj_df, aes(x=x,y=y)) +
      geom_tile(aes(fill = var)) +
      coord_equal() +
      scale_fill_gradientn( colours = as.vector(col.matrix), na.value="transparent") +
      guides(fill = guide_legend(reverse=T)) +
      theme(panel.grid.major = element_blank(), # remove grey plot backgound
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      theme(axis.title.x = element_text(face="bold", size=1),
            axis.text.x  = element_text(size=10)) +
      theme(axis.title.y = element_text(face="bold", size=1),
            axis.text.y  = element_text(size=10)) +
      theme(axis.line = element_line(colour="black", size=1))+
      theme(axis.ticks.length=unit(.15, "cm"))+
      coord_equal() 
   
    ## Add black frame around the area
    pp <- p + geom_path(data=e_small_poly_fortify, aes(x=long, y=lat, group=group), size=0.5, col="black") 
    
    pdf(paste0(pdf_OUT_SMALL_1, "/", i, "_SMALL_WGS_2022_TEST_COLOR", ".pdf"), width=25, height=10)
    plot(pp)
    dev.off()
    
    
    svg(paste0(pdf_OUT_SMALL_1, "/", i, "_SMALL_WGS_2022_TEST_COLOR", ".svg"))
    plot(pp)
    dev.off()
    rm(p, pp); gc()
    
  }
  

  

###---- Export bivariate color matrix as svg ----
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}
  
  
x11()
col.matrix<-colmat(nquantiles=10, upperleft="blue", upperright="yellow", bottomleft="green", bottomright="red", xlab="My x label", ylab="My y label")
  

