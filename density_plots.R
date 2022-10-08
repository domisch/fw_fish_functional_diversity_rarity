

###===================================================#
###--- Scatterplots and aggregation of SR and FD ----
###===================================================#

if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("hexbin")) { install.packages("hexbin", dependencies = T) ; library(hexbin)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}

path <- "/mnt/domisch/data/fw_fish"
setwd(path)
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))
path_FIGURES <- paste0(path, "/figures")
dir.create(path_FIGURES)
path_OUT <- paste0(path, "/traits_eco")


source(paste0(path, "/scripts/raster.as.data.table_coord.R"))

  
### Load data
SR <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/sprior_sum807_mex_masked_lakes_avg_gdal.tif"))
FDpg <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDpg_mex_masked_lakes_avg_gdal.tif"))
FDw <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDw_mex_masked_lakes_avg_gdal.tif"))
FR_rs <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_rangesize_mex_masked_lakes_avg_gdal.tif"))
FR_vol <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_volume_mex_masked_lakes_avg_gdal.tif"))
RS_med <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/med_RS_sprior807_mex_masked_lakes_avg_gdal.tif"))
RS_med_VOL <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/med_RS_volume_sprior807_mex_masked_lakes_avg_gdal.tif"))
RS_rarity <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_mex_masked_lakes_avg_gdal.tif"))
RS_rarity_volume <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_volume_mex_masked_lakes_avg_gdal.tif"))
median_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_median_ED_per_grid_mex_masked_lakes_avg_gdal.tif"))





### Elevation, stream order and flow
dem <- raster(paste0(path, "/layers_NA/dem_site_mex_masked.tif"))
str_order <- raster(paste0(path, "/layers_NA/str_order_lakes0_mex_masked.tif"))
flow <-  raster(paste0(path, "/layers_NA/flo1k_mean_mex_masked.tif"))
flow <- crop(flow, str_order)


### Categorize flow
# myval <- unique(values(flow))
myval <- values(flow)


min = cellStats(flow, min, na.rm=T)
max = cellStats(flow, max, na.rm=T)
diff <- max - min
std = cellStats(flow, sd, na.rm=T)

myval <- values(flow)
# names(myval) <- "myval"
# myval <- unique(myval$myval)


my_int_equal = seq(min, max, by = diff/10) # equal
flow_cut_equal <- cut(flow, breaks = my_int_equal, include.lowest = T)
names(flow_cut_equal) <- "flow_cut_equal"
# x11(); plot(aggregate(flow_cut_equal, fact=12, na.rm=T))

my_int_quant = quantile(myval, probs=seq(0, 1, by = 1/10), na.rm=T) # quantile
flow_cut_quant <- cut(flow, breaks = my_int_quant, include.lowest = T)
names(flow_cut_quant) <- "flow_cut_quant"
# x11(); plot(aggregate(flow_cut_quant, fact=12, na.rm=T))
# plot(aggregate(flow, fact=12, na.rm=T))



### log-transform
myval_log <- log(myval+1)
flow_log <- log(flow+1)


min_log = cellStats(flow_log, "min", na.rm=T)
max_log = cellStats(flow_log, "max", na.rm=T)
diff_log <- max_log - min_log
std_log = cellStats(flow_log, "sd", na.rm=T)


my_int_equal_log = seq(min_log, max_log, by = diff_log/10) # equal
flow_cut_equal_log <- cut(flow_log, breaks = my_int_equal_log, include.lowest = T)
names(flow_cut_equal_log) <- "flow_cut_equal_log"





### Load raster/domain template
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes_mex_masked.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, grid_id_long) # set key


### Align spatial extents
SR <- crop(SR, domain_raster)
FDpg<- crop(FDpg, domain_raster)
FDw<- crop(FDw, domain_raster)
FR_rs <- crop(FR_rs, domain_raster)
FR_vol <- crop(FR_vol, domain_raster)
RS_rarity <- crop(RS_rarity, domain_raster)
RS_rarity_volume <- crop(RS_rarity_volume, domain_raster)
RS_med <- crop(RS_med, domain_raster)
RS_med_VOL <- crop(RS_med_VOL, domain_raster)
median_ED <- crop(median_ED, domain_raster)

### Convert rasters to data.table
mystack <- stack(domain_raster, SR, FDpg, FDw, dem, str_order, flow,  FR_rs, FR_vol, RS_med, RS_med_VOL, RS_rarity, RS_rarity_volume, median_ED) #flow_cut_equal, flow_cut_equal_log, flow_cut_quant,
mystack_dt <- as.data.table.raster(mystack,  xy=T)
mystack_dt$seq_id <- seq.int(1:nrow(mystack_dt))

save(mystack_dt, file=paste0(path, "/stacked_outputs/lakes_avg/eco_datatable_domain_raster_SR_FDpg_FDw_dem_str_order_latitude.RData"))
gc()


### Load data
load(paste0(path, "/stacked_outputs/lakes_avg/eco_datatable_domain_raster_SR_FDpg_FDw_dem_str_order_latitude.RData"))


### Specify bins and colours
cr <- colorRampPalette(c("white", 'blue','orange'))


### Add density legend
# https://stackoverflow.com/questions/14271584/r-legend-for-color-density-scatterplot-produced-using-smoothscatter#comment11130419_8899096
add_density_legend <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}




###----Scatterplots eco -----
### Cut latitude below 30, avoid data plotting across Mexico
mystack_dt$y_30 <- ifelse(mystack_dt$y > 30, mystack_dt$y, NA)



### Get broad scale averages (the white boxes in the figure)
nbins <- 10
### Check the data ranges
range(mystack_dt$dem_site, na.rm=T) 
range(mystack_dt$y_30, na.rm=T) 

range(mystack_dt$sprior_sum807, na.rm=T) #  0.1428571 98.0000000
mystack_dt$sprior_sum807_cut <- cut(mystack_dt$sprior_sum807, nbins,  include.lowest = T, labels=rep(seq(10,100,10)))
mystack_dt$eco_FDpg_cut <- cut(mystack_dt$eco_FDpg, nbins,  include.lowest = T)
mystack_dt$eco_FDw_cut <- cut(mystack_dt$eco_FDw, nbins,  include.lowest = T)


nbins = 10 # ~350m
mystack_dt$dem_site_cut <- cut(mystack_dt$dem_site, nbins, include.lowest = T, labels=rep(seq(350,3500,350)))
mystack_dt$y_30_cut <- cut(mystack_dt$y_30, nbins, include.lowest = T, labels=rep(seq(33,60,3)))


mystack_dt$dem_site_cut_mean <- aggregate(mystack_dt["dem_site"], #the data frame
                          by=list(cut(mystack_dt$dem_site, nbins, include.lowest = T, labels=rep(seq(350,3500,350))) ), #the bins (see below)
                          mean) #the aggregating function



nbins = 10 # ~350m
### Create bins
mystack_dt$dem_site_cut <- cut(mystack_dt$dem_site, nbins, include.lowest = T, labels=rep(seq(350,3500,350)))
mystack_dt$y_30_cut <- cut(mystack_dt$y_30, nbins, include.lowest = T, labels=rep(seq(33,60,3)))
### stream order: needs factors
mystack_dt$str_order_lakes0_fact <- as.factor(mystack_dt$str_order_lakes0)
### flow: log first
mystack_dt$flo1k_mean_log <- log(mystack_dt$flo1k_mean+1)
range(mystack_dt$flo1k_mean_log, na.rm=T)
mystack_dt$flo1k_mean_log_cut <- cut(mystack_dt$flo1k_mean_log, nbins, include.lowest = T, labels=rep(seq(1,10,1)))







### ELEVATION BINS

### SD - aggregate across bins
tmp_sprior_sum807 <- mystack_dt[, as.double(mean(sprior_sum807, na.rm = T)), by = dem_site_cut]
names(tmp_sprior_sum807) <- c("elevation", "species_richness")
tmp_sprior_sum807$elevation <- as.numeric(as.character(tmp_sprior_sum807$elevation))

p <- ggplot(tmp_sprior_sum807, aes(x=elevation, y=species_richness))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(x = 0, y = 0)
p <- p + scale_x_continuous(breaks = tmp_sprior_sum807$elevation)

svg(paste0(path_FIGURES, "/elevation_SR_bin_10.svg"))
# png(paste0(path_FIGURES, "/elevation_SR_bin_10.png"))
plot(p)
dev.off()



### FD - aggregate across bins
tmp_eco_FDw <- mystack_dt[, as.double(mean(eco_FDw, na.rm = T)), by = dem_site_cut]
names(tmp_eco_FDw) <- c("elevation", "eco_FDw")
tmp_eco_FDw$elevation <- as.numeric(as.character(tmp_eco_FDw$elevation))

p <- ggplot(tmp_eco_FDw, aes(x=elevation, y=eco_FDw))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(x = 0, y = 0)
p <- p + scale_x_continuous(breaks = tmp_eco_FDw$elevation)

svg(paste0(path_FIGURES, "/elevation_FDw_bin_10.svg"))
# png(paste0(path_FIGURES, "/elevation_FDw_bin_10.png"))
plot(p)
dev.off()


### FD ~ SD bin plot
tmp_eco_FDw <- mystack_dt[, as.double(mean(eco_FDw, na.rm = T)), by = sprior_sum807_cut]
names(tmp_eco_FDw) <- c("species_richness", "eco_FDw")
tmp_eco_FDw$species_richness <- as.numeric(as.character(tmp_eco_FDw$species_richness))

p <- ggplot(tmp_eco_FDw, aes(x=species_richness, y=eco_FDw))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(x = 0, y = 0)
p <- p + scale_x_continuous(breaks = tmp_eco_FDw$species_richness)

svg(paste0(path_FIGURES, "/species_richness_FDw_bin_10.svg"))
# png(paste0(path_FIGURES, "/species_richness_FDw_bin_10.png"))
plot(p)
dev.off()





### LATITUDE BINS

### Richness - aggregate across bins
tmp_sprior_sum807 <- mystack_dt[, as.double(mean(sprior_sum807, na.rm = T)), by = y_30_cut]
names(tmp_sprior_sum807) <- c("latitude", "species_richness")
tmp_sprior_sum807$latitude <- as.numeric(as.character(tmp_sprior_sum807$latitude))

p <- ggplot(tmp_sprior_sum807, aes(x=latitude, y=species_richness))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_sprior_sum807$latitude)

svg(paste0(path_FIGURES, "/latitude_SR_bin_10.svg"))
# png(paste0(path_FIGURES, "/latitude_SR_bin_10.png"))
plot(p)
dev.off()



### FD - aggregate across bins
tmp_eco_FDw <- mystack_dt[, as.double(mean(eco_FDw, na.rm = T)), by = y_30_cut]
names(tmp_eco_FDw) <- c("latitude", "eco_FDw")
tmp_eco_FDw$latitude <- as.numeric(as.character(tmp_eco_FDw$latitude))

p <- ggplot(tmp_eco_FDw, aes(x=latitude, y=eco_FDw))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_eco_FDw$latitude)

svg(paste0(path_FIGURES, "/latitude_FDw_bin_10.svg"))
# png(paste0(path_FIGURES, "/latitude_FDw_bin_10.png"))
plot(p)
dev.off()




### STREAM ORDER BINS

### Richness - aggregate across bins
tmp_sprior_sum807 <- mystack_dt[, as.double(mean(sprior_sum807, na.rm = T)), by = str_order_lakes0]
names(tmp_sprior_sum807) <- c("str_order_lakes0", "species_richness")
tmp_sprior_sum807$str_order_lakes0 <- as.numeric(as.character(tmp_sprior_sum807$str_order_lakes0))

p <- ggplot(tmp_sprior_sum807, aes(x=str_order_lakes0, y=species_richness))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_sprior_sum807$str_order_lakes0)

svg(paste0(path_FIGURES, "/str_order_lakes0_SR_bin_10.svg"))
# png(paste0(path_FIGURES, "/str_order_lakes0_SR_bin_10.png"))
plot(p)
dev.off()



### FD - aggregate across bins
tmp_eco_FDw <- mystack_dt[, as.double(mean(eco_FDw, na.rm = T)), by = str_order_lakes0]
names(tmp_eco_FDw) <- c("str_order_lakes0", "eco_FDw")
tmp_eco_FDw$str_order_lakes0 <- as.numeric(as.character(tmp_eco_FDw$str_order_lakes0))

p <- ggplot(tmp_eco_FDw, aes(x=str_order_lakes0, y=eco_FDw))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_eco_FDw$str_order_lakes0)

svg(paste0(path_FIGURES, "/str_order_lakes0_FDw_bin_10.svg"))
# png(paste0(path_FIGURES, "/str_order_lakes0_FDw_bin_10.png"))
plot(p)
dev.off()




### FLOW (log) BINS
### Richness - aggregate across bins
tmp_sprior_sum807 <- mystack_dt[, as.double(mean(sprior_sum807, na.rm = T)), by = flo1k_mean_log_cut]
names(tmp_sprior_sum807) <- c("log_flow", "species_richness")
tmp_sprior_sum807$log_flow <- as.numeric(as.character(tmp_sprior_sum807$log_flow))

p <- ggplot(tmp_sprior_sum807, aes(x=log_flow, y=species_richness))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_sprior_sum807$log_flow)

svg(paste0(path_FIGURES, "/log_flow_SR_bin_10.svg"))
# png(paste0(path_FIGURES, "/log_flow_SR_bin_10.png"))
plot(p)
dev.off()



### FDw - aggregate across bins
tmp_eco_FDw <- mystack_dt[, as.double(mean(eco_FDw, na.rm = T)), by = flo1k_mean_log_cut]
names(tmp_eco_FDw) <- c("log_flow", "eco_FDw")
tmp_eco_FDw$log_flow <- as.numeric(as.character(tmp_eco_FDw$log_flow))

p <- ggplot(tmp_eco_FDw, aes(x=log_flow, y=eco_FDw))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=15)
p <- p +  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
p <- p + 
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=15))

p <- p + expand_limits(y = 0)
p <- p + scale_x_continuous(breaks = tmp_eco_FDw$log_flow)

svg(paste0(path_FIGURES, "/log_flow_FDw_bin_10.svg"))
# png(paste0(path_FIGURES, "/log_flow_FDw_bin_10.png"))
plot(p)
dev.off()









###----- Run various density plots to explore patterns ----

### elevation
mystack_dt$dem_site_cut <- as.numeric(as.character(mystack_dt$dem_site_cut))


# svg(paste0(path_FIGURES, "/elevation_SR_bin_10.svg"))
png(paste0(path_FIGURES, "/elevation_SR_bin_10_TEST.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum814 ~ mystack_dt$dem_site_cut_mean, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Species richness", xlab="Elevation [m]", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend) #xlim=c(0, 3500)
dev.off()



# svg(paste0(path_FIGURES, "/elevation_FDw_bin_10.svg"))
png(paste0(path_FIGURES, "/elevation_FDw_bin_10.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$dem_site_cut, nrpoints = 0, nbin=250, colramp=cr,
              ylab="FDw", xlab="Elevation [m]", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend) #xlim=c(0, 3500)
dev.off()



### latitude
mystack_dt$y_30_cut <- as.numeric(as.character(mystack_dt$y_30_cut))

# richness
# svg(paste0(path_FIGURES, "/latitude_SR_bin_10.svg"))
png(paste0(path_FIGURES, "/latitude_SR_bin_10.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum814 ~ mystack_dt$y_30_cut, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Species richness", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend) #xlim=c(0, 3500)
dev.off()


# FDw
# svg(paste0(path_FIGURES, "/latitude_FDw_bin_10.svg"))
png(paste0(path_FIGURES, "/latitude_FDw_bin_10.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$y_30_cut, nrpoints = 0, nbin=250, colramp=cr,
              ylab="FDw", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend) #xlim=c(0, 3500)
dev.off()



### Run different combinations

# FDpg - SR
# pdf(paste0(path_FIGURES, "/col_legend.pdf"))
svg(paste0(path_FIGURES, "/eco_FDpg_SR.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_FDpg_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDpg ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Functional diversity (PG)", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()

# FDw - SR
svg(paste0(path_FIGURES, "/eco_FDw_SR.svg")) #, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_FDw_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250,  colramp=cr, 
              ylab="Functional diversity (wt)", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

# FDpg - FDw
svg(paste0(path_FIGURES, "/eco_FDpg_FDw.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_FDpg_FDw.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDpg ~ mystack_dt$eco_FDw, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional diversity (PG)", xlab="Functional diversity (wt)", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

# FDw - latitude
svg(paste0(path_FIGURES, "/eco_FDw_latitude.svg")) #, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_FDw_latitude.png")) 
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$y_30, nrpoints = 0, nbin=250,  colramp=cr, 
              ylab="Functional diversity (wt)", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



# richness - latitude
svg(paste0(path_FIGURES, "/SR_lat.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/SR_lat.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Species richness", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# LOG richness - latitude
svg(paste0(path_FIGURES, "/SRlog_lat.svg"))
# png(paste0(path_FIGURES, "/SRlog_lat.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(log(mystack_dt$sprior_sum807+1) ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="log species richness", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# range size - richness
svg(paste0(path_FIGURES, "/RS_SR.svg"))
# png(paste0(path_FIGURES, "/RS_SR.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter( mystack_dt$sprior_sum807 ~ mystack_dt$med_RS_sprior807 , nrpoints = 0, nbin=250, colramp=cr, 
               ylab="Species richness", xlab="Range size", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# log(range size) - richness
svg(paste0(path_FIGURES, "/logRS_SR.svg"))
# png(paste0(path_FIGURES, "/logRS_SR.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ log(mystack_dt$med_RS_sprior807+1), nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Species richness", xlab="LOG range size",  cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# range size (VOLUME)- richness
svg(paste0(path_FIGURES, "/RSvol_SR.svg"))
# png(paste0(path_FIGURES, "/RSvol_SR.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter( mystack_dt$sprior_sum807 ~ mystack_dt$med_RS_volume_sprior807 , nrpoints = 0, nbin=250, colramp=cr, 
               ylab="Species richness", xlab="Range size (volume)", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# log(range size) (VOLUME) - richness
svg(paste0(path_FIGURES, "/logRSvol_SR.svg"))
# png(paste0(path_FIGURES, "/logRSvol_SR.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ log(mystack_dt$med_RS_volume_sprior807+1), nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Species richness", xlab="LOG range size (volume)",  cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



# range size - range size volume
svg(paste0(path_FIGURES, "/RS_RSvol.svg"))
# png(paste0(path_FIGURES, "/RS_RSvol.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter( mystack_dt$med_RS_volume_sprior807 ~ mystack_dt$med_RS_sprior807 , nrpoints = 0, nbin=250, colramp=cr, 
               ylab="range size volume", xlab="Range size", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()





# range size rarity - richness
svg(paste0(path_FIGURES, "/eco_RSrarity_SR.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_RSrarity_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$average_range_size_rarity ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Range size rarity", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



# range size rarity (volume)- richness
svg(paste0(path_FIGURES, "/eco_RSrarityVOL_SR.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_RSrarityVOL_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$average_range_size_rarity_volume ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Range size rarity (VOL)", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()




# range size rarity - latitude
svg(paste0(path_FIGURES, "/eco_RSrarity_latitude.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_RSrarity_latitude.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$average_range_size_rarity ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Range size rarity", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

# range size rarity (volume) - latitude
svg(paste0(path_FIGURES, "/eco_RSrarityVOL_latitude.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_RSrarityVOL_latitude.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$average_range_size_rarity_volume ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Range size rarity", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# functional rarity (rangesize) - richness
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_SR.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_rangesize ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) rangesize", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# functional rarity (volume) - richness
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_vol_SR.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_vol_SR.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_volume ~ mystack_dt$sprior_sum807, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) volume", xlab="Species richness", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



# functional rarity (rangesize) - latitude
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_latitude.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_latitude.png"))#
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_rangesize ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) rangesize", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

# functional rarity (volume) - latitude
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_vol_latitude.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_vol_latitude.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_volume ~ mystack_dt$y_30, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) volume", xlab="Latitude", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



# functional rarity - range size rarity
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_RSrarity.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_RSrarity.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_rangesize ~ mystack_dt$average_range_size_rarity, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer)", xlab="Range size rarity", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()


# functional rarity (rangesize)- average ED
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_ED.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_ED.png"))
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_rangesize ~ mystack_dt$eco_median_ED_per_grid, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) rangesize", xlab="median ED", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

# functional rarity (volume)- average ED
svg(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_ED.svg"))#, width=20, height=20)
# png(paste0(path_FIGURES, "/eco_funrar_rosauer_rs_ED.png"))#
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_funrar_rosauer_volume ~ mystack_dt$eco_median_ED_per_grid, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional rarity (Rosauer) volume", xlab="median ED", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()




### median flow plots
# richness - flow
svg(paste0(path_FIGURES, "/richness_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/richness_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="Species richness", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()

# log richness - log flow
svg(paste0(path_FIGURES, "/log_richness_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/log_richness_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(log(mystack_dt$sprior_sum807+1) ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="log Species richness", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()



# FDw - flow
svg(paste0(path_FIGURES, "/FDw_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/FDw_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="FDw", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()




# dem_site plots

### dem site ~ richness
svg(paste0(path_FIGURES, "/SR_dem.svg"))
# png(paste0(path_FIGURES, "/SR_dem.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ mystack_dt$dem_site, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Species richness", xlab="Elevation [m]", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()

### dem site ~ FDpg
svg(paste0(path_FIGURES, "/eco_FDpg_dem.svg"))
# png(paste0(path_FIGURES, "/eco_FDpg_dem.png")) #, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDpg ~ mystack_dt$dem_site, nrpoints = 0, nbin=250,  colramp=cr, 
              ylab="Functional diversity (PG)", xlab="Elevation [m]", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()

### dem site ~ FDw
svg(paste0(path_FIGURES, "/eco_FDw_dem.svg"))
# png(paste0(path_FIGURES, "/eco_FDw_dem.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$dem_site, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional diversity (wt)", xlab="Elevation [m]", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
dev.off()



### stream_order plots

### stream_order ~ richness
svg(paste0(path_FIGURES, "/SR_StrOrd.svg"))
# png(paste0(path_FIGURES, "/SR_StrOrd.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$sprior_sum807 ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Species richness", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()


### stream_order ~ range size
svg(paste0(path_FIGURES, "/RS_StrOrd.svg"))
# png(paste0(path_FIGURES, "/RS_StrOrd.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$med_RS_sprior807 ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Range size", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()


### stream_order ~ range size (volume)
svg(paste0(path_FIGURES, "/RSvol_StrOrd.svg"))
# png(paste0(path_FIGURES, "/RSvol_StrOrd.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$med_RS_volume_sprior807 ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Range size (volume)", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()


### stream_order ~ elevation
svg(paste0(path_FIGURES, "/dem_StrOrd.svg"))
# png(paste0(path_FIGURES, "/dem_StrOrd.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$dem_site ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250, colramp=cr,
              ylab="Elevation (m)", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()





### stream_order ~ FDpg
svg(paste0(path_FIGURES, "/eco_FDpg_StrOrd.svg"))
# png(paste0(path_FIGURES, "/eco_FDpg_StrOrd.png")) #, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDpg ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250,  colramp=cr, 
              ylab="Functional diversity (PG)", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()

### stream_order ~ FDw
svg(paste0(path_FIGURES, "/eco_FDw_StrOrd.svg"))
# png(paste0(path_FIGURES, "/eco_FDw_StrOrd.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$eco_FDw ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250, colramp=cr, 
              ylab="Functional diversity (wt)", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()




### stream_prder ~ log(flow)
svg(paste0(path_FIGURES, "/logFlow_StrOrd.svg"))
# png(paste0(path_FIGURES, "/logFlow_StrOrd.png")) #, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(mystack_dt$"flo1k_mean" ~ mystack_dt$str_order_lakes0, nrpoints = 0, nbin=250,  colramp=cr, 
              ylab="log Flow [m3/sec]", xlab="Stream order", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), xaxt="n", postPlotHook = add_density_legend, useRaster=T)
axis(1, at=seq(0, 8, by=1), labels = T, cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0))
dev.off()



# log richness - log flow
svg(paste0(path_FIGURES, "/log_richness_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/log_richness_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(log(mystack_dt$sprior_sum807+1) ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="log Species richness", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()


# log range size - log flow
svg(paste0(path_FIGURES, "/log_RS_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/log_RS_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(log(mystack_dt$med_RS_sprior807+1) ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="log range size", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()


# log range size (volume)- log flow
svg(paste0(path_FIGURES, "/log_RSvol_flow_LOG.svg"))
# png(paste0(path_FIGURES, "/log_RvolS_flow_LOG.png"))#, width=20, height=20)
par(mfrow=c(1, 1), bty='n', mar=c(5.1, 4.1, 4.1, 2.1)+3.5)
smoothScatter(log(mystack_dt$med_RS_volume_sprior807+1) ~ log(mystack_dt$flo1k_mean+1), nrpoints = 0, nbin=250, colramp=cr,
              ylab="log range size (volume)", xlab="log flow m3/sec", cex=4,  bty='l', cex.axis=2, cex.lab=2, mgp = c(3, 1.2, 0), postPlotHook = add_density_legend, useRaster=T)
# axis(side=1, pos=0, lwd.ticks=0)
dev.off()


