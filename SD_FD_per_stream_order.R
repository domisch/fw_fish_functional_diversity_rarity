

###===================================================#
###----- Total unique SR and FD per stream order -----
###===================================================#



path <- "/mnt/domisch/data/fw_fish"
setwd(path)
getwd()

path_FIGURES <- paste0(path, "/figures")
dir.create(path_FIGURES)

"%ni%" <- Negate("%in%")


if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("hexbin")) { install.packages("hexbin", dependencies = T) ; library(hexbin)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}
rasterOptions(tmpdir=paste0(path, "/R_temp_delete"))


source(paste0(path, "/scripts/raster.as.data.table_coord.R"))


### Load data
SR <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/sprior_sum807_mex_masked_lakes_avg_gdal.tif"))
FDpg <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDpg_mex_masked_lakes_avg_gdal.tif"))
FDw <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_FDw_mex_masked_lakes_avg_gdal.tif"))
FR <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_rangesize_mex_masked_lakes_avg_gdal.tif"))
FR_vol <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/eco_funrar_rosauer_volume_mex_masked_lakes_avg_gdal.tif"))
RS_med <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/med_RS_sprior807_mex_masked_lakes_avg_gdal.tif"))
RS_med_vol <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/med_RS_volume_sprior807_mex_masked_lakes_avg_gdal.tif"))
RS_rarity <- raster(paste0(path,  "/stacked_outputs/lakes_avg_gdal/average_range_size_rarity_mex_masked_lakes_avg_gdal.tif"))
median_ED <- raster(paste0(path, "/stacked_outputs/lakes_avg_gdal/eco_median_ED_per_grid_mex_masked_lakes_avg_gdal.tif"))



### These need to cover the entire range!
dem <- raster(paste0(path, "/layers_NA/dem_site.tif"))
str_order <- raster(paste0(path, "/layers_NA/str_order_lakes0.tif"))
flow <- raster(paste0(path, "/FLO1k/mean_flow_1km_1960_2015.tif"))
flow <- crop(flow, str_order)

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
str_order <- extend(str_order, e, value=NA)
flow <- extend(flow, e, value=NA)


### Convert rasters to data.table
mystack <- stack(domain_raster, SR, FDpg, FDw, str_order, flow)
mystack_dt <- as.data.table.raster(mystack)
mystack_dt$seq_id <- seq.int(1:nrow(mystack_dt))



### Load large species-per-grid table
tmp_rs <- fread("stacked_outputs/sprior_807_probs_per_grid.txt")
setnames(tmp_rs, c("probs", "range_size", "grid_id_long", "species"))
setkey(tmp_rs, grid_id_long)

tmp <- merge(tmp_rs, mystack_dt, by="grid_id_long", all.y=T) # masked mex. for SR etc, not env. layers (=dem, flow etc)
rm(tmp_rs, mystack_dt); gc()


water_volume_by_spp_per_streamorder <- tmp[, lapply(.SD, sum,  na.rm=T), by= c("species", "str_order_lakes0"), .SDcols=c("mean_flow_1km_1960_2015") ]
water_volume_by_spp_per_streamorder <- water_volume_by_spp_per_streamorder[-c(1:10),]

### Make wide format
water_volume_by_spp_per_streamorder_wide <- reshape(water_volume_by_spp_per_streamorder, idvar = "species", timevar = "str_order_lakes0", direction = "wide")
write.csv(water_volume_by_spp_per_streamorder_wide, paste0(path, "/water_volume_by_spp_per_streamorder_wide_Nov21.csv"), row.names=F, quote=F)


write.csv(water_volume_by_spp_per_streamorder, paste0(path, "/water_volume_by_spp_per_streamorder_Nov21.csv"), row.names=F, quote=F)
save(tmp, file=paste0(path, "/water_volume_SR_FD_strOrd_flow_per_grid_big_table_Nov21.RData"))


load(paste0(path, "/water_volume_SR_FD_strOrd_flow_per_grid_big_table_Nov21.RData"))


### Get the species per stream order
spp_0 <- subset(tmp, str_order_lakes0==0)
spp_0 <- unique(spp_0$species)

spp_1 <- subset(tmp, str_order_lakes0==1)
spp_1 <- unique(spp_1$species)

spp_2 <- subset(tmp, str_order_lakes0==2)
spp_2 <- unique(spp_2$species)

spp_3 <- subset(tmp, str_order_lakes0==3)
spp_3 <- unique(spp_3$species)

spp_4 <- subset(tmp, str_order_lakes0==4)
spp_4 <- unique(spp_4$species)

spp_5 <- subset(tmp, str_order_lakes0==5)
spp_5 <- unique(spp_5$species)

spp_6 <- subset(tmp, str_order_lakes0==6)
spp_6 <- unique(spp_6$species)

spp_7 <- subset(tmp, str_order_lakes0==7)
spp_7 <- unique(spp_7$species)

spp_8 <- subset(tmp, str_order_lakes0==8)
spp_8 <- unique(spp_8$species)


### Plot against stream order
unique_spp_per_ord <- data.frame(str_ord=seq(1:8), 
                                 n_unique_spp = c(length(spp_1), 
                                                  length(spp_2), 
                                                  length(spp_3),
                                                  length(spp_4),
                                                  length(spp_5),
                                                  length(spp_6),
                                                  length(spp_7),
                                                  length(spp_8))) 



library(ggplot2)
p<-ggplot(data=unique_spp_per_ord, aes(x=str_ord, y=n_unique_spp)) +
  geom_bar(stat="identity")
p


### Get total area nd volume per stream order
total_area <- read.csv(paste0(path, "/estimated_volume_in_stream_orders.csv"))
total_area <- subset(total_area, str_ord !=0) # omit lakes
total_area <- merge(total_area, unique_spp_per_ord, by="str_ord")
total_area$str_ord <- as.factor(total_area$str_ord)


### Volume in million m3
total_area$volume_per_strord / 1000000


### Which species in #2 not not occur in #1, and which in #3 do not occur in #1+2 etc...

### --> these find oly the non-matching from the FIRST argument:
# "%ni%" <- Negate("%in%")
funct_not_in_first <- function(x, y) {x[x %ni% y] }

funct_not_in_all <- function(list.a, list.b ) {
  c(setdiff(list.b, list.a), setdiff(list.a, list.b))
} 


ord_1 <- spp_1
ord_1_to_0 <- funct_not_in_all(spp_1, spp_0) # 49 occur in lakes, not in ord_1

funct_not_in_first(spp_1, spp_0)


spp_12 <- unique(Reduce(c, list(spp_1, spp_2)))
spp_123 <- unique(Reduce(c, list(spp_1, spp_2, spp_3 )))
spp_1234 <- unique(Reduce(c, list(spp_1, spp_2, spp_3, spp_4 )))
spp_12345 <- unique(Reduce(c, list(spp_1, spp_2, spp_3, spp_4, spp_5 )))
spp_123456 <-unique(Reduce(c, list(spp_1, spp_2, spp_3, spp_4, spp_5, spp_6  )))
spp_1234567 <- unique(Reduce(c, list(spp_1, spp_2, spp_3, spp_4, spp_5, spp_6, spp_7  )))

ord_1_to_2 <- funct_not_in_first(spp_2, spp_1)
ord_12_to_3 <- funct_not_in_first(spp_3, spp_12 )
ord_123_to_4 <- funct_not_in_first(spp_4, spp_123)
ord_1234_to_5 <- funct_not_in_first(spp_5, spp_1234)
ord_12345_to_6 <- funct_not_in_first(spp_6, spp_12345)
ord_123456_to_7 <-funct_not_in_first(spp_7, spp_123456)
ord_1234567_to_8 <- funct_not_in_first(spp_8, spp_1234567)



### Plot against stream order
unique_spp_per_ord <- data.frame(str_ord=seq(1:8), 
                                n_unique_spp = c(length(ord_1), 
                                                length(ord_1_to_2), 
                                                length(ord_12_to_3),
                                                length(ord_123_to_4),
                                                length(ord_1234_to_5),
                                                length(ord_12345_to_6),
                                                length(ord_123456_to_7),
                                                length(ord_1234567_to_8))) 












### Calculate FD of those unique species per stream order


### Load libraries
if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("maptools")) { install.packages("maptools", dependencies = T) ; library(maptools)}
if (!require("data.table")) { install.packages("data.table", dependencies = T) ; library(data.table)}
if (!require("plyr")) { install.packages("plyr", dependencies = T) ; library(plyr)}


if (!require("devtools")) { install.packages("devtools", dependencies = T) ; library(devtools)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = T) ; library(doParallel)}
if (!require("foreach")) { install.packages("foreach", dependencies = T) ; library(foreach)}
if (!require("FD")) { install.packages("FD", dependencies = T) ; library(FD)}
if (!require("cluster")) { install.packages("cluster", dependencies = T) ; library(cluster)}
if (!require("vegan")) { install.packages("vegan", dependencies = T) ; library(vegan)}



# install_github("ibartomeus/fundiv")
require(fundiv)
# zip -r  IN.zip IN
path <- "/mnt/domisch/data/fw_fish"
source(paste0(path, "/scripts/Xtree.R"))


### Load trait data
load(paste0(path, "/traits_eco/fish_traits_and_tree_with_weights.RData")) 
trait_weights <- c(rep(1, 7)) # eco


### Get species
results <- list()
myloop <- c("spp_0", "spp_1", "spp_2", "spp_3", "spp_4", "spp_5", "spp_6", "spp_7", "spp_8")

for (i in 1:length(myloop) ) {
  cat("Running species set", (myloop[i]), "\n")
  
  tmp_site <- get(myloop[i])


tmp_site <- data.frame(species=tmp_site, probs=1)
tmp_site <- na.omit(tmp_site)

library(tidyr)
tmp_site <- spread(data = tmp_site, 
             key = species,
             value = probs)


### Columns (=species) need to be in the same order as the rows in "trait_final"
species_order <- row.names(traits_final)

### Not all species modelled in each grid - add empty columns
### Which species do not occur in the subset?
"%ni%" <- Negate("%in%")
missing <- species_order[species_order %ni% names(tmp_site) ]

if(length(missing) >0) {
### Create dummy data (missing species)
tmp_dt <- data.table(matrix(ncol = length(missing), nrow = nrow(tmp_site)))
names(tmp_dt) <- missing

### Add the dummy data
tmp_site <- cbind(tmp_site, tmp_dt)
}
### Set column order (species names) as in traits
setcolorder(tmp_site, as.character(species_order)) # data.table

### Replace all NA with zero (probabilities)
tmp_site[is.na(tmp_site)] <- 0

### Fuction needs rownames
tmp_site <- as.data.frame(tmp_site)
row.names(tmp_site) <- tmp_site$V1



FD_per_stream_order <-  FD_dendro(S = traits_final[-1], A = tmp_site, w=trait_weights,
                                 Distance.method= "gower",  Cluster.method = "average", ord = "podani",
                                 Weigthedby = "abundance")

FD_per_stream_order$str_ord <- i-1 # stream order starts from 0 (lakes)
  
results[[i]] <- FD_per_stream_order
rm(tmp_site, FD_per_stream_order)
graphics.off()

}


results <- do.call(rbind, results)
results
write.csv(results, paste0(path, "/traits_eco/total_FD_per_stream_order.csv"), quote = F, row.names = F)


### Plot against volume and area (number pixels) per order

### Get total area nd volume per stream order
total_area <- read.csv(paste0(path, "/estimated_volume_in_stream_orders.csv"))
total_area <- merge(total_area, results, by="str_ord")
total_area$str_ord <- as.factor(total_area$str_ord)
write.csv(total_area, paste0(path, "/traits_eco/total_FD_per_stream_order.csv"), quote = F, row.names = F)



### Total area per stream order vs. FD
x11()
p <- ggplot(total_area, aes((log(number_pixels)), FDw )) + geom_point(size=8) # col="#a29b93"
p <- p +  theme(panel.grid.major = element_blank(), # remove grey plot backgound
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                text = element_text(size=20))
p  <- p + theme(axis.text.x = element_text(face="bold", color="black", size=20),
                axis.text.y = element_text(face="bold", color="black", size=20))
p <- p + xlab("log number 1km pixels")+ylab("total FD")





x11()
p <- ggplot(total_area, aes((log(volume_per_strord)), (FDw) )) + geom_point(size=8) # col="#a29b93"
p <- p +  theme(panel.grid.major = element_blank(), # remove grey plot backgound
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                text = element_text(size=20))
p  <- p + theme(axis.text.x = element_text(face="bold", color="black", size=20),
                axis.text.y = element_text(face="bold", color="black", size=20))
p <- p + xlab("log volume")+ylab("total FD")
p



total_area$volume_per_strord_log <- round(log(total_area$volume_per_strord+1),1)


p <- ggplot(total_area, aes(x=reorder(volume_per_strord_log, str_ord), y=FDw ))
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
p <- p + xlab("log volume")+ylab("total FD")

svg(paste0(path, "/total_FD_volume.svg"))
# png(paste0(path, "/total_FD_volume.png"))
plot(p)
dev.off()



### area
p <- ggplot(total_area, aes(x=reorder(log(number_pixels), -str_ord), y=FDw ))
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
p <- p + xlab("log area")+ylab("total FD")
p


svg(paste0(path, "/total_FD_area.svg"))
# png(paste0(path, "/total_FD_area.png"))
plot(p)
dev.off()



log(total_area$number_pixels)


### Area vs. water volume per stream order
x11()
p <- ggplot(total_area, aes(x=str_ord, y=log(number_pixels*1000000) ))
p <- p + geom_point(color="black", pch=22, fill="#ADD8E6", stroke = 2, size=12)
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
p <- p + xlab("stream order")+ylab("area log m2")
p

# 
svg(paste0(path, "/total_area_network_log.svg"))
# png(paste0(path, "/total_area_network_log.png"))
plot(p)
dev.off()




x11()
p <- ggplot(total_area, aes(x=str_ord, y=log(volume_per_strord)  ))
p <- p + geom_point(color="black", pch=24, fill="#ADD8E6", stroke = 2, size=12)
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
p <- p + xlab("stream order")+ylab("volume log m3")
p

svg(paste0(path, "/total_volume_network_log.svg"))
# png(paste0(path, "/total_volume_network_log.png"))
plot(p)
dev.off()







### Plot SR

### Volume
p <- ggplot(total_area, aes(x=reorder(volume_per_strord_log, str_ord), y=n_sp ))
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
p <- p + xlab("log volume")+ylab("total SR")
# p <- p + scale_x_continuous(breaks = pretty(as.numeric(total_area$volume_per_strord_log), n = 4))


svg(paste0(path, "/total_SR_volume.svg"))
# png(paste0(path, "/total_SR_volume.png"))
plot(p)
dev.off()




### area
p <- ggplot(total_area, aes(x=reorder(log(number_pixels), -str_ord), y=n_sp ))
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
p <- p + xlab("log area")+ylab("total SR")
p
# p <- p + scale_x_continuous(breaks = pretty(as.numeric(total_area$volume_per_strord_log), n = 4))

svg(paste0(path, "/total_SR_area.svg"))
# png(paste0(path, "/total_SR_area.png"))
plot(p)
dev.off()





### number of pixels
p <- ggplot(total_area, aes(x=reorder(log(number_pixels), -str_ord), y=n_sp ))
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
p <- p + xlab("log area")+ylab("total SR")
p
# p <- p + scale_x_continuous(breaks = pretty(as.numeric(total_area$volume_per_strord_log), n = 4))

svg(paste0(path, "/total_pixel_area.svg"))
# png(paste0(path, "/total_SR_area.png"))
plot(p)
dev.off()


