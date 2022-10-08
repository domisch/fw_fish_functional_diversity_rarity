

###---- Get distinctness (ED) as raster ----
library(plyr)

### Merge directly to big table
eco_ED <- read.csv(paste0(path, "/traits_eco/master_table_spp_fam_comm_ED_RS_volume_residuals.csv"))
tmp <- merge(tmp_rs, eco_ED[c("sequence_original", "ED")], by.x="species", by.y="sequence_original", all.x=T)

### Check if duplicates
tmp <- tmp[!duplicated(tmp),]


### Median ED per grid cell
med_ED_perGrid <- tmp[, as.double(median(ED)),by = grid_id_long]
names(med_ED_perGrid) <- c("grid_id_long", "med_ED")



### Load raster/domain template
load(paste0(path, "/global_layers/additional_layers/grid_id_long_with_lakes.RData"))
domain_raster <- grid_id_long  ; rm(grid_id_long); gc()
domain_cells <- as.data.table.raster(domain_raster)
names(domain_cells) <- "grid_id_long"
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key

### Assign the predicted values to the raster 
med_ED_NAwide <- merge(domain_cells, med_ED_perGrid, by="grid_id_long", all.x=T)
med_ED_NAwide <- arrange(med_ED_NAwide, seq_id) 
med_ED_NAwide_r <- setValues(domain_raster, med_ED_NAwide$med_ED)


### Mask Mexico
mask_mexico <- shapefile(paste0(path, "/shape_mexico/mask_mexico2.shp"))
### Border to the US kep original, rest simplified to speed up

med_ED_NAwide_r_mex <- mask(med_ED_NAwide_r, mask_mexico, inverse=T, updatevalue=NA)
writeRaster(med_ED_NAwide_r_mex, paste0(path, "/stacked_outputs/eco_median_ED_per_grid_mex_masked.tif"), overwrite=T)


