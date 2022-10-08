

### Get all possible freshwater range maps. Start with range map data from Page & Burr
library(maptools)
library(foreign)
library(dismo)
library(maps)
library(sp)
library(rgeos)
library(rgdal)
library(stringr)
library(rgeos)
library(Hmisc)
library(knitr)
library(data.table)

### Set path and options
setwd("/lustre/scratch/client/fas/sbsc/fw_fish/")
`%ni%` <- Negate(`%in%`)
options(stringsAsFactors = FALSE)


pageburr <- read.csv("all_species_NA_PageBurr.csv", h=F) # corrected for typos
names(pageburr) <- "species"
pageburr$species <- as.character(pageburr$species)


### Page & Burr contains 826 rangemaps in the shapefile (the book has more species but no maps, but these contain also species complexes / not yet resolved hybridizations, which we do not account in this study. 
### Supplement the additional ones from NatureServe
natureserve_species <- as.data.frame(dir("NatureServe/Shapefiles", pattern="dbf"))
names(natureserve_species) <- "species"
natureserve_species$species <- gsub("\\.dbf", "", natureserve_species$species)
natureserve_species$species <- gsub("_", " ", natureserve_species$species)


### Check names in EOL and Catalog of Fishes
library(taxize)
sources <- gnr_datasources()
knitr::kable(sources)
eol <- sources$id[sources$title == 'EOL']
cof <- sources$id[sources$title == 'Catalog of Fishes']


results <- list()

for (i in pageburr$species) {
### Paste loop progress
cat("Searching name for", i, "\n")
### Retrieve best match from all data bases
tmp <- try(gnr_resolve(names = i,
                      resolve_once = T,
                      stripauthority = T,
                      best_match_only = T, 
                      data_source_ids=eol), #cof
      silent=T)
### Save only if no error
if ( grepl("Error",tmp) == FALSE) { results[[i]] <- tmp }

}
### Export to fix small issues (remove remaining authority etc.)
checked_names <- do.call(rbind, results)
write.csv(checked_names, "pageburr_EOL.csv", row.names = F, quote = F) # fixed in excel
checked_names <- read.csv("pageburr_EOL.csv", h=T)
check_synonyms <- checked_names$matched_name




### Run in a loop
library(taxize)
checked_names <- read.csv("pageburr_EOL.csv", h=T)
check_synonyms <- checked_names$matched_name
results <- list()
for (i in check_synonyms) {
### Paste loop progress
cat("Searching synonyms for", i, "\n")
tmp <- try(synonyms(i, db="itis"),
silent=T)
### Save only if no error
if ( grepl("NA",tmp) == FALSE) { results[i] <- tmp } # only 1 square-bracket in this case!
}



check_synonyms_out <- rbindlist(results, fill=T,  use.names=T, idcol=T)
head(check_synonyms_out)
str(check_synonyms_out)
write.csv(check_synonyms_out, "pageburr_EOL_synonyms.csv", row.names = F, quote = F) 
### Added those species again that were dropped in synonyms() above
check_synonyms_out <- read.csv("pageburr_EOL_synonyms.csv", h=T) 

fix(check_synonyms_out)


### Check which original or synonum species names are not in NatureServe? Get these as addional species
"%ni%"  <- Negate("%in%") 
### Check in original names
additional_species <- natureserve_species$species [which(natureserve_species$species %ni% check_synonyms_out$.id) ] # 109 found

### Check if the can be found in synonyms
additional_species <- additional_species [which(additional_species %ni% check_synonyms_out$syn_name) ] # 16 found
###--> a total of 93 candidates could be added from NatureServe
### Check which ones of these are actual species level data
write.table(additional_species, "potential_additional_species_in_NatureServe.csv", row.names=F, col.names=F)

### Fixed in Excel. Check if these potential species can be found in Page & Burr after correction
additional_species <- read.csv("potential_additional_species_in_NatureServe.csv", h=T)
additional_species <- as.data.frame(na.omit(additional_species)) # 56 remain
additional_species <- additional_species [which(additional_species$corrected_species %ni% check_synonyms_out$.id), ] # 4 dropped
additional_species <- additional_species [which(additional_species$corrected_species %ni% check_synonyms_out$syn_name), ] # 0 dropped

### For the 53 new species, check for synonyms:
results <- list()
for (i in additional_species$corrected_species) {
### Paste loop progress
cat("Searching synonyms for", i, "\n")
tmp <- try(synonyms(i, db="itis"),
silent=T)
### Save only if no error
if ( grepl("NA",tmp) == FALSE) { results[i] <- tmp } 
}

library(data.table)
check_synonyms_out <- rbindlist(results, fill=T,  use.names=T, idcol=T)
head(check_synonyms_out)
str(check_synonyms_out)
write.csv(check_synonyms_out, "natureserve_56additional_species_with_synonyms.csv", row.names = F, quote = F) 
### Any missing species in the search?
additional_species[which(additional_species$corrected_species %ni% check_synonyms_out$.id),]
### --> add these 4 species in the table





### Get a table with all possible species names+synonyms, and the correct name as a new masterlist
### Use this table to go back to the shapefiles as well

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
masterlist <- read.csv("pageburr_natureserve_with_synonyms_lookup_table.csv", h=T) # 869 unique species
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

### use ".id" as the names in the shapefiles, and to scan for records
### Use "syn_name"  name to scan for additional records
### rangemap_source: where can the rangemap be found?
### raw_name_in_shp: species name in the rangemap shapefile




tmp1<- as.data.frame(check_synonyms_out$.id)
names(tmp1) <- "species"
tmp2<- as.data.frame(check_synonyms_out$syn_name)
names(tmp2) <- "species"

masterlist <- rbind(tmp1, tmp2); rm(tmp1, tmp2)
masterlist <- as.data.frame(masterlist[!is.na(masterlist)]) 
names(masterlist) <- "species"

### Remove only Genus -entries
masterlist = as.data.frame(masterlist[str_count(masterlist$species, " ") >= 1,])
names(masterlist) <- "species"
masterlist<- as.data.frame(masterlist[!duplicated(masterlist$species),])
names(masterlist) <- "species"
masterlist <- droplevels(masterlist) # clean up possible old names
write.csv(masterlist, "masterlist_PageBurr_Natureserve_species_synonyms.csv", row.names=F, quote=F)
# use this list to scan for species in all databases


### Create a template table with 
### species in shapefile | original name and possible synonyms
check_synonyms_out <- read.csv("pageburr_natureserv_with_synonyms.csv", h=T) 
orig <- tmp1<- as.data.frame(check_synonyms_out$.id)

### Use Page&Burr as main list, and supplement with NatureServe (only USA, Page&Burr has larger geographic extent)
### Combine both lists
masterlist <- rbind(pageburr, natureserve_species)
### Remove duplicates
masterlist <- masterlist[!duplicated(masterlist), ]

### Check synonyms
### Load raw species lists again:
`%ni%` <- Negate(`%in%`)
options(stringsAsFactors = FALSE)

pageburr <- read.csv("all_species_NA_PageBurr.csv", h=F) # corrected
names(pageburr) <- "species"
pageburr$species <- as.character(pageburr$species)


### Page & Burr have 909 species but only 826 rangemaps in the shapefile. Get the missing ones from NatureServe
natureserve_species <- as.data.frame(dir("NatureServe/fishShapefiles_old/Shapefiles", pattern="dbf"))
names(natureserve_species) <- "species"
natureserve_species$species <- gsub("\\.dbf", "", natureserve_species$species)
natureserve_species$species <- gsub("_", " ", natureserve_species$species)

### Which species do not overlap?
"%ni%"  <- Negate("%in%") 
additional_species <- natureserve_species$species [which(natureserve_species$species %ni% pageburr$species )] # 111 species are extra

###=============================================================================================#
### Read the table with all original species, and the synonyms for each
masterlist <- read.csv("original_corrected_synonyms_sp._removed.csv", h=T)
masterlist <- data.frame(lapply(masterlist, as.character), stringsAsFactors=FALSE)
# original = merged from Page&Burr and natureserve
# EOL_corrected = original names cheecked against EOL
# search_these = the EOL names, manually removed a few authorities, no "sp."" etc. Use this as the masterlist to scan for records
# type= is there a synonym
###=============================================================================================#

additional_species_syn <- as.data.frame(additional_species [which(additional_species %in% masterlist$search_these )]) # 66 are synonyms
names(additional_species_syn) <- "species"
### Get the original name for these
additional_species_syn <- merge(additional_species_syn, subset(masterlist, select=c(original, search_these)),
                                by.x="species", 
                                by.y="search_these", 
                                all.x=T)

additional_species_syn <-  data.frame(lapply(additional_species_syn, as.character), stringsAsFactors=FALSE)


### Correct the table: if "original" also present in shapefiles ("original" in masterlist, then remove)
additional_species_syn$double_entry <- ifelse(additional_species_syn$species == additional_species_syn$original, 1, 0)
additional_species_syn <- additional_species_syn[additional_species_syn$double_entry==0,]
### --> these species from NatureServe are NOT needed as they are synonyms from Page&Burr
### Remove from the "additional_species" table
additional_species <- as.data.frame(additional_species [which(additional_species %ni% additional_species_syn$species )]) 
names(additional_species) <- "species" 
### --> these are potentical candidates for additional species
### Check which ones are in the masterlist table (cleaned = there are no "sp." etc...)
additional_species <- data.frame(lapply(additional_species, as.character), stringsAsFactors=FALSE) # make as.charcter


additional_species_extra <- additional_species [which(additional_species$species %ni% masterlist$original ),]
additional_species_extra <- additional_species [which(additional_species$species %ni% masterlist$search_these ),]
as.data.frame(additional_species_extra) # check





###----------- GBIF dataset ------------------

library(foreach)
library(doParallel)
library(data.table)

genspp <- read.table("PageBurr_natureserve_only_unique_genus.txt", h=F)
genspp <- as.character(genspp$V1)


### Make cluster object
cl <- makePSOCKcluster(15) 
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers


#---- First run with download=F to check, then remove all zero-record-species and run again with download=T
gbif_results <- foreach(i=genspp, .final=c, .packages = c("dismo"), .verbose=T) %dopar% { # 
tmp <- gbif(genus=i, species='*', geo=F, removeZeros=F, ext=extent(-188, -51, 5, 85), download=F) # * wildcard used to get all spp. within each genus
if(length(tmp)>0) { return(tmp) } # return only if not empty
}
stopCluster(cl)

gbif_download <- as.data.frame(do.call(rbind, gbif_results))
gbif_download <- cbind(genspp, gbif_download) # attach species names
gbif_download <- gbif_download[gbif_download[2]>0,] # remove all zero-record-species
gbif_download <- as.character(gbif_download$genspp)

#---- Run again with download=T

### Make cluster object
cl <- makePSOCKcluster(15) 
registerDoParallel(cl) # register parallel backend
getDoParWorkers() # show number of workers

#---- Run again with download=T
gbif_results <- foreach(i=gbif_download, .final=c, .packages = c("dismo"), .verbose=T) %dopar% { # 
tmp <- try(gbif(genus=i, species='*', geo=F, removeZeros=F, ext=extent(-188, -51, 5, 85), download=T), 
           silent=T) # * wildcard used to get all spp. within each genus
# ### Save only if no error
if ( grepl("Error",tmp) == TRUE) { return( data.frame(error=1)) } else { return(tmp) }
}
stopCluster(cl)


gbif_results2 <- rbindlist(gbif_results, fill=T) # some species have additional columns?
save(gbif_results2, file="gbif_raw_genus_download.RData")




###=============================================================================================#
masterlist <- read.csv("original_corrected_synonyms_sp._removed.csv", h=T)

# original = merged from Page&Burr and natureserve
# EOL_corrected = original names cheecked against EOL
# search_these = the EOL names, manually removed a few authorities, no sp. etc. Use this as the masterlist to scan for records
# type= is there a synonym
###=============================================================================================#




###---------- Clean GBIF data ---------

### Load raw download data
load("GBIF/gbif_raw_genus_download.RData") # gbif_results2
str(gbif_results2)

gbif_results2$myid <- seq(1:nrow(gbif_results2))

### Get all records that have coordinates
gbif_results2_coord <- subset(gbif_results2, !is.na(gbif_results2$lat) | !is.na(gbif_results2$lon))
summary(gbif_results2_coord$lat)
summary(gbif_results2_coord$lon)

### Export those records without coordinates
"%ni%" <- Negate("%in%") # create a "not in" -function
gbif_results2_locality <- which(gbif_results2$myid %ni% gbif_results2_coord$myid) # get those ID's that are without coordinates
gbif_results2_locality <- gbif_results2[gbif_results2_locality,] 
gbif_results2_locality[50:70,15:25]
sort(unique(gbif_results2_locality$country))

### Filter by country
gbif_results2_locality_short <- subset(gbif_results2_locality, 
                                       gbif_results2_locality$country == "United States" | 
                                        gbif_results2_locality$country == "canada" |
                                         gbif_results2_locality$country == "Mexico")


gbif_results2_locality_short <- subset(gbif_results2_locality_short, 
                                       select=c(cloc, country, adm1, adm2, waterBody, species, year))

### Add waterbody to locality
gbif_results2_locality_short$locality <-  do.call(paste, c(gbif_results2_locality_short[c("cloc", "waterBody")], sep=", "))
gbif_geolocate <- cbind(gbif_results2_locality_short["locality"],
                        gbif_results2_locality_short["country"],
                        gbif_results2_locality_short["adm1"],
                        gbif_results2_locality_short["adm2"])

### Add columns
gbif_geolocate[,"latitude"] <- NA
gbif_geolocate[,"longitude"] <- NA
gbif_geolocate[,"correction_status"] <- NA
gbif_geolocate[,"precision"] <- NA
gbif_geolocate[,"error_polygon"] <- NA
gbif_geolocate[,"multiple_results"] <- NA
gbif_geolocate[,"uncertainty"] <- NA

gbif_geolocate <- cbind(gbif_geolocate, 
                        gbif_results2_locality_short["species"], 
                        gbif_results2_locality_short["year"])
gbif_geolocate$newid <- do.call(paste0, c(gbif_geolocate[c("locality", "country", "adm1", "adm2")]))
gbif_geolocate <- gbif_geolocate[!duplicated(gbif_geolocate),] # remove duplicates
write.csv(gbif_geolocate, "GBIF/gbif_geolocate.csv", row.names=F) 




### Check taxonomy using taxize
names <- gbif_geolocate$species
names <- names[!duplicated(names)]
names <- names[-2] # remove NA
### clean up and free memory
rm(gbif_results2, gbif_results2_coord, gbif_results2_locality, gbif_results2_locality_short); gc()

library(taxize)
sources <- gnr_datasources()
knitr::kable(sources)
eol <- sources$id[sources$title == 'EOL']
cof <- sources$id[sources$title == 'Catalog of Fishes']

results <- list()

for (i in names) {
### Paste loop progress
cat("Searching name for", i, "\n")
### Retrieve best match from all data bases
tmp <- try(gnr_resolve(names = i,
                      resolve_once = T,
                      stripauthority = T,
                      best_match_only = T, 
                      data_source_ids=eol), #cof
      silent=T)
### Save only if no error
if ( grepl("Error",tmp) == FALSE && length(tmp)>0 ) { results[[i]] <- tmp } # or | length(tmp)>0
}



### Export to fix small issues (remove remaining authority etc.)
checked_names <- do.call(rbind, results)


### Which names were not saved due to errors?
"%ni%" <- Negate("%in%")
names[names %ni% checked_names$submitted_name] #"Diplesium blennioides"  "Poxomis nigromaculatus", both not in catalogue of fishes

write.csv(checked_names, "GBIF/checked_names_for_geolocate_species_EOL.csv", row.names = F, quote = F) # fixed in excel
checked_names <- read.csv("GBIF/checked_names_for_geolocate_species_EOL.csv", h=T)

### Load output file from GeoLocate
gbif_geolocate <- read.csv("GBIF/gbif_geolocate_output.csv", h=T) 

### Merge corrected names to species records and clean
gbif_geolocate <- merge(gbif_geolocate, checked_names, by.x="species", by.y="submitted_name")
### Remove all records without lat/long
gbif_geolocate <- subset(gbif_geolocate, gbif_geolocate$latitude!=0)
gbif_geolocate <- subset(gbif_geolocate, gbif_geolocate$longitude!=0)

### Save file
save(gbif_geolocate, file="GBIF/gbif_geolocate_lat_long_species_cleaned.RData") # ready, use "matched_name"
load("GBIF/gbif_geolocate_lat_long_species_cleaned.RData") 



### Clip all records outside North America
min_lat <- 5
max_lat <- 85
min_lon <- -188
max_lon <- -51
  

### Identify the points within the bbox (==1, otherwise zero)
gbif_results2_coord$points_in_NA <- ifelse(gbif_results2_coord$lat >= min_lat & gbif_results2_coord$lat <= max_lat 
                                  & gbif_results2_coord$lon >= min_lon & gbif_results2_coord$lon <= max_lon, 1, 0)    
  
### Drop the points outside the bbox
gbif_results2_coord <- subset(gbif_results2_coord, gbif_results2_coord$points_in_NA == 1 )
gbif_results2_coord <- gbif_results2_coord[-182] # remove bounding-box identifier..



### The lon, lat, species, year, coord_uncertainty, abundance, number_visits
gbif_results2_coord <- gbif_results2_coord[c("lon", "lat", "species", "year", "coordinateUncertaintyInMeters")]
colnames(gbif_results2_coord)[3:5] <- c("species", "year", "coord_uncertainty")



### Make spatial
gbif_results2_coord <- as.data.frame(gbif_results2_coord)
coordinates(gbif_results2_coord) <- c("lon", "lat")
projection(gbif_results2_coord) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"


### LAEA in meters (the initial hydro1k prj)
dist_proj =  CRS("+proj=laea  +lat_0=45
                 +lon_0=-100
                 +x_0=0
                 +y_0=0")



gbif_results2_coord_laea <- spTransform(gbif_results2_coord,dist_proj)
# NA_shp_laea <- spTransform(NA_shp,dist_proj)


# load("GBIF/NA_fish_GBIF.Rdata") 
NA_shp_buff_50km <- readShapePoly("/lustre/scratch/client/fas/sbsc/fw_fish/shapes/basins_NA/NA_dissolved/na_dissolved_buff_50km.shp", verbose=T, proj4string=CRS("+proj=longlat +datum=WGS84"))
NA_shp_laea_buff <- spTransform(NA_shp_buff_50km,dist_proj)

gbif_results2_coord$dist_km <- gDistance(gbif_results2_coord_laea, NA_shp_laea_buff, byid=T)[1,]


## that adds 'distance' (in meters) from each point to the polygon, check out the values:
x11(); hist(gbif_final$dist_km/1000,xlab="km to border North America")


### Drop points that have a distance value >0 (outside the 50km buffer...)
gbif_results2_coord_clip  <- gbif_results2_coord[!gbif_results2_coord$dist_km > 0,]

### Map in LAEA
x11(); plot(spTransform(NA_shp,dist_proj))
points(spTransform(gbif_final_clip,dist_proj), col = "blue")


### Map in WGS
x11(); plot(NA_shp)
points(gbif_final_clip, col = "blue")

### Make GBIF (clipped to NA border) a df for later on...
gbif_NA_df <- as.data.frame(gbif_results2_coord_clip)
gbif_NA_df <- gbif_NA_df [-6] # drop column for used for clipping
names(gbif_NA_df)[1:2] <- c("lon", "lat")
head(gbif_NA_df)

# other method for identifying points within polygons..but RAM-intensive!  
# http://www.nceas.ucsb.edu/scicomp/usecases/point-in-polygon
# gbif_in_NA_shp <- !is.na(over(na_fish_gbif_826, as(NA_shp, "SpatialPolygons")))
save(gbif_NA_df, file="GBIF/gbif_NA_clipped_lon_lat_all_years.RData")
load("GBIF/gbif_NA_clipped_lon_lat_all_years.RData")


### Check names
names <- as.character(unique(gbif_NA_df$species))
names <- names[-1]

library(taxize)
sources <- gnr_datasources()
knitr::kable(sources)
eol <- sources$id[sources$title == 'EOL']
cof <- sources$id[sources$title == 'Catalog of Fishes']

results <- list()

for (i in names) {
### Paste loop progress
cat("Searching name for", i, "\n")
### Retrieve best match from all data bases
tmp <- try(gnr_resolve(names = i,
                      resolve_once = T,
                      stripauthority = T,
                      best_match_only = T, 
                      data_source_ids=eol), #cof
      silent=T)
### Save only if no error
if ( grepl("Error",tmp) == FALSE && length(tmp)>0 ) { results[[i]] <- tmp } # or | length(tmp)>0
}



### Export to fix small issues (remove remaining authority etc.)
checked_names <- do.call(rbind, results)

### Which names were not saved due to errors?
"%ni%" <- Negate("%in%")
names[names %ni% checked_names$submitted_name] #"Diplesium blennioides"  "Poxomis nigromaculatus", both not in catalogue of fishes

write.csv(checked_names, "GBIF/checked_names_with_lat_long_species_EOL.csv", row.names = F, quote = F) # fixed in excel
checked_names <- read.csv("GBIF/checked_names_with_lat_long_species_EOL.csv", h=T)

### Merge corrected names to species records and clean
gbif_NA_df <- merge(gbif_NA_df, checked_names, by.x="species", by.y="submitted_name")
save(gbif_NA_df, file="GBIF/gbif_NA_clipped_lon_lat_all_years_names_checked.RData") # use "matched_name"
# gbif data ready




###----------- USGS BioData ------------------
# https://aquatic.biodata.usgs.gov

### Add the BioData dataset, new download on 1Jul2016
biodata_fish <- read.csv("E:/xubuntu_shared/freshwater/NA/tables/species_download/BioData/20160630.1341.FishCount.csv", h=T)

### get the coordinates
biodata_sites <- read.csv("E:/xubuntu_shared/freshwater/NA/tables/species_download/BioData/20160630.1341.SiteInfo.csv", h=T)

# for each "SiteNumber" there is the "BioDataTaxonName", first merge the data:
biodata_merge_fish <- merge(biodata_sites, biodata_fish, by="SiteNumber")
head(biodata_merge_fish)


### Remove rows without coords
biodata_merge_fish <- biodata_merge_fish[!is.na(biodata_merge_fish$Longitude_dd),]
biodata_merge_fish <- biodata_merge_fish[!is.na(biodata_merge_fish$Latitude_dd),]

### Remove duplicates
biodata_merge_fish <- biodata_merge_fish[!duplicated(biodata_merge_fish),]

str(biodata_merge_fish)
### How many species and sites?
length(unique(biodata_merge_fish$Species)) # 497 (old: 488 species)
length(unique(biodata_merge_fish$SiteNumber)) # 1258 (old: 1206 sites)


### Subset the data for the follwoing projections
biodata_NAD27 <- subset(biodata_merge_fish, CoordinateDatum == "NAD27")
biodata_WGS84 <- subset(biodata_merge_fish, CoordinateDatum == "NAD83") # same a s WGS84
# biodata_WGS84 <- subset(biodata_merge_fish, CoordinateDatum == "WGS84")
biodata_OLDHI <- subset(biodata_merge_fish, CoordinateDatum == "OLDHI") # Hawaii

### Make these a spatial object, define prj, project to wgs and merge to one singe set
coordinates(biodata_NAD27) <- c("Longitude_dd", "Latitude_dd")
# coordinates(biodata_NAD83) <- c("Longitude_dd", "Latitude_dd") # same as WGS84
# coordinates(biodata_OLDHI) <- c("Longitude_dd", "Latitude_dd") # not needed

### Define projection
proj4string(biodata_NAD27) <- "+proj=longlat +datum=NAD27"
### Project to WGS84
biodata_NAD27_to_WGS84 <- spTransform(biodata_NAD27, CRS("+proj=longlat +datum=WGS84"))
### Prepare names and merge...
biodata_NAD27_to_WGS84 <- as.data.frame(biodata_NAD27_to_WGS84)


biodata_merged <- rbind(biodata_NAD27_to_WGS84, biodata_WGS84) # drop 2 locations from Hawaii


### Pick only coords, species, year, coord-uncertainty, 
biodata_merged_small <- biodata_merged[c("Longitude_dd", "Latitude_dd", "BioDataTaxonName", "CollectionYear")]

### Rename columns
names(biodata_merged_small) <- c("lon", "lat", "species", "year")

### Add column with coord_uncertainty
biodata_merged_small$coord_uncertainty <- NA

### Add column with abundance
biodata_merged_small$abundance <- biodata_merged$Abundance

### Add column with number of visits / samples
biodata_merged_small$number_visits <- NA

head(biodata_merged_small)

### Make these a spatial object, define prj
coordinates(biodata_merged_small) <- c("lon", "lat")
### Define projection
### http://spatialreference.org/
### http://www.remotesensing.org/geotiff/proj_list/
projection(biodata_merged_small) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

x11(); plot(NA_shp); points(biodata_merged_small, col="blue")

biodata_merged <- as.data.frame(biodata_merged_small)
save(biodata_merged, file="BioData/biodata_merged_species_all_years.RData")
# BioData fish ready





###----------- fishhabitat dataset ---------------------- # 2016: no new data available
# National fish habitat 
# http://ecosystems.usgs.gov/fishhabitat/ # download
# http://fishhabitat.org/
### Species records are also linked to NHDPlus (and match the river neywork..)

### Load geographic locations with lat/lon
fishhabitat_loc <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/National_Fish_Habitat/Sample_Descriptions_Sampled_Fish_export.dbf", as.is = F)


### Load the fish species found at those locations
fishhabitat_spp <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/species_download/National_Fish_Habitat/Sampled_Fish_export.dbf", as.is = F)

### Merge
intersect(names(fishhabitat_loc), names(fishhabitat_spp)) # "OBJECTID" "NFHAP_ID"

fishhabitat <- merge(fishhabitat_loc, fishhabitat_spp, "NFHAP_ID")

### Remove obsolete columns and remove duplicates
fishhabitat <- fishhabitat[c("LONG", "LAT", "SCIENTIFIC", "SAMP_DATE")]
fishhabitat$SAMP_DATE <- substring(fishhabitat$SAMP_DATE, 1, 4) # get first 4 numbers = year (1990-06-18)

### Rename columns
colnames(fishhabitat) <- c("lon", "lat", "species", "year")
fishhabitat <- fishhabitat[!duplicated(fishhabitat),]   # remove duplicate records

### Add coord_uncertainty, abundance and number of visists (have no data...)
fishhabitat$coord_uncertainty <- NA
fishhabitat$abundance <- NA
fishhabitat$number_visits <- NA


### Make a spatial object
fishhabitat_sp <- fishhabitat
coordinates(fishhabitat_sp) <- c("lon", "lat")

### Scan for species in masterlist
# fishhabitat_species <- merge(fishhabitat, masterlist, by="species")
# unique(fishhabitat$species)
save(fishhabitat, file="National_Fish_Habitat/fishhabitat_species_all_years.RData")

### Plot
x11(); plot(NA_shp); points(fishhabitat_species[c("lon", "lat")], col = "blue")
### fishhabitat data ready




###----- StreamNet data for North West North America -----
streamnet <- read.csv("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/StreamNet/lon_lat_species_yearCSV.csv", h=T)

colnames(streamnet) <- c("lon", "lat", "BarrierID", "species", "year")
streamnet <- streamnet[!is.na(streamnet$lon),]   
streamnet <- streamnet[!is.na(streamnet$lat),]   
length(unique(streamnet$species)) # 74294 rows with 43 species
streamnet <- streamnet[-3] # remove barrier column


### Add Coord.unvertainty, abundance and visists
streamnet$coord_uncertainty <- NA
streamnet$abundance <- NA
streamnet$number_visits <- NA

### Clean species names. 
### Remove: Cottus spp. and Oncorhynchus mykiss subsp.
unique(streamnet$species)
streamnet <- streamnet[streamnet$species != "Cottus spp.",]
streamnet <- streamnet[streamnet$species != "Oncorhynchus mykiss subsp.",]

### Make spatial object
streamnet_sp <- streamnet
coordinates(streamnet_sp) <- c("lon", "lat")

### Plot
x11(); plot(NA_shp); points(streamnet_sp, col = "blue", pch = 16, cex = 0.5)
x11(); plot(streamnet_sp, 
            pch = 16, cex = 0.5, 
#             xlim = c(-140, -100), ylim = c(40, 50), 
            col = "blue")
plot(NA_shp, add = T)

### Scan for species in masterlist
# streamnet_species <- merge(streamnet, masterlist, by="species")
# unique(fishhabitat$species)
save(streamnet, file="StreamNet/streamnet_species_all_years.RData")
### StreamNet data ready



###----- EPA Mercury Fish data set -------
epa_sites <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/EPA_Mercury_Fish/hgfish/HgData.dbf", as.is=F) 

epa_species <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/EPA_Mercury_Fish/hgfish/Species.dbf", as.is=F) 

### Get the species names
epa_species$latin <- do.call(paste, c(epa_species[c("GENUS", "SPECIES")], sep = " "))

# temp_epa <- merge(fish_check, epa_species, by.x = "Species", by.y = "latin")
# temp_epa <- merge(masterlist, epa_species, by.x = "species", by.y = "latin")
# epa_species <- temp_epa[c("species", "SPECIES_ID")] ; rm(temp_epa)

### Merge 
intersect(names(epa_species), names(epa_sites)) # "SPECIES_ID"
epa_spp_sites <- merge(epa_species, epa_sites, by = "SPECIES_ID")

### Reduce to lon, lat, species
epa_fish <- epa_spp_sites[c("LONGITUDE", "LATITUDE", "latin", "DATE", "NUMBER_IN_")]
colnames(epa_fish) <- c("lon", "lat", "species", "DATE", "abundance")


### Clean the date
epa_fish$DATE <- substring(epa_fish$DATE, 1, 4) # get first 4 numbers = year (1990-06-18)
colnames(epa_fish)[4] <- "year" 
colnames(epa_fish)[3] <- "species"

### Add column and re-order
epa_fish$coord_uncertainty <- NA
epa_fish$number_visits <- NA
epa_fish <- cbind(epa_fish[1:4], epa_fish[6], epa_fish[5], epa_fish[7])


### Clean and check
epa_fish <- epa_fish[!is.na(epa_fish$lon),]
epa_fish <- epa_fish[!is.na(epa_fish$lat),]
length(unique(epa_fish$species)) # 193

### Make spatial
epa_fish_sp <- epa_fish
coordinates(epa_fish_sp) <- c("lon", "lat")


### Plot
x11(); plot(NA_shp); points(epa_fish_sp, col = "blue", pch = 16, cex = 0.5)


save(epa_fish, file="EPA_Mercury_Fish/epa_fish_species_all_years.RData")
### EPA Mercury data ready




###------ Canada British Columbia ----- 

### Read coordinates with species codes (2 sets)
sites1_bc <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/Canada_British_Columbia/streams/lfdpbcgz_WGS.dbf", as.is=T)
sites2_bc <- read.dbf("/lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/Canada_British_Columbia/streams/lfdsbcgz_WGS.dbf", as.is=T)

### Reduce df
sites1_bc <- sites1_bc[c("POINT_X", "POINT_Y", "SPECIES_CO", "OFFSET2", "OFFSET1", "REFS_AND_D")]
sites2_bc <- sites2_bc[c("POINT_X", "POINT_Y", "SPECIES_CO", "OFFSET2", "OFFSET1", "REFS_AND_D")]

sites_bc <- rbind(sites1_bc, sites2_bc)

### Read species codes (used the older version of codes..) new version also downloaded
spp_codes_bc <- read.csv("E:/xubuntu_shared/freshwater/NA/tables/species_download/Canada_British_Columbia/alt/fish_codes_csvR.csv", h=T)

### Merge by codes
fish_bc <- merge(sites_bc, spp_codes_bc, by.x = "SPECIES_CO", by.y = "CODE")

### Clean data
### OFFSET1 and 2 have both zero-values
fish_bc <- fish_bc[c("POINT_X", "POINT_Y", "LATIN.NAME", "OFFSET2", "REFS_AND_D")]
colnames(fish_bc) <- c("lon", "lat", "species", "coord_uncertainty", "year")

### Clean the year
temp <-  sub("^.* ", "", fish_bc$year)
# This regular expression matches the beginning of the string (^), any character (.) repeated zero or more times (*), and underscore (_). 
# The ? makes the match "lazy" so that it only matches are far as the first underscore. That match is replaced with just an underscore. See ?regex for more details and references 

### Get year and attach to colums
temp <- substring(temp, 8, 11) # get characters from first 8 to 11 = year 
temp[temp=="r>"] <- NA # samples without a year 
fish_bc$year <-temp

### Add columns and clean data
fish_bc$abundance <- NA
fish_bc$number_visits <- NA

fish_bc <- fish_bc[!is.na(fish_bc$year),] # missing years?
fish_bc <- cbind(fish_bc[1:3], fish_bc[5], fish_bc[4], fish_bc[6:7])

# fish_bc_species <- merge(fish_bc, masterlist, by="species")
dim(fish_bc)


### Make spatial object
fish_bc_sp <- fish_bc
coordinates(fish_bc_sp) <- c("lon", "lat")

### Plot
x11(); plot(NA_shp); points(fish_bc_sp, col = "blue", pch = 16, cex = 0.5)
save(fish_bc, file="Canada_British_Columbia/fish_bc_species_all_years.RData")
### Canada British Columbia data ready






###--------- fishnet2 -------------------------
load("E:/Yale_docs/fishnet2_all_data_9Jun15/final/fishnet2_Jun15_taxonomy_checked.RData")
head(fishnet2_Jun15_taxonomy_checked)
fishnet <- fishnet2_Jun15_taxonomy_checked
### Restrict to North America
min_lat <- 5
max_lat <- 85
min_lon <- -188
max_lon <- -51

summary(fishnet)  

### Identify the points within the bbox (==1, otherwise zero)
fishnet$points_in_NA <- ifelse(fishnet$Latitude >= min_lat & fishnet$Latitude <= max_lat 
                         & fishnet$Longitude >= min_lon & fishnet$Longitude <= max_lon, 1, 0)    
  
### Drop the points outside the bbox
fishnet_final <- subset(fishnet, fishnet$points_in_NA == 1 )

names(fishnet_final)[2] <- "species"

### Check for species in the masterlist
fishnet_final <- merge(fishnet_final, masterlist, by="species")
dim(fishnet_final)

### Plot the data on a map
x11(); plot(fishnet_final[c("Longitude", "Latitude")]) # lat/long
plot(wrld, add=T)

### Make spatial
coordinates(fishnet_final) <- c("Longitude", "Latitude")
projection(fishnet_final) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"


### LAEA in meters (the initial hydro1k prj)
dist_proj =  CRS("+proj=laea  +lat_0=45
                 +lon_0=-100
                 +x_0=0
                 +y_0=0")

fishnet_final_laea <- spTransform(fishnet_final,dist_proj)

NA_shp_buff_50km <- readShapePoly("/lustre/scratch/client/fas/sbsc/fw_fish/shapes/basins_NA/NA_dissolved/na_dissolved_buff_50km.shp", verbose=T, proj4string=CRS("+proj=longlat +datum=WGS84"))
NA_shp_laea_buff <- spTransform(NA_shp_buff_50km,dist_proj)
x11(); plot(NA_shp_laea_buff)

fishnet_final$dist_km <- gDistance(fishnet_final_laea, NA_shp_laea_buff, byid=T)[1,]

## that adds 'distance' (in meters) from each point to the polygon, check out the values:
x11(); hist(fishnet_final$dist_km/1000,xlab="km to border North America")

### Drop points that have a distance value >0 (outside the 50km buffer...)
fishnet_final_clip  <- fishnet_final[!fishnet_final$dist_km > 0,] # 600 points are dropped

### Map in LAEA
x11(); plot(spTransform(NA_shp,dist_proj))
points(spTransform(fishnet_final_clip,dist_proj), col = "blue")


### Map in WGS
x11(); plot(NA_shp)
points(fishnet_final_clip, col = "blue")

### Make fishnet (clipped to NA border) a df for later on...
fishnet_NA_df <- as.data.frame(fishnet_final_clip)
fishnet_NA_df <- fishnet_NA_df [-6] # drop column for used for clipping
names(fishnet_NA_df)[1:2] <- c("lon", "lat")
head(fishnet_NA_df)

save(fishnet_NA_df, file="fishnet2/fishnet_species_NA_clipped_lon_lat_all_years.RData")
# fishnet2 data ready



###------- GeoLocate ------
geolocated <- read.csv("GeoLocate/upload_georef_desktop_added_species_date.csv", h=T)
geolocated <- subset(geolocated, select=c(species, longitude, latitude, year, precision, uncertainty))

### Check names with taxize
names <- as.character(geolocated$species)
names <- names[!duplicated(names)]
length(unique(names))


### Check names
library(taxize)
sources <- gnr_datasources()
knitr::kable(sources)
eol <- sources$id[sources$title == 'EOL']
cof <- sources$id[sources$title == 'Catalog of Fishes']

results <- list()
for (i in names) {
### Paste loop progress
cat("Searching name for", i, "\n")
### Retrieve best match from all data bases
tmp <- try(gnr_resolve(names = i,
                      resolve_once = T,
                      stripauthority = T,
                      best_match_only = T, 
                      data_source_ids=eol), #cof
      silent=T)
### Save only if no error
if ( length(tmp) > 0  ) { results[[i]] <- tmp }
# if ( grepl("Error",tmp) == FALSE) { results[[i]] <- tmp }

}
### Export to fix small issues (remove remaining authority etc.)
checked_names <- do.call(rbind, results)
write.csv(checked_names, "GeoLocate/geolocate_EOL.csv", row.names = F, quote = F)
geolocate_new <- read.csv("GeoLocate/geolocate_EOL_manually_corrected.csv", h=T)


### Merge the names back to raw data
geolocated <- merge(geolocated, geolocate_new, by.x="species", by.y="submitted_name")
geolocated <- geolocated[!duplicated(geolocated), ]
save(geolocated, file="GeoLocate/geolocated_lat_long_all_species_years.RData") # ready, use "matched_name"







### Snap to stream network -15km distance
### Download the Java-Tool from phycoweb.net
download.file("http://www.phycoweb.net/software/rasterGIS/moveCoordinatesToClosestDataPixel103.jar", 
paste(getwd(), "snap_points_to_streams_lakes/moveCoordinatesToClosestDataPixel103.jar", sep="/"), mode = "wb")


### Load and crop stream network again (=latest version)
mask <- raster("/lustre/scratch/client/fas/sbsc/fw_fish/layers_NA/dem_range.tif")
extent(mask)

### Write as ascii
# open .asc in Notepad++ and change "-9999.000000000000000" to "-9999"
writeRaster(mask, "snap_points_to_streams_lakes/raster_mask.asc", NAflag=-9999, overwrite=TRUE)



### Run Java tool: You may need to set the "path" variable in the system settings, 
### see https://www.java.com/en/download/help/path.xml

# system("cmd /c  java -version") # check if Java is installed.  
# system("cmd /c  java -jar moveCoordinatesToClosestDataPixel103.jar") # see options and flags

#    -i   input coordinates file (csv)
#    -r   raster used to determine which pixels have data (esri ascii format)
#    -o   output coordinates file (csv)
# 
# optional parameters
#    -md  maximum distance that new coordinates are allowed to be from original coordinates (in km)

###--- Snapping tolerance of 1 km ----
### Mask 1
system("cmd /c  java -jar /lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/snap_points_to_streams_lakes/moveCoordinatesToClosestDataPixel103.jar  -i /lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/snap_points_to_streams_lakes/points_for_snap.csv   -r  /lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/snap_points_to_streams_lakes/raster_mask.asc    -o  /lustre/scratch/client/fas/sbsc/fw_fish/tables/species_download/snap_points_to_streams_lakes/points_snapped.csv  -md 1")

### Very slow, though there are virtually no points outside the 50km buffer if North America.. 
### Run on Litoria, ~30min for 2Mio points
mkdir snap_points
cd snap_points
unzip points_for_snap.zip

java -jar moveCoordinatesToClosestDataPixel103.jar  -i points_for_snap.csv   -r  raster_mask.asc    -o  points_snapped_1km.csv  -md 1

### Read the raw points
na_fish <- read.csv("snap_points_to_streams_lakes/points_for_snap.csv", h=T)
### Read snapped points
na_fish_1km <- read.csv("snap_points_to_streams_lakes/points_snapped_1km.csv", h=T)

### Check number of species in each snapped set
length(unique(na_fish$species)) 
length(unique(na_fish_1km$species))


### Remove the old coordinate columns
pts_snapped <- subset(pts_snapped, select=-c(old_longitude, old_latitude))


### Which points were removed?
"%ni%" <- Negate("%in%") # create a "not in" -function
rows_removed <- which(na_fish$id %ni% na_fish_1km$id) # get those ID's that were not moved to the stream grids
pts_removed_1km <- na_fish[rows_removed,] # subset the raw SpatialPointsDataFrame


### get the percentage of retained points
kable(data.frame(proportion_retained=c(
"1 km" = 1 - nrow(pts_removed_1km) / nrow(na_fish),

x11()
plot(c(
1 - nrow(pts_removed_1km) / nrow(na_fish),
ylim=c(0, 1), ylab="% retained", xlab="distance set", pch=16)


### Add the removed points again to dataframe. Due to a bug in the Java-tool, some points are discarded even when they lie on a stream grid cell. 
# using raster::extract(), these points can be used anyway (and the points falling outside the stream network will then be finally discarded).

na_fish_1km <- subset(na_fish_1km, select=-c(old_longitude, old_latitude))

### Re-order columns
source("/lustre/scratch/client/fas/sbsc/fw_fish/move_column.R")
pts_removed_1km <- pts_removed_1km[moveme(names(pts_removed_1km), "longitude after id; latitude last")]
na_fish_1km <- rbind(na_fish_1km, pts_removed_1km)

### Remove duplicates
na_fish_1km <- na_fish_1km[!duplicated(subset(na_fish_1km, select=-c(syn_name))) ,]

### Filter by year: allow 16 years in both directions?
na_fish_1km$year_num <- as.numeric(as.character(na_fish_1km$year))


### Allow +-16 years to include as much data as possible
na_fish_1km_year <- subset(na_fish_1km, na_fish_1km$year_num >=1934 & na_fish_1km$year_num <= 2016 )


sort(unique(na_fish_1km$year))
x11(); hist(na_fish_1km_year$year_num, main="1 km")


### Check how many species are in the datasets
length(unique(na_fish_1km_year$species)) # 847 in each


### Export the data
save(na_fish_1km_year, 
     file="snapping_1934-2016.RData")
load("all_snapping_1934-2016.RData")


### This tables include:
# - points snapped to stream network and lakes (latest version as on earthenv.org)
# - non-snapped points added again, as can be used --> bug in Java-tool if point in the middle of cell?
# - species and synonyms merged
# - "raw_name_in_shp" column matches the expert range map polygons


### Merge the shapefiles from both sources using the masterlist
fish_poly <- readShapePoly("rangemaps/na_fish.shp")
proj4string(fish_poly) <- "+proj=longlat +ellps=WGS84"


### Get the species names of NatureServe
add_species <- unique(subset(masterlist, rangemap_source=="NatureServe")$raw_name_in_shp)
add_species <- gsub(" ", "_", add_species) # add underscore

### Read these shapefiles and merge into one file
filenames <- list.files("NatureServe/new_July2016/Shapefiles", pattern=".shp$", full.names=TRUE)
### Get only the paths of the additional species
tmp <- gsub("NatureServe/new_July2016/Shapefiles/", "",  filenames)
tmp <- gsub("\\.shp", "", tmp)
ids <- which(tmp %in% add_species)
filenames <- filenames[ids]

### Read files
# read in all shps, and prepend shapefile name to IDs
natureserve_shapes <- lapply(filenames, function(x){
  shp <- readShapePoly(x)
  shp <- spChFIDs(shp, paste0(x, '_', sapply(slot(shp, "polygons"), slot, "ID")))
  shp$latin <-  gsub("\\.shp", "", gsub("NatureServe/new_July2016/Shapefiles/", "", x))
  shp
})

# rbind to a single object
natureserve_shapes <- do.call(rbind, as.list(natureserve_shapes))
proj4string(natureserve_shapes) <- "+proj=longlat +ellps=WGS84"

### Crop Hawaii
NA_shp_buff_50km <- readShapePoly("/lustre/scratch/client/fas/sbsc/fw_fish/shapes/basins_NA/NA_dissolved/na_dissolved_buff_50km.shp", verbose=T, proj4string=CRS("+proj=longlat +datum=WGS84"))
 


x11()
plot(NA_shp_buff_50km)
plot(natureserve_shapes, add=T, col="blue")
crop_extent_USA<- drawExtent()
natureserve_shapes <- crop(natureserve_shapes, crop_extent_USA)

x11()
plot(natureserve_shapes)

### access through natureserve_shapes$species
writeOGR(natureserve_shapes, "rangemaps/natureserve_additional_shapes.shp", driver="ESRI Shapefile", layer="natureserve_additional_shapes.shp")
# cat(showWKT(proj4string(natureserve_shapes)),file="rangemaps/natureserve_additional_shapes.prj") 


### Merge both polygon files
fish_poly$source <- "PageBurr"
natureserve_shapes$source <- "NatureServe"

tmp <- subset(fish_poly, select=c(latin, source))
tmp2 <- subset(natureserve_shapes, select=c(latin, source))

fish_poly_pageburr_natureserv <- rbind(tmp, tmp2)
fish_poly_pageburr_natureserv$latin <- gsub(" ", "_", fish_poly_pageburr_natureserv$latin) 
writeOGR(fish_poly_pageburr_natureserv, "rangemaps/fish_poly_pageburr_natureserv.shp", driver="ESRI Shapefile", layer="fish_poly_pageburr_natureserv.shp")

# spp <- gsub(" ", "_", unique(fish_poly_pageburr_natureserv$latin)) # 873 species
write.table(spp, "species_sequence_873.txt", row.names = F, col.names = F, quote = F)

### Prepare data for each species on Omega
- shapefiles "genus_species"
- imported as "genus_species" (args[6]), subsetting shapefile
- subsetting point data using "raw_name_in_shp"
- 4 species in Hawai were dropped
- total of 822 species have >4 records

mean(gsub(" ", "_", masterlist$raw_name_in_shp) %in% spp)
masterlist[which(gsub(" ", "_", masterlist$raw_name_in_shp) %ni% spp),]
# spp[which(spp %ni% gsub(" ", "_", masterlist$raw_name_in_shp))]
# 4 species missing were dropped wjen removing Hawaii --> all 4 species occur only here 


### Check how many unique occurrences in each distance data set
r <- raster("snap_points_to_streams_lakes/raster_mask.asc")


coordinates(na_fish_1km_year) <- c("longitude", "latitude")

# rmp <- na_fish_1km_year[1:10,]
# pts_1km <- extract(r, rmp, df=T, small=T, cellnumbers=T)

pts_1km <- extract(r, na_fish_1km_year, sp=T, small=T, cellnumbers=T)

### Check unique cellnumbers for each species
length(unique(na_fish_1km_year$raw_name_in_shp)) # 847 species
### 1km snapping
pts_1km_agg <- plyr::count(pts_1km@data, c('raw_name_in_shp','cells'))
pts_1km_agg <- aggregate(pts_1km_agg, by=list(pts_1km_agg$raw_name_in_shp), function(x) { length(unique(x))})
nrow(subset(pts_1km_agg, pts_1km_agg$cells>=3)) # 836, 819

pts_1km_agg <- pts_1km[c("raw_name_in_shp", "cells")]
pts_1km_agg <- aggregate(raw_name_in_shp ~ cells, pts_1km_agg, function(x) { length(unique(x))})
head(pts_1km_agg)








###------ Merge the data sets ----- 

load("GBIF/gbif_NA_clipped_lon_lat_all_years.RData")
load("BioData/biodata_merged_species_all_years.RData")
load("National_Fish_Habitat/fishhabitat_species_all_years.R")
load("StreamNet/streamnet_species_all_years.R")
load("EPA_Mercury_Fish/epa_fish_species_all_years.RData")
load("Canada_British_Columbia/fish_bc_species_all_years.RData")
load("fishnet2/fishnet_species_NA_clipped_lon_lat_all_years.RData")



### Include the data source as an attribute
gbif_NA_df$data_source <- "GBIF"
biodata_merged$data_source <- "BioData"
fishhabitat_species$data_source <- "fishhabitat"
streamnet_species$data_source <- "StreamNet"
epa_fish$data_source <- "EPA_Mercury"
fish_bc_species$data_source <- "FISS_BC"
# fish_ont_final$data_source <- "Ontario_sport"
fishnet_NA_df$data_source <- "fishnet2"


### Correct fishnet data
names(fishnet_NA_df)[1] <- "species"
names(fishnet_NA_df)[10] <- "lon"
names(fishnet_NA_df)[9] <- "lat"
names(fishnet_NA_df)[21] <- "year"
fishnet_NA_df$abundance <- NA
### Toss columns
fishnet_NA_df <- subset(fishnet_NA_df, select=c(species, lon, lat, year, abundance, data_source))



###Convert "abundance" as.numeric (for rbindlist...)
gbif_NA_df$abundance <- NA
gbif_NA_df$abundance <- as.numeric(gbif_NA_df$abundance)
biodata_merged$abundance <- as.numeric(biodata_merged$abundance)
fishhabitat_species$abundance <- as.numeric(fishhabitat_species$abundance)
streamnet_species$abundance <- as.numeric(streamnet_species$abundance)
epa_fish$abundance <- as.numeric(epa_fish$abundance)
fish_bc_species$abundance <- as.numeric(fish_bc_species$abundance)
# fish_ont_final$abundance <- as.numeric(fish_ont_final$abundance)
fishnet_NA_df$abundance <- as.numeric(fishnet_NA_df$abundance)


library(data.table)
### Merge files, do not use rbind! rbindlist takes only few seconds and doesn't crash!
all_fish_NA <- rbindlist(list(gbif_NA_df,       # GBIF 
                     biodata_merged,            # BioData USGS
                     fishhabitat_species,       # FishHabitat
                     streamnet_species,         # StreamNet North-West US
                     epa_fish,                  # EPA Mercury data US
                     fish_bc_species,           # Canada British Columbia
                     fishnet_NA_df), 
                     fill=T, use.names=T)       # fishnet2 museum (freshwater) records


all_fish_NA <- as.data.frame(all_fish_NA) # 1,852,425  # 1,649,852 obs., possible duplicates
length(unique(all_fish_NA$species)) # 949 species

all_fish_NA <- merge(all_fish_NA, masterlist, by="species") 
length(unique(all_fish_NA$species))  # 807 species

### Which ones are missing?
"%ni%"  <- Negate("%in%") 
missing <- masterlist$species[which(masterlist$species %ni% all_fish_NA$species)]  
write.table(missing, "species_without_records.txt", row.names = F, col.names = F, quote=F)

save.image("NA_fish_point_processing.RData")

source("/lustre/scratch/client/fas/sbsc/fw_fish/clean_gbif_data.R")

### Export corrected names for visual check
write.csv(spp_names_corrected, "all_NA_fish_spp_names_corrected.csv")

### Load corrected species names...
spp_names_corrected <- read.csv("all_NA_fish_spp_names_corrected_edit_import.csv", h=T)

### .. and use output "spp_names_corrected" to merge to initial df
all_fish_NA_final <- merge(all_fish_NA, spp_names_corrected, by.x="species", by.y="species_raw")

summary(all_fish_NA_final)

### Remove NA's in the corrected species column.. (e.g. due to spp. and only genus information)
all_fish_NA_final <- all_fish_NA_final[!is.na(all_fish_NA_final$species_corrected),]

length(unique(all_fish_NA_final$species_corrected)) # 1539
str(fish_check) # 1094 , more species in list than actally in check-list due to sub-species and groups (cf.)


### Clean data
### Use only corrected names and delete old ones
all_fish_NA_final <- cbind(all_fish_NA_final[9], all_fish_NA_final[2:8])
colnames(all_fish_NA_final)[1] <- "species"


### Remove NAs
### Delete all entries without a year
all_fish_NA_final <- all_fish_NA_final[!is.na(all_fish_NA_final$year),]
### ..without species name
all_fish_NA_final <- all_fish_NA_final[!is.na(all_fish_NA_final$species),]
### ..withoyt coordinates
all_fish_NA_final <- all_fish_NA_final[!is.na(all_fish_NA_final$lon),]
length(unique(all_fish_NA_final$species)) # 1527 species
dim(all_fish_NA_final) # 1.669.613 records 

### Remove duplicates
all_fish_NA_final_dupl_del <- all_fish_NA_final[!duplicated(all_fish_NA_final),]
dim(all_fish_NA_final_dupl_del) # 1.316.598 records


### Set filter for the sampling year =  1950 (start of worldclim...)
all_fish_NA_final_dupl_del_year <- subset(all_fish_NA_final_dupl_del, all_fish_NA_final_dupl_del$year >= 1950)
dim(all_fish_NA_final_dupl_del_year) # 1.267.239

### Export
save(all_fish_NA_final_dupl_del_year, file = "all_fish_NA_final_dupl_del_since_1950.RData")

### Make spatial object
all_fish_NA_sp <- all_fish_NA_final_dupl_del_year
coordinates(all_fish_NA_sp) <- c("longitude", "latitude")
proj4string(all_fish_NA_sp) <- "+proj=longlat +datum=WGS84"

### Plot
x11(); plot(NA_shp); 
points(subset(all_fish_NA_sp, data_source == "GBIF"), col = "blue", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "fishhabitat"), col = "green2", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "BioData"), col = "red", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "StreamNet"), col = "gray", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "EPA_Mercury"), col = "yellow", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "FISS_BC"), col = "black", cex = 0.1)
# points(subset(all_fish_NA_sp, data_source == "Ontario_sport"), col = "purple", cex = 0.1)
points(subset(all_fish_NA_sp, data_source == "fishnet2"), col = "brown", cex = 0.1)

### Add legend
legend(-169,44, # coordinates, get via locator(1)
       title="Fish records", 
       cex=1, pch=16, # size of text and symbol code
       col=c("blue", # GBIF
             "red", # fishhabitat
             "green2", # BioData
             "gray",  # StreamNet
             "yellow", # EPA_Mercury
             "black", # FISS_BC
#               "purple", # Ontario
             "purple"), # fishnet2
       legend=c("GBIF", 
                "fishhabitat", 
                "BioData", 
                "StreamNet", 
                "EPA_Mercury",
                "FISS_BC",
#                 "Ontario_Sport", # text in same sequence
                "fishnet2"),
       
       ncol=1,  # number of columns?
       border = F, # show border of symbols?
       bty="n") # turn off putting a box around legend

