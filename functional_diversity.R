


### Calculate fish functional diversity from traits
if (!require("raster")) { install.packages("raster", dependencies = T) ; library(raster)}
if (!require("devtools")) { install.packages("devtools", dependencies = T) ; library(devtools)}
if (!require("ape")) { install.packages("ape", dependencies = T) ; library(ape)}
if (!require("caper")) { install.packages("caper", dependencies = T) ; library(caper)}
if (!require("taxize")) { install.packages("taxize", dependencies = T) ; library(taxize)}
if (!require("rfishbase")) { install.packages("rfishbase", dependencies = T) ; library(rfishbase)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = T) ; library(dplyr)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}

path <- "D:/projects/fw_fish"
setwd(path)
source("D:/codes/fw_fish/Xtree.R")

### Set output path
path_OUT <- paste0(path, "/traits_eco")
dir.create(path_OUT)


### Define functions
options(stringsAsFactors = FALSE)
"%ni%" <- Negate("%in%")

### Merge dataframes with different column names, keep only those that match
rbind.match.columns <- function(input1, input2) {
  n.input1 <- ncol(input1)
  n.input2 <- ncol(input2)
  
  if (n.input2 < n.input1) {
    TF.names <- which(names(input2) %in% names(input1))
    column.names <- names(input2[, TF.names])
  } else {
    TF.names <- which(names(input1) %in% names(input2))
    column.names <- names(input1[, TF.names])
  }
  
  return(rbind(input1[, column.names], input2[, column.names]))
}



### Calculate tree, stores the tree as an object
FDtree_customized <- 
  function (S, w = NA, Distance.method = "gower", ord = c("podani", "metric"), 
            Cluster.method = c(ward = "ward", single = "single", complete = "complete", 
                               UPGMA = "average", UPGMC = "centroid", WPGMC = "median", WPGMA = "mcquitty")) 
  {
    if (is.na(w)[1]) {
      w <- rep(1, ncol(S))
    }
    if (Distance.method == "gower") {
      D <- gowdis(S, w = w, ord = ord) # FD package
    }
    else {
      if (stand.x == TRUE) {
        S2 <- scale(S, center = TRUE, scale = TRUE)
        D <- dist(S2, method = Distance.method) # stats
      }
      else {
        D <- dist(S, method = Distance.method)
      }
    }
    tree <<- hclust(D, method = Cluster.method) # create output object 
    # plot(tree) # sami. commented to avoid automatic plottingf
    xtree <<- Xtree(tree)
    c_distance <- cor(D, cophenetic(tree))
    print(paste("The quality of the dendrogram is", round(c_distance, 
                                                         2)))
    xtree
  }


### Load the "fw_fish_FD_dendro" -function, uses apply() and foreach()
source("D:/codes/github/na_fish_models/na_fish_all/fw_fish_FD_dendro.R")




###---- Read data ----
traits <- read.csv(paste0(path, "/fishtrait_info_merged.csv"), h=T)
species <- read.csv(paste0(path, "/sequence_trait_lookup.csv"), h=T)


###---- Impute the missing values from from genus-level data
traits_firstcols <- subset(traits, select=c(my_species, genus))

traits_filled <- traits
was_filled <- list()

for (i in 1:nrow(traits)) {
  
  ### Subset each row
  tmp <- traits[i,]
  cat("Checking genus level traits for", tmp$my_species,"\n" )
  ### If numeric
  nums <- sapply(tmp, is.numeric)
  num_traits <- tmp[ , nums]
  num_traits <- cbind(traits_firstcols[1,], num_traits)
  # head(num_traits)
  
  ### If logical
  logs <- sapply(tmp, is.logical)
  logi_traits <- tmp[ , logs]
  logi_traits <- cbind(traits_firstcols[1,], logi_traits)
  # head(logi_traits)
  
  ### If categorical
  facts <- sapply(tmp, is.factor)
  facts_traits <- tmp[ , facts]
  facts_traits <- cbind(traits_firstcols[1,], facts_traits)
  # head(facts_traits)
  
  ### Get the genus that shoudl be aggregated
  tmp_genus <- num_traits$genus
  
  ### Subset the trait tabke for  this genus
  traits_num_sub <- subset(traits, genus==tmp_genus)
  traits_logi_sub <- subset(traits, genus==tmp_genus)
  traits_facts_sub <- subset(traits, genus==tmp_genus)
  
  ### Get only numeric, logical of factor columns
  traits_num_sub <- traits_num_sub[names(num_traits)]
  traits_logi_sub <- traits_logi_sub[names(logi_traits)]
  traits_facts_sub <- traits_facts_sub[names(facts_traits)]
  
  
  ### Numeric: Aggregate by the given genus
  num_agg <- aggregate(traits_num_sub, by=list(traits_num_sub$genus), mean, na.rm=TRUE, na.action=NULL)
  num_agg$Group.1 <- as.character(num_agg$Group.1)
  num_agg <- subset(num_agg, select=-c(my_species, genus))
  names(num_agg)[1] <- "genus"
  # str(num_agg)
  
  ### Logical: Aggregagte by genus
  logi_agg <- aggregate(traits_logi_sub, by=list(traits_logi_sub$genus), modal, na.rm=TRUE, na.action=NULL)
  logi_agg$Group.1 <- as.character(logi_agg$Group.1)
  logi_agg <- subset(logi_agg, select=-c(my_species, genus))
  names(logi_agg)[1] <- "genus"
  # str(logi_agg)
  
  ### Factor: Aggregagte by genus
  facts_agg <- aggregate(traits_facts_sub, by=list(traits_facts_sub$genus), modal, na.rm=TRUE, na.action=NULL)
  facts_agg$Group.1 <- as.character(facts_agg$Group.1)
  facts_agg <- subset(facts_agg, select=-c(my_species, genus))
  names(facts_agg)[1] <- "genus"
  # str(facts_agg)
  
  ### Put all together
  genus_traits <- cbind(num_agg, subset(logi_agg, select=-c(genus)), subset(facts_agg, select=-c(genus)))
  
  ### Get the columns that are NA in the origial data, i.e. which one needs to be filled?
  fill_these <- colnames(tmp[-c(1:14)])[colSums(is.na(tmp[-c(1:14)])) > 0] # <-- specific to fishtraits.org!
  
  ### Fill each trait
  if (length(fill_these) > 0) {
    
    for (myname in fill_these)  {
      # print(myname)
      traits_filled[i, myname] <- genus_traits[myname]
    } 
    ### Write extra column in table to keep track
    was_filled[i] <- "filled"
  } else {
    was_filled[i] <- "skipped"
  }

}

tmp <- as.data.frame(do.call(rbind, was_filled))
names(tmp) <- "traits_filled"
traits_filled <- cbind(traits_filled, tmp)

all_traits <- traits_filled
traits <- traits_filled


### Add species that have only modelled distributions, but no traits
species_genus <- subset(species, average01==1)
genus_list <- unique(species_genus$use_traits_from)
genus_list # 24 genera


traits_firstcols <- subset(traits, select=c(my_species, genus))


### If numeric
nums <- sapply(traits, is.numeric)
num_traits <- traits[ , nums]
num_traits <- cbind(traits_firstcols, num_traits)
head(num_traits)


### Aggregate by genus
num_genus <- aggregate(num_traits, by=list(num_traits$genus), mean, na.rm=TRUE, na.action=NULL)
num_genus$Group.1 <- as.character(num_genus$Group.1)
num_genus <- subset(num_genus, select=-c(my_species, genus))
names(num_genus)[1] <- "genus"
str(num_genus)


### If logical
logs <- sapply(traits, is.logical)
logi_traits <- traits[ , logs]
logi_traits <- cbind(traits_firstcols, logi_traits)
head(logi_traits)

### Aggregagte by genus
logi_genus <- aggregate(logi_traits, by=list(logi_traits$genus), modal, na.rm=TRUE, na.action=NULL)
logi_genus$Group.1 <- as.character(logi_genus$Group.1)
logi_genus <- subset(logi_genus, select=-c(my_species, genus))
names(logi_genus)[1] <- "genus"
str(logi_genus)


### Put genus data together
genus_traits <- cbind(num_genus, subset(logi_genus, select=-c(genus)))
names(genus_traits)

### select only those genera that are listed as missing species
genus_traits_sub <- subset(genus_traits, genus %in% genus_list)
genus_list[genus_list %ni% genus_traits_sub$genus] # check, none missing


### Check names in original traits
intersect(names(traits), names(genus_traits))


### Get species with traits
tmp <- subset(species, average01==0)
tmp <- merge(tmp, traits, by.x="use_traits_from", by.y="my_species")

### Merge the averaged traits to those species that rely on genus data
tmp2 <- subset(species, average01==1)
tmp2 <- merge(tmp2, genus_traits, by.x="use_traits_from", by.y="genus")

### Put all together
all_traits <- rbind.match.columns(tmp, tmp2)
dim(all_traits)



###---- Select traits ----
str(all_traits)

### Feeding - each gets 1/11 weight
traits_trophic <- 
  subset(all_traits, select=c(sequence_original, # the modelled species
                              nonfeed, # Adults do not feed [binary]
                              benthic, # Benthic feeder [binary]
                              surwcol, # Surface or water column feeder [binary]
                              algphyto, # Algae or phytoplankton [binary]
                              macvascu, # Any part of macrophytes and vascular plants [binary]
                              detritus, # Detritus or unidentifiable vegetative matter [binary]
                              invlvfsh, # Aquatic and terrestrial invertebrates [binary]
                              fshcrcrb, # Larger fishes, crayfishes, crabs, frogs, etc. [binary]
                              blood, # For parasitic lampreys that feed mainly on blood [binary]
                              eggs, # Eggs of fishes, frogs, etc. [binary]
                              other)) # Other different to above [binary] 



### Reproductive Ecology, 1/25 weight (?)
traits_reproductive <- 
  subset(all_traits, select=c(sequence_original, # the modelled species
                              a_1_1, # Nonguarders; Open substratum spawners; Pelagophils [binary]
                              a_1_2, # Nonguarders; Open substratum spawners; Litho-pelagophils.
                              a_1_3a, # Nonguarders; Open substratum spawners; Lithophils (rock-gravel).
                              a_1_3b, # Nonguarders; Open substratum spawners; Lithophils (gravel-sand).
                              a_1_3c, # Nonguarders; Open substratum spawners; Lithophils (silt-mud).
                              a_1_4, # Nonguarders; Open substratum spawners; Phyto-lithophils.
                              a_1_5, # Nonguarders; Open substratum spawners; Phytophils.
                              a_1_6, # Nonguarders; Open substratum spawners; Psammophils.
                              a_2_3a, # Nonguarders; Brood hiders; Lithophils (rock-gravel).
                              a_2_3b, # Nonguarders; Brood hiders; Lithophils (gravel-sand).
                              a_2_3c, # Nonguarders; Brood hiders; Lithophils (mud).
                              a_2_4a, # Nonguarders; Brood hiders; Speleophils (rock cavity).
                              a_2_4c, # Brood hiders; Speleophils
                              b_1_3a, # Guarders; Substratum choosers; Lithophils.
                              b_1_4, # Guarders; Substratum choosers; Phytophils
                              b_2_2, # Guarders; Nest spawners; Polyphils
                              b_2_3a, # Guarders; Nest spawners; Lithophils (rock-gravel).
                              b_2_3b, # Guarders; Nest spawners; Lithophils (gravel-sand).
                              b_2_4, # Guarders; Nest spawners; Ariadnophils
                              b_2_5, # Guarders; Nest spawners; Phytophils
                              b_2_6, # Guarders; Nest spawners; Psammophils
                              b_2_7a, # Guarders; Nest spawners; Speleophils (rock cavity/roof).
                              b_2_7b, # Guarders; Nest spawners; Speleophils
                              b_2_7c, # Guarders; Nest spawners; Speleophils (cavity generalist
                              c1_3_4_c2_4)) # A lumping of all bearers, substrate-indifferent


### Life fistory, each has weight of 1
traits_life_history <- 
  subset(all_traits, select=c(sequence_original, # the modelled species
                              maxtl, # Maximum total length in cm [continuous]
                              matuage, # age at maturity in years of females [continuous]
                              longevity, # Longevity in years based on life in the wild [continuous]
                              fecundity, # Maximum reported fecundity [count]
                              serial, # Serial or batch spawner [binary]
                              season, # length of the spawning season [continuous]
                              potanadr)) #Potamodromous/anadromous, sign. movement regarding spawning [binary] 
        

### Merge traits and specify weights

# ### eco_behav
# traits_final <- cbind(traits_trophic[1], # only species names
#                       traits_trophic[-1],
#                       traits_reproductive[-1],
#                       traits_life_history[-1])
# trait_weights <- c(
#                    rep(1/11,11), # traits_trophic
#                    rep(1/25,25), # traits_reproductive
#                    rep(1,7))     #  traits_life_history


# ### behav
# traits_final <- cbind(traits_trophic, traits_reproductive[-1])
# trait_weights <- c(rep(1/11, 11), rep(1/25, 25))

### eco
traits_final <- cbind(traits_life_history)
trait_weights <- c(rep(1,7))





### Transform logical values to binary 0-1 (but required as numeric)
traits_final_uniqueness <- traits_final

cols <- sapply(traits_final, is.logical)
traits_final[,cols] <- lapply(traits_final[,cols], as.numeric)
head(traits_final)
str(traits_final)

### Data for uniqueness (requires categorical binary data)
cols <- sapply(traits_final_uniqueness, is.logical)
traits_final_uniqueness[,cols] <- lapply(traits_final_uniqueness[,cols], as.factor)
head(traits_final_uniqueness)
str(traits_final_uniqueness)


### Write species name as row.names
row.names(traits_final) <- as.character(traits_final$sequence_original)
row.names(traits_final_uniqueness) <- as.character(traits_final_uniqueness$sequence_original)



###--- Create functional dedrogram ----

### Set output path
# path_OUT <- paste0(path, "/traits_eco_behav")
# path_OUT <- paste0(path, "/traits_behav")
path_OUT <- paste0(path, "/traits_eco")
write.csv(traits_final, paste0(path_OUT, "/fish_traits_eco.csv"), quote = F, row.names = F )

### Calculate cophenetic scores
ex1 <- FDtree_customized(S = traits_final[-1], w = trait_weights,
              Distance.method = "gower", ord = "podani", Cluster.method = "average")

# cophenetic only life history = 0.89
# cophenetic all = 0.87
# cophenetic only trophic and reproductive = 0.74
# cophenetic trophic, reproductive, life history = 0.87


### Get tree
fw_fish_tree <- as.phylo(tree)

### Save trait data
save(traits_final, fw_fish_tree, file=paste0(path_OUT, "/fish_traits_and_tree_with_weights.RData"))
save(traits_final_uniqueness, file=paste0(path_OUT, "/traits_final_uniqueness.RData"))



if (!require("ape")) { install.packages("ape", dependencies = T) ; library(ape)}
if (!require("Biostrings")) { install.packages("Biostrings", dependencies = T) ; library(Biostrings)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = T) ; library(ggplot2)}
if (!require("tibble")) { install.packages("tibble", dependencies = T) ; library(tibble)}
if (!require("tidyr")) { install.packages("tidyr", dependencies = T) ; library(tidyr)}
if (!require("treeio")) { install.packages("treeio", dependencies = T) ; library(treeio)}
library(ggtree)



### Plot tree with genus as labels
groupInfo <- split(fw_fish_tree$tip.label, gsub("_\\w+", "", fw_fish_tree$tip.label))
fishtree <- groupOTU(fw_fish_tree, groupInfo)
fishtree$tip.label <- gsub("_", " ", fw_fish_tree$tip.label)
# x11(); ggtree(fishtree, aes(color=group), layout='fan') + geom_tiplab(size=0.6, aes(angle=angle)) #circular

x11(); ggtree(fishtree, aes(color=group), layout='fan') + geom_tiplab2(size=0.55)  #+  theme_tree("black")





###================================#
###--- Functional distictness ----
###================================#

path <- "D:/projects/fw_fish"
setwd(path)

### Set output path
# path_OUT <- paste0(path, "/traits_eco_behav")
path_OUT <- paste0(path, "/traits_eco")
# path_OUT <- paste0(path, "/traits_behav")


library(caper)
load(paste0(path_OUT, "/fish_traits_and_tree_with_weights.RData"))
load(paste0(path_OUT, "/traits_final_uniqueness.RData"))

fishCM <- clade.matrix(fw_fish_tree)
fishED <- ed.calc(fishCM) # uses fair proportion
save(fishED, file=paste0(path_OUT, "/fishED.RData"))







path <- "D:/projects/fw_fish"
setwd(path)

species <- read.csv(paste0(path, "/sequence_trait_lookup.csv"), h=T)
head(species)

### Add ED
species <- merge(species, fishED$spp, by.x="sequence_original", by.y="species")

### Clean
head(species)
species <- subset(species, select=-c(othernames))

### Export
write.csv(species, paste0(path_OUT, "/master_table_spp_fam_comm_ED.csv"), row.names=F)



### Add range size for each species
rs <- read.csv(paste0(path, "/rangesize_all_spp.csv"), h=T) 
species <- merge(species, rs, by.x="sequence_original", by.y="species", all.x=T)
head(species)


### Add volume for each species
volume <- read.csv(paste0(path, "/estimated_volume_per_species.csv"), h=T)
species <- merge(species, volume, by.x="sequence_original", by.y="species", all.x=T)
head(species)
str(species)

### Export
write.csv(species, paste0(path_OUT, "/master_table_spp_fam_comm_ED_RS_volume_Aug21.csv"), row.names=F)




###--- Create quadrant plot of range size and ED -----

library(ggplot2)
path <- "D:/projects/fw_fish/traits_from_Aug21"
setwd(path)

### Set output path
# path_OUT <- paste0(path, "/traits_eco_behav")
path_OUT <- paste0(path, "/traits_eco")
# path_OUT <- paste0(path, "/traits_behav")


### Get percentiles
df <- read.csv(paste0(path_OUT, "/master_table_spp_fam_comm_ED_RS_volume_Aug21.csv"))

missing <- read.csv(paste0(path, "/species_outside_domain.csv"))
"%ni%" <- Negate("%in%")
df <- df[df$sequence_original %ni% missing$species  ,]
dim(df)


### Transfrom km2 to m2
df$rs_sprior_m2 <- df$rs_sprior * 1000000

### Get percentiles
ED_10 <- quantile(df$ED, probs=0.9, na.rm=T) # 10% max ED
RS_10 <- quantile(df$rs_sprior_m2, probs=0.1, na.rm=T) # 10% smallest range
VOL_10 <-  quantile(df$volume, probs=0.1, na.rm=T) # 10% smallest volume

### 50% --> position where to put the quadrants / lines in the plot
RS_50 <- quantile(df$rs_sprior_m2, probs=0.5, na.rm=T) # 10% smallest range
VOL_50 <-  quantile(df$volume, probs=0.5, na.rm=T) # 10% smallest volume 



### Create the subsets
### Top 10% smallest range?
df$top10_RS <- ifelse(df$rs_sprior_m2 <= RS_10, 1, 0)
sum(df$top10_RS)

### Top 10% highest ED?
df$top10_ED <- ifelse(df$ED >= ED_10, 1, 0)
sum(df$top10_ED)

### Top 10% smallest volume
df$top10_VOL <- ifelse(df$volume <= VOL_10, 1, 0)
sum(df$top10_VOL)

### Top 10% smallest range AND top 10% highest ED?
df$top10_RS_ED <- ifelse(df$rs_sprior_m2 <= RS_10 &  df$ED >= ED_10, 1, 0)
sum(df$top10_RS_ED)

df$top10_VOL_ED <- ifelse(df$volume <= VOL_10 &  df$ED >= ED_10, 1, 0)
sum(df$top10_VOL_ED)

### Write out
df_sub_RS_ED <- df[df$top10_RS_ED==1, ]
df_sub_VOL_ED <- df[df$top10_VOL_ED==1, ]


summary(df_sub_VOL_ED)


write.csv(df_sub_RS_ED, paste0(path_OUT, "/eco_top10_RS_ED.csv"), row.names = F, quote = F)
write.csv(df_sub_VOL_ED, paste0(path_OUT, "/eco_top10_ED_VOL.csv"), row.names = F, quote = F)

