
###===========================================================#
###---- Load model evaluation metrics for the 769 species ----
###===========================================================#


### Set variables for the cluster:
if(Sys.info()[["sysname"]]=="Linux") path="/mnt/domisch/data/fw_fish"
if(Sys.info()[["sysname"]]=="Linux") n_cores=2

setwd(path)
path=getwd()
library(maptools)
library(raster)
dir.create("/data/domisch/fw_fish/R_temp_delete")
rasterOptions(tmpdir="/data/domisch/fw_fish/R_temp_delete")
S <- read.table("sequence_trait_lookup.txt", h=F)


### Create output lists
AUC_offset <- list()
AUC_glm <- list()

sens_offset <- list()
sens_glm <- list()

TSS_offset <- list()
TSS_glm <- list()


for (i in S[[1]]) {
# remove underscore
# i <- gsub("_", " ", i) # get original species name without underscore..
#   print(i)
  if(file.exists(paste0(getwd(), "/species_folders/", i, "/spatial_priors/eval_best_model.RData"))) {
    
   load(paste0(getwd(), "/species_folders/", i, "/spatial_priors/eval_best_model.RData")) # object "eval_bestPred"
   load(paste0(getwd(), "/species_folders/", i, "/spatial_priors/eval_maxent.RData"))  # object  "eval_glm"

   AUC_offset[[i]] <- eval_bestPred$eval@auc # AUC
   AUC_glm[[i]] <-  eval_maxent$eval@auc
  
   sens_offset[[i]] <- eval_bestPred$thres$value[4] # sensitivity
   sens_glm[[i]] <- eval_maxent$thres$value[4]
   
   TSS_offset[[i]] <- eval_bestPred$thres$value[1] # TSS
   TSS_glm[[i]] <- eval_maxent$thres$value[1] 
   
  ### Show progress
  cat(gsub(" ", "_", i), "done", "\n")


  }
}


### Merge results
AUC_offset <- do.call(rbind, AUC_offset); dim(AUC_offset)
AUC_glm <- do.call(rbind, AUC_glm)

sens_offset <- do.call(rbind, sens_offset)
sens_glm <- do.call(rbind, sens_glm)

TSS_offset <- do.call(rbind, TSS_offset)
TSS_glm <- do.call(rbind, TSS_glm)

full_table_all <- cbind(AUC_offset, AUC_glm, sens_offset, sens_glm, TSS_offset, TSS_glm)
colnames(full_table_all) <- cbind("AUC_offset", "AUC_glm", "sens_offset", "sens_glm", "TSS_offset", "TSS_glm")

full_table_all <- cbind(S, full_table_all)

write.csv(full_table_all, "full_table_eval_metrics769_modelled_within_domain.csv", row.names = F, quote = F)
summary(full_table_all)


### Create summary data frame
summary_all <- data.frame(
           metric = c("AUC_offset", "AUC_glm", 
                     "sens_offset", "sens_glm", 
                     "TSS_offset", "TSS_glm"),
           average=c(mean(AUC_offset), mean(AUC_glm), 
                     mean(sens_offset), mean(sens_glm), 
                     mean(TSS_offset), mean(TSS_glm)),
           
           SD=c(sd(AUC_offset), sd(AUC_glm), 
                     sd(sens_offset), sd(sens_glm), 
                     sd(TSS_offset), sd(TSS_glm))
                         )

library(knitr)
kable(summary_all)

write.table(summary_all, "summary_eval_metrics769.txt", row.names=F, col.names=T, quote=F)


### How often was maxent the best model?
model_type <- list()

for (i in S[[1]]) {
# remove underscore
# i <- gsub("_", " ", i) # get original species name without underscore..
#   print(i)
  if(file.exists(paste0(getwd(), "/species_folders/", i, "/spatial_priors/best_model_type.txt"))) {
    
   best_model_type <- read.table(paste0(getwd(), "/species_folders/", i, "/spatial_priors/best_model_type.txt"), h=F) 
   maxent_model_type <- read.table(paste0(getwd(), "/species_folders/", i, "/spatial_priors/maxent_model_type.txt"), h=F) 

   if (as.character(best_model_type$V1) ==  as.character(maxent_model_type$V1)) {
     model_type[[i]] <- as.character(maxent_model_type$V1)
       ### Show progress
  cat(gsub(" ", "_", i), "has identical models", "\n")
   } #if
  } # if
} # for loop

model_type <- do.call(rbind, model_type); dim(model_type)
write.table(rownames(model_type), "identical_model_type769.txt", row.names=F, col.names=F, quote=F)



### Create boxplots
full_table_all <- read.csv(paste0(path, "/full_table_eval_metrics769_modelled_within_domain.csv"), h=T)

usePackage <- function(p){
  if (!is.element(p, installed.packages()[,1])) install.packages(p, dep = TRUE) 
  library(p, character.only = TRUE)
}


usePackage("reshape2")
auc_melt <- melt(full_table_all[c("AUC_offset", "AUC_glm")])
sens_melt <- melt(full_table_all[c("sens_offset", "sens_glm")])
tss_melt <- melt(full_table_all[c("TSS_offset", "TSS_glm")])


path_FIGURES="/mnt/domisch/data/fw_fish/figures"



usePackage("ggplot2")
### AUC
p <- ggplot(auc_melt, aes(factor(variable), value)) + 
     geom_boxplot(aes(fill = factor(variable)), lwd=3, outlier.size=5) + 
  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     # panel.border = element_blank(),
     axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank(), # remove title background
     strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=30),
           axis.text.x  = element_text(angle=0, vjust=0.5, size=30)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=30),
           axis.text.y  = element_text(angle=0, vjust=0.5, size=30)) +
  xlab("Model type") + ylab("AUC")+                           # axis title
  scale_x_discrete(labels=c("GLM with offset", "GLM")) +     # set tick labels
  theme(axis.ticks = element_line(size = 2)) +                # tick size
  theme(axis.ticks.length = unit(.3, "cm")) +
  theme(legend.position="none")                               # remove legend
# x11(); p

svg(paste0(path_FIGURES, "/auc.svg"))
p
dev.off()


### Sensitivity
p <- ggplot(sens_melt, aes(factor(variable), value)) + 
     geom_boxplot(aes(fill = factor(variable)), lwd=3, outlier.size=5) + 
  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     # panel.border = element_blank(),
     axis.line = element_line(colour = "black"))  +
  theme(strip.background = element_blank(), # remove title background
     strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=30),
           axis.text.x  = element_text(angle=0, vjust=0.5, size=30)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=30),
           axis.text.y  = element_text(angle=0, vjust=0.5, size=30)) +
  xlab("Model type") + ylab("Sensitivity")+                           # axis title
  scale_x_discrete(labels=c("GLM with offset", "GLM")) +     # set tick labels
  theme(axis.ticks = element_line(size = 2)) +                # tick size
  theme(axis.ticks.length = unit(.3, "cm")) +
  theme(legend.position="none")                               # remove legend
# x11(); p
svg(paste0(path_FIGURES, "/sens.svg"))
p
dev.off()


### TSS
p <- ggplot(tss_melt, aes(factor(variable), value)) + 
     geom_boxplot(aes(fill = factor(variable)), lwd=3, outlier.size=5) + 
  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     # panel.border = element_blank(),
     axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank(), # remove title background
     strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=30),
           axis.text.x  = element_text(angle=0, vjust=0.5, size=30)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=30),
           axis.text.y  = element_text(angle=0, vjust=0.5, size=30)) +
  xlab("Model type") + ylab("TSS")+                           # axis title
  scale_x_discrete(labels=c("GLM with offset", "GLM")) +     # set tick labels
  theme(axis.ticks = element_line(size = 2)) +                # tick size
  theme(axis.ticks.length = unit(.3, "cm")) +
  theme(legend.position="none")                               # remove legend
# x11(); p
svg(paste0(path_FIGURES, "/tss.svg"))
p
dev.off()


graphics.off()



### Plot all matrics for sprior in one figure
full_table_all$sens_offset <- full_table_all$sens_offset /100
all_melt <- melt(full_table_all[c("AUC_offset", "TSS_offset", "sens_offset")])

# full_table_all$sens_glm <- full_table_all$sens_glm /100
# all_melt <- melt(full_table_all[c("AUC_glm", "TSS_glm", "sens_glm")])


library(ggplot2)
### all
p <- ggplot(all_melt, aes(factor(variable), value)) + 
  geom_boxplot(aes(fill = factor(variable)), lwd=2, outlier.size=3) +
  # geom_boxplot(aes(fill = '#56B4E9'), lwd=2, outlier.size=3) + 
  
  # geom_jitter(size=0.3) +
  theme(panel.grid.major = element_blank(), # remove grey plot backgound
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(strip.background = element_blank(), # remove title background
        strip.text.x = element_blank()) +  # remove title
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=20),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=20)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=20),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=20)) +
  xlab("Evaluation metric") + ylab("Score")+                           # axis title
  scale_x_discrete(labels=c("AUC", "TSS", "SE")) +     # set tick labels
  theme(axis.ticks = element_line(size = 2)) +                # tick size
  theme(axis.ticks.length = unit(.3, "cm")) +
  theme(legend.position="none")                               # remove legend
p <- p + scale_fill_manual(values=c("#ADD8E6", "#ADD8E6", "#ADD8E6"))
p <- p + expand_limits(y = 0)
x11(); p

path_FIGURES <- paste0(path, "/figures")

dir.create(path_FIGURES)
svg(paste0(path_FIGURES, "/auc_tss_sens.svg"))
plot(p)
dev.off()


svg(paste0(path, "/auc_tss_sens_maxent.svg"))
plot(p)
dev.off()

graphics.off()


