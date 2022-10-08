

cat("------ Initiate -----\n")
library(raster)
library(ggplot2)
library(rasterVis)
library(bossMaps)
library(magrittr)
library(tidyr)
library(dplyr)
library(geosphere)
library(devtools)
library(knitr)
library(doParallel)
library(arm)
library(BIOMOD)
library(dismo)
library(SDMTools)



### Set variables for the cluster:
if(Sys.info()[["sysname"]]=="Linux") path="/lustre/scratch/client/fas/sbsc/na_fish" 
n_cores=2
setwd(path)


if(Sys.info()[["sysname"]]=="Linux") {
args <- commandArgs()
print(args)
spp <- args[6]
# spp <- gsub("_", " ", args[6]) # get original species name without underscore..
}


### Create folder for species
path_spp <- paste0(path, "/species_folders/", spp)
dir.create(path_spp)

### Create and set the temporary folder for the raster-package (will be deleted after each species) 
dir.create(paste0(path_spp, "/R_temp/"))

### Delete the files in the R-tmp folder
unlink(paste0(path_spp, "/R_temp/*.gr*"), recursive=T, force = T)

### Create directories
dir.create(paste0(path_spp,"/spatial_priors/priors/"))
dir.create(paste0(path_spp,"/spatial_priors/maps/"))
# setwd(paste0(path_spp, "/spatial_priors"))


### Set up the cluster
cl <- makePSOCKcluster(n_cores, outfile="")
registerDoParallel(cl)
getDoParWorkers() 


### Set raster-package options depending on the platform
rasterOptions(tmpdir=paste0(path_spp, "/R_temp/"), progress="text")
# rasterOptions(chunksize=1000,maxmemory=1000, tmpdir=paste0(path_spp, "/R_temp/"), progress="text")
if(Sys.info()[["sysname"]]=="Linux") rasterOptions(tmpdir=paste0(path_spp, "/R_temp/"))

### Create custom functions
col <- rev(rainbow(100, start = 0, end = 0.7)) 

### Edited version or rangeOffset-function, without rmRaster()
rangeOffset_edit <- function (rdist, dists = NULL, parms, returnRaster = T, normalize = T, 
  verbose = T, writeRaster = F, writeMetadata = T, ...) 
{
  require(raster)
  if (is.null(names(parms))) 
    error(paste0("parms must be a named vector"))
  if (is.null(dists)) 
    dists = freq(rdist, digits = 0, useNA = "no")
  dists = na.omit(dists)
  if (verbose) 
    writeLines(paste0("Running Optimization"))
  rp = sum(dists[, "count"][dists[, "value"] <= 0])/sum(dists[, 
    "count"])
  upper1 = as.numeric(parms["prob"]/rp)
  lower1 = as.numeric((1 - parms["prob"])/(1 - rp))
  parms = c(parms, upper = upper1)
  fixparms = unlist(c(parms, buffer = 0))
  out.lower1 = nlm(f = minimize.fn, p = c(lower = lower1), 
    pname = "lower", fixparms = fixparms, dists = dists, 
    gradtol = 1e-10, iterlim = 200)
  fixparms = unlist(c(parms, lower = out.lower1$estimate))
  out.out = optimize(f = minimize.fn, interval = c(0, max(dists[, 
    "value"], na.rm = T)), pname = "buffer", fixparms = fixparms, 
    dists = dists)
  pout = c(fixparms, fitbuffer = round(out.out$minimum))
  pout["pinside"] = pinside(dists, pout)
  pout["pinside_fitbuffer"] = pinside(dists = dists, parms = c(pout, 
    buffer = as.numeric(pout["fitbuffer"])))
  if ((as.numeric(parms["prob"]) - as.numeric(pout["pinside"])) < 
    -0.01) {
    warning(paste0("You asked for ", 100 * pout["prob"], 
      "% inside the range, but the parameters", " in combination with your expert range and domain resulted in ", 
      round(100 * pout["pinside"], 2), "%. If the fitted % inside the range is less than the requested value,", 
      " you're telling the algorithm that the expert", 
      " map is very accurate (a relatively high % inside) but also that there's a slow decay in", 
      " relative occurrence rate outside the range, which isn't consistent with the", 
      " expert map being very accurate. Try adjusting the the prior probability inside", 
      " the expert map (make it lower), decay rate (make it faster), or skew (make it higher) if you are unhappy about this...", 
      " However, moving the boundary outwards by ", pout["fitbuffer"], 
      "km can produce ", " a probability 'inside' of ", 
      round(100 * pout["pinside_fitbuffer"], 2), "%.", 
      "If the fitted % inside the range is greater than the requested value, try a lower value of rate.", 
      "Due to the range geometry, a slow decay rate is inconsistent with a % inside the expert map."))
  }
  if (!returnRaster) 
    return(as.data.frame(t(pout)))
  if (returnRaster) {
    if (verbose) 
      writeLines(paste0("Calculating the prior over the domain"))
    prior = calc(rdist, function(x) logistic(x, parms = pout))
    if (verbose) 
      writeLines(paste0("Normalizing Prior"))
    if (normalize) {
      nprior = normalize(prior)
    }
    if (!normalize) 
      nprior = prior
    pmeta = c(pout, prange = rp, ncells = ncell(nprior), 
      nNAcells = cellStats(is.na(nprior), sum), normalize = normalize)
    m <- list(parms = pmeta, pnames = names(pmeta))
    metadata(nprior) <- m
    # lapply(list(prior), function(x) rmRaster(x)) # not working
    if (verbose) 
      writeLines(paste0("Saving normalized prior with the following parameters:\n", 
        paste(names(pout), "=", format(unlist(pout), 
          digits = 2), collapse = "\n")))
    if (writeRaster) {
      result = writeRaster(nprior, ...)
      metadata(result) = m
      # rmRaster(nprior)
      if (writeMetadata) {
        write.csv(m, file = paste0(result@file@name, 
          ".csv"))
      }
      return(result)
    }
    if (!writeRaster) 
      return(nprior)
  }
}


### Calculate sensitivity, specificity and cut-off (sens=spec)
get_eval_metrics <- function(r) 
{
  temp <- as.data.frame(raster::extract(r, vdata, sp=T))[1:2]
  temp <- temp[!is.na(rev(temp)[1]),]
  temp <- as.data.frame(cbind(seq(1:nrow(temp)), temp$presences, rev(temp)[1])) # add ID
  temp[2] <- as.numeric(as.character(temp[,2])) # make obs. numeric
  # temp[3] <- round((temp[,3]*1000),5) # scale between 0 and 1000
  # temp[3] <- temp[,3]*1000 # scale between 0 and 1000
  colnames(temp) <- c("ID", "obs", "fit")
  # dev <- calc.deviance(obs=temp$obs, pred=temp$fit, family="binomial", calc.mean = TRUE)
  # return(optimal.cutpoints(X="fit", status="obs",tag.healthy=1, methods=c("ValueSe", "ValueSp", "SpEqualSe", "MaxSpSe"), data=temp, direction = c(">"))) ; #return(dev); #rm(temp)
  # use only Youden --> = TSS
  cut_1 <- as.data.frame(t(as.data.frame(KappaRepet(temp$obs, temp$fit, TSS=T)))) # old BIOMOD version
  names(cut_1) <- "value"
  cut_1$value2 <- NA
  cut_2 <- as.data.frame(t(as.data.frame(optim.thresh(temp$obs, temp$fit)))) # SDMTools
  if (ncol(cut_2)==1) { cut_2$value2 <- NA } # add 2nd dummy column just in case to keep function running
  names(cut_2) <- c("value", "value2")
  p <- subset(temp, temp$obs>=1)$fit
  a <- subset(temp, temp$obs==0)$fit
  eval <- evaluate(p=p, a=a) # dismo
  cut_3 <- as.data.frame(t(threshold(eval)))
  names(cut_3) <- "value"
  cut_3$value2 <- NA
  out <- rbind(cut_1, cut_2, cut_3)

  # SDMTools:
  # min.occurence.prediction: is the minimum prediction for the occurrence (presence) records
  # mean.occurence.prediction: is the mean prediction for the occurrence (presence) records
  # '10.percent.omission': is the threshold value or range in values that excludes approx. 10 percent of the occurrence records
  # 'sensitivity=specificity': is the threshold value or range in values where sensitivity is equal to sensitivity
  # 'max.sensitivity+specificity' :is the threshold value or range in values that maximizes sensitivity plus specificity
  # maxKappa: is the threshold value or range in values with the maximum Kappa statistic
  # max.prop.correct: is the threshold value or range in values with the maximum proportion of presence and absence records correctly identified
  # min.ROC.plot.distance: is the threshold value or range in values where the ROC curve is closest to point (0,1) (or perfect fit)
  ### dismo:
  # kappa: the threshold at which kappa is highest ("max kappa")
  # spec_sens: the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
  # no_omission: the highest threshold at which there is no omission
  # prevalence: modeled prevalence is closest to observed prevalence
  # equal_sens_spec: equal sensitivity and specificity
  # sensitivty: fixed (specified) sensitivity
  return(list(thres=out, eval=eval)) # eval contains AUC
}




###========================#
###----- Load layers ------
###========================#
cat("----- Load layers ------\n")
load(paste0(path, "/global_layers/layers_global.RData"))
env <- layers_global[[c("dem_range", "hydro_avg_01", "hydro_avg_07", "hydro_avg_12", "hydro_avg_15", "lc_wavg_01", "lc_wavg_02", "lc_wavg_03", "lc_wavg_07", "geo_wsum_20", "geo_wsum_42", "geo_wsum_47", "geo_wsum_48", "geo_wsum_56", "geo_wsum_76", "geo_wsum_80", "soil_avg_02")]]

### Load and prepare environmental data for spatial-prior-models
rdist=raster(paste0(path_spp, "/dist_rmap_corr.tif"))
names(rdist) <- "rangeDist"


### Crop near-global layers to the rdist (=domain= rangemap + 5degrees), mask and scale 
rescale_fun <- function(x) {arm::rescale(x, "full")} # scale by 2 SD, then add binary untransformed layer
senv <- foreach(i=names(env), r = unstack(env), .combine=stack, .packages=c("raster", "arm"), .verbose=T) %dopar% {
  tmp <- crop(r, rdist)
  tmp <- mask(tmp, rdist)
  calc(tmp, rescale_fun)
}

names(senv) <- names(env)

### Add binary layer
senv <- addLayer(senv, mask(crop(layers_global[["lentic_lotic01"]], rdist), rdist))



### Check correlation across North America
### Crop near-global layers to the rdist (=domain= rangemap + 5degrees), mask and scale 
# e <- extent(-145, -52, 5, 60) # North America
# tmp <- foreach(i=names(env), r = unstack(env), .combine=stack, .packages=c("raster")) %dopar% {
#   tmp <- crop(r, e)
# }
# 
# names(tmp) <- names(env)
# cor <- layerStats(tmp,'pearson', na.rm=T)
# save(cor, file=paste0(path,"correlation_results_17variables.RData"))
# corr_matrix=cor$'pearson correlation coefficient'
# write.csv(corr_matrix, paste0(path, "corr_matrix.csv"))


### Add a ID-layer
id_long <- setValues(senv[[1]], cellsFromExtent(senv[[1]], extent(senv[[1]]), expand=FALSE))
id_long <- mask(id_long, senv[[1]])
names(id_long) <- "id_long"
senv <- addLayer(senv, id_long)


### Load the previously prepared data ("prepare_na_fish.R")
load(paste0(path_spp, "/", spp, ".RData"))

### Load point data, presences and absences
points <- na.omit(subset(data, select=c(x, y, presences)))
presences <- subset(points, points$presences>0)
absences <- subset(points, points$presences==0)

# rename presences to points and set projection
points <- presences
coordinates(points) <- c("x", "y")
proj4string(points) <- "+proj=longlat +ellps=WGS84"
coordinates(absences) <- c("x", "y")
proj4string(absences) <- "+proj=longlat +ellps=WGS84"


### Load rangemap
range <- readShapePoly(paste0(path_spp, "/", spp, "_rangemap.shp"))
proj4string(range) <- "+proj=longlat +ellps=WGS84"
### Create list object with points and range
spdata <- list(points, range)
names(spdata) <- c("points", "range")
### Check if any points are outside 1000km the range (should not be!)
# dsp=MOLclean(spdata,dropdist=1000,drop=T)



#' # Calculate Priors
cat("-------- Calculate Priors ------\n")
#' 
#' ## Define curves
## ------------------------------------------------------------------------
#== proportion of presences inside the expert map --------------------
# proportion of presences inside the expert map.
expert.r <- raster(paste0(path_spp, "/rangemap_r_corr_final.tif")) # corrected rangemap
values(expert.r) <- ifelse(values(expert.r) ==1, 1, NA)

# expert.r=rasterize(range,senv) #rasterized expert map
tmp1=raster::extract(expert.r,points)
(expert.accuracy=round(sum(!is.na(tmp1))/length(tmp1),2))
write.table(expert.accuracy, paste0(path_spp,"/spatial_priors/",spp,'_ratio.txt'), row.names=F, col.names=F)

#' ## calculate frequency table of distances
## ------------------------------------------------------------------------
dists=freq(rdist,useNA="no",digits=2)
# knitr::kable(head(dists))

#need to tell users to parallelize
#commented cuz already ran, and it's slow
# pdf(paste0(path_spp,'/spatial_priors/checkRates.pdf'),w=10,h=10)
c=checkRates(rdist, dists = dists, prob = seq(0.1, 1, len = 20),rate = exp(seq(log(0.01), log(100), len = 20)), skew = c(0.2), shift = 0, verbose = T, plot = F)
# dev.off()
# system(paste0('open ',paste0(path_spp,'/spatial_priors/checkRates.pdf')))

probs= c(expert.accuracy,expert.accuracy-.1,expert.accuracy-.2,expert.accuracy+.05)

vars=rbind(expand.grid(
  prob=probs,
  rate=c(.05,.1,.2,.4,1),
  skew=0.2,
  shift=0,
  stringsAsFactors=F),
  c(1e-6,0,1e-6,0))
vars=vars[vars$prob < 1 & vars$prob > 0,]  # toss errors in expert accuracy



#' ## Calculate all the curves
## ------------------------------------------------------------------------
uvars=unique(vars[,c('prob',"rate","skew","shift")])

x=seq(-150,500,len=500)

erd=do.call(rbind,
            lapply(1:nrow(uvars),function(i) {
              y=logistic(x,parms=unlist(c(lower=0,upper=1,uvars[i,])))  
              return(cbind.data.frame(group=i,c(uvars[i,]),x=x,y=y))
            })
)  


#' ## Visualize potential decay parameters
## ------------------------------------------------------------------------
# x11(20,10)
# ggplot(erd, aes(x=x,y=y,linetype=as.factor(skew),colour=as.factor(rate),group=group)) + 
#     geom_vline(aes(0),colour="red",xintercept=1)+
#     geom_line()+
#     ylab("Prior value (not normalized)")+
#     xlab("Distance to range edge (km)")

#' # Process priors
#' ## Calculate the priors
## ----warning=F, results="hide"-------------------------------------------
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
# registerDoParallel(7)


# file.remove(paste0(path_spp,"/spatial_priors/priors/"), list.files(paste0(path_spp,"/spatial_priors/priors/")))
# dir.create(paste0(path_spp,"/spatial_priors/priors/"))
foreach(i=1:nrow(vars),.options.multicore=mcoptions) %dopar% {
  ## calculate the expert range prior
  	fo0=paste0(path_spp,'/spatial_priors/priors')
  	if(!file.exists(fo0)) dir.create(fo0)
    fo=paste0(fo0,'/',spp,"_prior_",paste(vars[i,],collapse="_"),".tif")
    #if(file.exists(fo)) return(NULL)
    expert=rangeOffset_edit(rdist,dists=dists,
                       parms=c(unlist(vars[i,]),upper=1),
                       normalize=T,verbose=T,
                       writeRaster=T,filename=fo,overwrite=T,datatype="FLT4S")
}




#' ## Evaluate curves  
cat("----- Evaluate curves -----\n")
#' Calculate which curve parameter combinations were able to get the requested % inside 
## ------------------------------------------------------------------------
fs=list.files(paste0(path_spp,'/spatial_priors/priors'),pattern="prior.*tif$",full=T,recursive=T)
# fs=fs[grepl(species,fs)&!grepl("old",fs)]

res= foreach(f=fs,.combine=rbind_list,.packages=c("raster", "dplyr")) %dopar% {
    dt=raster(f)
    tmp <- read.csv(paste0(f, ".csv"), h=T)
    mt <- append(list(tmp$parms), list(tmp$pnames))
    names(mt) <- c("parms", "pnames")
    # mt=metadata(dt)
    names(mt$parms)=mt$pnames
    tres=data.frame(species=spp,t(mt$parms),file=f)
    return(tres)
    }


#' ## Stack the priors
## ------------------------------------------------------------------------
priors=stack(fs)

## build prior table from metadata
# priorf=foreach(i=1:nlayers(priors),.combine=rbind_list) %do% {
#   t1=metadata(priors[[i]])
#   t2=t1$parms
#   names(t2)=t1$pnames
#   return(data.frame(id=i,t(t2)))
# }

### Run direclty on the csv-file
priorf=foreach(f=fs, i=1:length(fs), .combine=rbind_list) %do% {
  t1 <- read.csv(paste0(f, ".csv"), h=T)
  t2=t1$parms
  names(t2)=t1$pnames
  return(data.frame(id=i,t(t2)))
}

names(priors)=paste0("prior_",priorf$prob,'_',priorf$rate,'_',priorf$skew)#basename(fs[wp])

  
#' ## Assemble modeling dataset
## ------------------------------------------------------------------------
## build single raster stack of all needed data (env and priors)
# senv <- dropLayer(senv, "id_long")
rdata=stack(senv,priors)




### ---- Subset data 70/30 % for validation -----------
"%ni%" <- Negate("%in%")
points$ID <- seq(1:nrow(points))
points_val <- points[sample.int(nrow(points), round(nrow(points) * 0.3, 0), replace=F), ]
points_fit <- points[points$ID %ni% points_val$ID, ]


## generate presence and non-detection datasets
pres=cbind.data.frame(presence=1,raster::extract(rdata,points_fit,df=T,ID=F))
abs=cbind.data.frame(presence=0,raster::extract(rdata,absences,df=T,ID=F))

dim(pres)
dim(abs)

### Get columns that are potentially entirely NA
drop_vars <- colnames(pres)[colSums(is.na(pres)) > 0]


fdata <- rbind(pres, abs)
fdata <- fdata[colSums(!is.na(fdata)) > 0] # Remove columns 

### Remove absences if in the same cell as presences
fdata <- fdata[which(abs(fdata$presence) == ave(fdata$presence, fdata$id_long, 
                            FUN=function(x) max(abs(x)))), ]

### Add the absences in the validation data set --> use this df, corrected for duplicate entries in cells with presences
### Get coordinates from random sample
fdata_absences <- subset(fdata, presence==0)
tmp <- xyFromCell(rdata[["id_long"]], fdata_absences$id_long, spatial=T)
tmp_data <- data.frame(presences=0, dummy=1:nrow(fdata_absences))[-2]
tmp <- SpatialPointsDataFrame(tmp, tmp_data)
proj4string(tmp) <- proj4string(points_val)

vdata <- rbind(points_val["presences"],tmp) # remove "ID", both spatial points, validate models (30% of presences)
vdata$presences <- ifelse(vdata$presences>=1, 1, 0) # truncate all presences to "1" 

### Remove redundant columns
fdata <- subset(fdata, select=-c(ID,id_long))

### Clean data frames
data <- droplevels(data)
fdata <- droplevels(fdata)
vdata <- droplevels(vdata)


### First remove all variables in the tables that have NA's (no coverage at all)
fdata[ is.na(fdata) ] <- NA
fdata <- fdata[,colSums(is.na(fdata))<nrow(fdata)]

### Remove all variables in the tables that are entirely zero (landcover and geology)
fdata <- fdata[, colSums(abs(fdata)) != 0]


cat("------ Variables names in the model: ------\n")
cat(names(fdata), sep="\n") # print one per line


#' ## Fit models
cat("------ Fit models ------\n")
## ---- warning=F----------------------------------------------------------
fdata$weight=1e-6
best.var=names(senv)
best.var <- intersect(best.var, names(fdata)) # do not use those rasters with zero values

ncell_water <- priorf$ncells[1] - priorf$nNAcells[1] #ncell() includes also NA
write.table(ncell_water, paste0(path_spp, "/spatial_priors/number_cells.txt"), row.names=F, col.names=F)
fdata$weight[fdata$presence == 0] = ncell(ncell_water)/sum(fdata$presence == 0)

## Set up models
formulas=paste("presence/weight ~ offset(log(",grep("prior",colnames(fdata),value=T),'))+', 
               paste0(best.var,collapse='+'))

mods=foreach(f=formulas) %dopar%{
  glm(as.formula(f), family=poisson(link=log),data=fdata,weights=weight, maxit=100)
}


# file.remove(paste0(pred.path,species,'/',list.files(paste0(pred.path,species))))
#' ### Calculate AIC
## ------------------------------------------------------------------------
priorf$AIC=round(unlist(lapply(mods,function(x) AIC(x))))
priorf=data.frame(priorf)
write.csv(priorf, paste0(path_spp, "/spatial_priors/",spp, "_priorf.csv"), row.names=F)


#' ## Make spatial predictions
## ------------------------------------------------------------------------

mcoptions <- list(preschedule=FALSE, set.seed=FALSE)

ptype="response"

psi=1:nrow(priorf)


unlink(paste0(path_spp, "/spatial_priors/maps/*tif*"), recursive=T, force = T)
unlink(paste0(path_spp, "/spatial_priors/maps/*aux*"), recursive=T, force = T)

### IMPORTANT: remove previously omitted variables from the stack
tmp <- subset(rdata, best.var) # get variables only
tmp2 <- subset(rdata, which(grepl("prior", names(rdata)))) # get priors only
rdata <- stack(tmp, tmp2)


cat("----- Predict models ------\n")
# bossMaps::normalize
foreach(i=psi,.options.multicore=mcoptions, .errorhandling="remove", .packages = c("raster", "bossMaps"))%dopar%{
  fo=paste0(path_spp,"/spatial_priors/maps/",spp,"_posterior_",priorf$prob[i],"_", priorf$rate[i], "_",priorf$skew[i],"_",priorf$shift[i])
  tmp_norm <- normalize(predict(rdata,mods[[i]],type=ptype), file=paste0(fo, "_norm.tif"), overwrite=T) # original
}



### Load normalized maps
### Cumulative output
psf_norm=list.files(paste0(path_spp,"/spatial_priors/maps/"),pattern="posterior.*.norm.tif",full=T)
ps_norm=stack(psf_norm)
names(ps_norm)=sub("prior","posterior",names(priors)[psi])





###------------------------------#
###--- Evaluate predictions -----
###------------------------------#


### Get the plain maxent prediction
maxentPred <- ps_norm[["posterior_1e.06_0_1e.06"]]

### Get evaluation metrics: a list with 2 elements
eval_maxent <- get_eval_metrics(maxentPred)
eval_maxent$thres
eval_maxent$eval

eval_maxent$eval@auc # AUC
eval_maxent$thres["TSS",] $value # TSS

save(eval_maxent, file=paste0(path_spp,"/spatial_priors/eval_maxent.RData"))
write.table(names(maxentPred), paste0(path_spp,"/spatial_priors/maxent_model_type.txt"), row.names=F, col.names=F)




###------------------------------------#
###--- Which model has lowest AIC? ----
###------------------------------------#
bestAIC <- priorf[which.min(priorf$AIC),]

### Subset the stack to get the particular prediction
bestPred <- ps_norm[[bestAIC$id]]
write.table(names(bestPred), paste0(path_spp,"/spatial_priors/best_model_type.txt"), row.names=F, col.names=F)

eval_bestPred <- get_eval_metrics(bestPred)
eval_bestPred$thres
eval_bestPred$eval

### AUC and TSS
eval_bestPred$eval@auc # AUC
eval_bestPred$thres["TSS",] $value # TSS

save(eval_bestPred, file=paste0(path_spp,"/spatial_priors/eval_best_model.RData"))




###--------------------------------------------#
###--- Create binary presence-absence maps ----
###--------------------------------------------#


### Maxent output

### Use the best model to create a variety of binary maps
### Extract single cutoffs
cut_TSS <- eval_maxent$thres["CutOff",] $value
cut_min.occurence.pred <- eval_maxent$thres["min.occurence.prediction",] $value
cut_mean.occurence.pred <- eval_maxent$thres["mean.occurence.prediction",] $value
cut_X10.percent.omission <- eval_maxent$thres["X10.percent.omission",] $value
cut_sensitivity_specificity <- eval_maxent$thres["sensitivity=specificity",] $value
cut_max.sensitivity.specificity <- eval_maxent$thres["max.sensitivity.specificity",] $value
cut_maxKappa <- eval_maxent$thres["maxKappa",] $value
cut_max.prop.correct <- eval_maxent$thres["max.prop.correct",] $value
cut_min.ROC.plot.distance <- eval_maxent$thres["min.ROC.plot.distance",] $value
cut_kappa <- eval_maxent$thres["kappa",] $value
cut_spec_sens <- eval_maxent$thres["spec_sens",] $value
cut_no_omission <- eval_maxent$thres["no_omission",] $value
cut_prevalence <- eval_maxent$thres["prevalence",] $value
cut_equal_sens_spec <- eval_maxent$thres["equal_sens_spec",] $value
cut_sensitivity <- eval_maxent$thres["sensitivity",] $value


### Use a subset of thresholds
maxentPred_bin <- stack(
calc(maxentPred, fun=function(x) { ifelse(x >= cut_TSS, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_min.occurence.pred, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_mean.occurence.pred, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_kappa, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_spec_sens, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_prevalence, 1, 0) } ),
calc(maxentPred, fun=function(x) { ifelse(x >= cut_sensitivity, 1, 0) } )
)

### Fix the names
names(maxentPred_bin) <- c("TSS", 
                         "min_occurence_prediction",
                         "mean.occurence_prediction",
                         "max_kappa", 
                         "spec_sens", 
                         "prevalence", 
                         "sensitivity")

### Explanations
# TSS: minimizes the difference between sensitivity and specificity, i.e. maximizes TSS
# min.occurence.prediction: is the minimum prediction for the occurrence (presence) records
# mean.occurence.prediction: is the mean prediction for the occurrence (presence) records
# max.kappa: the threshold at which kappa is highest ("max kappa")
# spec_sens: the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest
# prevalence: modeled prevalence is closest to observed prevalence
# sens: fixed (specified) sensitivity

### Plot all thresholded maps
res2=1e7
pdf(paste0(path_spp,"/spatial_priors/maxent_thresholds_single.pdf"))
plot(maxentPred_bin, maxpixels=res2, col=col, zlim=c(0,1))
dev.off()

### Export single threshold layers
dir.create(paste0(path_spp, "/spatial_priors/binary_maxent"), recursive=T)
writeRaster(maxentPred_bin, paste0(path_spp, "/spatial_priors/binary_maxent/bin.tif"), bylayer=T, suffix=names(maxentPred_bin), datatype="INT1U", overwrite=T)


### Get sum, average and SD
maxentPred_bin_sum <- calc(maxentPred_bin, sum, na.rm=T)


### Plot
# x11(15,5); 
pdf(paste0(path_spp,"/spatial_priors/maxent_thresholds.pdf"))
plot(aggregate(maxentPred_bin_sum, fact=4, fun=max, na.rm=T), maxpixels=res2, col=col, zlim=c(0.1,7), main="sum");  # how many thresholds "predict" same output?
plot(range, border="black", lwd=1, add=T); 
points(points, pch=1, cex=1) 
dev.off()


### Export
writeRaster(maxentPred_bin_sum, paste0(path_spp, "/spatial_priors/binary_maxent/bin_sum.tif"), datatype="INT1U", overwrite=T)



### Use the best model to create a variety of binary maps
### Extract single cutoffs. Overwrite the previous maxent cutoffs
cut_TSS <- eval_bestPred$thres["CutOff",] $value
cut_min.occurence.pred <- eval_bestPred$thres["min.occurence.prediction",] $value
cut_mean.occurence.pred <- eval_bestPred$thres["mean.occurence.prediction",] $value
cut_X10.percent.omission <- eval_bestPred$thres["X10.percent.omission",] $value
cut_sensitivity_specificity <- eval_bestPred$thres["sensitivity=specificity",] $value
cut_max.sensitivity.specificity <- eval_bestPred$thres["max.sensitivity.specificity",] $value
cut_maxKappa <- eval_bestPred$thres["maxKappa",] $value
cut_max.prop.correct <- eval_bestPred$thres["max.prop.correct",] $value
cut_min.ROC.plot.distance <- eval_bestPred$thres["min.ROC.plot.distance",] $value
cut_kappa <- eval_bestPred$thres["kappa",] $value
cut_spec_sens <- eval_bestPred$thres["spec_sens",] $value
cut_no_omission <- eval_bestPred$thres["no_omission",] $value
cut_prevalence <- eval_bestPred$thres["prevalence",] $value
cut_equal_sens_spec <- eval_bestPred$thres["equal_sens_spec",] $value
cut_sensitivity <- eval_bestPred$thres["sensitivity",] $value


### Pick a subset of thresholds. There is no right or wrong, rather use a variety to estimate the spread
bestPred_bin <- stack(
calc(bestPred, fun=function(x) { ifelse(x >= cut_TSS, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_min.occurence.pred, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_mean.occurence.pred, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_kappa, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_spec_sens, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_prevalence, 1, 0) } ),
calc(bestPred, fun=function(x) { ifelse(x >= cut_sensitivity, 1, 0) } )
)

### Fix the names
names(bestPred_bin) <- c("TSS", 
                         "min.occurence.prediction",
                         "mean.occurence.prediction",
                         "max.kappa", 
                         "spec.sens", 
                         "prevalence", 
                         "sensitivity")

### Plot
### Plot all thresholded maps
pdf(paste0(path_spp,"/spatial_priors/best_model_thresholds_single.pdf"))
plot(bestPred_bin, maxpixels=res2, col=col, zlim=c(0,1))
dev.off()


## Export single threshold layers
dir.create(paste0(path_spp, "/spatial_priors/binary_best_model"), recursive=T)
writeRaster(bestPred_bin, paste0(path_spp, "/spatial_priors/binary_best_model/bin.tif"), bylayer=T, suffix=names(bestPred_bin), datatype="INT1U", overwrite=T)


### Get sum, average and SD
bestPred_bin_sum <- calc(bestPred_bin, sum, na.rm=T)


### Plot
# x11(15,5); 
pdf(paste0(path_spp,"/spatial_priors/best_model_thresholds.pdf"))
plot(aggregate(bestPred_bin_sum, fact=4, fun=max, na.rm=T), maxpixels=res2, col=col, zlim=c(0.1,7), main="sum");  # how many thresholds "predict" same output?
plot(range, border="black", lwd=1, add=T); 
points(points, pch=1, cex=1) 
dev.off()


### Export sum
writeRaster(bestPred_bin_sum, paste0(path_spp, "/spatial_priors/binary_best_model/bin_sum.tif"), datatype="INT1U", overwrite=T)








### Plotting....
cat("------ Plot raw outputs ------\n")
### Normalized output
pdf(paste0(path_spp,"/spatial_priors/",spp,'_norm.pdf'),w=10,h=10)
  theme_set(theme_bw())
	res2=1e7
	p_psn=gplot(aggregate(ps_norm, fact=4, fun=mean, na.rm=T),max=res2)+geom_raster(aes(fill=value))
	p_psn$data$id=as.factor(sub("prior|posterior","",p_psn$data$variable))
	p_psn$data$rate=priorf$rate[p_psn$data$id] 
	p_psn$data$prob=priorf$prob[p_psn$data$id]
	p_psn+#facet_wrap(~rate_prob)+
		facet_grid(prob~rate,labeller = label_both)+
		scale_fill_gradientn(colours=hcols(100,bias=.7),trans = "log",
												 name="Relative Occurrence Rate\np(x|Y=1)",na.value="transparent")+
	geom_polygon(aes(x=long,y=lat,group=group),
								 data=fortify(range),
								 fill="transparent",col="black",size=.2)+
	geom_text(aes(xmin(ps_norm)+.5*(xmax(ps_norm)-xmin(ps_norm)),ymin(ps_norm)+ .1*(ymax(ps_norm)-ymin(ps_norm)), 
									label=paste0(
										"AIC=",round(AIC)), 
										#ifelse(round(fitbuffer)==0,"",paste0("\nFitDist=", round(fitbuffer)," km"))
									group=NULL), data=priorf,hjust=1,size=2)+
	ylab("Latitude")+xlab("Longitude")+
			coord_equal() +
	  theme(panel.grid.major = element_blank(), # remove grey plot backgound
     panel.grid.minor = element_blank())
dev.off()





cat("------ Plot maxent model ------\n")
### Plot the model without any priors
pdf(paste0(path_spp,"/spatial_priors/maxent_model.pdf"),w=10,h=10)
# x11() # comment if exporting e.g. as pdf()
theme_set(theme_bw())
res2=1e7
	p_psn=gplot(aggregate(maxentPred, fact=4, fun=max, na.rm=T),max=res2)+geom_raster(aes(fill=value))
	p_psn+scale_fill_gradientn(colours=hcols(100,bias=.7),trans = "log",
												 name="Relative Occurrence Rate\np(x|Y=1)",na.value="transparent")+
	geom_polygon(aes(x=long,y=lat,group=group),
								 data=fortify(range),
								 fill="transparent",col="black",size=1)+
	geom_point(aes(x = x, y = y), 
							 data = as.data.frame(coordinates(points)),
							 col="grey20",size=2,pch=16,lwd=1)+
	ylab("Latitude")+xlab("Longitude")+
			coord_equal()
dev.off()





cat("------ Plot best model ------\n")
### Plot the best model
### Normalized output
pdf(paste0(path_spp,"/spatial_priors/best_model.pdf"),w=10,h=10)
theme_set(theme_bw())
res2=1e7
	p_psn=gplot(aggregate(bestPred, fact=4, fun=max, na.rm=T),max=res2)+geom_raster(aes(fill=value))
	p_psn+scale_fill_gradientn(colours=hcols(100,bias=.7),trans = "log",
												 name="Relative Occurrence Rate\np(x|Y=1)",na.value="transparent")+
	geom_polygon(aes(x=long,y=lat,group=group),
								 data=fortify(range),
								 fill="transparent",col="black",size=1)+
	geom_point(aes(x = x, y = y), 
							 data = as.data.frame(coordinates(points)),
							 col="grey20",size=2,pch=16,lwd=1)+
	ylab("Latitude")+xlab("Longitude")+
			coord_equal()
dev.off()





###----- Clean up -----
### Delete all the .grd-and .gri files
unlink(paste0(path_spp, "/spatial_priors/*.grd*"), recursive=T, force = T)
unlink(paste0(path_spp, "/spatial_priors/*.gri*"), recursive=T, force = T)

### Delete the files in the R-tmp folder
unlink(paste0(path_spp, "/R_temp/*.gr*"), recursive=T, force = T)
stopCluster(cl)



