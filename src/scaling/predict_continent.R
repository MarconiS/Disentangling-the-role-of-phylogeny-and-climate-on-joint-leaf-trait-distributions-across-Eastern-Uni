setwd("/Users/sergiomarconi/Documents/GitHub/traits_tradeoffs_East")
library(brms)
library(tidyverse)
#source("./src/phylo_ensamble.R")
#source("./src/utilities.R")


get_predictions <- function(x){
  library(brms)
  library(tidyverse)
  #source("./src/phylo_ensamble.R")
  #source("./src/utilities.R")
  
  bjtdm<-readRDS( paste("./models/ch2perc_gauss_species_noml.rds", sep=""))
  #tryCatch({
  
  #sp_cov = readRDS("./outdir/phylo_mat.rds")
  fit <- bjtdm$mod
  sp_cov =fit$data2$phylogeny #%>% as.matrix %>% data.frame
  scaling_y = bjtdm$scale[c(1:8),]
  
  responses_scale <- bjtdm$scale[1:8,]
  train_species <- unique(fit$data$taxonID)
  real_species <- x %>% select(treeID, taxonID, LAT, LON, statecd, plot, countycd, unitcd)
  x <- x %>% dplyr::filter(taxonID != "NA")
  x <- get_closest_species(x, nclose = 2, cov_ranef = sp_cov,
                           list_train = as.character(train_species))
  rm(bjtdm)
  #"NYSY"    "CRATEG"  "QUERCUS" "SALI"    "LARIX"   "HALES"   "JUGLA"  
  x$taxonID[x$taxonID == "NYSY"] = "NYBI"
  x$taxonID[x$taxonID == "CRATEG"] = "CRMO2"
  x$taxonID[x$taxonID == "QUERCUS"] = "QUBI"
  x$taxonID[x$taxonID == "SALI"] = "NYBI"
  x$taxonID[x$taxonID == "LARIX"] = "LALA"
  x$taxonID[x$taxonID == "HALES"] = "HACA3"
  x$taxonID[x$taxonID == "JUGLA"] = "JUCI"
  
  if(nrow(x)==1){
    x = rbind.data.frame(x, x)
  }
  red_fia <- predict(fit, newdata = x, allow_new_levels = T, nsamples = 200, robust = T)
  #saveRDS(red_fia, paste("./preds/2020_full/FIA", unique(x$sep_id),"section.rds",  sep="_"))
  
  predictions_dims <- list()
  for(tr in 1:dim(red_fia)[3]){
    tr_lb <- dimnames(red_fia)[[3]][tr]
    rescaled_tr <- red_fia[, , tr] * scaling_y$scale[[tr]] + responses_scale$mean[[tr]]
    rescaled_tr <- exp(rescaled_tr)
    predictions_dims[[tr_lb]] <- rescaled_tr
  }
  predictions_dims <- do.call(cbind.data.frame, predictions_dims)
  outlier_filter_vars <- colnames(predictions_dims)
  predictions_dims <- cbind(x, predictions_dims) #features
  
  #outlier_remove
  #predictions_dims[outlier_filter_vars] <- apply(predictions_dims[outlier_filter_vars], 2,  
  #  function(x)remove_outliers(x))
  
  #predictions_dims <- predictions_dims %>% group_by(treeID) %>%
  #  summarise_all(list(~if(is.numeric(.)) median(., na.rm=T) else first(.)))
  mincol <- which(colnames(predictions_dims) == "nitrogenPercent.Estimate")
  predictions_dims  <- predictions_dims %>% group_by(treeID) %>% 
    summarise_at(vars(colnames(predictions_dims)[mincol:ncol(predictions_dims)]), funs(weighted.mean(., w=weight)))
  
  
  predictions_dims <- left_join(real_species, predictions_dims, by = "treeID")
  summary(predictions_dims)
  write_csv(predictions_dims, paste("./outdir/2020/sp/FIA", unique(x$sep_id),"section.csv",  sep="_"))
  #},error=function(e){message(paste(unique(x$sep_id),"resulted in error!"))})
  return()
}


get_closest_species <- function(raw_dat, cov_ranef, list_train, nclose = 1){
  # raw_dat$taxonID[raw_dat$taxonID == "DIVIS"] = "DIVI5"
  # raw_dat$taxonID[raw_dat$taxonID == "QUERCUS"] = "QULA2"
  list_to_predict <- raw_dat %>% unique
  unknown_species <- list_to_predict[!(list_to_predict %in% list_train)]
  sample_cov_phylo <- cov_ranef[rownames(cov_ranef) %in% list_train, 
                                        colnames(cov_ranef) %in% unknown_species] #%>%
  #data.frame
  take_five <- list()
  for(sp in 1:length(raw_dat$taxonID)){
    if(raw_dat$taxonID[sp] %in% colnames(sample_cov_phylo)){
      weights <- tail(sort(sample_cov_phylo[,(raw_dat$taxonID[sp])]),nclose)
      taxID <- rownames(sample_cov_phylo)[which(sample_cov_phylo[,(raw_dat$taxonID[sp])] %in% weights)]
      foo <- data.frame(mefa:::rep.data.frame(raw_dat[sp,], length(taxID)), 
                        weight = sample_cov_phylo[taxID, raw_dat$taxonID[sp]])
      foo$taxonID <- taxID
      foo$weight <- softmax(foo$weight)
      take_five[[sp]] <- foo
    }else{
      take_five[[sp]] <- data.frame(raw_dat[sp,], weight = 1)
    }
  }
  res <- do.call(rbind.data.frame, take_five) 
  return(res)
}


get_data_pca_transforamtion <- function(joint_fia_clim) {
  #PCA by variable
  var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
  climate_pca = NULL
  climate_features = NULL
  for(ii in var){
    climate_pca[[ii]] <- joint_fia_clim %>%
      select(contains(ii)) %>%
      FactoMineR::PCA(ncp=1)
    #get pca transformed data
    climate_features[[ii]] = climate_pca[[ii]]$ind$coord
  }
  climate_features = do.call(cbind.data.frame, climate_features)
  colnames(climate_features) <- var
  climate_features <- data.frame(joint_fia_clim[["identifierID"]], climate_features)
  colnames(climate_features)[1] <- "individualID"
  return(list(climate_pca=climate_pca, climate_features=climate_features))
}


#main()
final_fia <- readr::read_csv("./indir/2020f_FIA_features.csv") 
#colnames(final_fia)[19] <- "DTM"
final_fia$snow_melt = 0
final_fia$treeID <- seq(1:dim(final_fia)[1])
#final_fia <- final_fia[complete.cases(final_fia),]
final_fia$sep_id <- paste(final_fia$statecd, final_fia$countycd, sep="_")
dt <- split.data.frame(final_fia, f = final_fia$sep_id)
for(x in 1:length(dt)){
  get_predictions(dt[[x]])
}
#x=dt[[1]]
# check_missing = list.files("./preds/2020_full_csv/")
# nms = names(dt) 
# nms = paste("FIA_",nms,"_section.csv", sep="")
# nms = nms %in% check_missing  
# #main 1132:1195
# dt = dt[-1132]
# 
# library(parallel)
# # Calculate the number of cores
# no_cores <- 64
# # Initiate cluster
# cl <- makeCluster(no_cores)
# fia_object <- parLapply(cl = cl, dt, get_predictions)
# stopCluster(cl)

