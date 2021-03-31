library(tidyverse)
library(brms)
get_predictions <- function(x){
  library(brms)
  library(tidyverse)
  #source("./src/phylo_ensamble.R")
  get_closest_species <- function(raw_dat, cov_ranef, list_train, nclose = 3){
    softmax <- function(x) {
      x <- x[!is.na(x)]
      
      exp(x) / sum(exp(x))
    }
    # raw_dat$taxonID[raw_dat$taxonID == "DIVIS"] = "DIVI5"
    # raw_dat$taxonID[raw_dat$taxonID == "QUERCUS"] = "QULA2"
    list_to_predict <- raw_dat$taxonID %>% unique
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
  x = data.frame(x)
  
  colnames(x)[5] = "DTM"
  fit<- readRDS( paste("./models/ch2perc_full.rds", sep=""))
  sp_cov = fit$mod$data2$phylogeny
  mod = fit$mod
  features_scale <- full_model$scale[9:18,]
  responses_scale <- full_model$scale[1:8,]

  train_species <- unique(mod$data$taxonID)
  real_species <- x %>% select(treeID, taxonID)
  #scale predictors
  feat_x = x %>% select(rownames(features_scale))
  for(tr in 1:ncol(feat_x)){
    scaled_x <- scale(feat_x[tr], scale= unlist(features_scale$scale[tr]),
      center = unlist(features_scale$mean[tr]))
    ind = colnames(x) == rownames(features_scale)[tr]
    x[ind] = scaled_x
  }  

  
  x <- get_closest_species(x, cov_ranef = sp_cov, 
                            list_train = train_species)
  red_fia <- predict(mod, newdata = x, allow_new_levels = T, nsamples = 200)
  #saveRDS(red_fia, paste("./outdir/predictions/2020/FIA", unique(x$sep_id),"section.rds",  sep="_"))
  
  predictions_dims <- list()
  for(tr in 1:dim(red_fia)[3]){
    tr_lb <- dimnames(red_fia)[[3]][tr]
    rescaled_tr <- red_fia[,, tr] * unlist(responses_scale$scale[tr]) + 
      unlist(responses_scale$mean[tr])
    rescaled_tr <- exp(rescaled_tr)
    predictions_dims[[tr_lb]] <- rescaled_tr
  }
  predictions_dims <- do.call(cbind.data.frame, predictions_dims) 
  predictions_dims <- cbind(x, predictions_dims) #features
  predictions_dims <- predictions_dims %>% group_by(treeID) %>% 
    summarise_all(list(~if(is.numeric(.)) mean(.) else first(.)))
  
  predictions_dims <- predictions_dims %>% 
    select(-one_of(c("taxonID", rownames(features_scale))))
  
  predictions_dims <- left_join(predictions_dims, real_species, by = "treeID")
  write_csv(predictions_dims, paste("./outdir/2020/full/FIA", unique(x$sep_id),"section.csv",  sep="_"))
  
  return()
}

# final_fia <- readr::read_csv("//blue/ewhite/s.marconi/traits_tradeoffs_East/indir/2020_FIA_features.csv") %>%
#   dplyr::select(-one_of(c("plot", "unitcd", "dbh", "hgt")))
final_fia = data.table::fread("./indir/2020f_features.csv")
taxa = data.table::fread("./indir/2020_FIA_features.csv", 
                         select = c("spid","taxonID", "species_sc"))
IDS = data.frame(as.character(final_fia$spid), stringsAsFactors = F)
colnames(IDS)="spid"
taxa = unique(taxa) %>% data.frame
taxa$spid = as.character(taxa$spid)
IDS = left_join(IDS, taxa)
final_fia = cbind.data.frame(final_fia, IDS[-1])

final_fia$treeID <- seq(1:dim(final_fia)[1])
final_fia <- final_fia[complete.cases(final_fia),]

final_fia$sep_id <- paste(final_fia$statecd, final_fia$countycd, sep="_")

dt <- split.data.frame(final_fia, f = final_fia$sep_id)

library(parallel)
# Calculate the number of cores
no_cores <- 64
# Initiate cluster
cl <- makeCluster(no_cores)
fia_object <- parLapply(cl = cl, dt, get_predictions)
stopCluster(cl)

# predictions_dims <- do.call(rbind.data.frame, fia_object) 
# write_csv(predictions_dims,  "./outdir/FIA_traits_full.csv")
#saveRDS(fia_object, "./outputs/FIA_traits_full.rds")
