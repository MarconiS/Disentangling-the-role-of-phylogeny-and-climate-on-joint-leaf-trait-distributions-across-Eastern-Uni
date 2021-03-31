# create a dataset with external traits data sources
#BIEN
library(reshape2)
library(tidyverse)
BIEN <- readr::read_csv("./indir/external_validation/bien_data.csv")
species_tweaked <- read_csv("indir/external_validation/species_.csv") %>% dplyr::select(SCI_NAME, taxonID)
colnames(species_tweaked)[1] <- "name_matched" 
BIEN = left_join(BIEN, species_tweaked)
BIEN_topo <- readr::read_csv("./indir/external_validation/bien_topography.csv")
BIEN = left_join(BIEN, BIEN_topo)
BIEN = BIEN[complete.cases(BIEN), ]
BIEN_plots = BIEN %>% select(latitude, longitude, aspect, slope, elevation) %>% unique

BIEN_traits = BIEN %>% select(trait_name, trait_value, latitude, longitude, taxonID)
BIEN_traits = BIEN_traits %>% unique
BIEN_traits = dcast(BIEN_traits, latitude + longitude + taxonID ~ trait_name, 
                    value.var="trait_value", fun.aggregate = median)

BIEN = left_join(BIEN_traits, BIEN_plots)
BIEN$individualID = paste("BIEN", 1:88, sep="_")
colnames(BIEN) = c("latitude", "longitude", "taxonID", "LMA",  "N", "P", "ect", "ope", "DTM", "individualID")
BIEN$C = NA
BIEN = BIEN %>% select(individualID, taxonID, latitude, longitude, ect, ope, DTM, N, C,  LMA)
BIEN$LMA = 1/BIEN$LMA *1000

full_dataset = read_csv("./indir/external_validation/external_dataset.csv")


#full_dataset %>% select(individualID, latitude, longitude) %>% write_csv("./indir/external_validation/ext_loc.csv")


#download climate variables
library(daymetr)

# download_daymet_batch(file_location = "./indir/external_validation/ext_loc.csv",
#                       start = 1995,
#                       end = 2015,
#                       internal = F,
#                       path = "./daymet")

ext_clim = melt_daymet(path ="./daymet/")
# load model and transform climate accordingly
clim_pca =  readRDS("./indir/2020climate_FIA.rds")

full_model <- readRDS("./models/ch2perc_gauss_combined.rds")
var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
#make PCA per each feature
climate_features <- predict_data_pca_transforamtion(ext_clim, 
                                                    clim_pca$climate_pca)
climate_features = climate_features$climate_features

full_dataset["ect"] = cos(full_dataset["ect"]* 0.0174533)
full_dataset["ope"] = sin(full_dataset["ope"]* 0.0174533)
full_dataset = inner_join(full_dataset, climate_features)
scaling_x = full_model$scale[9:18,]
scaling_y = full_model$scale[1:8,]

full_dataset = full_dataset %>% select(daylength, prec, rad, snow_melt, tmax, tmin, vp, DTM, ope, ect,
                        individualID, taxonID, latitude, longitude, N, C, LMA)

ext_norm <- sapply(1:10, function(x){
  scale(data.frame(full_dataset)[x], center = scaling_x$mean[x],
        scale = scaling_x$scale[x])})

full_dataset[1:10] = ext_norm

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

ext_dataset = get_closest_species(raw_dat = full_dataset, cov_ranef = full_model$mod$data2$phylogeny, 
                                   list_train = unique(full_model$mod$data$taxonID), nclose = 2)
ext_dataset[ext_dataset$taxonID=="BEPAP", "taxonID"] = "BEPA"
ext_dataset[ext_dataset$taxonID=="CAGL8", "taxonID"] = "CAOV2"
ext_dataset[ext_dataset$taxonID=="CATO", "taxonID"] = "CATO6"
ext_dataset[ext_dataset$taxonID=="NYSY", "taxonID"] = "NYBI"
ext_dataset[ext_dataset$taxonID=="QUGE", "taxonID"] = "QUGE2"
ext_dataset[ext_dataset$taxonID=="QUMO", "taxonID"] = "QUMO4"
ext_dataset = ext_dataset %>% filter(!taxonID %in% c("CASTA", "DIVIS", "GOLQ", "LARIX","MAGNO", "AEPA", "ROPS"))

predictions = predict(full_model$mod, newdata = ext_dataset, nsamples = 200)
no_taxa_predictions = predict(full_model$mod, newdata = full_dataset, allow_new_levels = T,  nsamples = 200)

predictions_dims <- list()
no_relations <- list()

for(tr in 1:dim(predictions)[3]){
  tr_lb <- dimnames(predictions)[[3]][tr]
  rescaled_tr <- predictions[, , tr] * scaling_y$scale[[tr]] + scaling_y$mean[[tr]]
  rescaled_tr <- exp(rescaled_tr)
  predictions_dims[[tr_lb]] <- rescaled_tr
  # no using species relations
  rescaled_tr <- no_taxa_predictions[, , tr] * scaling_y$scale[[tr]] + scaling_y$mean[[tr]]
  rescaled_tr <- exp(rescaled_tr)
  no_relations[[tr_lb]] <- rescaled_tr
  
}
predictions_dims <- do.call(cbind.data.frame, predictions_dims)
no_relations <- do.call(cbind.data.frame, no_relations)

outlier_filter_vars <- colnames(predictions_dims)
predictions_dims <- cbind(ext_dataset, predictions_dims) #features

#outlier_remove
#predictions_dims[outlier_filter_vars] <- apply(predictions_dims[outlier_filter_vars], 2,  
#  function(x)remove_outliers(x))

#predictions_dims <- predictions_dims %>% group_by(treeID) %>%
#  summarise_all(list(~if(is.numeric(.)) median(., na.rm=T) else first(.)))
observations = full_dataset %>% select(individualID, taxonID, N, C, LMA, latitude, longitude)
mincol <- which(colnames(predictions_dims) == "nitrogenPercent.Estimate")
predictions_dims  <- predictions_dims %>% group_by(individualID) %>% 
  summarise_at(vars(colnames(predictions_dims)[mincol:ncol(predictions_dims)]), funs(weighted.mean(., w=weight)))
predictions_dims = inner_join(observations, predictions_dims) %>% unique

predictions_dims = predictions_dims %>% filter(taxonID != "ROPS")
predictions_dims$dataset = substr(predictions_dims$individualID, 1, 4) 
predictions_dims$dataset[!predictions_dims$dataset == "BIEN"]="Dimensions"

no_relations = cbind(observations, no_relations)
#get known and unknown species
predictions_dims$un_known = "unknown"
predictions_dims$un_known[predictions_dims$taxonID %in%  unique(full_model$mod$data$taxonID)] = "known"
no_relations$un_known = "unknown"
no_relations$un_known[no_relations$taxonID %in%  unique(full_model$mod$data$taxonID)] = "known"
no_relations$dataset = "both"
#predictions_dims = predictions_dims[predictions_dims$latitude < 30 | predictions_dims$latitude > 35,]
#no_relations = no_relations[no_relations$latitude < 30 | no_relations$latitude > 35,]

N = predictions_dims %>% select(N, nitrogenPercent.Estimate, nitrogenPercent.Q2.5, nitrogenPercent.Q97.5, 
                                taxonID, dataset, un_known) %>% filter(!is.na(N)) 
  
N_coverage = sum(N$N < N$nitrogenPercent.Q97.5 & N$N > N$nitrogenPercent.Q2.5)/sum(!is.na(N$N))
ggplot(N, aes(x = N, y = nitrogenPercent.Estimate)) +
  geom_point(aes(color=(dataset), shape=factor(un_known)), size=3) + 
  scale_shape_manual(values=c(21, 13),
                    name = "species")+
  geom_abline() + 
  theme_bw() +
  coord_flip()+
  ylim(0,4)+xlim(0.,4)
lm(N~nitrogenPercent.Estimate, data = N) %>% summary
nophylo = N %>% filter(un_known == "unknown")
rmse(N~nitrogenPercent.Estimate, data = nophylo) %>% summary
lm(N~nitrogenPercent.Estimate, data = nophylo) %>% summary
noN_coverage = sum(nophylo$N < nophylo$nitrogenPercent.Q97.5 & nophylo$N > nophylo$nitrogenPercent.Q2.5)/
  sum(!is.na(nophylo$N))


C = predictions_dims %>% select(C, carbonPercent.Estimate, carbonPercent.Q2.5, carbonPercent.Q97.5, 
                                taxonID, dataset, un_known) %>% filter(!is.na(C)) 
C_coverage = sum(C$C <= C$carbonPercent.Q97.5 & C$C >= C$carbonPercent.Q2.5)/sum(!is.na(C$C))
ggplot(C, aes(x = C, y = carbonPercent.Estimate)) +
  geom_point(aes(shape=factor(un_known)), color='cyan3', size=3) + 
  scale_shape_manual(values=c(21, 13),
                     name = "species")+
  geom_abline() + 
  theme_bw() +
  coord_flip()+
  ylim(40,58)+xlim(40,58) +
scale_alpha_manual(values=c(1, 0)) 
lm(C~carbonPercent.Estimate, data = C) %>% summary
nophylo = C %>% filter(un_known == "unknown")
lm(C~carbonPercent.Estimate, data = nophylo) %>% summary
noC_coverage = sum(nophylo$C <= nophylo$carbonPercent.Q97.5 & nophylo$C >= nophylo$carbonPercent.Q2.5)/sum(!is.na(nophylo$C))



LMA = predictions_dims %>% select(LMA, leafMassPerArea.Estimate, leafMassPerArea.Q2.5, leafMassPerArea.Q97.5, 
                                  taxonID, dataset, un_known) %>% filter(!is.na(LMA)) 
LMA_coverage = sum(LMA$LMA <= LMA$leafMassPerArea.Q97.5 & LMA$LMA >= LMA$leafMassPerArea.Q2.5)/sum(!is.na(LMA$LMA))

ggplot(LMA, aes(x = LMA, y = leafMassPerArea.Estimate)) +
  geom_point(aes(color=factor(dataset), shape=factor(un_known)), size=3) + 
  scale_shape_manual(values=c(21, 13),
                     name = "species")+
  geom_abline() + 
  theme_bw() +
  coord_flip()+
  ylim(35,360)+xlim(35,360) +
scale_alpha_manual(values=c(1, 0)) 
lm(LMA~leafMassPerArea.Estimate, data = LMA) %>% summary
nophylo = LMA %>% filter(un_known == "unknown")
lm(LMA~leafMassPerArea.Estimate, data = nophylo) %>% summary
no_coverage = sum(nophylo$LMA <= nophylo$leafMassPerArea.Q97.5 & nophylo$LMA >= nophylo$leafMassPerArea.Q2.5)/
  sum(!is.na(nophylo$LMA))

#summary plot of unknowns
full_md = predictions_dims %>% select(individualID, taxonID, N, C, LMA, 
                                      nitrogenPercent.Estimate, leafMassPerArea.Estimate,carbonPercent.Estimate)

no_corr_md = no_relations %>% select(individualID, taxonID, N, C, LMA, 
                                      nitrogenPercent.Estimate, leafMassPerArea.Estimate,carbonPercent.Estimate)
full_md["model"] = "phylogeny"
no_corr_md["model"] = "null"

phylo_null = rbind.data.frame(full_md, no_corr_md)
phylo_null = phylo_null %>% filter(!taxonID %in% unique(full_model$mod$data$taxonID))
ave_phylo = phylo_null %>% group_by(taxonID, model) %>% summarize_all(mean)
ggplot(ave_phylo, aes(x = C, y = carbonPercent.Estimate, color = model))+
  geom_abline(intercept = 0,slope = 1) + coord_flip()+
  geom_point() + theme_bw() +  stat_ellipse()

ggplot(ave_phylo, aes(x = N, y = nitrogenPercent.Estimate, color = model))+
  geom_abline(intercept = 0,slope = 1) +  coord_flip()+
  geom_point() + theme_bw() +  stat_ellipse()

ggplot(ave_phylo, aes(x = LMA, y = leafMassPerArea.Estimate, color = model))+
  geom_abline(intercept = 0,slope = 1) +  coord_flip()+
  geom_point() + theme_bw() +  stat_ellipse()


RMSE = function(m, o){
  sqrt(mean((m - o)^2, na.rm = T))
}
untrained_full = full_md %>% filter(!taxonID %in% unique(full_model$mod$data$taxonID))
untrained_nophylo = no_corr_md  %>% filter(!taxonID %in% unique(full_model$mod$data$taxonID))

rmses = rbind.data.frame(c(RMSE(untrained_full$nitrogenPercent.Estimate, untrained_full$N) / sd(untrained_full$nitrogenPercent.Estimate), "phylo", "N%"), #
c(RMSE(untrained_nophylo$nitrogenPercent.Estimate, untrained_nophylo$N)/ sd(untrained_full$nitrogenPercent.Estimate), "no", "N%"),#

c(RMSE(untrained_full$carbonPercent.Estimate, untrained_full$C)/ sd(untrained_full$carbonPercent.Estimate),"phylo","C%"), #
c(RMSE(untrained_nophylo$carbonPercent.Estimate, untrained_nophylo$C)/ sd(untrained_full$carbonPercent.Estimate), "no","C%"),# 

c(RMSE(untrained_full$leafMassPerArea.Estimate, untrained_full$LMA)/ sd(untrained_full$leafMassPerArea.Estimate),"phylo", "LMA"), #
c(RMSE(untrained_nophylo$leafMassPerArea.Estimate, untrained_nophylo$LMA)/ sd(untrained_full$leafMassPerArea.Estimate), "no", "LMA")) #

colnames(rmses) = c("RMSE", "model", "Trait")
ggplot(rmses, aes(y=as.numeric(RMSE), fill = model, x = Trait))+geom_bar(stat="identity", position=position_dodge()) + theme_bw()
