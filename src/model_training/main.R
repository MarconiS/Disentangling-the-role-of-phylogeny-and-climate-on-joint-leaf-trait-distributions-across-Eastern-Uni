#load libraries and functions first
library(tidyverse)
library(brms)
source("./src/model_training/get_climate_pca_transformation.R")
source("./src/model_training/my_bayes_r2.R")
base_data <- readr::read_csv("./indir/final_field_dataset.csv")
corrected_data = readr::read_csv("./indir/pigments.csv") 
corrected_data = corrected_data %>%select(individualID, chlA_perc_dw, chlB_perc_dw, carot_perc_dw)
base_data = inner_join(base_data, corrected_data, by="individualID")
base_data = base_data %>% filter(!is.na(chlA_perc_dw))
base_data= base_data %>%filter(!individualID %in% c("NEON.PLA.D08.LENO.00079", "NEON.PLA.D02.SERC.13196","NEON.PLA.D17.SOAP.08298", "NEON.PLA.D08.TALL.02155"))
base_data = base_data %>% filter(H > 2)

base_data = base_data %>% group_by(individualID) %>% slice(1)
base_data %>% summary

climate_FIA <- readRDS("./indir/2020climate_FIA.rds")
base_data$taxonID <- factor(base_data$taxonID)
#base_data = left_join(base_data, chloro_with_mass, by=c("extractChlAConc", "extractChlBConc", "carotConc_DW")) 
#base_data = base_data %>% filter(!is.na(chlAConc_DW), !is.na(chlBConc_DW), !is.na(carotConc_DW))
set.seed(1987)
train_data <- base_data %>% 
  dplyr::select(individualID, taxonID, siteID) %>%
  group_by(taxonID, siteID) %>%
  sample_frac(0.8)
traits_names <- c("nitrogenPercent",  "carbonPercent",  "ligninPercent",  "cellulosePercent", 
                  "leafMassPerArea", "chlA_perc_dw", "chlB_perc_dw", "carot_perc_dw")
local_environment <- c("ope", "ect", "DTM")


#log transform, plot responses
leaf_traits_data <- base_data[colnames(base_data) %in% traits_names] 
leaf_traits_data = leaf_traits_data %>%
  mutate_all(log)# %>%
leaf_traits_data <-  cbind.data.frame(base_data[["individualID"]], leaf_traits_data)
colnames(leaf_traits_data)[1] <- "individualID"


# load climate trends
if(file.exists("./indir/neon_climate.csv")){
  climate_trends <- readr::read_csv("./indir/neon_climate.csv")
}else{
  # download ckimate from Daymet
  ids = base_data %>% select(individualID, plotLatitude, plotLongitude)
  readr::write_csv(ids, './DMT_retriever/neon_coords.csv')
  
  download_daymet_batch(file_location = './DMT_retriever/neon_coords.csv',
                        start = 1995,
                        end = 2015,
                        internal = F,
                        path = './DMT_retriever/daymet_neon/')
  climate_trends = melt_daymet(path = "./DMT_retriever/daymet_neon/", 
                               outfile = "./indir/neon_climate.csv")
}
climate_trends <- climate_trends %>%
  dplyr::filter(individualID %in% base_data[["individualID"]])
climate <- predict_data_pca_transforamtion(climate_trends, 
                                           climate_FIA$climate_pca)
climate_features = climate$climate_features

#get phylogeny from opentree of life
#append list of species from FIA and the species in the dataset
species_list <- base_data %>%
  dplyr::select(taxonID, scientificName)
sp_tmp <- lapply(1:nrow(base_data),
                 function(x)str_split(base_data$scientificName[x], pattern = " ")[[1]][1:2])
sp_tmp <- do.call(rbind, sp_tmp)
species_list = cbind.data.frame(base_data$taxonID, paste(sp_tmp[,1], sp_tmp[,2]))
colnames(species_list) <- c("taxonID","scientificName")

#remove ACRUR typo in NEON data
species_list = species_list %>% filter(taxonID != "ACRUR")
sp_list2 <- readr::read_csv("./indir/phylogeny/full_list.csv") %>%
  dplyr::select(taxonID, species_sc) 
colnames(sp_list2) <- c("taxonID", "scientificName")
sp_list2 = sp_list2 %>% filter(!(scientificName %in% unique(species_list$scientificName)))
# make sure that each species counted only once, then calculate cophenetic distance 
species_list <- rbind.data.frame(species_list, sp_list2) %>%
  unique 
species_list <-  get_cophenetic_distance(species_list, sc_name = "scientificName")
species_list = readRDS("./indir/phylogeny/2020_phylo.rds")
#create 
species_features <- base_data %>% 
  dplyr::select(individualID, scientificName)
species_features[["scientificName"]] <- word(species_features[["scientificName"]], 1,2)
species_features[["scientificName"]] <- 
  str_replace(species_features[["scientificName"]], " spp.", "")
species_features[["scientificName"]] <- 
  str_replace(species_features[["scientificName"]], " sp.", "")
species_features[["scientificName"]] <- tolower(species_features[["scientificName"]])
species_features <-left_join(species_features, species_list$resolved_names, 
                             by = c("scientificName" = "search_string")) %>%
  select(individualID, taxonID)
colnames(species_features) <- c("individualID", "taxonID")
species_features = species_features[complete.cases(species_features), ]

site_features <- base_data[colnames(base_data) %in% local_environment] 
site_features <-  cbind.data.frame(base_data[["individualID"]], site_features)
colnames(site_features)[1] <- "individualID"





full_dataset <- base_data %>% select(individualID, siteID, taxonID) %>%
  ungroup() %>%
  dplyr::select(-one_of("taxonID")) %>%
  inner_join(leaf_traits_data) %>%
  inner_join(climate_features) %>%
  inner_join(site_features) %>%
  inner_join(species_features) %>%
  unique

#scale angular data using the cosine in order to guarantee that close angle are closer
full_dataset["ect"] = cos(full_dataset["ect"]* 0.0174533)
full_dataset["ope"] = sin(full_dataset["ope"]* 0.0174533)


# clean taxonID data that can't be traced/ typos
full_dataset<- full_dataset[!is.na(full_dataset[["taxonID"]]),]
full_dataset[!full_dataset[["taxonID"]] 
             %in% colnames(species_list$cov_taxa_eff),"taxonID"] %>% 
  unique %>% 
  print
full_dataset[full_dataset[["taxonID"]]=="ABIES","taxonID"] = "ABBA"
full_dataset[full_dataset[["taxonID"]]=="ACSA3","taxonID"] = "ACSA2"
full_dataset[full_dataset[["taxonID"]]=="BETUL","taxonID"] = "BELE"
full_dataset[full_dataset[["taxonID"]]=="BOURR","taxonID"] = "BOSU2"
full_dataset[full_dataset[["taxonID"]]=="FRAXI","taxonID"] = "FRAM2"
full_dataset[full_dataset[["taxonID"]]=="HALES","taxonID"] = "HACA3"
full_dataset[full_dataset[["taxonID"]]=="NYSY","taxonID"] = "NYBI"
full_dataset[full_dataset[["taxonID"]]=="OXYDE","taxonID"] = "OXAR"
full_dataset[full_dataset[["taxonID"]]=="PINUS","taxonID"] = "PIEC2"
full_dataset[full_dataset[["taxonID"]]=="SASSA","taxonID"] = "SAAL5"

#remove desert point (eventually all western sites?)
full_dataset = full_dataset %>% filter(siteID != "JORN")


set.seed(1987)
#train/test split
train_data <- full_dataset %>% 
  filter(individualID %in% train_data[["individualID"]])

train_data = train_data %>% filter(siteID != "MLBS")
test_data <- full_dataset %>% 
  filter(!individualID %in% train_data[["individualID"]])
test_data = test_data %>% filter(siteID != "MLBS")

colnames(train_data)
train_data=train_data[complete.cases(train_data),]
test_data <- test_data[complete.cases(test_data),]

#scale both responses and features
ind <- sapply(train_data, is.numeric)
train_data[ind] <- lapply(train_data[ind], scale)
cls=1:ncol(test_data)
tst_scaled <- lapply(cls[as.vector(ind)], function(x){
  scale(test_data[x], center = attr(train_data[[x]],"scaled:center"),
        scale = attr(train_data[[x]],"scaled:scale"))})
tst_scaled <- do.call(cbind.data.frame, tst_scaled)
test_data <- cbind.data.frame(test_data[!ind], tst_scaled)
phylogeny = as.matrix(species_list$cov_taxa_ef)
set.seed(1987)
phylogeny = Matrix::nearPD(phylogeny, corr = T, keepDiag =T)

#define formula and build the model
#+ (1| gr(taxonID, cov=phylogeny)),
list_formulas <- paste(paste("mvbind("
                             #responses
                       , paste(colnames(leaf_traits_data)[-1], collapse = " , "), ") ~ ", sep = "")
                       #climate variables
                       , paste("s(",colnames(climate_features)[-c(1,2,5)],  ")", collapse = " + ", sep = "")  
                       , " + "
                       , paste("s(",local_environment,")", collapse = " + ", sep = "") 
                       , " + "
                       , "(1| gr(taxonID, cov=phylogeny))"
                       #, collapse = "+"
)
fit <- brm(list_formulas, data = train_data, seed = 1987
           #, set_rescor(F)
           , prior = prior(horseshoe()), 
           , family=gaussian(), cores=2, chains = 2
           , iter = 2000
           , data2 = list(phylogeny = (phylogeny$mat))
)
# bR2 <-  retrotransformed_r2(fit, sd_y = scaling_ys$scale, 
#                             mean_y = scaling_ys$center, newdata = train_data)
#print(bR2)
scaling_ys =  lapply(traits_names, function(x){
  list(center = attr(train_data[[x]],"scaled:center"),
       scale = attr(train_data[[x]],"scaled:scale"))})
scaling_ys = do.call(rbind.data.frame, scaling_ys)
scaling_ys$tr_name = traits_names
bR2_retro <- retrotransformed_r2(fit, sd_y = scaling_ys$scale, 
                                    mean_y = scaling_ys$center, newdata = test_data)
print(bR2_retro)
#get scale of responses and features for predictions
scale = list()
for(ii in colnames(train_data[ind])){
  scale[[ii]] = cbind(ii, attributes(train_data[[ii]])[2], attributes(train_data[[ii]])[3])
}
scale = do.call(rbind.data.frame, scale)
colnames(scale) <- c("trname", "mean", "scale")
parameters_summary <- VarCorr(fit)


#
#this probably goes in another function
#
predictions_full = predict(fit, newdata = test_data)
predictions_dims <- list()
for(tr in 1:dim(predictions_full)[3]){
  tr_lb <- dimnames(predictions_full)[[3]][tr]
  rescaled_tr <- predictions_full[,1, tr] * unlist(scale$scale[tr]) + unlist(scale$mean[tr])
  rescaled_tr <- exp(rescaled_tr)
  predictions_dims[[tr_lb]] <- rescaled_tr
}
predictions_dims <- do.call(cbind.data.frame, predictions_dims) 
#predictions_dims <- cbind(test_data["individualID"], predictions_dims) #features
parameters_summary <- VarCorr(fit)
#get the LES 
obj = list(mod = fit, R2 = bR2_retro, phylo = species_list, 
           params =parameters_summary, scale_ys = scaling_ys, 
           scale = scale, env_pca=list(climate))#, soil))
saveRDS(obj, "./models/perc_gauss_combined.rds")
