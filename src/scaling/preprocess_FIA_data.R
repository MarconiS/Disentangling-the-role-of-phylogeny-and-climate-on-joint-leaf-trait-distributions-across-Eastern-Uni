library(tidyverse)
joint_fia_clim <- data.table::fread("./indir/climate_plts_east_fia.csv")#, select = c("LAT", "LON", "daylength_9"))
geolocations = data.table::fread("./indir/2020_features.csv", select = c("statecd",
                                      "plot", "countycd", "unitcd", "spid","LAT","LON", "elevation","identifierID"))
joint_fia_clim = left_join(geolocations, joint_fia_clim)
joint_fia_clim$elevation <- joint_fia_clim$elevation/3.281

climate_FIA <- readRDS("./indir/2020climate_FIA.rds")
#climate_features = climate_FIA$climate_features
climate_features = NULL
var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
features = data.frame(joint_fia_clim[,10:93])

# # # apply PCA on FIA data
library(FactoMineR)
climate_features = NULL
climate_trends = NULL
for(ii in var){
  climate_trends[[ii]] <- features %>%
    dplyr::select(contains(ii))
  climate_FIA$climate_pca[[ii]]$var$coord = data.frame((climate_FIA$climate_pca[[ii]]$var$coord))
  foo = FactoMineR::predict.PCA(climate_FIA$climate_pca[[ii]], climate_trends[[ii]])
  climate_features[[ii]]  <- foo["cos2"]
}
rm(climate_trends)
climate_features = do.call(cbind.data.frame, climate_features)
colnames(climate_features) <- var
#climate_features = climate_FIA$climate_features
#update dataset with crunched climate features
joint_fia_clim = cbind.data.frame(joint_fia_clim[,c(1:9)], climate_features)
rm(climate_features)


# ppend terrain data to  climate, and save into fixed values. Reuse fixed_features.csv 
# any time you want to apply on a different model. Scaling will depend on the model fitted
topo_data = data.table::fread("/Users/sergiomarconi/Documents/Data/Surveys/FIA/FIA.csv")

#colnames(climate_features)[9] = "DTM"
colnames(topo_data) <- c("statecd", "unitcd", "countycd", "plot", "LAT", "LON","invyr", 
                         "subplot", "ope", "ect", "wdepth")

topo_data = topo_data %>% group_by(statecd, unitcd, countycd, plot) %>%
  filter(invyr < 2020) %>%
  filter(invyr > 2014) %>%
  summarize_all(mean)
topo_data = topo_data %>% select(-one_of("LAT", "LON", "subplot", "invyr"))

#rescale slope and aspect on a non angular scale
topo_data["ect"] = cos(topo_data["ect"]* 0.0174533)
topo_data["ope"] = sin(topo_data["ope"]* 0.0174533)
topo_data["individualID"] = paste(topo_data$statecd, topo_data$unitcd, topo_data$countycd, topo_data$plot, sep="_")
colnames(joint_fia_clim)[9] = "individualID"
features_fia = inner_join(joint_fia_clim, topo_data)
write_csv(features_fia, "./indir/2020f_features.csv")

#use model's scaling to rescale 
features_fia = readr::read_csv("./indir/2020f_features.csv")
bjtdm = readRDS( paste("./models/ch2perc_full.rds", sep=""))
colnames(features_fia)[c(10:16, 8,17, 18)]
colnames(features_fia)[8] =  "DTM"

scaling_x = bjtdm$scale[9:18,]
scaling_y = bjtdm$scale[1:8,]
features_fia = features_fia[complete.cases(features_fia),]

fia = features_fia %>% select(rownames(scaling_x))#[c(12:18, 9,19, 20)]

fia <- sapply(1:ncol(fia), function(x){
  scale(data.frame(fia)[x], center = unlist(scaling_x$mean[x]),
        scale = unlist(scaling_x$scale[x]))})

colnames(fia) = rownames(scaling_x)

fia = cbind.data.frame(features_fia[,-c(10:16, 8,17, 18)], fia)

#fia_scaled = cbind.data.frame(features_fia[,-c(11:17,9,20,21)], fia_scaled)

dict_ls<- readr::read_csv("~/Documents/Data/Surveys/FIA/species_list_fia.csv")
colnames(dict_ls)[4] <-  "spid"
joint_fia_clim <- left_join(fia, dict_ls)
joint_fia_clim  = joint_fia_clim  %>% select(-one_of("ott_id"))
write_csv(joint_fia_clim, "./indir/2020f_FIA_features.csv")
