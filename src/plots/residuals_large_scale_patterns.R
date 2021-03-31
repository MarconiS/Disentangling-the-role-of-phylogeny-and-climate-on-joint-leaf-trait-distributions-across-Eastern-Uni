#figure 3
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = 50, returnclass = "sf")
full = data.table::fread("./outdir/2020//FIA_2020_no_mlbs_full.csv")
full = final_plot_preds
sp = data.table::fread("./outdir/2020//FIA_2020_no_mlbs_sp.csv")
clim = data.table::fread("./outdir/2020//FIA_2020_no_mlbs_env.csv")
#clim = readr::read_csv("./outdir/outdir/FIA_full_map.csv") %>% select(statecd,  plot, countycd, unitcd, LON, LAT)
#full = joint_fia_clim
full=full[complete.cases(full),] %>% data.frame
per_mass_data = full %>% select(contains(c("Estimate", "LON", "LAT")))

mlm2 <- lm(log10(nitrogenPercent.Estimate) ~ log10(leafMassPerArea.Estimate), data = per_mass_data)
test_res = residuals(mlm2)
per_mass_data$nitrogenPercent.Estimate = test_res
#world <- sf::read_sf("/Volumes/Stele/eco-us-shp/eco_us.shp")
#world = sf::st_transform(world, 4326)
#full[c("nitrogenPercent.Estimate")]=(test_res)
full$individualID = paste(full$statecd, full$unitcd, full$countycd, full$plot, sep="_")
full$nitrogenPercent.Est.Error = full$nitrogenPercent.Q97.5 - full$nitrogenPercent.Q2.5 
full$carbonPercent.Est.Error = full$carbonPercent.Q97.5 - full$carbonPercent.Q2.5 
full$ligninPercent.Est.Error = full$ligninPercent.Q97.5 - full$ligninPercent.Q2.5 
full$cellulosePercent.Est.Error = full$cellulosePercent.Q97.5 - full$cellulosePercent.Q2.5 
full$carotpercdw.Est.Error = full$carotpercdw.Q97.5 - full$carotpercdw.Q2.5 
full$chlApercdw.Est.Error = full$chlApercdw.Q97.5 - full$chlApercdw.Q2.5 
full$chlBpercdw.Est.Error = full$chlBpercdw.Q97.5 - full$chlBpercdw.Q2.5 
full$leafMassPerArea.Est.Error = full$leafMassPerArea.Q97.5 - full$leafMassPerArea.Q2.5 
full= full[complete.cases(full), ]
full = full %>% group_by(individualID)%>% summarize_if(is.numeric, mean)
for(tr in colnames(full)[c(10,14,18,22,26,30,34,38)]){
  plt = ggplot(data = world) +  ylim(25, 50)+ xlim(-96,-67) +
    geom_point(data = full, aes_string(x = "LON", y = "LAT", 
                                       fill = paste("cut(",tr, 
                                                    ", quantile(",tr, ", probs = seq(0, 1, 0.1), na.rm = FALSE))")),  
               size = 0.8, shape = 23, stroke = 0, alpha = 0.5) +  scale_fill_viridis_d(name = tr) + 
    theme_bw()+  theme(legend.position="bottom", panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) + geom_sf(aes(),alpha = 0, show.legend=FALSE)
  
  ggsave(plot=plt, filename = paste("./outdir/maps/", tr, ".png", sep=""))  
}
# 
# full_perc["individualID"] = paste(full_perc$statecd, full_perc$unitcd, full_perc$countycd, full_perc$plot, sep="_")
# full_perc = full_perc %>% group_by(individualID)%>% summarize_if(is.numeric, mean)
# ggplot(data = world) + geom_sf(aes(alpha = 0))+
#   ylim(25, 50)+ xlim(-96,-67)+
#   geom_point(data = full_perc, aes(x = LON, y = LAT, 
#                               fill = cut(full_perc$leafMassPerArea.Estimate,quantile(full_perc$leafMassPerArea.Estimate))),  
#              size = 0.7, shape = 23, stroke = 0) +  scale_fill_viridis_d()+ theme_bw()+  theme(legend.position="bottom") 
# foo = data.frame(log10(full_perc["leafMassPerArea.Estimate"]), test_res)
# colnames(foo)= c("LMA", "N", "C", "lignin", "cellulose")
# check = cor(foo)
# corrplot::corrplot(check, type = "lower")
