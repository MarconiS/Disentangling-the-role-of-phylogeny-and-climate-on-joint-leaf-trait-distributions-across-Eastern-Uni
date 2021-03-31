library(tidyverse)
library(sf)
library(raster)
library(sp)

fia_locations = readr::read_csv("~/Documents/Data/wordclim_fia.csv")
mods = c("./outdir/2020//FIA_2020_no_mlbs_full.csv", 
         "./outdir/2020//FIA_2020_no_mlbs_env.csv", 
         "./outdir/2020//FIA_2020_no_mlbs_sp.csv")

mt = c("blue", "orange", "lightgreen")
lbls = rbind(c("(a)", "(b)", "(c)"),
             c("(d)", "(e)", "(f)"),
             c("(g)", "(h)", "(i)"))
for(pt in  1:3){
  full = data.table::fread(mods[pt]) %>% data.frame
  full = full[complete.cases(full),]
  #nitrogen decorrelated
  Nml = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), full)
  N_res = full
  N_res[9]= residuals(Nml)
  #carbon decorrelated
  Nml = lm(log10(carbonPercent.Estimate)~log10(leafMassPerArea.Estimate), full)
  N_res[11]= residuals(Nml)
  
  
  N_res = N_res[complete.cases(N_res),]
  N_perc = N_res%>% group_by(LAT, LON) %>% summarize_if(is.numeric, mean)
  N_perc$cat_LAT = N_perc$LAT %>% as.integer %>% factor
  N_sp_perc = N_perc %>% group_by(cat_LAT) %>% summarize_if(is.numeric, mean)
  N_lat = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=LAT)) + geom_point(alpha=1) +
    geom_smooth(method="loess", color = mt[pt]) + theme_bw()
  
  C_lat = ggplot(N_sp_perc, aes(y=carbonPercent.Estimate, x=LAT)) + geom_point(alpha=1) +
    geom_smooth(method="loess", color = mt[pt]) + theme_bw()
  
  LMA_lat = ggplot(N_sp_perc, aes(y=leafMassPerArea.Estimate, x=LAT)) + geom_point(alpha=1) +
    geom_smooth(method="loess", color = mt[pt]) + theme_bw()
  # #link with elevation
  # categories = arules::discretize(N_perc$elevation, method = "interval", breaks = 30)
  # N_perc$elev_c = categories
  # N_sp_perc = N_perc %>% group_by(elev_c) %>% summarize_if(is.numeric, mean)
  # N_ele = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=elevation)) + geom_point(alpha=1) + 
  #   geom_smooth(method="loess", color = mt[pt]) + theme_bw()
  # #link with tmean
  # N_perc$bio1 = N_perc$bio1/10
  # #N_perc = N_perc %>% filter(bio1
  # categories = arules::discretize(N_perc$bio1, method = "interval", breaks = 30)
  # N_perc$tave = categories
  # N_sp_perc = N_perc %>% group_by(tave) %>% summarize_if(is.numeric, mean)
  # N_tave = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=bio1)) + geom_point(alpha=1) + 
  #   geom_smooth(method="loess",  color = mt[pt]) + theme_bw() + theme(legend.position = "none")
  
  ggarrange(N_lat, C_lat, LMA_lat, labels = lbls[pt,], ncol = 3)
  
}