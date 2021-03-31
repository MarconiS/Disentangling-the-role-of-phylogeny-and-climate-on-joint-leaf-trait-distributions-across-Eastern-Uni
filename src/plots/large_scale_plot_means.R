library(tidyverse)
library(sf)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.15, .85), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}
library(raster)
library(sp)


#getData('CMIP5', var='tmin', res=10, rcp=85, model='AC', year=70)
# subcontinental dataset
full = data.table::fread("./outdir/2020//FIA_2020_no_mlbs_full.csv")
fia_locations = readr::read_csv("/Users/sergiomarconi/Documents/GitHub/traits_tradeoffs_East/TOS_retriever/out/wordclim_fia.csv")
# r <- raster::getData("worldclim",var="bio",res=2.5)
# fia_locations = full[3:8] %>% unique
# fia_locations =  st_as_sf(fia_lcoations, coords=c("LON", "LAT"), crs=4326) 
# Tave = raster::extract(r, fia_locations)
# elevation =   raster::getData("alt",country="USA")
# elevation = raster::extract(elevation[[1]], fia_locations)
# fia_locations = cbind.data.frame(fia_locations, Tave, elevation)
# fia_locations=fia_locations %>% dplyr::select(statecd, plot, countycd, unitcd, bio1, bio2, bio3, bio4, 
#                                               bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13,
#                                               bio14, bio15, bio16, bio17, bio18, bio19, elevation) %>% unique
# rm(r, Tave, elevation)
# clim = readr::read_csv("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/2020_FIA_features.csv") %>%
#   select(statecd,  plot, countycd, unitcd, tmax,  tmin,elevation, ope) %>% unique
full = full[complete.cases(full),]
N_perc = full%>% group_by(LAT, LON) %>% summarize_if(is.numeric, mean)
N_perc$cat_LAT = N_perc$LAT %>% as.integer %>% factor
N_sp_perc = N_perc %>% group_by(cat_LAT) %>% summarize_if(is.numeric, mean)
N_lat = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=LAT)) + geom_point(alpha=1) + geom_smooth(method="loess") + theme_bw()

#link with elevation
N_perc = left_join(N_perc, fia_locations, by=c("statecd",  "plot", "countycd", "unitcd")) %>% unique
categories = arules::discretize(N_perc$elevation, method = "interval", breaks = 30)
N_perc$elev_c = categories
N_sp_perc = N_perc %>% group_by(elev_c) %>% summarize_if(is.numeric, mean)
N_ele = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=elevation)) + geom_point(alpha=1) + geom_smooth(method="loess") + theme_bw()
#link with tmean
N_perc$bio1 = N_perc$bio1/10
#N_perc = N_perc %>% filter(bio1 > 0.4)
categories = arules::discretize(N_perc$bio1, method = "interval", breaks = 30)
N_perc$tave = categories
N_sp_perc = N_perc %>% group_by(tave) %>% summarize_if(is.numeric, mean)
N_tave = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=bio1)) + geom_point(alpha=1) + geom_smooth(method="loess") + theme_bw() 

library(ggpubr)
ggarrange(N_lat, N_ele, N_tave, labels = c("(a)", "(b)", "(c)"), ncol = 3)
#nitrogen decorrelated
Nml = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), N_perc)
N_res = N_perc
N_res[8]= residuals(Nml)
#carbon decorrelated
Nml = lm(log10(carbonPercent.Estimate)~log10(leafMassPerArea.Estimate), N_perc)
N_res[9]= residuals(Nml)
N_res$cat_LAT = N_res$LAT %>% as.integer %>% factor
N_sp_perc = N_res %>% group_by(cat_LAT) %>% summarize_if(is.numeric, mean)
N_lat = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=LAT)) + geom_point(alpha=1) + 
  geom_smooth(method="loess", color="orange") + theme_bw() + theme(legend.position = "none")

#link with elevation
categories = arules::discretize(N_res$elevation, method = "interval", breaks = 30)
N_res$elev_c = categories
N_sp_perc = N_res %>% group_by(elev_c) %>% summarize_if(is.numeric, mean)
N_ele = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=elevation)) + geom_point(alpha=1) + 
geom_smooth(method="loess", color="orange") + theme_bw() + theme(legend.position = "none")

#link with tmean
categories = arules::discretize(N_res$bio1, method = "interval", breaks = 30)
N_res$tave = categories
N_sp_perc = N_res %>% group_by(tave) %>% summarize_if(is.numeric, mean)
N_tave = ggplot(N_sp_perc, aes(y=nitrogenPercent.Estimate, x=bio1)) + geom_point(alpha=1) + 
 geom_smooth(method="loess",  color = "orange") + theme_bw() + theme(legend.position = "none")

ggarrange(N_lat, N_ele, N_tave, labels = c("(d)", "(e)", "(f)"), ncol = 3)
