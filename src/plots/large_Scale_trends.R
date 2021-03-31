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
full <- readr::read_csv("./outdir/2020/FIA_2020_no_mlbs_full.csv")
full = full %>% dplyr::filter(leafMassPerArea.Estimate < 600)
r <- raster::getData("worldclim",var="bio",res=2.5)
elevation = data.table::fread("./indir/2020f_features.csv",) %>% dplyr::select(statecd, plot, countycd, unitcd,elevation) %>%unique

fia_locations = full[3:8] %>% unique
fia_locations = left_join(fia_locations, elevation) %>% group_by(statecd, plot, countycd, unitcd) %>% summarize_all(mean)
fia_locations =  st_as_sf(fia_locations, coords=c("LON", "LAT"), crs=4326) 
Tave = raster::extract(r, fia_locations)
fia_locations = cbind.data.frame(fia_locations, Tave)
fia_locations=fia_locations %>% dplyr::select(statecd, plot, countycd, unitcd, bio1, bio2, bio3, bio4, 
                               bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13,
                               bio14, bio15, bio16, bio17, bio18, bio19, elevation) %>% unique
rm(r, Tave, elevation)
# clim = readr::read_csv("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/2020_FIA_features.csv") %>%
#   select(statecd,  plot, countycd, unitcd, tmax,  tmin,elevation, ope) %>% unique
N_perc = full[c(1:9,23)]
N_perc$cat_LAT = N_perc$LAT %>% as.integer %>% factor
N_sp_perc = N_perc %>% group_by(taxonID, cat_LAT) %>% summarize_if(is.numeric, mean)
plot(N_sp_perc$cat_LAT, N_sp_perc$nitrogenPercent.Estimate)

#link with elevation
library(arules)
N_perc = full[c(1:9,17)]
N_perc = left_join(N_perc, fia_locations, by=c("statecd",  "plot", "countycd", "unitcd")) %>% unique
categories = discretize(N_perc$elevation, method = "interval", breaks = 30)
N_perc$elev_c = categories
N_sp_perc = N_perc %>% group_by(taxonID, elev_c) %>% summarize_if(is.numeric, mean)
plot(N_sp_perc$elev_c, N_sp_perc$nitrogenPercent.Estimate)
#link with tmean
N_perc$bio1 = N_perc$bio1/10
categories = discretize(N_perc$bio1, method = "interval", breaks = 30)
N_perc$tave = categories
N_sp_perc = N_perc %>% group_by(taxonID, tave) %>% summarize_if(is.numeric, mean)
plot(N_sp_perc$tave, N_sp_perc$nitrogenPercent.Estimate)



Nml = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), N_perc)
N_res = N_perc
N_res[9]= residuals(Nml)
N_res$cat_LAT = N_res$LAT %>% as.integer %>% factor
N_sp_perc = N_res %>% group_by(taxonID, cat_LAT) %>% summarize_if(is.numeric, mean)
summary(mgcv::gam(nitrogenPercent.Estimate~cat_LAT, data=N_res))
plot(N_sp_perc$cat_LAT, N_sp_perc$nitrogenPercent.Estimate)

#link with elevation
categories = discretize(N_res$elevation, method = "interval", breaks = 30)
N_res$elev_c = categories
N_sp_perc = N_res %>% group_by(taxonID, elev_c) %>% summarize_if(is.numeric, mean)
summary(mgcv::gam(nitrogenPercent.Estimate~elev_c, data=N_res))
plot(N_sp_perc$elev_c, N_sp_perc$nitrogenPercent.Estimate)
#link with tmean
#N_res$bio1 = N_res$bio1/10
categories = discretize(N_res$bio1, method = "interval", breaks = 30)
N_res$tave = categories
N_sp_perc = N_res %>% group_by(taxonID, tave) %>% summarize_if(is.numeric, mean)
summary(mgcv::gam(nitrogenPercent.Estimate~tave, data=N_res))
plot(N_sp_perc$tave, N_sp_perc$nitrogenPercent.Estimate)
