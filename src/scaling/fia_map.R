library("ggspatial")
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library(readr)
library(sf)
library(sp)
library(tidyverse)
FIA_full <- read_csv("./outdir/fia_predictions/FIA_full_map.csv")
FIA_full <- FIA_full[complete.cases(FIA_full),]
FIA_full <- sf::st_as_sf(FIA_full, coords = c("LON", "LAT"), crs = 4326)
sf::write_sf(FIA_full, "./outdir/fia_predictions/Final_full.shp")

#species
FIA_sp <- read_csv("./preds/FIA_2020ssp_map.csv")
FIA_sp <- sf::st_as_sf(FIA_sp, coords = c("LON", "LAT"), crs = 4326)
sf::write_sf(FIA_sp, "./outdir/fia_predictions/Final_species.shp")

#climate
FIA_clim <- read_csv("./outdir/fia_predictions/FIA_climfinal_map.csv")
FIA_clim[11:82] <- exp(FIA_clim[11:82])
FIA_clim <- sf::st_as_sf(FIA_clim, coords = c("LON", "LAT"), crs = 4326)
sf::write_sf(FIA_clim, "./outdir/fia_predictions/Final_climate.shp")

full_pe<- FIA_full %>% dplyr::select(contains("Estimate"))
sp_pe<- FIA_sp %>% dplyr::select(contains(("Estimate")))

cl_pe <- FIA_clim %>% dplyr::select(contains("Estimate"))

tr_nm <-  c("LMA", "lignin", "Cellulose", "[P]", "[K]", "[Ca]", "[ChlA]",  "[ChlB]",
                                              "[Carotenoids]",  "[Mg]",  "[S]",  "[Mn]",  "[Fe]",  "[Cu]",  "[B]",  "[Zn]",
                                              "N", "C", "geometry")
colnames(full_pe) <-  colnames(sp_pe) <- colnames(cl_pe) <- tr_nm 
common_sp <- sf::st_join(full_pe, sp_pe, left = T)
common_clim <- sf::st_join(full_pe, cl_pe, left = T)

ii = 17
dt_sp <- dt_cl <- common_sp["geometry"]
for(ii in 1:18){
  delta_sp <- common_sp %>% mutate(delta = .[[ii]] - .[[ii + 19]]) %>%
    select(delta)
  dt_sp[tr_nm[ii]] <- delta_sp$delta
  
  delta_cl <- common_clim %>% mutate(delta = .[[ii]] - .[[ii + 19]]) %>%
    select(delta)
  dt_cl[tr_nm[ii]] <- delta_cl$delta
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_full.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  full_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xaxs = "i", yaxs = "i",
         xlim = c(-100, -67),  breaks = quantile(full_pe[[ii]], probs = seq(0.02,0.98, length.out = 15), na.rm = T), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1))
  dev.off()
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_delta_species.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  delta_sp %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(delta_sp[[1]], probs = seq(0.02,0.98, length.out = 15), na.rm = T), 
         pal = colorspace::diverging_hcl) 
  dev.off()
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_delta_clim.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  delta_cl %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(delta_cl[[1]], probs = seq(0.02,0.98, length.out = 15), na.rm = T), 
         pal = colorspace::diverging_hcl) 
  dev.off()
  
  #simulations having only climate or phylogeny
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_species.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  sp_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(full_pe[[ii]], probs = seq(0.02,0.98, length.out = 15), na.rm = T), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1)) 
  dev.off()
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_clim.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  cl_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(full_pe[[ii]], probs = seq(0.02,0.98, length.out = 15), na.rm = T), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1)) 
  dev.off()
}

sf::write_sf(dt_sp, "./outdir/fia_predictions/delta_sp.shp")
sf::write_sf(dt_cl, "./outdir/fia_predictions/delta_clim.shp")

#SE error
full_pe<- FIA_full %>% dplyr::select(contains("Error"))
sp_pe<- FIA_sp %>% dplyr::select(contains("Error"))
cl_pe<- FIA_clim %>% dplyr::select(contains("Error"))

colnames(full_pe) <-  colnames(sp_pe) <-  colnames(cl_pe)  <- c("LMA", "lignin", "Cellulose", "[P]", "[K]", "[Ca]", "[ChlA]",  "[ChlB]",
                                            "[Carotenoids]",  "[Mg]",  "[S]",  "[Mn]",  "[Fe]",  "[Cu]",  "[B]",  "[Zn]",
                                            "N", "C", "geometry")


for(ii in 1:18){
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_full_SE.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  full_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xaxs = "i", yaxs = "i",
         xlim = c(-100, -67),  breaks = quantile(full_pe[[ii]], probs = seq(0.02,0.98, length.out = 15)), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1))
  dev.off()
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_species_SE.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  sp_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(sp_pe[[ii]], probs = seq(0.02,0.98, length.out = 15)), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1)) 
  dev.off()
  
  filenm = paste("./outdir/maps/", colnames(full_pe)[ii], "_clim_SE.tif", sep="")
  tiff(filenm, width = 550, height = 550)
  cl_pe[ii] %>% 
    plot(cex = 0.2, pch = 19,  axes = F, key.pos = 1, xlim = c(-100, -67),  xaxs = "i", yaxs = "i",
         breaks = quantile(cl_pe[[ii]], probs = seq(0.02,0.98, length.out = 15)), 
         pal = colorspace::sequential_hcl(14, "blues",rev = T, power = 1)) 
  dev.off()
}


