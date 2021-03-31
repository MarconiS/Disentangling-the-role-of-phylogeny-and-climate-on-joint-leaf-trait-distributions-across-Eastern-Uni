library(sf)
clim = data.table::fread("indir/2020f_FIA_features.csv")

standard_climate = clim %>% data.frame
standard_climate[10:19] = scale(standard_climate[10:19]) 
standard_climate = sf::st_as_sf(standard_climate, coords = c("LON", "LAT"), crs = 4326)
ecoregions = sf::read_sf("/Users/sergiomarconi/Documents/Data/Metadata/eco_us.shp")
ecoregions = sf::st_transform(ecoregions, st_crs(standard_climate))
standard_climate = sf::st_join(standard_climate, ecoregions)
standard_climatesd_clim = standard_climate %>% dplyr::select(ECO_US_ID, prec,rad,tmax,tmin,vp,DTM,ope,ect) %>% group_by(ECO_US_ID) %>%  
  summarise_if(is.numeric, ~ sd(.x, na.rm = TRUE))
full_sf = sf::st_as_sf(full,coords=c("LON", "LAT"), crs=4326)
full_sf = sf::st_join(full_sf, ecoregions)
full_sf = full_sf %>% data.frame %>% dplyr::select("PROVINCE", "ECO_US_ID",
                                                 "leafMassPerArea.Estimate", "nitrogenPercent.Estimate",
                                                 "carbonPercent.Estimate", "ligninPercent.Estimate", 
                                                 "cellulosePercent.Estimate", "chlApercdw.Estimate", 
                                                 "chlBpercdw.Estimate", "carotpercdw.Estimate")
full_sf = full_sf[complete.cases(full_sf),]
for(ii in 4:10){
  tmp = full_sf[,c(ii, 2)]
  colnames(tmp) = c("y", "lma")
  md = lm(y ~ lma, data = log10(tmp))
  full_sf[,ii] = (residuals(md))
}

sd_traits = full_sf %>% filter(leafMassPerArea.Estimate < 600) %>% data.frame %>% dplyr::select(ECO_US_ID, 
                              nitrogenPercent.Estimate, chlApercdw.Estimate, chlBpercdw.Estimate, carbonPercent.Estimate,
                              carotpercdw.Estimate, ligninPercent.Estimate, 
                               cellulosePercent.Estimate, leafMassPerArea.Estimate) %>% group_by(ECO_US_ID) %>%
  summarise_if(is.numeric, ~ sd(.x, na.rm = TRUE))

standard_climate = standard_climate %>% data.frame() %>% select(c("ECO_US_ID", "PROVINCE")) %>% unique
sd_traits = left_join(sd_traits, standard_climate) %>% unique
ave_sd_clim = standard_climatesd_clim %>% data.frame %>% select(prec,rad,tmax,tmin,vp,DTM,ope,ect) %>% apply( 1, mean)
check = inner_join(sd_traits, cbind.data.frame(standard_climatesd_clim["ECO_US_ID"], ave_sd_clim ))
check = check  %>% data.frame %>% select(PROVINCE, ave_sd_clim, nitrogenPercent.Estimate, 
                         carotpercdw.Estimate,  chlApercdw.Estimate, chlBpercdw.Estimate,
                            leafMassPerArea.Estimate, carbonPercent.Estimate,ligninPercent.Estimate,cellulosePercent.Estimate)

colnames(check) = c("Ecoprovince", "Heterogenity", "Nm",  "Carot", "ChlA", "ChlB", "LMA", "Cm",  "lign", "cell")

check = gather(check, trait, variation, Nm:cell)
check = check[complete.cases(check),]
check$Ecoprovince = gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',check$Ecoprovince,perl = TRUE)
ggplot(check,aes(x=Heterogenity, y= variation, color=Ecoprovince))+ geom_point()+facet_wrap(.~trait, scale="free", nrow = 2)+
  theme_bw() + theme(legend.position = "none")
ggplot(check,aes(x=Heterogenity, y= variation, color=Ecoprovince))+ geom_point(alpha=0)+facet_wrap(.~trait, scale="free", nrow = 2)+
  theme_bw() + theme(legend.position = "bottom") + geom_smooth(method = "lm", se = F)
#version with both in 1 plot
ggplot(check,aes(x=Heterogenity, y= variation, color=Ecoprovince))+ geom_point(alpha=0.5, color="grey50", size=0.7)+facet_wrap(.~trait, scale="free", nrow = 2)+
  theme_bw() + theme(legend.position = "bottom") + geom_smooth(method = "lm", se = F)
# colnames(check) = c("Nm", "Cm",  "Carot", "lign", "cell", "LMA", 
#                     "Prec", "Rad", "Tmax", "Tmin", "Vp", "Elev", "Slope", "Aspect")
