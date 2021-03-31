library(sf)
library(tidyverse)

out_rm  <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.05, .95), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}


eco_us <- sf::read_sf("/Volumes/Data/Chapter2/NA_CEC_Eco_Level3/NA_CEC_Eco_Level3.shp")

FIA_full <- sf::read_sf("~/Downloads/rstudio-export/Final_full_model_means.shp")#"/Volumes/Data/Chapter2/Maps/median/FIA_lfull_map.shp")
FIA_sp <- sf::read_sf("~/Downloads/rstudio-export/delta_sp.shp")#"/Volumes/Data/Chapter2/Maps/median/FIA_lfull_map.shp")
FIA_clim <- sf::read_sf("~/Downloads/rstudio-export/delta_clim.shp")#"/Volumes/Data/Chapter2/Maps/median/FIA_lfull_map.shp")

east_us <- sf::read_sf("./Maps/Ecoregions/Study_area.shp")
eco_us <- st_transform(eco_us, crs = 4326)
east_us <- st_transform(east_us, crs = 4326)
new_fia <- st_join(FIA_full, eco_us, join = st_intersects)
new_sp <- st_join(FIA_sp, eco_us, join = st_intersects)
new_cl <- st_join(FIA_clim, eco_us, join = st_intersects)


sf::write_sf(new_fia, "/Volumes/Data/Chapter2/Maps/median/FIA_delta_sp_ER.shp")
sf::write_sf(new_sp, "/Volumes/Data/Chapter2/Maps/median/FIA_delta_sp_ER.shp")
sf::write_sf(new_cl, "/Volumes/Data/Chapter2/Maps/median/FIA_delta_sp_ER.shp")


#test_region <- new_fia %>% as.data.frame() %>% select(NA_L3NAME, contains("_E")) 
test_region <- new_fia %>% as.data.frame() %>% select(NA_L3NAME, X.N)#, X.P., X.Zn.) 
sp_region <- new_sp %>% as.data.frame() %>% select(NA_L3NAME, X.N)#, X.P., X.Zn.) 
cl_region <- new_cl %>% as.data.frame() %>% select(NA_L3NAME, X.N)#, X.P., X.Zn.) 

delta <-  st_join(new_cl, new_sp, join = st_intersects)

dt <- delta %>% as.data.frame() %>% select(NA_L3NAME.x, contains("X.S")) 
dt <- new_fia %>% as.data.frame() %>% select(NA_L3NAME, lfMPA_E, flrPhC_E, ntP_E, flIC_E)
dt[-1] <- scale(log(dt[-1]))
dt <- gather(dt, trait, value, 2:5, factor_key=TRUE)

test_region[-1] <- apply(test_region[-1], 2, out_rm)
test_region <- test_region %>% select(-one_of(s_error))
#test_region[-1] <- scale(test_region[-1])
test_region <- gather(test_region, trait, value, 2, factor_key=TRUE) #lfMPA_E:crP_E
# ggplot(test_region, aes(x = value, fill = NA_L3NAME)) + geom_density(alpha=.3) + 
#   facet_wrap( ~ trait, ncol=6, scales="free") #+ xlim(-5, 5)

dt_c <- dt %>% filter(NA_L3NAME == "Southeastern Plains")
ggplot(dt_c, aes(x = value, fill = trait)) + geom_density(adjust = 3) + 
  facet_wrap( ~ trait, ncol=2, scales="free") + theme_bw() + xlim(-3.5, 3) +
  scale_fill_manual(values=c("coral", "cyan3", "antiquewhite", "green"))
  # scale_fill_manual(values=c("coral", "coral1", "antiquewhite", 
  #                            "cyan", "cyan1", "grey55", "green",
  #                            "green1", "green2", "grey56", "cyan2", "indianred1",
  #                            "grey57", "grey58", "grey59", "indianred2", "cyan3", "cyan4")) #+ xlim(-5, 5)

ggsave("./outdir/trends/lma_n.png")

