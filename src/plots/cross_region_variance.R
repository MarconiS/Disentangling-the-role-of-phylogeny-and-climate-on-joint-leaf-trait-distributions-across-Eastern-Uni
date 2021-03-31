#discontinued


sf_full = sf::st_as_sf(full, coords = c("LON", "LAT"), crs = 4326)
ecoregions = sf::read_sf("/Users/sergiomarconi/Documents/Data/Metadata/eastern.shp")
ecoregions = st_transform(ecoregions, sf::st_crs(sf_full))

sf_full = st_join(sf_full, clim)
ave_rg = sf_full %>% group_by(PROVINCE) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
sd_rg = sf_full %>% group_by(PROVINCE) %>%
  summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))

sd_tot = sf_full %>%
  summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))

sd_inter = ave_rg %>%
  summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))
sd_intra = sd_rg %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
compare_var = rbind.data.frame(sd_inter, sd_intra, sd_tot) %>% 
  select(nitrogenPercent.Estimate, carbonPercent.Estimate, extractChlAConc.Estimate,
         extractChlBConc.Estimate, extractCarotConc.Estimate, ligninPercent.Estimate, 
         cellulosePercent.Estimate, leafMassPerArea.Estimate)
compare_var$grp = c("cross_region", "within_region", "continent")
compare_var = reshape2::melt(compare_var, id.vars = c("grp", "geometry"))
ggplot(compare_var, aes(x = variable, y=value, fill = grp)) + theme_bw()+
  geom_bar(position="dodge", stat="identity", color="black") + facet_wrap(variable ~.,ncol = 4, nrow = 2, scale= "free")
