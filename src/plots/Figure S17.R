norm =  readr::read_csv("./outdir/2020/FIA_2020_no_mlbs_full.csv")
norm = sf::st_as_sf(norm, coords=c("LON", "LAT"), crs = 4326)
ecoregions = sf::read_sf("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/atlantic/atlantic_outer.shp")
ecoregions = sf::st_transform(ecoregions, st_crs(norm))
norm = sf::st_join(norm, ecoregions)
norm = norm %>% filter(!is.na(ECO_US_))
yep= norm #%>% #group_by(statecd, countycd, unitcd, plot) %>%
  filter(taxonID == "QULA3") 
yep$N_LMA = yep$nitrogenPercent.Estimate/yep$leafMassPerArea.Estimate
crds = st_coordinates(yep)
yep = cbind.data.frame(crds, yep)
yep = yep %>% group_by(Y, taxonID) %>% summarize_if(is.numeric, mean)
# ggplot(data = world) +
#   xlim(-95,-65)+ ylim(24,50)+
#   theme_bw() +
#   geom_point(data = yep, aes(x = LON, y = LAT, 
#                              fill = cut(yep$N_LMA,quantile(yep$N_LMA))), alpha = 0.6, 
#              size = 0.7, shape = 23, stroke = 0) +  scale_fill_viridis_d()+
#   geom_sf(alpha = 0)+
#   geom_point(data = ACSA, aes(x = LON_site, y = LAT_site, color = "red", alpha = 0.9), 
#              size = 1, shape =4, stroke = 1) 
# #scale_colour_gradient("red")

model_95 = quantile(yep$N_LMA,probs=c(.025,.975))
#plot(yep$LAT, yep$N_LMA)
nLm=ggplot(yep, aes(x = Y, y = N_LMA))+ geom_density2d() + geom_smooth(method = "loess") + 
  theme_bw()+ facet_wrap(.~taxonID, scales = "free")
ggsave(nLm, "n_lma_ratio.png")
