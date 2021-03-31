library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = 50, returnclass = "sf")
us_eco = sf::read_sf("/Users/sergiomarconi/Documents/Data/Metadata/eastern.shp")
transferability = readr::read_csv("./indir/external_validation/external_dataset.csv")
neon_data = base_data %>% select(siteID, latitude, longitude) %>% group_by(siteID) %>% summarize_all(mean)
neon_data = neon_data %>% filter(siteID !="MLBS")
fia = data.table::fread("./outdir/2020/FIA_2020_no_dl_full.csv")
fia = fia %>% select(LAT, LON) %>% unique %>% data.frame
try = read_csv("./indir/TRY_species_sample.csv")
try$SpeciesName %>% table %>% 
  sort(decreasing = T) %>% head(40)

try = try %>% filter(SpeciesName %in% c("Fagus grandifolia", "Acer rubrum", "Abies balsamea"))
  
ggplot(data = world) +
  ylim(10,70)+ xlim(-160,-50)+
  theme_bw() +
  geom_point(data = fia, aes(x = LON, y = LAT), color = "lightblue", size = 0.2, shape = 23) +
  geom_sf(alpha = 0)+
  geom_sf(data=us_eco, alpha=0)+
  geom_point(data = try, aes(x = LON_site, y = LAT_site), color = "red", 
             size = 1, shape =4, stroke = 1) + 
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank())+
  geom_point(data = neon_data, aes(x = longitude, y = latitude), color = "blue",pch=21) + 
  geom_point(data = transferability, aes(x = latitude, y = longitude), color = "red", 
             size = 2, shape =4, stroke = 1)+ theme(legend.position="bottom") #+   
#scale_colour_gradient("red")

