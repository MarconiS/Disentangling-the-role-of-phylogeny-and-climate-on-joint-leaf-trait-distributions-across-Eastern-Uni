library(tidyverse)
library(sf)# subcontinental dataset
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = 50, returnclass = "sf")
data <- data.table::fread("./outdir/2020/FIA_2020_no_mlbs_full.csv")
yep= data %>% group_by(statecd, countycd, unitcd, plot) %>%
  summarize_if(is.numeric, mean)
yep = yep[complete.cases(yep), ]
#yep = sf::st_as_sf(yep, coords = c("LAT", "LON"), crs = 4326)

ggplot(data = world) +
  xlim(-96,-65)+ ylim(24,50)+
  theme_bw() +
  geom_point(data = yep, aes(x = LON, y = LAT, 
                              fill = cut(carotpercdw.Estimate,
                                         quantile(carotpercdw.Estimate, probs = seq(0, 1, 0.1)))), alpha = 0.6, 
             size = 0.7, shape = 23, stroke = 0) +  scale_fill_viridis_d()+
  geom_sf(alpha = 0)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.position="bottom") #+   


neon = neon %>% group_by(siteID) %>% top_n(1)
#plot FAGR points on map
ggplot(data = world) +
  theme_bw() +
  geom_point(data = neon, aes(x = longitude, y = latitude), 
             color = "blue") +  geom_sf(aes(alpha = 0))

