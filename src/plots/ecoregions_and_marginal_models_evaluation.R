library(tidyverse)
library(sf)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
# load data
world <- ne_countries(scale = "medium", returnclass = "sf")

# subcontinental dataset
full <- fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_full.csv")
species <- fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_sp.csv")
environment <- fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_env.csv")
selected_ecoregions = sf::read_sf("/Volumes/Data/Metadata/na_cec_eco_l4/us_eco_l4_no_st.shp")
ecoregions = sf::read_sf("/Volumes/Data/Metadata/na_cec_eco_l4/us_eco_l4_no_st.shp")
sampled_eco = c("Western Allegheny Plateau", #non influent
                "Ridge and Valley", #non influent
                "Western Corn Belt Plains", # env -
                "Ozark Highlands" # env +
)
selected_ecoregions = selected_ecoregions %>% filter(US_L3NAME %in% sampled_eco)
ecoregions = ecoregions %>% group_by(US_L3NAME) %>% summarize_if(is.numeric, mean)


selected_ecoregions = selected_ecoregions %>% group_by(US_L3NAME) %>% summarize_if(is.numeric, mean)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(0, 1), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}
#remove outliers
vl_full <- 
  apply(full[,9:ncol(full), with=FALSE], 2,remove_outliers) %>% data.frame

vl_sp <- 
  apply(species[,9:ncol(species), with=FALSE], 2,remove_outliers) %>% data.frame

vl_env <- 
  apply(environment[,9:ncol(environment), with=FALSE], 2,remove_outliers) %>% data.frame


#get divergence from environment
delta_env = vl_full - vl_env
delta_sp = vl_full - vl_sp



#reassamble
vl_full = cbind.data.frame(full[,1:8], vl_full)
vl_sp = cbind.data.frame(species[,1:8], vl_sp)
vl_env = cbind.data.frame(environment[,1:8], vl_env)
delta_env = cbind.data.frame(full[,1:8], delta_env)
delta_sp = cbind.data.frame(full[,1:8], delta_sp)

vl_full=vl_full[complete.cases(vl_full),]
vl_sp=vl_sp[complete.cases(vl_sp),]
vl_env=vl_env[complete.cases(vl_env),]

#get plot averages
vl_full = vl_full %>% group_by(LAT, LON) %>% summarize_if(is.numeric, funs(mean(., na.rm=T)))
vl_env = vl_env %>% group_by(LAT, LON) %>% summarize_if(is.numeric, funs(mean(., na.rm=T)))
vl_sp = vl_sp %>% group_by(LAT, LON) %>% summarize_if(is.numeric, funs(mean(., na.rm=T)))
delta_env = delta_env %>% group_by(LAT, LON) %>% summarize_if(is.numeric, funs(mean(., na.rm=T)))
delta_sp = delta_sp %>% group_by(LAT, LON) %>% summarize_if(is.numeric, funs(mean(., na.rm=T)))

col.pal <- colorRampPalette(c("khaki", "yellow", "#FF7F00", "red", "#7F0000"))
#col.pal <- colorRampPalette(c("white", "darkblue"))

fill.colors <- col.pal(64)
vl_sp = vl_sp[complete.cases(vl_sp),]
vl_env = vl_env[complete.cases(vl_env),]
vl_full = vl_full[complete.cases(vl_full),]

#get the deltas
library(wesanderson)
delta_env = delta_env[complete.cases(delta_env), ]
delta_sp = delta_sp[complete.cases(delta_sp), ]

pal <- wes_palette("Zissou1", 100, type = "continuous")
ggplot(data = world) +
  theme_bw() +  ylim(25, 50)+ xlim(-96,-67)+
  geom_point(data = delta_sp, aes(x = LON, y = LAT, 
                                  fill = nitrogenPercent.Estimate), 
             size = 0.8, shape = 23, stroke = 0) +theme(legend.position="bottom") +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(limits = c(-1.5,1.5), midpoint = 0, low = "springgreen4", mid = "white",
                       high = "royalblue1", space = "Lab" )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_sf(data = ecoregions, alpha = 0, size = 0.05)+
  geom_sf(data = selected_ecoregions, alpha = 0, color="magenta")+
  geom_sf(alpha = 0)

ggplot(data = world) +
  theme_bw() +  ylim(25, 50)+ xlim(-96,-67)+
  geom_point(data = delta_env, aes(x = LON, y = LAT, 
                                   fill = nitrogenPercent.Estimate), 
             size = 0.8, shape = 23, stroke = 0) +theme(legend.position="bottom") +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(limits = c(-1.5,1.5), midpoint = 0, low = "orange", mid = "white",
                       high = "royalblue1", space = "Lab" )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_sf(data = ecoregions, alpha = 0, size = 0.05)+
  geom_sf(data = selected_ecoregions, alpha = 0, color="magenta")+
  geom_sf(alpha = 0)

# get distributions per ecoregion by appending ecoregions to dataframe
ecoregions = st_transform(ecoregions, crs = 4326)
full = sf::st_as_sf(full, coords = c("LON", "LAT"), crs = 4326)
full = sf::st_join(full, ecoregions)

species = sf::st_as_sf(species, coords = c("LON", "LAT"), crs = 4326)
species = sf::st_join(species, ecoregions)

environment = sf::st_as_sf(environment, coords = c("LON", "LAT"), crs = 4326)
environment = sf::st_join(environment, ecoregions)

delta_intra = sf::st_as_sf(delta_sp, coords = c("LON", "LAT"), crs = 4326)
delta_intra = sf::st_join(delta_intra, ecoregions) %>% data.frame 
delta_intra$modelID = "Species"
delta_inter = sf::st_as_sf(delta_env, coords = c("LON", "LAT"), crs = 4326)
delta_inter = sf::st_join(delta_inter, ecoregions) %>% data.frame 
delta_inter$modelID = "Environment"

delta_ = rbind.data.frame(delta_intra,  delta_inter)
traits = delta_inter %>% select(contains("Estimate")) %>% colnames

delta_$PID = gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',delta_$US_L3NAME,perl = TRUE)
cbind(delta_$PID, delta_$US_L3NAM) %>% unique
for(tr in traits){
  hsits_pl=ggplot(delta_, aes_string(x = tr, fill = "modelID"))+ geom_histogram(alpha=0.5, position = "identity") + 
    theme_bw() + facet_wrap( .~ PID, scales="free") +  
    scale_fill_manual(values=c("#E69F00", "#2C830F")) + 
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_vline(xintercept = 0)
  ggsave(paste("./plots/delta_dist_",tr, ".png", sep=""), hsits_pl, width = 2*4.96, height = 2* 3.86)
}
pairwise = ttest_sp = ttest_env = list()
for(eco in unique(delta_inter$US_L3NAME)){
  if(!is.na(eco)){
    single_t_env = delta_inter %>% filter(US_L3NAME == eco) %>% ungroup %>% select(traits)
    single_t_sp = delta_intra %>% filter(US_L3NAME == eco) %>% ungroup %>% select(traits)
    
    pairwise[[eco]] = lapply(1:8, function(x)t.test(single_t_env[,x], single_t_sp[,x]))
    ttest_sp[[eco]] = lapply(1:8, function(x)t.test(single_t_env[,x]))
    ttest_env[[eco]] = lapply(1:8, function(x)t.test(single_t_sp[,x]))
  }
}
pairwise = lapply(1:length(pairwise), function(x) pairwise[[x]][[1]]$p.value) %>% unlist %>% cbind(names(pairwise))
ttest_sp = lapply(1:length(ttest_sp), function(x) ttest_sp[[x]][[1]]$p.value) %>% unlist %>% cbind(names(ttest_sp))
ttest_env = lapply(1:length(ttest_env), function(x) ttest_env[[x]][[1]]$p.value) %>% unlist %>% cbind(names(ttest_env))

sum(as.numeric(pairwise[,1]) <0.0001) /dim(pairwise)[1]
