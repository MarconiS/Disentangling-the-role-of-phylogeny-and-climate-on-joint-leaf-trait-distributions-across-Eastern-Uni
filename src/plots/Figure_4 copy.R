library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
# subcontinental dataset
outdir="/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/"
full <- readr::read_csv(paste(outdir, "/FIA_2020_no_mlbs_full.csv", sep=""))
species <- readr::read_csv(paste(outdir, "/FIA_2020_no_mlbs_sp.csv", sep=""))
environment <- readr::read_csv(paste(outdir, "/FIA_2020_no_mlbs_env.csv", sep=""))
ecoregions = sf::read_sf(paste(outdir, "/us_eco_l3/us_eco_l3.shp", sep=""))
#all_regions = sf::read_sf(paste(outdir, "/FIA_2020_no_mlbs_env.csv", sep=""))




#remove outliers
# vl_full <- 
#   apply(full[,9:ncol(full), with=FALSE], 2,remove_outliers) %>% data.frame
# 
# vl_sp <- 
#   apply(species[,9:ncol(species), with=FALSE], 2,remove_outliers) %>% data.frame
# 
# vl_env <- 
#   apply(environment[,9:ncol(environment), with=FALSE], 2,remove_outliers) %>% data.frame


#get divergence from environment
delta_env = full[,9:ncol(full), with=FALSE] - 
  environment[,9:ncol(environment), with=FALSE]
delta_sp =  full[,9:ncol(full), with=FALSE] - 
  species[,9:ncol(species), with=FALSE]



#reassamble
#vl_full = cbind.data.frame(full[,1:8], vl_full)
#vl_sp = cbind.data.frame(species[,1:8], vl_sp)
#vl_env = cbind.data.frame(environment[,1:8], vl_env)
delta_env = cbind.data.frame(full[,1:8], delta_env)
delta_sp = cbind.data.frame(full[,1:8], delta_sp)

vl_full=full[complete.cases(full),]
vl_sp=species[complete.cases(species),]
vl_env=environment[complete.cases(environment),]

#decouple 
mlm2 = lm(cbind(log10(nitrogenPercent.Estimate), log10(carbonPercent.Estimate), 
                 log10(ligninPercent.Estimate), log10(cellulosePercent.Estimate)) ~ 
             log10(leafMassPerArea.Estimate), data = (vl_full))
test_res = residuals(mlm2)
vl_full[c("nitrogenPercent.Estimate", "carbonPercent.Estimate" , "ligninPercent.Estimate" , "cellulosePercent.Estimate")]=(test_res)

mlm2 = lm(cbind(log10(nitrogenPercent.Estimate), log10(carbonPercent.Estimate), 
                 log10(ligninPercent.Estimate), log10(cellulosePercent.Estimate)) ~ 
             log10(leafMassPerArea.Estimate), data = (vl_sp))
test_res = residuals(mlm2)
vl_sp[c("nitrogenPercent.Estimate", "carbonPercent.Estimate" , "ligninPercent.Estimate" , "cellulosePercent.Estimate")]=(test_res)

mlm2 = lm(cbind(log10(nitrogenPercent.Estimate), log10(carbonPercent.Estimate), 
                 log10(ligninPercent.Estimate), log10(cellulosePercent.Estimate)) ~ 
             log10(leafMassPerArea.Estimate), data = (vl_env))
test_res = residuals(mlm2)
vl_env[c("nitrogenPercent.Estimate", "carbonPercent.Estimate" , "ligninPercent.Estimate" , "cellulosePercent.Estimate")]=(test_res)


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

traits = c("nitrogenPercent.Estimate", "extractCarotConc.Estimate", "extractChlAConc.Estimate", "extractChlBConc.Estimate",
           "ligninPercent.Estimate","carbonPercent.Estimate", "cellulosePercent.Estimate", "leafMassPerArea.Estimate")
for(tr in traits){
  foo <- ggplot(data = world) +
    theme_bw() +   ylim(25, 50)+ xlim(-96,-67)+
    geom_point(data = vl_sp, aes_string(x = "LON", y = "LAT", fill = tr), 
               alpha = 0.7, size = 1.2, shape = 23, stroke = 0) +
    theme(legend.position="bottom") +
    scale_fill_viridis_c()+
    #scale_fill_gradient(low = "white", high = "blue")+
    #scale_fill_gradientn(colours=fill.colors)+#, 
    #breaks = quantile(vl_sp$nitrogenPercent.Estimate, probs = seq(0,1,0.05)))+
    geom_sf(alpha = 0)
  ggsave(paste("./outdir/maps/species_", tr, ".png",  sep=""), foo)
}

#get the deltas
library(wesanderson)
delta_env = delta_env[complete.cases(delta_env), ]
delta_sp = delta_sp[complete.cases(delta_sp), ]
world <- ne_countries(scale = "medium", returnclass = "sf")

pal <- wes_palette("Zissou1", 100, type = "continuous")
delta_Ne <- ggplot(data = world) +
  theme_bw() +  ylim(25, 50)+ xlim(-96,-67)+
  geom_point(data = delta_sp, aes(x = LON, y = LAT, 
                                   fill = nitrogenPercent.Estimate), 
             size = 0.8, shape = 23, stroke = 0) +theme(legend.position="bottom") +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(midpoint = 0, low = "springgreen4", mid = "white",
                        high = "royalblue1", space = "Lab" )+
  geom_sf(alpha = 0)

foo = delta_env %>% filter(nitrogenPercent.Estimate>-2)
delta_Ns <- ggplot(data = world) +
  theme_bw() +   ylim(25, 50)+ xlim(-96,-67)+
  geom_point(data = foo, aes(x = LON, y = LAT, fill = nitrogenPercent.Estimate), 
             size = 1.5, shape = 23, stroke = 0) +theme(legend.position="bottom") +
  scale_fill_gradient2(midpoint = 0, low = "orangered3", mid = "white",
                       high = "royalblue1", space = "Lab" )+
  geom_sf(alpha = 0)
# get distributions per ecoregion by appending ecoregions to dataframe
all_regions = st_transform(all_regions, crs = 4326)
full = sf::st_as_sf(full, coords = c("LON", "LAT"), crs = 4326)
full = sf::st_join(full, all_regions)

species = sf::st_as_sf(species, coords = c("LON", "LAT"), crs = 4326)
species = sf::st_join(species, all_regions)

environment = sf::st_as_sf(environment, coords = c("LON", "LAT"), crs = 4326)
environment = sf::st_join(environment, all_regions)

delta_intra = sf::st_as_sf(delta_sp, coords = c("LON", "LAT"), crs = 4326)
delta_intra = sf::st_join(delta_intra, all_regions)

delta_inter = sf::st_as_sf(delta_env, coords = c("LON", "LAT"), crs = 4326)
delta_inter = sf::st_join(delta_inter, all_regions)

full$SECTION %>% unique
reg1 = full%>%filter(SECTION == "Adirondack Highlands Section")

#plot LMA and N per ecoregion
eco_dist = full %>% select(nitrogenPercent.Estimate, leafMassPerArea.Estimate, SECTION)
eco_delta = delta_intra %>% select(nitrogenPercent.Estimate, leafMassPerArea.Estimate, SECTION)
env_delta = delta_inter %>% select(nitrogenPercent.Estimate, leafMassPerArea.Estimate, SECTION)

#scale to make it more comparable in plot
eco_dist$nitrogenPercent.Estimate = scale(eco_dist$nitrogenPercent.Estimate)
eco_dist$leafMassPerArea.Estimate = scale(eco_dist$leafMassPerArea.Estimate)

eco_dist = gather(eco_dist, trait, value, nitrogenPercent.Estimate:leafMassPerArea.Estimate)
eco_dist = eco_dist %>% filter(!is.na(SECTION))
eco_delta = eco_delta %>% filter(!is.na(SECTION))
env_delta = env_delta %>% filter(!is.na(SECTION))

hsits_pl = ggplot(eco_dist, aes(x = value, fill = trait))+ geom_histogram(position="identity", alpha=0.5) + 
  theme_bw() + facet_wrap( .~ SECTION, scales="free_y") + xlim(-2.5,5)
ggsave("./eco_dist_NLMA.png", hsits_pl)
hsits_pl = ggplot(eco_delta, aes(x = nitrogenPercent.Estimate))+ geom_histogram() + 
  theme_bw() + facet_wrap( .~ SECTION, scales="free_y") + xlim(-1,1) + geom_vline(xintercept = 0)
ggsave("./sp_dist_N.png", hsits_pl)

hsits_pl = ggplot(eco_delta, aes(x = leafMassPerArea.Estimate))+ geom_histogram() + 
  theme_bw() + facet_wrap( .~ SECTION, scales="free_y") + xlim(-100,100) + geom_vline(xintercept = 0)
ggsave("./sp_dist_LMA.png", hsits_pl)

hsits_pl = ggplot(env_delta, aes(x = nitrogenPercent.Estimate))+ geom_histogram() + 
  theme_bw() + facet_wrap( .~ SECTION, scales="free_y") + xlim(-1,1) + geom_vline(xintercept = 0)
ggsave("./env_dist_N.png", hsits_pl)
#"Florida Coastal Lowlands (Western) Section")
#"Mississippi Alluvial Basin Section"
#Northern Great Lakes Section
#Southern Superior Uplands Section
reg1 = full%>%filter(SECTION == "Green, Taconic, Berkshire Mountains Section")
reg1sp = species%>%filter(SECTION == "Green, Taconic, Berkshire Mountains Section")
reg1env = environment%>%filter(SECTION == "Green, Taconic, Berkshire Mountains Section")

SP = data.frame(c(reg1$nitrogenPercent.Estimate, reg1sp$nitrogenPercent.Estimate), 
                        c(rep("Full", nrow(reg1)), rep("Species", nrow(reg1))))
ENV = data.frame(c(reg1$nitrogenPercent.Estimate, reg1env$nitrogenPercent.Estimate), 
                      c(rep("Full", nrow(reg1)), rep("Environment", nrow(reg1))))
colnames(SP) = colnames(ENV) = c("N(%)", "Simulation")
ggplot(data = SP, aes(x = `N(%)`, fill = Simulation)) + geom_density(alpha=.5, adjust = 2) +
  xlim(0.5, 3.5) + ylim(0,1.7)+
  scale_discrete_manual(aesthetics = "fill", values = c("royalblue1", "orange")) + theme_bw()

ggplot(data = ENV, aes(x = `N(%)`, fill = Simulation)) + geom_density(alpha=.5, adjust = 2) + 
  xlim(0.5, 3.5) + ylim(0,1.7)+
  scale_discrete_manual(aesthetics = "fill", values = c("springgreen4", "royalblue1")) + theme_bw()

