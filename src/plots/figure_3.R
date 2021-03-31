#figure 3
library(tidyverse)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = 50, returnclass = "sf")
data <- data.table::fread("./outdir/2020/FIA_2020_no_mlbs_full.csv")
dt_sp <- data.table::fread("./outdir/2020/FIA_2020_no_mlbs_sp.csv")
NEtry = read_csv("./indir/TRY_species_sample.csv")
neon <- data.table::fread("/Users/sergiomarconi/Documents/Data/Traits/CFC/final_dataset_with_itc_0620.csv") %>% 
  filter(taxonID %in% c("ACRU")) %>%
  filter(itcLongitude > -100) 
yep= data %>% #group_by(statecd, countycd, unitcd, plot) %>%
  filter(taxonID == "FAGR") 

sp= dt_sp %>% #group_by(statecd, countycd, unitcd, plot) %>%
  filter(taxonID == "ACRU") 

NEtry$SpeciesName %>% table %>% 
  sort(decreasing = T) %>% head(40)

TRY = NEtry %>% filter(SpeciesName == "Acer rubrum")
TRY$OrigValueStr = as.numeric(TRY$OrigValueStr) 
N = TRY %>% filter(TraitName == "Leaf nitrogen (N) content per leaf dry mass") %>% 
  filter( OriglName != "N_senesced_leaf") %>%
  #filter( OriglName == "Leaf %N") %>%
  filter(OrigValueStr < 6) #%>%

#yep = full %>% filter(taxonID %in% c("PIEN"))
yep = yep[complete.cases(yep), ]
#yep  = quantile(yep$nitrogenPercent,probs=c(.025,.975))
TRY %>% colnames()

TRYplot = TRY %>% group_by(LON_site,LAT_site) %>% summarize_if(is.numeric, mean)
colnames(TRYplot)[c(1,2,9)] = c("LON", "LAT", "nitrogenPercent.Estimate")
TRYplot$dataset = "TRY"
cplot = yep %>% group_by(LAT, LON) %>%  summarize_if(is.numeric, mean)
TRYplot = full_join(TRYplot, cplot)
TRYplot$dataset[is.na(TRYplot$dataset)] = "Md"
ggplot(TRYplot, aes(x = LAT, y = LON, color = nitrogenPercent.Estimate, shape = dataset))+geom_point(size = 1.5, alpha=0.5)


ggplot(data = world) +
  xlim(-96,-65)+ ylim(24,50)+
  theme_bw() +
  geom_point(data = yep, aes(x = LON, y = LAT, 
                             fill = cut(yep$nitrogenPercent.Estimate,quantile(yep$nitrogenPercent.Estimate))), alpha = 0.6, 
              size = 0.7, shape = 23, stroke = 0) +  scale_fill_viridis_d()+
  geom_sf(alpha = 0)+
  geom_point(data = N, aes(x = LON_site, y = LAT_site, color = "red", alpha = 0.9), 
             size = 1, shape =4, stroke = 1) + theme(
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()
             )+
  geom_point(data = neon, aes(x = itcLongitude, y = itcLatitude, color = "red", alpha = 0.9), 
           size = 1, shape =4, stroke = 1) +theme(legend.position="bottom") #+   
  #scale_colour_gradient("red")

neon_95 = quantile(neon$nitrogenPercent,probs=c(0.005,.995))
model_95 = quantile(yep$nitrogenPercent.Estimate,probs=c(0.005,.995))
try_95 = quantile(N$OrigValueStr,probs=c(0.005,.995))
sp_95 =  quantile(sp$nitrogenPercent.Estimate,probs=c(0.005,.995))


model_range = cbind.data.frame(round(min(yep$nitrogenPercent.Estimate),2), round(mean(yep$nitrogenPercent.Estimate), 2),
                round(max(yep$nitrogenPercent.Estimate),2), "Model")
sp_range = cbind.data.frame(round(min(sp$nitrogenPercent.Estimate),2), round(mean(sp$nitrogenPercent.Estimate), 2),
                                          round(max(sp$nitrogenPercent.Estimate),2), "Model")
try_range = cbind.data.frame(round(min(N$OrigValueStr), 2), round(mean(N$OrigValueStr), 2), 
              round(max(N$OrigValueStr), 2), "TRY")
neon_range = cbind.data.frame(round(min(neon$nitrogenPercent),2), round(mean(neon$nitrogenPercent), 2),
               round(max(neon$nitrogenPercent),2), "NEON")

model_range = cbind.data.frame(round(min(model_95),2), round(mean(yep$nitrogenPercent.Estimate), 2),
                               round(max(model_95),2), "Combined-Model")
try_range = cbind.data.frame(round(min(try_95),2), round(mean(N$OrigValueStr), 2), 
                             round(max(try_95),2), "TRY")
neon_range = cbind.data.frame(round(min(neon_95),2), round(mean(neon$nitrogenPercent), 2),
                              round(max(neon_95),2), "NEON")
sp_range = cbind.data.frame(round(min(sp_95),2), round(mean(sp$nitrogenPercent.Estimate), 2),
                               round(max(sp_95),2), "Species-Only-Model")



colnames(sp_range)<- colnames(model_range)<- colnames(try_range)<- colnames(neon_range)<- c("min", "mean", "max", "dataset")

ranges = rbind.data.frame(sp_range, model_range, neon_range, try_range)

ggplot(ranges, aes(x = factor(dataset, levels = c("Species-Only-Model", "Combined-Model", "NEON", "TRY")), 
                   y = mean, ymin = min, ymax = max)) + 
  coord_flip()+ theme_bw()+  
  scale_y_continuous(breaks = seq(0.5, 3.9, 0.3))+
  geom_linerange() + 
  geom_pointrange()


model_range = c(round(min(yep$ligninPercent.Estimate),2), round(mean(yep$ligninPercent.Estimate), 2),
                round(max(yep$ligninPercent.Estimate),2), "Model")
try_range = c(round(min(FAGR_l$OrigValueStr)), round(mean(FAGR_l$OrigValueStr), 2), 
              round(max(FAGR_l$OrigValueStr),2), "TRY")


ranges = rbind.data.frame(model_range, try_range)
ranges[1:3] <- as.numeric(ranges[1:3])
colnames(ranges)<- c("min", "mean", "max", "dataset")
ggplot(ranges, aes(x = dataset, y = mean, ymin = min, ymax = max)) + 
  coord_flip()+ theme_bw()+
  geom_linerange() + 
  geom_pointrange() +
  scale_x_continuous(breaks = seq(1, 3.5, 0.3))
  

df <- data.frame(
  trt = factor(c("Model", "TRY")),
  resp = c(round(mean(yep$ligninPercent.Estimate), 2), round(mean(FAGR_l$OrigValueStr), 2)),
  upper = c(round(max(yep$ligninPercent.Estimate), 2), round(max(FAGR_l$OrigValueStr),2)),
  lower = c(round(min(yep$ligninPercent.Estimate),2), round(min(FAGR_l$OrigValueStr), 2))
)

p <- ggplot(df, aes(trt, resp))
p + geom_linerange(aes(ymin = lower, ymax = upper)) +geom_point() +theme_bw()+coord_flip()

df <- data.frame(
  trt = factor(c("Model", "TRY")),
  resp = c(round(mean(yep$nitrogenPercent.Estimate), 2), round(mean(FAGR_n$OrigValueStr), 2)),
  upper = c(round(max(yep$nitrogenPercent.Estimate), 2), round(max(FAGR_n$OrigValueStr),2)),
  lower = c(round(min(yep$nitrogenPercent.Estimate),2), round(min(FAGR_n$OrigValueStr), 2))
)

p <- ggplot(df, aes(trt, resp))
p + geom_linerange(aes(ymin = lower, ymax = upper)) +geom_point() +theme_bw()+coord_flip()
