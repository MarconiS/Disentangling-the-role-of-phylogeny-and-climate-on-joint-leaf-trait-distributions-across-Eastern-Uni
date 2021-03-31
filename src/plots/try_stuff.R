NEtry = read_csv("./indir/TRY_species_sample.csv")
NEtry$SpeciesName %>% table %>% 
  sort(decreasing = T) %>% head(40)

sp = NEtry %>% filter(SpeciesName == "Acer rubrum")
sp$OrigValueStr = as.numeric(sp$OrigValueStr) 
N = sp %>% filter(TraitName == "Leaf nitrogen (N) content per leaf dry mass") %>% 
  filter( OriglName != "N_senesced_leaf") %>%
  #filter( OriglName == "Leaf %N") %>%
  filter(OrigValueStr < 6) #%>%
C = sp %>% filter(TraitName == "Leaf carbon (C) content per leaf dry mass") %>% 
  mutate(OrigValueStr = OrigValueStr/10) #%>%
LMA = sp %>% filter(TraitName == "Leaf area per leaf dry mass (SLA or 1/LMA): undefined if petiole and rachis are in- or excluded") %>% 
  mutate(OrigValueStr = OrigValueStr) #%>%
LMA$OrigValueStr = 1 / LMA$OrigValueStr *10000
idLMA = LMA %>% group_by(LON_site, LAT_site) %>% summarise_if(is.numeric, mean)
idN = N %>% group_by(LON_site, LAT_site) %>% summarise_if(is.numeric, mean)
idC = C %>% group_by(LON_site, LAT_site) %>% summarise_if(is.numeric, mean)
try_obs = inner_join(idLMA, idN, by=c("LON_site", "LAT_site"))
plot(x=log(try_obs$OrigValueStr.x), y=log(try_obs$OrigValueStr.y))
try_obs = inner_join(idC, idN, by=c("LON_site", "LAT_site"))
plot(x=log(try_obs$OrigValueStr.x), y=log(try_obs$OrigValueStr.y))
try_obs = inner_join(idC, idLMA, by=c("LON_site", "LAT_site"))
plot(x=log(try_obs$OrigValueStr.x), y=log(try_obs$OrigValueStr.y))

