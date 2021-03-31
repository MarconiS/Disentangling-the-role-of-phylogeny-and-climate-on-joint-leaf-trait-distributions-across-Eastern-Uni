full = readr::read_csv("./outdir/2020/FIA_2020_no_mlbs_full.csv")
env = readr::read_csv("./outdir/2020/FIA_2020_no_mlbs_env.csv")
phylo = readr::read_csv("./outdir/2020/FIA_2020_no_mlbs_sp.csv")
clim = readr::read_csv("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/2020_FIA_features.csv")
clim = sf::st_as_sf(clim, coords = c("LON", "LAT"), crs = 4326)
eco = sf::read_sf("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/eco-us-shp/eco_us.shp")
eco = sf::st_transform(eco, crs = sf::st_crs(clim))
ecoregions = sf::read_sf("/Users/sergiomarconi/Documents/Data/Data_products/Chapter3_product/atlantic/atlantic_outer.shp")
sum95 = function(x){
  x = x[x < quantile(x, 0.95), ]
  x = x[x > quantile(x, 0.05), ]
  sum(x)
}
#average values per plot
# full[-c(1:4)] = full[-c(1:4)] %>% group_by(statecd,unitcd,countycd, plot) %>% summarize_if(is.numeric, sum)
# env = env %>% group_by(statecd,unitcd,countycd, plot) %>% summarize_if(is.numeric, sum)
# phylo = phylo %>% group_by(statecd,unitcd,countycd, plot) %>% summarize_if(is.numeric, sum)

library(sf)
library(tidyverse)
foo = full %>% dplyr::select(statecd,unitcd, countycd,  plot)
foo = apply(foo, 1, function(x) paste(x, collapse = "_"))
full["identifierID"] = foo
env$identifierID = full$identifierID
phylo$identifierID = full$identifierID

# foo = env %>% select(statecd, countycd, unitcd, plot)
# foo = apply(foo, 1, function(x) paste(x, collapse = "_"))
# env["indentifierID"] = foo
# 
# foo = phylo %>% select(statecd, countycd, unitcd, plot)
# foo = apply(foo, 1, function(x) paste(x, collapse = "_"))
# phylo["indentifierID"] = foo


full_ave  = full %>% group_by(identifierID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
# 
# env_ave = env %>% group_by(identifierID) %>% 
#   summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
# 
# phylo_ave   = phylo %>% group_by(identifierID) %>%
#   summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

denv = (full[9:24] - env[9:24]) 
dphylo  = (full[9:24] - phylo[9:24])

denv = cbind.data.frame( full[25], full[1:8], denv)
dphylo = cbind.data.frame(full[25], full[1:8], dphylo)

#full[c("nitrogenPercent.Estimate", "carbonPercent.Estimate" , "ligninPercent.Estimate" , "cellulosePercent.Estimate")]=(test_res)


denv_plot = denv %>% group_by(identifierID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

dphylo_plot = dphylo %>% group_by(identifierID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
clim = sf::st_join(clim, eco) 
clim=clim %>% dplyr::select(identifierID, daylength,prec,rad, snow_melt,
                            tmax,tmin,vp,elevation,ope,ect, ECO_US_ID, ECOCODE, DOMAIN, PROVINCE) %>%unique
clim = clim %>% group_by(identifierID) %>% slice(1)
denv_plot = left_join(denv_plot, clim[-c(2:9)], by = "identifierID")
dphylo_plot = left_join(dphylo_plot, clim[-c(2:9)],  by = "identifierID")


# denv = st_as_sf(denv, coords = c("LON", "LAT"), crs =4326)
# dphylo = st_as_sf(dphylo, coords = c("LON", "LAT"), crs =4326)
full_sf = st_as_sf(full_ave, coords = c("LON", "LAT"), crs =4326)
denv_sf = st_as_sf(denv, coords = c("LON", "LAT"), crs =4326)
dephylo_sf = st_as_sf(dphylo, coords = c("LON", "LAT"), crs =4326)

# topo = st_as_sf(topo, coords = c("LON", "LAT"), crs =4326)
ecoregions = st_transform(ecoregions,  crs = 4326)
full_sf = st_join(full_sf, eco)

sd_clim = clim %>% group_by(ECOCODE) %>%  
  summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))

sd_traits = full_sf %>% group_by(ECOCODE) %>%
  summarise(across(where(is.numeric), ~ sd(.x, na.rm = TRUE)))


denv_sf = st_join(denv_sf, ecoregions)
dephylo_sf = st_join(dephylo_sf, ecoregions)

denv_sf$LAT = st_coordinates(denv_sf)[,2]
dephylo_sf$LAT = st_coordinates(dephylo_sf)[,2]
full_sf$LAT = st_coordinates(full_sf)[,2]

#plot distribution of traits across ecoprovinces
check = full_sf %>% data.frame %>% dplyr::select("PROVINCE", 
                         "leafMassPerArea.Estimate", "nitrogenPercent.Estimate",
                         "carbonPercent.Estimate", "ligninPercent.Estimate", 
                         "cellulosePercent.Estimate", "chlApercdw.Estimate", 
                         "chlBpercdw.Estimate", "carotpercdw.Estimate")


colnames(check) = c("PROVINCE", "LMA", "rsN", "rsC", "rslignin", "rsCell", "carot", "chlA", "chlB")
check=check[complete.cases(check),]
check$PROVINCE = gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',check$PROVINCE,perl = TRUE)
for(ii in 3:9){
  tmp = check[,c(ii, 2)]
  colnames(tmp) = c("y", "lma")
  md = lm(y ~ lma, data = log10(tmp))
  check[,ii] = (residuals(md))
}

#difference in traits distribution?
p_val = pair_nm = dunn_tr=list()
for(ii in 2:9){
  dunn_tst = conover.test::conover.test(check[,ii], check[,1],method = "bh")
  p_val[[ii]] = dunn_tst$P.adjusted
  pair_nm[[ii]] = dunn_tst$comparisons
  dunn_tr[[ii]] = rep(colnames(check)[ii], length(dunn_tst$P.adjusted))
}
dunn_res = cbind.data.frame(do.call(c, p_val), do.call(c, pair_nm), do.call(c, dunn_tr))
colnames(dunn_res) = c("pval", "couple", "trait")
dunn_res$signi = dunn_res$pval < 0.025
write_csv(dunn_res, "./outdir/outdir/provice_conover_significance.csv")
# # check if same results with conover test
# lma_diff = conover.test::conover.test(check$LMA, check$PROVINCE,method = "bh")

check = reshape2::melt(check)
provinces_order = c("EP", "OCPMFP","LMRFP", "SMFP", "OMFMP", "OBFMP", "PPTP", 
  "EBFCP", "EBFOP", "CABFCFMP","LMFP", "ANEMFCFAMP")
check = check %>% filter(!is.na(PROVINCE))
check$PROVINCE = factor(check$PROVINCE, levels = provinces_order)
ggplot(data = check, aes(x=PROVINCE, y=value))+geom_violin()+
  facet_grid(variable~., scale = "free")+theme_bw()+ theme(legend.position = "none",
                              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) 

full_sf = st_as_sf(full, coords = c("LON", "LAT"), crs =4326)
east_env = denv_sf %>% filter(!is.na(ECOCODE)) %>% group_by(statecd, unitcd, countycd) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

east_phy = dephylo_sf %>% filter(!is.na(ECOCODE))%>% group_by(statecd, unitcd, countycd) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

east_taxa = full %>% dplyr::select(statecd, unitcd, countycd, LAT, taxonID)
east_taxa = sf::st_join(full_sf, ecoregions)
east_taxa = east_taxa %>% filter(!is.na(PERIMETER))
eat_taxa = east_taxa %>%
  group_by(statecd, unitcd, countycd,taxonID) %>% 
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

#east_taxa = left_join(data.frame(eat_taxa)[-6], data.frame(east_taxa)[1:4]) %>% unique
eat_latLon= st_coordinates(sf::st_centroid(eat_taxa))
colnames(eat_latLon) = c("LON", "LAT")
atl_freq = data.frame(eat_taxa)
atl_freq = cbind.data.frame(atl_freq, eat_latLon)
atl_lat = data.frame(east_taxa) %>% group_by(statecd, unitcd, countycd) %>% 
  summarize_if(is.numeric, mean)
atl_freq = left_join(atl_freq, atl_lat) 
atl_freq$genusID = substr(atl_freq$taxonID, start = 1, stop = 2)
atl_freq$genusID[atl_freq$taxonID %in% c("PIPU", "PIRU")] = "PI2"
atl_freq$genusID[atl_freq$taxonID %in% c("JUCI", "JUNI")] = "JU2"
atl_freq$genusID[atl_freq$taxonID %in% c("LIST2")] = "LI2"

dat = atl_freq %>% filter(freq > 0.1) %>% group_by(LAT, genusID) %>% summarize_all(mean)
genuses_ok = dat$genusID %>% table %>% data.frame(stringsAsFactors = F) 
colnames(genuses_ok)[1]="genusID"
genuses_ok=genuses_ok[order(genuses_ok$Freq, decreasing = T),]
genuses_ok = as.character(genuses_ok[1:9,1])
dat = dat %>% filter(genusID %in% c(genuses_ok, "TA"))
# dat$group = "tropical"
# dat$group[dat$LAT >29] = "mixed_pines"
# dat$group[dat$LAT >35] = "mixed_broad"
# dat$group[dat$LAT >38] = "broad"
# dat$group[dat$LAT >42.5] = "mixed-cold"
# dat$group = factor(dat$group, levels = c("tropical", "mixed_pines", "mixed_broad", "broad", "mixed-cold"))
ggplot(dat, aes(x = LAT, y = freq, color = genusID))+geom_point(alpha=0.0) +
  geom_smooth(se=F, method="loess", size = 0.8) + xlim(25,48)+
  theme_bw() + #facet_wrap(.~group, ncol = 5, nrow = 1, scales = "free_x") + 
  theme(legend.position = "bottom",  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  ylim(0.00,0.6)#scale_color_manual(values=as.vector(polychrome())) + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F, size = 0.8)
#east_env$LAT = st_coordinates(east_env)[,2]
east_env$model = "environment"
#east_phy$LAT = st_coordinates(east_phy)[,2]
east_phy$model = "species"

#plot LMA divergence
plot_divergence = rbind.data.frame(east_env[c("leafMassPerArea.Estimate", "LAT", "model")],
                                   east_phy[c("leafMassPerArea.Estimate", "LAT", "model")])
plot_divergence = plot_divergence %>% data.frame %>% dplyr::select(-one_of("geometry"))
ggplot(plot_divergence, aes(x=LAT, y = leafMassPerArea.Estimate, color = model))+ 
  stat_density2d(alpha=1)+scale_color_manual(values = c('orange', "lightgreen")) +theme_bw() +
  theme(legend.position = "bottom",  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  + xlim(25,48)+
  geom_smooth(method="gam") +  geom_hline(yintercept = 0) + ylim(130,-130)

#plot N divergence
plot_divergence = rbind.data.frame(east_env[c("nitrogenPercent.Estimate", "LAT", "model")],
                                   east_phy[c("nitrogenPercent.Estimate", "LAT", "model")])
plot_divergence = plot_divergence %>% data.frame %>% dplyr::select(-one_of("geometry"))
ggplot(plot_divergence, aes(x=LAT, y = nitrogenPercent.Estimate, color = model))+ 
  stat_density2d(alpha=1)+scale_color_manual(values = c('orange', "lightgreen")) +theme_bw() +
  theme(legend.position = "bottom",  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  + xlim(25,48)+
  geom_smooth(method="gam") +  geom_hline(yintercept = 0) #+ ylim(130,-130)

#plot C divergence
plot_divergence = rbind.data.frame(east_env[c("carbonPercent.Estimate", "LAT", "model")],
                                   east_phy[c("carbonPercent.Estimate", "LAT", "model")])
plot_divergence = plot_divergence %>% data.frame %>% dplyr::select(-one_of("geometry"))
ggplot(plot_divergence, aes(x=LAT, y = carbonPercent.Estimate, color = model))+ 
  stat_density2d(alpha=1)+scale_color_manual(values = c('orange', "lightgreen")) +theme_bw() +
  theme(legend.position = "bottom",  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  + xlim(25,48)+
  geom_smooth(method="gam") +  geom_hline(yintercept = 0) #+ ylim(130,-130)

plot_divergence = rbind.data.frame(east_env[c("chlApercdw.Estimate", "LAT", "model")],
                                   east_phy[c("chlApercdw.Estimate", "LAT", "model")])
plot_divergence = plot_divergence %>% data.frame %>% dplyr::select(-one_of("geometry"))
ggplot(plot_divergence, aes(x=LAT, y = chlApercdw.Estimate, color = model))+ 
  stat_density2d(alpha=1)+scale_color_manual(values = c('orange', "lightgreen")) +theme_bw() +
  theme(legend.position = "bottom",  panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  + xlim(25,48)+
  geom_smooth(method="gam") +  geom_hline(yintercept = 0) #+ ylim(130,-130)


denv_p = denv %>% 
  group_by(identifierID) %>%
  summarize_if(is.numeric, mean)
dphi_p = dphylo %>% 
  group_by(identifierID) %>%
  summarize_if(is.numeric, mean)

ggplot(data = world) +
  xlim(-96,-65)+ ylim(24,50)+
  theme_bw() +
  geom_point(data = denv_p, aes(x = LON, y = LAT, fill = leafMassPerArea.Estimate), 
      size = 0.7, shape = 23, stroke = 0) +  
  geom_sf(alpha = 0)+
  scale_fill_gradient2(midpoint=0, low="steelblue", mid="white",
                        high="orange", space ="Lab")+
  geom_sf(data=eco, alpha=0)+
  geom_sf(data=ecoregions, alpha=0, color="red")+
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank())+
  theme(legend.position="bottom") #+   

ggplot(data = world) +
  xlim(-96,-65)+ ylim(24,50)+
  theme_bw() +
  geom_point(data = dphi_p, aes(x = LON, y = LAT, fill = leafMassPerArea.Estimate), 
             size = 0.7, shape = 23, stroke = 0) +  
  geom_sf(alpha = 0)+
  scale_fill_gradient2(midpoint=0, low="steelblue", mid="white",
                       high="green", space ="Lab")+
  geom_sf(data=eco, alpha=0)+
  geom_sf(data=ecoregions, alpha=0, color="red")+
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank())+
  theme(legend.position="bottom") #+   
# 
# write_csv(denv, "~/Documents/Data/Chapter3/delta_env.csv")
# write_csv(dphylo, "~/Documents/Data/Chapter3/delta_dphylo.csv")
# 
# den_plot = data.frame(denv) %>% group_by(unitcd, statecd, countycd, plot) %>% summarize_if(is.numeric, meannona)
# dphylo_plot = dphylo %>% group_by(unitcd, statecd, countycd, plot) %>% summarize_if(is.numeric, meannona)
# 
# ave_env = denv_plot %>% group_by(statecd, unitcd) %>% summarize_if(is.numeric, mean)
# sd_ave = denv_plot %>% group_by(statecd, unitcd) %>% summarize_if(is.numeric, sd)
# ave_phy = dphylo_plot %>% group_by(statecd, unitcd) %>% summarize_if(is.numeric, mean)
# sd_phy = dphylo_plot %>% group_by(statecd, unitcd) %>% summarize_if(is.numeric, sd)
# 
# heteroenv = cbind.data.frame(ave_env[1:31], sd_ave[32:43])
# heterophy = cbind.data.frame(ave_phy[1:31], sd_phy[32:43])
# 
# heteropy = cbind.data.frame(sd_clim, sd_traits)
# check = heteropy %>% select(nitrogenPercent.Estimate, carbonPercent.Estimate,
#                          extractCarotConc.Estimate, ligninPercent.Estimate, 
#                          cellulosePercent.Estimate, leafMassPerArea.Estimate,
#                          prec,rad,tmax,tmin,vp,elevation,ope,ect) %>% data.frame
# 
# colnames(check) = c("Nm", "Cm",  "Carot", "lign", "cell", "LMA", 
#                     "Prec", "Rad", "Tmax", "Tmin", "Vp", "Elev", "Slope", "Aspect")
# 
# 
# #large to long
# #check = scale(check)
# check = reshape2::melt(data.frame(check), id.vars = c("Prec", "Rad", "Tmax", "Tmin", "Vp", "Elev", "Slope", "Aspect"))
# check = reshape2::melt(check, id.vars = c("variable", "value"), value.name = "feature")
# colnames(check)= c("trait", "trait_value", "env", "env_value")
# check = check[complete.cases(check),]
# 
# lm_eqn <- function(df){
#   m <- lm(trait_value ~ env_value, df);
#   eq <- substitute(~~italic(r)^2~"="~r2, 
#                    list(r2 = format(summary(m)$r.squared, digits = 2)))
#   as.character(as.expression(eq));
# }
# 
# 
# ggplot(check, aes(y = trait_value, x = env_value, color=))
# 
# 
# 
# 
# env_ch = check %>% filter(trait == "LMA")
# require(plyr)
# eq <- ddply(env_ch,.(env),lm_eqn)
# ggplot(env_ch, aes(x = env_value, y = trait_value))+  facet_wrap(env~., scales = "free")  +
#   geom_point(alpha=0.5, size=0.6) + ylim(0,60)+
#   geom_text(data=eq,aes(x = 0.3, y = 55,label=V1), parse = TRUE, inherit.aes=FALSE) + facet_wrap(env~.,  scales = "free")+
#   geom_smooth(method = "lm", col = "black") + theme_bw()
#   
# ggplot(check, aes(x = env_value, y = trait_value))+ # facet_wrap(trait~., scales = "free")  +
#   geom_point(alpha=0.2, size=0.3, col = "grey") +
#   geom_smooth(method = "lm") + theme_bw()




# ggplot((check), aes(x = env_value, y = (trait_value), color = trait))+  facet_wrap(env~., scales = "free")  +
#   geom_smooth(method = "lm") + ylim(-3,3)+
#   geom_point(alpha=0) 
# 
# ggplot(check, aes(x = env_value, y = trait_value, color = env))+  facet_wrap(trait~., scales = "free")  +
#   geom_smooth(method = "lm") +
#   geom_point(alpha=0) 