#marginal distribution of temperature
library(brms)
library(tidyverse)
library(sf)
model = readRDS("./models/ch2perc_gauss_full_noml.rds")


predictions = data.table::fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_full.csv")
taxa_per_plot = predictions %>% group_by(LAT, LON, taxonID) %>% summarize_if(is.numeric, mean)

# quantify intraspecies 98% ranges
range98_N = predictions %>% group_by(taxonID)%>%
  mutate(low_bound = quantile(nitrogenPercent.Estimate, 0.01),
         up_bound = quantile(nitrogenPercent.Estimate, 0.99)) %>%
  select(taxonID, low_bound, up_bound) %>% unique

eco = sf::read_sf("/Volumes/Data/Metadata/us_eco_l3/us_eco_l3.shp")
sp_presence = sf::st_as_sf(taxa_per_plot, coords = c("LON", "LAT"), crs = 4326)
eco = sf::st_transform(eco, sf::st_crs(sp_presence))
sp_presence = sf::st_join(sp_presence, eco)
sp_spread = sp_presence %>% group_by(taxonID, NA_L3NAME) %>% slice(1)
sp_spread = sp_spread$taxonID %>% table %>% data.frame
colnames(sp_spread) = c("taxonID", "n_regions")
range98_N = inner_join(range98_N, sp_spread)

#plot intraN vs n regions
ggplot(range98_N, aes(x = n_regions, y = range))+
  geom_point()+theme_bw()

#conditions <- make_conditions(model$mod, "taxonID")
# 
# pred_ = taxa_per_plot
# pred_nlma = lm(nitrogenPercent.Estimate~leafMassPerArea.Estimate, taxa_per_plot)
# taxa_per_plot[["nitrogenPercent.Estimate"]] = predict(pred_nlma, taxa_per_plot)

#append tmax
climate = data.table::fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/2020f_FIA_features.csv")
climate2 = data.table::fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/climate_plts_east_fia.csv")

climate = climate %>% group_by(LAT, LON) %>% summarize_if(is.numeric, mean)
climate[["identifierID"]] = paste(climate$statecd, climate$unitcd, climate$countycd, climate$plot, sep="_")
lat_lon = climate #%>% select(identifierID, LAT, LON)
climate2 = left_join(climate2, lat_lon)
climate2 = climate2 %>% group_by(LAT, LON) %>% summarize_if(is.numeric, mean)


p_value = p_value

taxa_per_plot = left_join(taxa_per_plot, climate)
taxa_per_plot = left_join(taxa_per_plot, climate2)


which_taxa_out = taxa_per_plot$taxonID %>% table %>% data.frame
which_taxa_out = which_taxa_out %>% filter(Freq <30)
tx_ = taxa_per_plot %>% filter(!taxonID %in% which_taxa_out$.)
# classic plot :
tx_$max_tmax = apply(tx_[,86:97], 1, max)
tx_$mean_tmax = apply(tx_[,86:97], 1, mean)
tx_$min_tmax = apply(tx_[,86:97], 1, min)

pred_ = tx_
pred_nlma = lm(nitrogenPercent.Estimate~leafMassPerArea.Estimate, tx_)
tx_[["N_not_LMA"]] = predict(pred_nlma, tx_)

traits = colnames(tx_)[seq(9,24,by=2)]
features = colnames(tx_)[28:34]
#calculate p-value for each species
for(tr in traits){
  for(env_var in features){
    tx_ = taxa_per_plot %>% filter(!taxonID %in% which_taxa_out$.)
    # classic plot :
    tr_cols = tx_ %>% select(contains(env_var))
    tx_[["mean_tr"]] = apply(tr_cols, 1, mean)
    
    p_tmax = corr_tmax=pdlma_tmax = corrdlma_tmax  = list()
    for(taxa in unique(tx_$taxonID)){
      tmp = tx_ %>% filter(taxonID  == taxa)
      lm_for = as.formula(paste(tr, "~ mean_tr",sep = ""))
      linmod = summary(lm(lm_for, data = tmp))
      
      p_tmax[[taxa]] = linmod$coefficients[2,4]
      corr_tmax[[taxa]] = linmod$r.squared 
    }
    p_tmax = do.call(rbind, p_tmax)
    corr_tmax = do.call(rbind, corr_tmax)
    #get significant species
    p_tmax = cbind.data.frame(rownames(p_tmax), p_tmax, corr_tmax)
    colnames(p_tmax)[1]="taxonID"
    tx_ = left_join(tx_, p_tmax)
    tx_.sig = tx_ %>% filter(p_tmax < 0.001)
    #plot relation for ith traits env combination
    p<-ggplot()
    for(i in unique(tx_$taxonID)){ 
      tmp = tx_ %>% ungroup %>% filter(taxonID == i)
      tmp = tmp %>% 
        filter(mean_tr > quantile(tmp$mean_tr, 0.01, na.rm=T)) %>%
        filter(mean_tr < quantile(tmp$mean_tr, 0.99,na.rm = T))
      cor_taxa = tmp$corr_tmax %>% unique
      p<-p+stat_smooth(data = tmp, geom='line',method = 'loess',  se = T, 
                       aes_string(x="mean_tr", y=tr),
                       alpha = cor_taxa) 
    }
    p + theme_bw()+# + xlim(15, 32) + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
    ggsave(paste("./plots/taxa_effects/", tr, env_var, ".png", sep = ""), 
           width = 5, height = 5, units="cm")
  }
  
}

p_tmax = do.call(rbind, p_tmax)
corr_tmax = do.call(rbind, corr_tmax)
pdlma_tmax = do.call(rbind, pdlma_tmax)
corrdlma_tmax = do.call(rbind, corrdlma_tmax)

sum(p_tmax > 0.001)/nrow(p_tmax)
sum(pdlma_tmax > 0.001)/nrow(pdlma_tmax)

p_tmax = cbind.data.frame(rownames(p_tmax), p_tmax, corr_tmax, pdlma_tmax, corrdlma_tmax)
colnames(p_tmax)[1]="taxonID"
tx_ = left_join(tx_, p_tmax)
tx_.sig = tx_ %>% filter(p_tmax < 0.001)
# range01 <- function(x){(x-min(x))/(max(x)-min(x))}
# 
# ggplot(tx_.sig, aes(x=tmax, y=nitrogenPercent.Estimate, group=taxonID)) +
#   #geom_point(size=0.1) +
#   stat_smooth(stat="smooth",method = "loess", se=F,
#             #linetype ="dashed",
#             aes(alpha = range01(abs(corr_tmax)))) + theme_bw()+
#   theme(legend.position="none") 

df = tx_.sig %>% ungroup %>% select(taxonID, tmax, max_tmax, min_tmax, mean_tmax,
                                    nitrogenPercent.Estimate, p_tmax, corr_tmax,
                                    N_not_LMA, pdlma_tmax, corrdlma_tmax)

plotAllLayers<-function(df){
  p<-ggplot()
  df = foo
  #df$corr_tmax =  range01(df$corr_tmax)
  for(i in unique(df$taxonID)){ 
    tmp = df %>% filter(taxonID == i)
    cor_taxa = tmp$corr_tmax %>% unique
    p<-p+stat_smooth(data = tmp, geom='line',method = 'loess',  se = T, 
                     aes(x=mean_tmax, y=nitrogenPercent.Estimate),
                     alpha = cor_taxa) 
  }
  p + theme_bw() + xlim(15, 32) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    stat_smooth(data = df, geom='line',method = 'loess',  se = T, 
                aes(x=mean_tmax, y=nitrogenPercent.Estimate), color="orange")
  return(p)
}

plotAllLayers(df)
# marginal = ce$nitrogenPercent.nitrogenPercent_tmax %>% data.frame
# coeffs_for_taxa = tx_.sig %>% ungroup %>% select(taxonID, p_tmax, corr_tmax) %>% unique
# marginal = inner_join(marginal, coeffs_for_taxa)
# colnames(marginal)[20] = "nitrogenPercent.Estimate"
# plotAllLayers(marginal)
# ggplot(marginal, aes(x = effect1__, y = nitrogenPercent.Estimate, color = cond__, alpha = 0.3)) + 
#   geom_point() + theme(legend.position = "none")

# 
# p2 <- ggMarginal(p, type="density", margMapping = aes(colour = taxonID))
# ggMarginal(p = p, marginalGroup = list(colourGroup = TRUE, 
#                                        colourAlpha = .4, fillGroup = FALSE, fillAlpha = NA))
# 
# ThaiEdu_New %>%
#   data_grid(SEX, PPED) %>%
#   add_fitted_draws(Bayes_Model_Binary) %>%
#   ggplot(aes(x = .value, y = interaction(SEX, PPED))) +
#   stat_pointintervalh(.width = c(.68, .95)) +
#   coord_flip() +
#   xlab("predicted probability") +
#   scale_x_continuous(breaks = seq(0, 0.24, 0.02))
# 
# 
# 
# conditional_effects(model$mod,
#                     conditions = conditions,
#                     effects = "tmax",
#                     nsamples = 200) %>% 
#   plot(points = T,
#        point_args = c(alpha = 1/2, size = 1),
#        line_args = c(colour = "black"))
# 
# OSFDP_17ha_4_28_2020$sp %>% table
# 
full = data.table::fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_full.csv")
neon_field = data.table::fread("~/Documents//Data/Inventories/NEON/final_field_dataset.csv")  
taxa_n_plots = full %>%
  group_by(taxonID, LAT, LON) %>%slice(1)%>%
  ungroup %>%
  select(taxonID) %>% table %>% data.frame

colnames(taxa_n_plots) = c("taxonID", "plots_detected")
#calculate range of intraspecies variation in %N
N_perc_range = full %>% group_by(taxonID) %>%
  mutate(nmin = quantile(nitrogenPercent.Estimate, 0.005), 
         nmedian = median(nitrogenPercent.Estimate), 
         nmax = quantile(nitrogenPercent.Estimate, 0.995),
         taxa_frequency = n()) %>%
  select(taxonID, nmin, nmedian, nmax, taxa_frequency) %>% unique

field_n = neon_field %>%
  mutate(nmin = min(nitrogenPercent), 
         nmedian = median(nitrogenPercent), 
         nmax = max(nitrogenPercent)) %>%
  select(nmin, nmedian, nmax) %>% unique

ggplot(N_perc_range, aes(x = factor(taxonID), y = nmedian)) + 
  geom_pointrange(aes(ymin = nmin, ymax = nmax)) + coord_flip() + 
  geom_hline(yintercept = field_n$nmin) + 
  geom_hline(yintercept = field_n$nmedian) + 
  geom_hline(yintercept = field_n$nmax) 

N_perc_range = N_perc_range %>% group_by(taxonID) %>%
  mutate(nrange = nmax - nmin)

N_perc_range$fraction = N_perc_range$nrange / (field_n$nmax - field_n$nmin)

N_perc_range = inner_join(N_perc_range, taxa_n_plots)
hist(N_perc_range$fraction, breaks = 30)
