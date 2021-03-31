library("rnaturalearth")
library("rnaturalearthdata")
library(ggpubr)
world <- ne_countries(scale = "medium", returnclass = "sf")
trait = "leafMassPerArea.Estimate"
full <- data.table::fread("./outdir/FIA_2020_no_mlbs_full.csv")
taxa = readr::read_csv("/Users/sergiomarconi/Documents/Data/Metadata/list_of_species.csv")
ecoregions = sf::read_sf("/Users/sergiomarconi/Documents/Data/eco-us-shp/eco_us.shp")

colnames(taxa)[1]="taxonID"
full = left_join(full,taxa )
full = full %>% filter(pft %in% c("d", "e", "nn")) %>% filter(leafMassPerArea.Estimate < 600)
full = full[complete.cases(full),]
eco_full = sf::st_as_sf(full, coords = c("LON", "LAT"), crs = 4326)

rs_corr = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["nitrogenPercent.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(carbonPercent.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["carbonPercent.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(ligninPercent.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["ligninPercent.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(cellulosePercent.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["cellulosePercent.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(chlApercdw.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["chlApercdw.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(chlBpercdw.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["chlBpercdw.Estimate"] = residuals(rs_corr)
rs_corr = lm(log10(carotpercdw.Estimate)~log10(leafMassPerArea.Estimate), data=(full))
full["carotpercdw.Estimate"] = residuals(rs_corr)

ecoregions = sf::st_transform(ecoregions, sf::st_crs(eco_full))
provinces = ecoregions[8] %>% unique
eco_full = sf::st_join(eco_full, ecoregions)
dat = full %>% filter(pft == "EBF")
tmp = full %>% filter(leafMassPerArea.Est.Error < quantile(leafMassPerArea.Est.Error, 0.975),
                      leafMassPerArea.Est.Error > quantile(leafMassPerArea.Est.Error, 0.25))
ggplot(data = world) +
  xlim(-95, -65)+ ylim(25,50)+
  theme_bw() +
  geom_point(data = tmp, aes(x = LON, y = LAT,
            fill = leafMassPerArea.Est.Error),
             alpha = 0.9,  size = 0.8, shape = 23, stroke = 0) +  
            scale_fill_viridis_c()+ 
            theme(legend.position = "bottom") +
  geom_sf(alpha = 0)

ggplot(data = ecoregions) +
  #xlim(-95, -65)+ ylim(25,50)+
  theme_bw() +
  geom_sf(aes(fill = PROVINCE))+   theme(legend.position = "bottom") 

dat = full %>% data.frame %>% dplyr::select(pft, LAT, LON, nitrogenPercent.Estimate, carbonPercent.Estimate, chlApercdw.Estimate,
                                     chlBpercdw.Estimate, carotpercdw.Estimate, leafMassPerArea.Estimate, ligninPercent.Estimate, cellulosePercent.Estimate)
dat$type = "Broadleaf"
dat$type[dat$pft=="e"] = "Needleleaf"
#dat = reshape2::melt(dat, id.vars = c("PFT", "LAT", "LON"))
#colnames(dat) = c("PFT", "LAT", "LON", "Trait", "Estimate")
N = ggplot(data = dat, aes(x= LAT, y = nitrogenPercent.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
C = ggplot(data = dat, aes(x= LAT, y = carbonPercent.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
lig = ggplot(data = dat, aes(x= LAT, y = ligninPercent.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
cell = ggplot(data = dat, aes(x= LAT, y = cellulosePercent.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
LMA = ggplot(data = dat, aes(x= LAT, y = leafMassPerArea.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
chla = ggplot(data = dat, aes(x= LAT, y = chlApercdw.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
chlb = ggplot(data = dat, aes(x= LAT, y = chlBpercdw.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")
carot = ggplot(data = dat, aes(x= LAT, y = carotpercdw.Estimate,  color = type))+ geom_density_2d(alpha=0.7)+geom_smooth(se = F, linetype="dashed")  +
  theme_bw()+scale_color_jco() + theme(legend.position = "none")

ggarrange(N, chla, chlb, carot, LMA, C, lig, cell + rremove("x.text"), 
          labels = c("A", "B", "C", "D",  "E",  "F", "G", "H"),
          ncol = 4, nrow = 2)
  

standard_sd =  sd_clim %>% select(ECOCODE, prec,rad ,   tmax,   tmin,    vp, elevation, ope, ect)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
standard_sd = data.frame(standard_sd)
standard_sd = standard_sd %>% select(-one_of("geometry"))
standard_sd[-1]= apply(standard_sd[-1], 2, range01)
standard_sd$mean = apply(standard_sd[-c(1,10:11)], 1, mean)

colnames(standard_sd) = c("ECOCODE", "sd_prec","_sdrad","sd_tmax","sd_tmin","sd_vp","sd_elevation", "sd_ope","sd_ect","sd_mean")
eco_full = left_join(eco_full[1:32], standard_sd)
#get regions neighboring mississipi 2032, 2053,
#ggplot(mississipi, aes(y= leafMassPerArea.Estimate, x = factor(ECO_US_)))+ geom_boxplot() + coord_flip() + theme_bw()
eco_full = sf::st_as_sf(full, coords = c("LON", "LAT"), crs = 4326)
eco_full = sf::st_join(eco_full, ecoregions)
mississipi = eco_full %>% filter(ECO_US_ %in% c(2039, 1956, 2948, 2121, 2074, 2015, 2065, 2105,  1878, 1937, 1829))

mississipi = mississipi %>% select(ECO_US_,SECTION,PROVINCE, taxonID, leafMassPerArea.Estimate, 
                                   nitrogenPercent.Estimate, carbonPercent.Estimate, ligninPercent.Estimate, 
                                   cellulosePercent.Estimate,carotpercdw.Estimate, chlApercdw.Estimate, chlBpercdw.Estimate)
dat = data.frame(mississipi) %>% select(-one_of("geometry"))
colnames(dat) = c("ECO_US_","SECTION","PROVINCE", "taxonID", "LMA", "N", "C", "Lignin", "Cellulose", "Carotenoids", "ChlA", "ChlB")
dat[dat$SECTION=="Mississippi Alluvial Basin Section", "SECTION"] = "Mississipi"
dat[dat$SECTION != "Mississipi", "SECTION"]="Neighbors"
species_composition = dat %>% group_by(SECTION) %>% select(taxonID) %>% table 
tot_points = dat %>% group_by(SECTION) %>% select(taxonID) %>% tally
sp_rel_freq = lapply(1:2, function(x)species_composition[x,]/tot_points$n[x])
sp_rel_freq = do.call(rbind.data.frame, sp_rel_freq)
colnames(sp_rel_freq) = colnames(species_composition)
sp_rel_freq["SECTION"] = tot_points$SECTION
sp_rel_freq = reshape2::melt(sp_rel_freq)
which_sp = sp_rel_freq %>% filter(value > 0.01) %>% select(variable)%>% unique
sp_rel_freq2 = sp_rel_freq %>% filter(variable %in% unlist(which_sp))
ggplot(data = sp_rel_freq2, aes(x = (variable), y = value, color = SECTION,  group=SECTION)) + 
  geom_point(aes(alpha=0.5), position=position_dodge(width = 0.30)) + 
   theme_bw() + #stat_summary(fun = sum, geom="line") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

library(ggsci)
mississipi = eco_full %>% filter(ECO_US_ %in% c(2039, 1956, 2948, 2121, 2074, 2015, 2065, 2105,  1878, 1937, 1829))
mississipi= mississipi %>% data.frame %>% group_by(ECO_US_,SECTION,PROVINCE, statecd, plot, countycd, unitcd) %>%
  summarize_if(is.numeric, mean)
dat = mississipi %>% ungroup %>% select(ECO_US_,SECTION,PROVINCE, leafMassPerArea.Estimate, 
                                        nitrogenPercent.Estimate, carbonPercent.Estimate, ligninPercent.Estimate, 
                                        cellulosePercent.Estimate,carotpercdw.Estimate, chlApercdw.Estimate, chlBpercdw.Estimate)

colnames(dat) = c("ECO_US_","SECTION","PROVINCE",  "LMA", "N", "C", "Lignin", "Cellulose", "Carotenoids", "ChlA", "ChlB")

dat[dat$SECTION=="Mississippi Alluvial Basin Section", "SECTION"] = "Mississipi"
dat[dat$SECTION != "Mississipi", "SECTION"]="Neighbors"

dat$ECO_US_ = factor(dat$ECO_US_)
dat$PROVINCE = gsub('\\b(\\pL)\\pL{2,}|.','\\U\\1',dat$PROVINCE,perl = TRUE)
dat$PROVINCE = factor(dat$PROVINCE)
dat$SECTION = factor(dat$SECTION)

p_val = pair_nm = kwt =dunn_tr=list()

for(ii in 4:11){
  dunn_tst = conover.test::conover.test(unlist(dat[ii]), unlist(dat[3]),method = "bh")
  kwt[[ii]]=kruskal.test(unlist(dat[ii]), unlist(dat[2]),method = "bh") 
  p_val[[ii]] = dunn_tst$P.adjusted
  pair_nm[[ii]] = dunn_tst$comparisons
  dunn_tr[[ii]] = rep(colnames(dat)[ii], length(dunn_tst$P.adjusted))
}
dunn_res = cbind.data.frame(do.call(c, p_val), do.call(c, pair_nm), do.call(c, dunn_tr))
colnames(dunn_res) = c("pval", "couple", "trait")
dunn_res$signi = dunn_res$pval < 0.025
dd = dunn_res[!dunn_res$signi,]

dat = reshape2::melt(dat)
dat2 = dat %>% filter(variable %in% c("LMA", "N", "C", "Carotenoids"))

ggplot(data = dat, aes(y = (SECTION), x = value, color = SECTION)) + 
  facet_wrap(variable ~.,ncol=4, nrow=2,scale = "free")+ geom_violin() + 
  #geom_boxplot(position=position_dodge(width = 1))
  theme_bw() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8))

# 
# dat = eco_full %>% filter(! ECO_US_ %in% c(2577))   
# dat = dat %>% 
#   select(ECOCODE, PROVINCE,taxonID, leafMassPerArea.Estimate, nitrogenPercent.Estimate, 
#          carbonPercent.Estimate, carotpercdw.Estimate, chlApercdw.Estimate, chlBpercdw.Estimate, sd_mean, mean)
# dat = dat %>% group_by(PROVINCE) %>% mutate(sd_mean_in_group = range01(sd_mean))
# dat = data.frame(dat) %>% select(-one_of("geometry"))
# dat$ECOCODE = factor(dat$ECOCODE)
# dat = reshape2::melt(dat, measure.vars = 4:7)
# ggplot(data = dat, aes(x = (ECOCODE), y = value, color = PROVINCE, fill = sd_mean_in_group)) + 
#   facet_wrap(variable ~.,ncol=1, nrow=4,scale = "free")+
#   geom_boxplot(position=position_dodge(width = 1)) + theme_bw() + 
#   theme(axis.text.x = element_blank(), legend.position = "none")+
#   scale_color_d3(palette = c("category20"))+
#   theme(axis.text.x = element_text(angle = 45,hjust = 0.8))
# 
# vars = dat %>% group_by(ECOCODE) %>% 
#   summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
# get_variance = dat %>% group_by(ECOCODE) %>% 
#   summarise_all(funs(if(is.numeric(.)) sd(., na.rm = TRUE) else first(.)))
# 
# dat = reshape2::melt(get_variance, measure.vars = 4:7)
# ggplot(data = dat, aes(x = (mean), y = value, color = PROVINCE)) + 
#   facet_wrap(variable ~.,ncol=1, nrow=4,scale = "free")+
#   geom_point(position=position_dodge(width = 1), alpha=1) + theme_bw() + 
#   #geom_smooth(method = "lm", se=F)+
#   theme(axis.text.x = element_blank(), legend.position = "none")+
#   scale_color_d3(palette = c("category20"))+xlim(0,1)+
#   theme(axis.text.x = element_text(angle = 45,hjust = 0.8))
