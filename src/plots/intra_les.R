library(tidyverse)
full <- fread("/Volumes/Data/Derived_surveys/Eastern_US_traits/East_US_traits_estimates/FIA_2020_no_mlbs_full.csv")
pfts = readr::read_csv("/Volumes/Data/Metadata/list_of_species.csv")
pfts_ = pfts$taxonID[pfts$pft %in% c("d","e", "nn")]
full_nores = full %>%
  #filter(leafMassPerArea.Estimate < 600, nitrogenPercent.Estimate < 6)  %>% 
  filter(taxonID %in% pfts_)
##plot LES
inter_NLMA = full_nores %>% group_by(taxonID) %>% summarize_if(is.numeric, mean)
inter_cor = inter_NLMA %>% select(leafMassPerArea.Estimate, nitrogenPercent.Estimate, carbonPercent.Estimate,ligninPercent.Estimate,
                                  cellulosePercent.Estimate,  chlApercdw.Estimate,  chlBpercdw.Estimate, carotpercdw.Estimate)


ggplot(data=full_nores, aes(x = carbonPercent.Estimate,y = nitrogenPercent.Estimate))+ theme_bw() + 
  geom_density2d_filled(size=0.5)+ xlim(43,53)+ylim(0,3.5)+
  #geom_hex(alpha=0.5, size=0.5)+
  geom_smooth(method="lm")#+ ylim(0,5)+ xlim(35,60)


md = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), data = inter_NLMA)
#calculate log-log relationships: N:LMA
obs_range = data.frame(range(inter_NLMA$leafMassPerArea.Estimate))
colnames(obs_range) = "leafMassPerArea.Estimate"
pred_range = data.frame(predict.lm(md, newdata = obs_range, interval="predict"))
inter_NLMA = cbind.data.frame("Inter_species", t(log10(obs_range)), t(pred_range$fit), t(pred_range$lwr), t(pred_range$upr))
colnames(inter_NLMA) = c("taxonID", "x1", "x2", "y1", "y2", "lwr1",  "lwr2", "upr1", "upr2")
inter_cor = cor(inter_cor)
inter_slope = md$coefficients[2]
taxa_list = full_nores$taxonID %>% unique
intra_NLMA = slope =  intra_corr = res_corr=lmds = list()
for(tx in taxa_list){
  tx_full = full_nores %>% filter(taxonID == tx)
  if(nrow(tx_full) >2){
    md = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    lmds[[tx]] = md
    tx_full = tx_full %>% select(leafMassPerArea.Estimate, nitrogenPercent.Estimate, carbonPercent.Estimate,ligninPercent.Estimate,
                                 cellulosePercent.Estimate,  chlApercdw.Estimate,  chlBpercdw.Estimate, carotpercdw.Estimate)
    tx_res = tx_full
    #calculate log-log relationships: N:LMA
    obs_range = data.frame(range(tx_full$leafMassPerArea.Estimate))
    colnames(obs_range) = "leafMassPerArea.Estimate"
    pred_range = predict.lm(md, newdata = obs_range, interval="predict") %>% data.frame
    tmp = cbind.data.frame(tx, t(log10(obs_range)), t(pred_range$fit), t(pred_range$lwr), t(pred_range$upr))
    colnames(tmp) =  c("taxonID", "x1", "x2", "y1", "y2", "lwr1",  "lwr2", "upr1", "upr2")
    intra_NLMA[[tx]] = tmp
    slope[[tx]] = md$coefficients[2]
    md_res = lm(log10(nitrogenPercent.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["nitrogenPercent.Estimate"] = residuals(md_res)
    md_res = lm(log10(carbonPercent.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["carbonPercent.Estimate"] = residuals(md_res)
    md_res = lm(log10(ligninPercent.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["ligninPercent.Estimate"] = residuals(md_res)
    md_res = lm(log10(cellulosePercent.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["cellulosePercent.Estimate"] = residuals(md_res)
    md_res = lm(log10(chlApercdw.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["chlApercdw.Estimate"] = residuals(md_res)
    md_res = lm(log10(chlBpercdw.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)
    tx_res["chlBpercdw.Estimate"] = residuals(md_res)
    md_res = lm(log10(carotpercdw.Estimate)~log10(leafMassPerArea.Estimate), data = tx_full)    
    tx_res["carotpercdw.Estimate"] = residuals(md_res)
    res_corr[[tx]] = cor((tx_res))
    intra_corr[[tx]] = cor((tx_full))
  }
}
taxa = names(slope)
slope = do.call(rbind.data.frame, slope) 
slope = data.frame(taxa,slope)
slope["pft"] = "DBF"
slope$pft[slope$taxa %in% pfts$taxonID[pfts$pft =="e"]] = "ENF"
slope$pft[slope$taxa %in% pfts$taxonID[pfts$pft =="nn"]] = "EBF"
colnames(slope)[2] = "intra_species"
intra_NLMA = do.call(rbind.data.frame, intra_NLMA)
intra_NLMA["pft"] = "DBF"
intra_NLMA$pft[intra_NLMA$taxonID %in% pfts$taxonID[pfts$pft =="e"]] = "ENF"
intra_NLMA$pft[intra_NLMA$taxonID %in% pfts$taxonID[pfts$pft =="nn"]] = "EBF"
ave_intra = (apply(intra_NLMA[-c(1,10)], 2, mean)) %>% t %>% data.frame
colnames(ave_intra) =  c("x1", "x2", "y1", "y2", "lwr1",  "lwr2", "upr1", "upr2")
#plot intra and inter_species 
library(ggpubr)
library(ggsci)
intra_NLMA$pft = factor(intra_NLMA$pft, levels = c("DBF", "ENF", "EBF"))
intra_NLMA$alpha = 1
intra_NLMA$alpha[intra_NLMA$pft =="DBF" ]=0.9
ggplot() +
  geom_segment(data = intra_NLMA, aes(x = x1, y = y1, xend = x2, yend = y2, colour=pft, alpha=as.integer(pft))) + theme_bw()+
  scale_color_jco()+  scale_alpha(range = c(0.2, 1))+
  #scale_color_manual(values = c("yellow", "lightblue", "lightgrey"))+
   theme(legend.position = "bottom") + 
  geom_segment(data = inter_NLMA, aes(x = x1, y = y1, xend = x2, yend = y2), size = 1.2)


osnas_globe=-0.62
osnas_intra=-0.18

intra_slope = mean(slope$intra_species)

ggplot(slope, aes(x=intra_species, fill=pft))+geom_histogram(alpha=0.5,colour = "grey32")+theme_bw() + scale_fill_jco()+
  geom_vline(xintercept = inter_slope, color = "orange", size = 1) +# xlim(-1,1.1)+
  geom_vline(xintercept = intra_slope, color = "orange", linetype="dashed",size = 1) + 
  geom_vline(xintercept = osnas_globe, color = "black",size = 1)+
  geom_vline(xintercept = osnas_intra, color = "black", linetype="dashed",size = 1) + theme(legend.position = "below")

ggplot(slope, aes(x=intra_species, fill=pft))+geom_histogram(alpha=0.5, colour = "grey32")+theme_bw() + scale_fill_jco()+
  geom_vline(xintercept = inter_slope, color = "green", size = 1) +
  geom_vline(xintercept = intra_slope, color = "green", linetype="dashed", size = 1) +
  geom_vline(xintercept = osnas_globe, color = "black", size = 1)+
  geom_vline(xintercept = osnas_intra, color = "black", linetype="dashed", size = 1)+
  geom_vline(xintercept = mean(slope$intra_species[slope$pft=="DBF"]), color = "#00AFBB", linetype="dashed", size = 1) +
  geom_vline(xintercept = mean(slope$intra_species[slope$pft=="EBF"]), color = "#E7B800", linetype="dashed", size = 1) +
  geom_vline(xintercept = mean(slope$intra_species[slope$pft=="ENF"]), color = "#FC4E07", linetype="dashed", size = 1) 

hist(slope$les, breaks = 40, col = "darkorange")
abline(v = inter_slope, add=T, col = "darkgreen",lwd=3)


mean_intra = Reduce("+", res_corr) / length(res_corr)
mean_intra= apply(simplify2array(intra_corr), 1:2, mean)

colnames(mean_intra) = rownames(mean_intra) = c("LMA", "N", "C", "lignin", "cellulose", "ChlA", "ChlB", "Carotenoids")
colnames(inter_cor) = rownames(inter_cor) = c("LMA", "N", "C", "lignin", "cellulose", "ChlA", "ChlB", "Carotenoids")
corrplot::corrplot(inter_cor, method = "pie", diag = F, addCoefasPercent = T, type = "lower")
corrplot::corrplot(mean_intra, method = "pie", diag = F, addCoefasPercent = T, type = "lower")


environmental_features = readr::read_csv("./TOS_retriever/out/wordclim_fia.csv")
full_nores = left_join(full_nores, environmental_features)
intra_env = r2s = lmds = list()
for(tx in intra_NLMA$taxonID){
  tx_full = full_nores %>% filter(taxonID == tx)
  tx_full = tx_full[complete.cases(tx_full),]
  if(nrow(tx_full) >10){
    md = lm(log10(leafMassPerArea.Estimate)~log10(nitrogenPercent.Estimate), data = tx_full)
    tx_full["nitrogenPercent.Estimate"]= residuals(md)
    tx_full = tx_full %>% select(nitrogenPercent.Estimate, bio1, elevation, LAT, LON, taxonID)
    tx_full = tx_full %>% group_by(LAT, LON) %>% summarize_all(mean)
    #lm(Sale~poly(Year,2,raw=TRUE))
    #tave=gam(nitrogenPercent.Estimate~bio1, data = tx_full)
    #tave=lm(nitrogenPercent.Estimate~poly(bio1,2), data = tx_full)
    tave=lm(nitrogenPercent.Estimate~bio1  + I(bio1^2), data = tx_full)
    
    plot(tx_full$bio1, tx_full$nitrogenPercent.Estimate)
    
    r2tave = summary(tave)
    #ele=lm(nitrogenPercent.Estimate~poly(elevation,2), data = tx_full)
    #ele=gam(nitrogenPercent.Estimate~elevation, data = tx_full)
    ele=lm(nitrogenPercent.Estimate~elevation + I(elevation^2), data = tx_full)
    
    r2ele = summary(ele)
    intra_env[[tx]] = cbind.data.frame(tx_full, predict(tave), predict(ele))
    r2s[[tx]] = c(r2tave$r.sq, r2ele$r.sq)
  }
}

taxa = names(r2s)
r2s = do.call(rbind.data.frame, r2s) 
r2s = data.frame(taxa,r2s)
colnames(r2s) = c("taxa", "Tave", "Elev")
r2s["pft"] = "DBF"
r2s$pft[r2s$taxa %in% pfts$taxonID[pfts$pft =="e"]] = "ENF"
r2s$pft[r2s$taxa %in% pfts$taxonID[pfts$pft =="nn"]] = "EBF"
intra_env = do.call(rbind.data.frame, intra_env)
intra_env["pft"] = "DBF"
intra_env$pft[intra_env$taxonID %in% pfts$taxonID[pfts$pft =="e"]] = "ENF"
intra_env$pft[intra_env$taxonID %in% pfts$taxonID[pfts$pft =="nn"]] = "EBF"

mean(r2s$Tave)
mean(r2s$Elev)
ggplot(data = r2s, aes(x=Tave, fill=pft, alpha=0.5))+ geom_histogram(colour="black") + theme_bw() + scale_fill_jco() + 
  geom_vline(xintercept = mean(r2s$Tave), color = "black", linetype="dashed", size = 1) +   theme(legend.position = "none")

ggplot(data = r2s, aes(x=Elev, fill=pft, color="black", alpha=0.5))+ geom_histogram(colour="black") + theme_bw() + scale_fill_jco()+ 
  geom_vline(xintercept = mean(r2s$Elev), color = "black",linetype="dashed", size = 1)+  theme(legend.position = "none")

colnames(intra_env) = c("N", "Tave","elevation", "LAT", "taxonID", "p_tave", "p_ele", "PFT")
ggplot(data = intra_env, aes(x = Tave, y = p_tave, color = taxonID)) + geom_point(size=0.3)+facet_grid(PFT~., scales = "free")+
  theme(legend.position = "none")




#trying to understand what drives the N:LMA patterns within species
N_coeff = list()
LMA_coeffs = list()
library(mgcv)
for(tx in taxa_list){
  tx_full = full_nores %>% filter(taxonID == tx)
  tx_full = tx_full[complete.cases(tx_full),]
  if(nrow(tx_full) >10){
    
  #temperature
  t_n=gam(nitrogenPercent.Estimate~bio1, data = tx_full)
  t_lma=gam(leafMassPerArea.Estimate~bio1, data = tx_full)
  t_n = summary(t_n)
  t_lma = summary(t_lma)
  rt_n = c(t_n$p.coeff[2], t_n$r.sq)
  rt_lma = c(t_lma$p.coeff[2], t_lma$r.sq)
  #precipitation
  t_n=gam(nitrogenPercent.Estimate~bio16, data = tx_full)
  t_lma=gam(leafMassPerArea.Estimate~bio16, data = tx_full)
  t_n = summary(t_n)
  t_lma = summary(t_lma)
  rp_n = c(t_n$p.coeff[2], t_n$r.sq)
  rp_lma = c(t_lma$p.coeff[2], t_lma$r.sq)
  #elevation
  t_n=gam(nitrogenPercent.Estimate~elevation, data = tx_full)
  t_lma=gam(leafMassPerArea.Estimate~elevation, data = tx_full)
  t_n = summary(t_n)
  t_lma = summary(t_lma)
  re_n = c(t_n$p.coeff[2], t_n$r.sq)
  re_lma = c(t_lma$p.coeff[2], t_lma$r.sq)
  
  #append coeffs
  tmp = c(tx, rt_n, rp_n, re_n)
  names(tmp) = c("taxonID", "aveT", "r2_aveT", "aveP", "r2_aveP", "aveEle", "r2_aveEle" )
  N_coeff[[tx]] = tmp
  tmp = c(tx, rt_lma, rp_lma, re_lma)
  names(tmp) = c("taxonID", "aveT", "r2_aveT", "aveP", "r2_aveP", "aveEle", "r2_aveEle" )
  LMA_coeffs[[tx]] = tmp
}
}
lma_c = do.call(rbind.data.frame, LMA_coeffs)
n_c = do.call(rbind.data.frame, N_coeff)
colnames(n_c)=colnames(lma_c)=c("taxonID", "aveT", "r2_aveT", "aveP", "r2_aveP", "aveEle", "r2_aveEle" )

fofo = inner_join(lma_c, n_c, by="taxonID")
fofo_t = fofo %>% filter(r2_aveT.x > 0.05, r2_aveT.y>0.05)
plot(fofo_t$aveT.x, fofo_t$aveT.y)

fofo_p = fofo %>% filter(r2_aveP.x > 0.05, r2_aveP.y>0.05)
plot(fofo_p$aveP.x, fofo_p$aveP.y)
