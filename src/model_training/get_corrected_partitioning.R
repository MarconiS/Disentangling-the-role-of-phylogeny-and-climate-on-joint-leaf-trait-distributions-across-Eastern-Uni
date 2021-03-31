bj1 <- readRDS("./outdir/bjsTDM_full_log3d_full.rds")
set.seed(1987)
R2full <- retrotransformed_r2(bj1$model, 
                    newdata = bj1$tst_set, 
                    sd_y = bj1$norm_params$tr_dat$sd, 
                    mean_y = bj1$norm_params$tr_dat$mu, 
                    summary = F)

bj2 <- readRDS("./outdir/bjsTDM_clim_3d.rds")
set.seed(1987)
R2clim <- retrotransformed_r2(bj2$model, 
                              newdata = bj2$tst_set, 
                              sd_y = bj2$norm_params$tr_dat$sd, 
                              mean_y = bj2$norm_params$tr_dat$mu
                              , summary = F)

bjstdm <- readRDS("./outdir/bjsTDM_sp.rds")
set.seed(1987)
R2sp <- retrotransformed_r2(bjstdm$model, 
                              newdata = bjstdm$tst_set, 
                              sd_y = bjstdm$norm_params$tr_dat$sd, 
                              mean_y = bjstdm$norm_params$tr_dat$mu
                            , summary = F)



single_R2s = list()
for(ii in 1:nrow(R2sp)){
  sR2 <- cbind.data.frame(t(R2full[ii,]), t(R2clim[ii,]), t(R2sp[ii,]))
  sR2 <- t(sR2) %>% data.frame
  # sR2$Features <- c(rep("Full model", 18), rep("Environmental", 18), rep("Evolutionary", 18))
  # tr_names <- c("LMA", "lignin", "Cellulose", "P", "K", "Ca", "ChlA",  "ChlB",
  #               "Carotenoids",  "Mg",  "S",  "Mn",  "Fe",  "Cu",  "B",  "Zn",
  #               "N", "C")
  # sR2$Trait <- rep(tr_names, 3)
  single_R2s[[ii]] <- sR2
}
single_R2s = do.call(cbind.data.frame, single_R2s)
colnames(single_R2s) <- c("N", "C", "ChlA", "ChlB", "Carotenoids", "LA", "LMA", "lignin", "Cellulose")
single_R2s$Features <- c(rep("Full model", 4), rep("Environmental", 4), rep("Evolutionary", 4))
# tr_names <- c("LMA", "lignin", "Cellulose", "P", "K", "Ca", "ChlA",  "ChlB",
#                "Carotenoids",  "Mg",  "S",  "Mn",  "Fe",  "Cu",  "B",  "Zn",
#                "N", "C")
single_R2s = gather(single_R2s, Trait, Value, N:Cellulose, factor_key=TRUE)
readr::write_csv(single_R2s, "./outdir/single_extractions_R2.csv")


ggplot(data=single_R2s, aes(x=Trait, y=Value, fill=Features)) +
  geom_bar(stat="identity", position=position_dodge(), color="black")+
  theme_bw() + theme(legend.position="bottom") + 
  scale_fill_manual(values = c("orange", "olivedrab3", "steelblue2")) 
#c("FFB266", "009900", "004C99"))

all_methodsR2 <- rbind.data.frame(R2full, R2clim, R2sp)
all_methodsR2$Features <- c(rep("Full model", 18), rep("Environmental", 18), rep("Evolutionary", 18))
tr_names <- c("LMA", "lignin", "Cellulose", "P", "K", "Ca", "ChlA",  "ChlB",
              "Carotenoids",  "Mg",  "S",  "Mn",  "Fe",  "Cu",  "B",  "Zn",
              "N", "C")
all_methodsR2$Trait <- rep(tr_names, 3)
write_csv(all_methodsR2, "./outdir/var_partitioning.csv")

#coverage
bj1_predictions <- predict(bj1$model, bj1$tst_set)
bj2_predictions <- predict(bj2$model, bj2$tst_set)
bj3_predictions <- predict(bj3$model, bj3$tst_set)

coverage_bj1 <- coverage_bj2 <- coverage_bj3 <- list()
predictions_bj1 <- predictions_bj2 <- predictions_bj3 <- list()
for(tr in 1:dim(bj2_predictions)[3]){
  tr_lb <- dimnames(bj2_predictions)[[3]][tr]
  #coverage can be in saled data
  coverage_bj1[[tr]] <- bj1_predictions[,3,tr_lb] < bj1$tst_set[tr_lb] &
    bj1_predictions[,4,tr_lb] > bj1$tst_set[tr_lb]
  coverage_bj1[[tr]] <- sum(coverage_bj1[[tr]]) / nrow(bj1$tst_set) 
  
  coverage_bj2[[tr]] <- bj2_predictions[,3,tr_lb] <= bj2$tst_set[tr_lb] &
    bj2_predictions[,4,tr_lb] >= bj2$tst_set[tr_lb]
  coverage_bj2[[tr]] <- sum(coverage_bj2[[tr]]) / nrow(bj2$tst_set) 
  
  coverage_bj3[[tr]] <- bj3_predictions[,3,tr_lb] <= bj3$tst_set[tr_lb] &
    bj3_predictions[,4,tr_lb] >= bj3$tst_set[tr_lb]
  coverage_bj3[[tr]] <- sum(coverage_bj3[[tr]]) / nrow(bj3$tst_set) 
  
  r1_tr <- bj2_predictions[,, tr] * bj1$norm_params$tr_dat$sd[tr] + 
    bj1$norm_params$tr_dat$mu[tr]
  r2_tr <- bj2_predictions[,, tr] * bj2$norm_params$tr_dat$sd[tr] + 
    bj2$norm_params$tr_dat$mu[tr]
  r3_tr <- bj2_predictions[,, tr] * bj3$norm_params$tr_dat$sd[tr] + 
    bj3$norm_params$tr_dat$mu[tr]
  predictions_bj1[[tr_lb]] <- exp(r1_tr)
  predictions_bj2[[tr_lb]] <- exp(r2_tr)
  predictions_bj3[[tr_lb]] <- exp(r3_tr)
  
}

  
coverage_bj1 <- do.call(rbind.data.frame, coverage_bj1)
coverage_bj2 <- do.call(rbind.data.frame, coverage_bj2)
coverage_bj3 <- do.call(rbind.data.frame, coverage_bj3)
names(coverage_bj1) <- names(coverage_bj2) <- names(coverage_bj3) <- "Coverage"
coverage <- rbind.data.frame(coverage_bj1, coverage_bj2, coverage_bj3)
coverage$Features <- c(rep("Full model", 18), rep("Environmental", 18), 
                       rep("Evolutionary", 18))
coverage$Trait <- rep(tr_names, 3)

ggplot(coverage, aes(y = Coverage, x = Trait, color = Features)) + 
  geom_hline(yintercept = 0.95,  colour= "grey76") + ylim(0.8,1) +
  geom_hline(yintercept = mean(coverage_bj1$Coverage), linetype="dashed", colour= "grey25") +
  geom_hline(yintercept = mean(coverage_bj2$Coverage), linetype="dashed", colour= "coral") +
  geom_hline(yintercept = mean(coverage_bj3$Coverage), linetype="dashed", colour= "cyan3") +
  scale_color_manual(values=c("coral", "cyan3", "grey25"))+
  theme_bw()+ geom_point(position=position_dodge(width=0.5)) 
