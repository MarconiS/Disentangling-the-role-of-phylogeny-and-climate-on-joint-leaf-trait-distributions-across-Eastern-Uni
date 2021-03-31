library(tidyverse)
library(brms)
full = readRDS("./models/ch2perc_gauss_full_noml.rds")

pop_effs = fixef(full$mod)
#pop_effs = cbind.data.frame(rownames(pop_effs), pop_effs)
tr_fe = strsplit(rownames(pop_effs), split = "_")
tr_fe = do.call(rbind.data.frame,tr_fe )[1:2]
colnames(tr_fe) = c("Trait", "Driver")

foo = cbind.data.frame(tr_fe, pop_effs)
#foo2 = dplyr::filter(foo, grepl("Intercept",Driver))
foo = dplyr::filter(foo, !grepl("Intercept",Driver))
foo$Driver=rep(c("Precpt", "Radtn", "Tmax", "Tmin", "VP", "slope", "aspect", "elevation"),8)
foo$Trait = c(rep("N[%]", 8), rep("C[%]", 8), rep("lignin[%]", 8), rep("cellulose[%]", 8), 
              rep("LMA", 8), rep("chlorophyllA[%]", 8), rep("chlorophyllB[%]", 8), rep("carotenoids[%]", 8))
ggplot(foo, aes(x = Driver)) + geom_bar(aes(weight = Estimate),color="black", fill="dodgerblue3") +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(vars(Trait), nrow = 2, scales = "free")+
  theme_bw() + #ylim(0,1) + 
  theme(legend.position="bottom",  axis.text.x = element_text(angle = 45, hjust = 1))
coord_flip() 