#get species weighted averages
library(modEvA)
require(venneuler)
library(tidyverse)
sp = readRDS("./models/ch2perc_gauss_species_noml.rds")
env = readRDS("./models/ch2perc_gauss_env_noml.rds")
full = readRDS("./models/ch2perc_gauss_full_noml.rds")
result = rbind.data.frame(env$R2, sp$R2, full$R2)
result$Features = c(rep("Environmental",8), rep("Phylogenetic",8), rep("Combined",8))
result$Trait = rep(c("N", "C", "lignin", "cellulose", "LMA", "chlA", "chlB", "carot"),3)
source("./src/model_training/myvarPart.R")
part <- list()
err <-list()
for(tr in unique(result$Trait)){
  dat <- result %>% filter(Trait == tr)
  part[[tr]]  <- myvarPart(A = unlist(dat[dat$Features=="Environmental","Estimate"]),
                           B = unlist(dat[dat$Features=="Phylogenetic","Estimate"]),
                  AB = unlist(dat[dat$Features=="Combined","Estimate"]), 
                  A.name = "Environmental", B.name = "Phylogenetic",  
                  plot = TRUE, plot.digits = 3, cex.names = 3, cex.values = 1.2, main = tr, cex.main = 2, 
                  plot.unexpl = TRUE)
}
ff <- do.call(cbind.data.frame, part)
colnames(ff)<- unique(result$Trait)
ff <- t(ff) %>% data.frame()
colnames(ff) <- c("Environmental", "Phylogenetic", "Joint", "Unexplained")

ff <- ff[order(ff$Phylogenetic, decreasing = F) , ]
ff <- ff %>% select(Environmental, Joint, Phylogenetic, Unexplained)
ff$Trait = c("ChlB [%]", "ChlA [%]", "Carotenoids [%]", "N [%]",
             "Cellulose [%]", "C [%]", "Lignin [%]", "LMA")
ff$Trait <- factor(ff$Trait , levels = ff$Trait )

readr::write_csv(ff, "./outdir/var_part_fit.csv")
ff_nounex = ff %>% select(-one_of("Unexplained")) %>%
  gather(Feature, R2, Environmental:Phylogenetic, factor_key=TRUE)

ggplot(ff_nounex, aes(Trait)) + geom_bar(aes(weight = R2, fill = Feature)) +
  #geom_bar(aes(weight = R2, fill = Feature)) +
  theme_bw() + ylim(0,1) + 
  scale_fill_manual(values=c("#E69F00", "#999999", "#BFDF3D", "#56B4E9")) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  #geom_abline(slope = 0, intercept = 0.4)

result$Trait <- factor(result$Trait, levels = rownames(ff))
result$Features = factor(result$Features, levels=c("Environmental", "Combined", "Phylogenetic"))
result%>% #filter(Features == "Combined") %>%
  #slice(match(rownames(ff), Trait)) %>% #position=position_dodge(.9
ggplot(aes(x= Trait, y = Estimate, color = Features)) + geom_point(position=position_dodge(.5)) +
  theme_bw() + ylim(0,1) +     theme(legend.position = "bottom") +
  geom_errorbar(aes(ymin=Q2.5, ymax=Q97.5), width=.2, position=position_dodge(.5))+
  scale_color_manual(values=c("#E69F00","#56B4E9", "#BFDF3D")) 
  coord_flip() + geom_abline(slope = 0, intercept = 0.4) 

