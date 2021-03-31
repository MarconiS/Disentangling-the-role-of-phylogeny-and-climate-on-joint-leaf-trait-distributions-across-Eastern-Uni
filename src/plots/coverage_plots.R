#figure S1
library(tidyverse)
library(brms)
full = readRDS("./models/ch2perc_gauss_full_noml.rds")
sp = readRDS("./models/ch2perc_gauss_species_noml.rds")
env = readRDS("./models/ch2perc_gauss_env_noml.rds")


p_full = predict(full$mod, newdata = test_data)
p_sp = predict(sp$mod, newdata = test_data)
p_env = predict(env$mod, newdata = test_data)
y = test_data[c(4:11)]

q25 = p_full[,3,]
q975 = p_full[,4,]
lb_full = q25 < y
ub_full = q975 > y
full_coverage = lb_full & ub_full
full_coverage = apply(full_coverage, 2, sum)
full_coverage = full_coverage / nrow(lb_full)
unlist(full_coverage)[1:8] %>% as.numeric %>% mean
full_coverage$method = "combined"
pcf = mean(as.numeric(unlist(full_coverage)[1:8]))

q25 = p_sp[,3,]
q975 = p_sp[,4,]
lb_sp = q25 < y
ub_sp = q975 > y
sp_coverage = lb_sp & ub_sp
sp_coverage = apply(sp_coverage, 2, sum)
sp_coverage = sp_coverage / nrow(lb_sp)

sp_coverage$method = "species"
pcs = mean(as.numeric(unlist(sp_coverage)[1:8]))

q25 = p_env[,3,]
q975 = p_env[,4,]
lb_env = q25 < y
ub_env = q975 > y
env_coverage = lb_env & ub_env
env_coverage = apply(env_coverage, 2, sum)
env_coverage = env_coverage / nrow(lb_env)
env_coverage$method = "environment"
pce = mean(as.numeric(unlist(env_coverage)[1:8]))


names(full_coverage) = names(sp_coverage) = names(env_coverage) = c("N", "C", "ChlA", "ChlB", "Carot", "Lignin", "Cellulose", "LMA", "method")
final_coverage = rbind.data.frame(full_coverage, sp_coverage, env_coverage)
final_coverage = reshape2::melt(final_coverage, id.vars = "method")

colnames(final_coverage) = c("Method", "Trait", "Coverage")
ggplot(final_coverage, aes(x = Trait, y = Coverage, fill = Method)) + geom_point(colour="black",pch=21, size=3, position = "jitter") + 
  scale_fill_manual(values=c("#56B4E9", "#E69F00","lightgreen")) +
  geom_hline(yintercept=pcf, linetype="dashed", color = "#56B4E9") + 
  geom_hline(yintercept=pcs, color = "lightgreen", linetype="dashed") + 
  geom_hline(yintercept=pce, color = "#E69F00", linetype="dashed")+ ylim(0.9, 1)+
  geom_hline(yintercept=0.95) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))
  
