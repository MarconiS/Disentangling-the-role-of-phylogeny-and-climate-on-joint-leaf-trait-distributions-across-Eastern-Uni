FIA_full_map2 <- FIA_full_map2 %>% select(-one_of("yr", "plot", "treeid", "unitcd", "species_sc"))
dat <- rbind.data.frame(FIA_full_map, FIA_full_map2)
dat <- dat %>% group_by(LAT, LON) %>%       
  summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else first(.)))
write_csv(dat, "./outdir/outdir/temporary_full_map.csv")


summary <- readr::read_csv("./Summary_evaluation.csv")
varPartR2 <- function(dat){
  names(dat) <- c("A", "B", "AB")
  b_only = dat["AB"] - dat["A"] 
  a_only = dat["AB"] - dat["B"] 
  shared = dat["B"]  - b_only
  return(list(A = a_only, B = b_only, AB = shared))
}

library(dplyr)
partin <- apply(summary[-1],1, varPartR2) 
partin <- do.call(rbind.data.frame, partin) 
colnames(partin) <- c("Species", "Environment", "Shared")
partin$Element <- summary$Element
partin$Unexplained <- apply(partin[1:3], 1, sum)
partin$Unexplained <- 1 - partin$Unexplained
partin <- partin[order(partin$Species, decreasing = T) , ]
partin$Element <- factor(partin$Element, levels = partin$Element)

library(ggplot2)
library(tidyverse)
data_long <- gather(partin, Effect, Variance,  
      c("Environment", "Unexplained", "Shared", "Species"), factor_key=TRUE)
ggplot(data_long, aes(Element)) + geom_bar(aes(weight = Variance, fill = Effect)) +
  theme_bw() + ylim(0,1) + 
  scale_fill_manual(values=c("#E69F00", "#999999", "#BFDF3D", "#56B4E9")) +
  coord_flip() 

ggplot(data_long, aes(x = reorder(Element, -perc), y = perc)) + geom_bar(stat = "identity")
