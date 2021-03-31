library(BIEN)
library(tidyverse)
lma <- BIEN_trait_trait("leaf area per leaf dry mass", all.taxonomy = TRUE, political.boundaries = TRUE) %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched)
leaf_c <- BIEN_trait_trait("leaf carbon content per leaf dry mass", all.taxonomy = TRUE, political.boundaries = TRUE) %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched)#
leaf_n <- BIEN_trait_trait("leaf nitrogen content per leaf dry mass", all.taxonomy = TRUE, political.boundaries = TRUE) %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched)
leaf_p <- BIEN_trait_trait("leaf phosphorus content per leaf dry mass", all.taxonomy = TRUE, political.boundaries = TRUE) %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched)


all_bien <- list(lma, leaf_c, leaf_n, leaf_p)
saveRDS(all_bien, "./indir/BIEN_dat.rds")

all_bien <- readRDS("./indir/BIEN_dat.rds")
lma <- all_bien[[1]] %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched) %>%
  mutate(trait_value = replace(trait_value, values = as.numeric(trait_value)^(-1) * 1000 ))
leaf_c <- all_bien[[2]] %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched) %>%
  mutate(trait_value = replace(trait_value, values = as.numeric(trait_value)/10 ))

leaf_n <- all_bien[[3]] %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched) %>%
  mutate(trait_value = replace(trait_value, values = as.numeric(trait_value)/10))

leaf_p <- all_bien[[4]] %>%
  select(trait_name, trait_value, latitude, longitude, elevation_m, country, state_province, name_matched)%>%
  mutate(trait_value = replace(trait_value, values = as.numeric(trait_value)/10))


bien_dataset <- rbind.data.frame(lma, leaf_c, leaf_n, leaf_p) %>%
  filter(country =="United States")

leaf_p %>%  filter(country =="United States") %>% select(state_province) %>% unique

bien_dataset$trait_value <- as.numeric(bien_dataset$trait_value)
bien_dataset <- bien_dataset[complete.cases(bien_dataset),]

bien_dataset %>% filter(trait_name =="leaf phosphorus content per leaf dry mass") %>% select(state_province) %>% unique
bien_dataset %>% filter(trait_name =="leaf carbon content per leaf dry mass") %>% select(state_province) %>% unique
bien_dataset %>% filter(trait_name =="leaf nitrogen content per leaf dry mass") %>% select(state_province) %>% unique
bien_dataset %>% filter(trait_name =="leaf area per leaf dry mass") %>% select(state_province) %>% unique


write_csv(bien_dataset, "./outdir/bien_data.csv")
