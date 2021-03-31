retrotransform_features<- function(dat, scale = T, logtr = T, ft_nm){
  traits <- readr::read_csv("./TOS_retriever/out/utm_dataset.csv") %>% #siteID.x, taxonID.x,
    dplyr::select(c(individualID,ft_nm)) %>%
    group_by(individualID) %>%
    summarise_all(mean) %>%
    unique
  
  
  traits <- traits[complete.cases(traits), ft_nm] 
  if(logtr == T){
    traits[!colnames(traits) %in% c("individualID", "d15N", "d13C")] <- traits %>% 
      dplyr::select(-one_of("d15N", "d13C")) %>% log()
  }
  mean_traits <-  apply(traits,2, mean, na.rm = T)
  sd_traits <-  apply(traits,2, sd, na.rm = T)
  aa <- dat
  if(!is.na(dim(aa)[3])){
    for(jj in 1:dim(dat)[3]){
      aa[,-2,jj] <- dat[,-2,jj] * sd_traits[jj] + mean_traits[jj]
      #aa[,2,jj] <- dat[,2,jj] * sd_traits[jj] + mean_traits[jj]
      # aa[,3,jj] <- dat[,3,jj] * sd_traits[jj] + mean_traits[jj]
      # aa[,4,jj] <- dat[,4,jj] * sd_traits[jj] + mean_traits[jj]
    }
  }else{
    for(jj in 1:dim(aa)[2]){
      aa[,jj] <- dat[,jj] * sd_traits[jj] + mean_traits[jj]
    }
  }
  if(logtr == T){
    if(is.na(dim(aa)[3])){
      aa[!colnames(aa) %in% c("d15N", "d13C")] <- aa %>% dplyr::select(-one_of("d15N", "d13C")) %>% exp()
    }else{
      aa[,,1:15] <- exp(aa[,,1:15])
    }
  }
  return(aa)
}

