get_data_pca_transforamtion <- function(climate_trends) {
  #PCA by variable
  var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
  climate_pca = NULL
  climate_features = NULL
  for(ii in var){
    climate_pca[[ii]] <- climate_trends %>%
      select(contains(ii)) %>%
      FactoMineR::PCA(ncp=1)
    #get pca transformed data
    climate_features[[ii]] = climate_pca[[ii]]$ind$coord
  }
  climate_features = do.call(cbind.data.frame, climate_features)
  colnames(climate_features) <- var
  climate_features <- data.frame(climate_trends[["individualID"]], climate_features)
  colnames(climate_features)[1] <- "individualID"
  return(list(climate_pca=climate_pca, climate_features=climate_features))
}

predict_data_pca_transforamtion <- function(climate_trends, pca_mod) {
  library(FactoMineR)
  #PCA by variable
  var <- c("daylength", "prec", "rad", "snow_melt", "tmax", "tmin", "vp")
  climate_features = NULL
  # climate_trends[-1] = scale(climate_trends[-1], center = scaling_vars$center[-c(1:7)],
  #                           scale = scaling_vars$scale[-c(1:7)])
  climate_pca = NULL
  for(ii in var){
    climate_features[[ii]] <- climate_trends %>%
      select(contains(ii)) 
    pca_mod[[ii]]$var$coord = data.frame((pca_mod[[ii]]$var$coord))
    climate_pca[[ii]] <-  predict(pca_mod[[ii]], (climate_features[[ii]]))
    #get pca transformed data
    climate_features[[ii]] = climate_pca[[ii]]["cos2"]
  }
  climate_features = do.call(cbind.data.frame, climate_features)
  colnames(climate_features) <- var
  climate_features <- data.frame(climate_trends[1], climate_features)
  colnames(climate_features)[1] <- "individualID"
  return(list(climate_pca=climate_pca, climate_features=climate_features))
}
