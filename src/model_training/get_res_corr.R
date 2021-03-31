get_res_corr <- function(fit){
  
  # given a brmsfit with joint distribution (residual correlations = T), 
  # return the responsees correlation structure
  
  # arguments :
  #fit = brmsfit object
  
  #results:
  #corr_structure
  
  library(tidyverse)
  
  #get residuals 
  res_pred <- residuals(fit)
  #calculate response covariance matrix
  corr_mat <- VarCorr(fit)
  res_corr <- corr_mat$residual__$cor
  #get significant populatin level effects
  pop_effects_corr <- vcov(fit)
  pop_effects_corr[abs(pop_effects_corr) < 0.001] = NA
  
  #get most important features for each trait
  main_valuable_covariance_population <- list()
  tr = rownames(corr_mat$residual__$sd)
  for(ii in tr){
    main_valuable_covariance_population[[ii]] <- pop_effects_corr %>% 
      data.frame %>%  
      slice(which(rownames(pop_effects_corr) == paste(ii, "_Intercept", sep = ""))) %>%
      sort(decreasing = T) %>% t %>% head(30)
  }
  
  #get residuals correlation matrix
  res_corr <- corr_mat$residual__$cor
  res_corr_sd <- corr_mat$residual__$sd
  resid_pred = sigma_mu = sigma_sd = list()
  for(ii in 1: dim(corr_mat$residual__$cov)[1]){
    sigma_mu[[ii]] = corr_mat$residual__$cov[,1,ii]
    sigma_sd[[ii]] = corr_mat$residual__$cov[,2,ii]
    resid_pred[[ii]] = res_pred[,1,ii]
  }
  sigma_mu <- do.call(cbind.data.frame, sigma_mu)
  sigma_sd <- do.call(cbind.data.frame, sigma_sd)
  resid_pred <- do.call(cbind.data.frame, resid_pred)
  colnames(sigma_mu) = rownames(corr_mat$residual__$sd)
  rownames(sigma_mu) = rownames(corr_mat$residual__$sd)
  colnames(sigma_sd) = rownames(corr_mat$residual__$sd)
  rownames(sigma_sd) = rownames(corr_mat$residual__$sd)
  colnames(sigma_mu) = rownames(corr_mat$residual__$sd)
  colnames(resid_pred) = rownames(corr_mat$residual__$sd)
  
  res_cor = (cor(resid_pred))
  corstr <- Hmisc::rcorr(as.matrix(resid_pred))
  fixed_corr_filter <- corstr$P <0.01
  fixed_corr_filter[is.na(fixed_corr_filter)] = T
  fixed_cov_str = cor(resid_pred) * fixed_corr_filter 
  colnames(res_cor)= rownames(res_cor) = c("N%", "C%", "ChlA", "ChlB", "Carot", "lignin%", "cell%", "LMA")
  corrplot::corrplot(res_cor, type="lower", order="hclust", tl.srt=45,# addCoef.col = "black",
                      sig.level = 0.01, diag = F)
  plot(hclust(dist(res_cor)))
  results <- list(fixed_cov_str, pop_effects_corr, res_cor, resid_pred, sigma_mu, sigma_sd)
  return(results)
}
