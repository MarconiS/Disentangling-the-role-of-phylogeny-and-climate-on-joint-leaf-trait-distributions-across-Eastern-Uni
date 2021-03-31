get_pca_effects <- function(bjtdm, nComp){
  library(plotly)  
  res_pca <- bjtdm$norm_params$clim$weights
  aa <- fviz_contrib(res_pca, choice = "var", axes = 1, top = 81)
  
  #get importance of variables er components
  var_portion <- rownames(bjtdm$norm_params$clim$weights$rotation) %>% data.frame
  for(ii in 1:nComp){
    var_portion = cbind.data.frame(var_portion, fviz_contrib(res_pca, choice = "var", axes = 1, top = 81)$data["contrib"])
  }
  colnames(var_portion) <- c("features", paste("PC", 1:nComp, sep="_"))
  #aa <-fviz_cos2(res_pca, choice = "var", axes = 1)
  
  #get importance of each components per trait from the model
  vc <- VarCorr(bjtdm$model)
  
  cfs <- summary(bjtdm$model) 
  cfs <- cfs$splines %>% data.frame
  cfs$cf <- rownames(cfs)
  cfs$signmificance <- (cfs$l.95..CI * cfs$u.95..CI) >0
  cfs[!cfs$signmificance, "Estimate"] = 0
  
  
  #multiply rotation per spline coefficient
  tr = "leafMassPerArea"
  
  clim_effects = NULL
  ff <- list()
  for(tr in colnames(bjtdm$y_hat_tst[-1])){
    #get spline terms
    spline_terms <- cfs %>% filter(grepl(paste(tr, "sClim", sep="_"),cf))
    
    #multiply spline terms by PC rotation
    foo <- as.matrix(var_portion[-1]) %*% as.matrix(spline_terms[1]) %>% data.frame
    #foo <- as.matrix(spline_terms[1:4]) %>% data.frame
    
    ff[[tr]] <- as.matrix(spline_terms[1:4])
    foo$trait <- tr
    #add to data
    foo["feature"] = rownames(foo)
    clim_effects <- rbind.data.frame(clim_effects, foo)
    
  }
  #clim_effects$feature = rownames(clim_effects)
  
  ft_grp <- strsplit(clim_effects$feature, split = "_") 
  clim_effects["ft_grp"] <- do.call(rbind.data.frame, ft_grp)[1]
  clim_effects["MOY"] <- do.call(rbind.data.frame, ft_grp)[2]
  ggplot(clim_effects, aes(y = MOY, x = Estimate, color=trait)) + geom_point() + 
    facet_grid(ft_grp ~ ., scales = "free", space = "free") +
    theme(strip.text.y = element_text(angle = 0))
  
  spidePlots <- list()
  #spider plot
  for(tr in colnames(bjtdm$y_hat_tst[-1])){
    dat <- clim_effects %>% filter(trait == tr) %>%
      select(Estimate) %>% t  #, , Estimate
    
    feat <- clim_effects %>% filter(trait == tr) %>%
      select(feature) %>% unlist %>% as.character
    #feat <- paste("PC", 1:8, sep="_")
    dat <- as.vector(dat)
    #
    p <- plot_ly(
      type = 'scatterpolar',
      r = dat,
      theta = feat,
      fill = 'toself'
    ) %>%
      layout(
        polar = list(
          radialaxis = list(
            visible = T,
            range = c(min(dat),max(dat))
          )
        ),
        showlegend = F
      )
    spidePlots[[tr]] <- print(p)
  }
  
  
  species <- vc$taxonID
  
}

get_species_influence <- function(bjtdm){
  
}
