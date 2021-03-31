get_epsg_from_utm <- function(utm){
  utm <-  substr(utm,1,nchar(utm)-1)
  epsg <- paste("326", utm, sep="")
  return(epsg)
}

get_lat_long <- function(field_data, epsg = 32617){
  library(daymetr)
  library(rgdal)
  new_dat <- NULL
  for(NeonSites in unique(field_data$siteID)){
    dat <- field_data[field_data$siteID %in% NeonSites,] 
    epsg <- get_epsg_from_utm(unique(dat$utmZone))
    dat <- dat[complete.cases(dat$UTM_E), ]
    coordinates(dat) <- c("UTM_E", "UTM_N")
    proj4string(dat) <- CRS(paste("+init=epsg:", epsg, sep ="")) 
    CRS.new <- CRS("+init=epsg:4326")
    dat <- spTransform(dat, CRS.new)
    coords_dat <- dat@coords
    new_dat <- rbind(new_dat, cbind(dat@data ,dat@coords))
  }
  return(new_dat)
}

download_point_daymet <- function(listSites = NULL, path = "./TOS_retriever/out/field_data.csv", 
                                  outpath = './DMT_retriever/daymet_ts/'){
  require(tidyverse)
  field_data <- readr::read_csv(path) %>% unique
  colnames(field_data) 
  listSites <- field_data$siteID %>% unique
  if(is.null(listSites)){
    listSites = unique(field_data$siteID)
  }
  if(is.null(field_data$UTM_N) || is.na(field_data$UTM_N)){
    daymet_coords <- cbind(field_data$individualID,field_data$decimalLatitude, field_data$decimalLongitude) %>% unique
  }else{
    new_dat <- get_lat_long(field_data)
    daymet_coords <- cbind(field_data$individualID,new_dat$UTM_N, new_dat$UTM_E) %>% unique
  }
  readr::write_csv(data.frame(daymet_coords), './DMT_retriever/points_coords.csv', col_names=F)
  # ssurgo_coords <- cbind(new_dat$individualID, new_dat$siteID, new_dat$UTM_N, new_dat$UTM_E) %>% unique
  # readr::write_csv(data.frame(ssurgo_coords), './AOP_retriever/indir/tos_coords_full.csv', col_names=F)
  
  library(daymetr)
  
  download_daymet_batch(file_location = './DMT_retriever/points_coords.csv',
                        start = 1995,
                        end = 2015,
                        internal = F,
                        path = outpath)
  #ls_dat <- list.files("./DMT_retriever/daymet_ts//", pattern = ".csv")
  #climate_features <- melt_daymet(ls_dat)
  #save(df, file= paste(outpath, "Daymet_traits_preds.RData", sep="/"))
}

melt_daymet <- function(path = "//orange/ewhite/s.marconi/FIA_climate/Daymet_ts/"){
  library(readr)
  library(daymetr)
  library(tidyverse)
  library(lubridate)
  dataset = data.frame(matrix(NA, ncol = 9, nrow = 0))
  colnames(dataset) <- c("month","ts_daylength", "ts_prec","ts_rad","ts_melt","ts_tmax","ts_tmin","ts_vp",
                         "identifierID")
  ls_dat = list.files(path, pattern = "csv")
  for(ii in ls_dat){
    id_clim <- read.csv(paste(path,ii, sep="/"), skip=7)
    colnames(id_clim) <- c("year", "month", "daylength", "prec", "srad", "snow_melt", "tmax", "tmin", "vp")
    id_clim$month <- as.Date(id_clim$month-1, origin = "1995-01-01") %>%
      month
    point_features <- id_clim %>%
      dplyr::group_by(month, year) %>%
      dplyr::summarise(daylength = mean(daylength), prec = sum(prec), 
                       srad = mean(srad), snow_melt = sum(snow_melt), 
                       tmax = max(tmax), tmin = min(tmin), vp = mean(vp))
    
    point_features <- point_features %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(daylength = mean(daylength), prec = mean(prec), 
                       rad = mean(srad), snow_melt = mean(snow_melt), 
                       tmax = mean(tmax), tmin = mean(tmin), vp = mean(vp))
    point_features$identifierID <- gsub('.{14}$', '', ii)
    dataset = rbind(dataset, point_features)
  }
  #write_csv(dataset, "./DMT_retriever/climate_features.csv")
  library(data.table) ## v >= 1.9.6
  clim_dat <- dcast(melt(dataset, id.vars=c("identifierID", "month")), 
                    identifierID~variable+month)
  readr::write_csv(clim_dat, "./DMT_retriever/bien_climate_features.csv")
  return(clim_dat)
}

get_train_test_climate <- function(tr_lb, tst_lb, dat, nComp, ft_nm, pcaT = T, in_pt = "./DMT_retriever/climate_features.csv"){
  library(dplyr)
  mean_nona <- function(x)mean(x, na.rm=T)
  feat_dat <- readr::read_csv(in_pt) %>% 
    group_by(identifierID) %>% summarize_all(mean_nona) %>%
    select(-one_of("snow_melt_7",   "snow_melt_8",   "snow_melt_9"))  
  clim_dat <- feat_dat %>% 
    filter(identifierID %in% tr_lb)  %>% unique
  
  sum(is.na(clim_dat))
  x_mean<- clim_dat %>% dplyr::select(-one_of("identifierID")) %>% apply(2, mean)
  x_sd <- clim_dat %>% dplyr::select(-one_of("identifierID")) %>% apply(2, sd)
  Xt = clim_dat %>% dplyr::select(-one_of("identifierID")) %>% 
    scale(center = x_mean, scale = x_sd)
  
  Yt = dat %>%
    filter(individualID %in% tr_lb)  %>% 
    dplyr::select(c("individualID", ft_nm)) %>%
    group_by(individualID) %>%
    summarize_all(mean_nona)
  
  library(plsdepot)  
  
  
  if(pcaT ==T){
    x_pca =  prcomp(Xt, scale. = T)
    x_scores = x_pca$x
    modpca <- x_pca
  }else{
    x_weights <- plsreg2(Xt, Yt[ft_nm], nComp)
    x_scores <- as.matrix(scale(Xt)) %*% x_weights$mod.wgs
  }
  colnames(x_scores)[1:nComp] <- paste("Clim", seq(1:nComp), sep=".")
  pls_tr_data = data.frame(Yt, x_scores[, 1:nComp])
  
  #test data
  X_tst <- feat_dat %>% 
    filter(identifierID %in% tst_lb) %>%
    dplyr::select(-one_of("identifierID"))
  
  sum(is.na(X_tst))
  
  
  if(pcaT ==T){
    X_tst <- scale(X_tst, center = x_mean,
                          scale =  x_sd)
    x_tst_scores = predict(x_pca, X_tst)[,1:nComp]  %>%
      data.frame
  } else{
    x_tst_scores <- scale(X_tst, center = x_mean, 
                          scale =  x_sd) %*% x_weights$mod.wgs %>%
      data.frame
  }
  colnames(x_tst_scores) <- colnames(x_scores)[1:nComp]  
  x_tst_scores$individualID <- tst_lb
  
  scores = list(train = pls_tr_data, test = x_tst_scores, center = x_mean, scale = x_sd, weights =modpca)
  return(scores)
}

get_climate <- function(path, outpath){
  download_point_daymet(path = path, outpath = outpath)
  melt_daymet(path = outpath)
}
