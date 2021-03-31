remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.1, .9), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}

drop_nodata <- function(x){
  if(length(x[x<=0])>0){
    x[which(x<=0)] <- NA
  }
  return((x))
}
fappend_fia_predictions <- function(fld = "./outdir/2020/sp/", stat = "mean"){
  
  # load data, make average by plot, and appen
  library(tidyverse)
  #unzip(paste("./outdir/predictions_",fld,".zip", sep=""))
  pt = paste(fld, sep="")
  dat.list <- list.files(pt)
  
  #simple version (NO NA remove)
  final_plot_preds <- readr::read_csv(paste(pt,  dat.list[1], sep = "/")) 
  
  # final_plot_preds[9:ncol(final_plot_preds)] <- apply(final_plot_preds[9:ncol(final_plot_preds)], 2,
  #                                                        remove_outliers)
  # 
  final_plot_preds[,c("LAT", "LON")] <-  sapply(final_plot_preds[,c("LAT", "LON")],as.character) 
  # final_plot_preds <- final_plot_preds %>% group_by(plot) %>%
  #   summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
  # 
  final_plot_preds = final_plot_preds[c(1:10, 13:14,17:18,21:22,25:26,29:30,33:34,37:38)]
  summary(final_plot_preds)
  missing_groups <- list()
  tkn = 0
  for(ii in dat.list[-1]){
    tryCatch({
      foo <-  readr::read_csv(paste(pt,  ii, sep = "/")) 
      foo[,c("LAT", "LON")] <- sapply(foo[,c("LAT", "LON")],as.character) 
      
      #foo[9:ncol(foo)] <- apply(foo[9:ncol(foo)], 2, remove_outliers)
      # foo <- foo %>% group_by(plot) %>%
      #   summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.)))
      foo = foo[c(1:10, 13:14,17:18,21:22,25:26,29:30,33:34,37:38)]
      
      final_plot_preds <- rbind.data.frame(final_plot_preds, foo)
    }, error = function(e) {
      tkn= tkn + 1
      missing_groups[[tkn]] <- which(dat.list == ii)
    })
  }
  final_plot_preds = final_plot_preds[complete.cases(final_plot_preds),]
  write_csv(final_plot_preds, paste("./outdir/2020//FIA_2020_no_mlbs_sp.csv", sep=""))
  
}

