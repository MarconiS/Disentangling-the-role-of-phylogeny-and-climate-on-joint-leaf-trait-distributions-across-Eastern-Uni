remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.05, .95), na.rm = na.rm, ...)
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
# load data, make average by plot, and append
library(tidyverse)
#dat.list <- list.files("./outdir/pred_species_csv/")
pt = "./outdir/2020/full_no_dl//"

dat.list <- list.files(pt)
#dat.list <- dat.list[68:length(dat.list)]

#simple version (NO NA remove)
final_plot_preds <- readr::read_csv(paste(pt,  dat.list[1], sep = "")) 
#summary(final_plot_preds)
# final_plot_preds[11:ncol(final_plot_preds)-1] <- apply(final_plot_preds[11:ncol(final_plot_preds)-1], 2, 
#                                                        remove_outliers)
#summary(final_plot_preds)
# final_plot_preds[13:ncol(final_plot_preds)-1] <- apply(final_plot_preds[13:ncol(final_plot_preds)-1], 2,
#                                                        remove_outliers)
#final_plot_preds <- apply(final_plot_preds, 2,  function(x)remove_outliers(x))
# final_plot_preds <- final_plot_preds %>% group_by(plot) %>%
#   summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else first(.)))
# 
# summary(final_plot_preds)
#for(ii in dat.list[-c(1, 362, 736, 1024)]){

missing_groups <- list()
tkn = 0
for(ii in dat.list[-1]){ #518 - 538
  tryCatch({
    foo <-  readr::read_csv(paste(pt,  ii, sep = "")) 
    # foo[13:ncol(foo)-1] <- apply(foo[13:ncol(foo)-1], 2, remove_outliers)
    # foo <- foo %>% group_by(plot) %>%
    #   summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else first(.)))
    final_plot_preds <- rbind.data.frame(final_plot_preds, foo)
  }, error = function(e) {
    tkn= tkn + 1
    missing_groups[[tkn]] <- which(dat.list == ii)
  })
}
write_csv(final_plot_preds, "./outdir/2020/FIA_full_nodl_map.csv")




# #simple version (NO NA remove)
# final_plot_preds <- readr::read_csv(paste("~/Documents/manuscript_ongoing/FIA_predictions/pred_full_csv/", 
#                                           dat.list[1], sep = ""))
# final_plot_preds <- final_plot_preds %>% group_by(LAT, LON) %>%
#   summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else first(.)))
# for(ii in dat.list[752:1151]){
#   foo <-  readr::read_csv(paste("~/Documents/manuscript_ongoing/FIA_predictions/pred_full_csv/", 
#                                 ii, sep = ""))
#   foo <- foo %>% group_by(LAT, LON) %>%
#     summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else first(.)))
#   final_plot_preds <- rbind.data.frame(final_plot_preds, foo)
# }
# write_csv(final_plot_preds, "./outdir/outdir/FIA_full_map.csv")
