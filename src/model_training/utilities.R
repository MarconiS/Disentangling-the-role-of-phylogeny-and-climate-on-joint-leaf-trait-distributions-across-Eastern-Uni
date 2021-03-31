library(outliers)
# OutVals = apply(y_obs[-1], 2, function(x)head(sort(boxplot(x, plot=FALSE)$out, decreasing = T)))

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.05, .95), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

get_mod_r2 <- function(pred, obs){
  #1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs, na.rm=T))^2)
} 

# OutVals = apply(y_obs[-1], 2,  function(x)remove_outliers(x))
# summary(OutVals)


plot_spectra<-function(plt_dat){
  #plot reflectances
  plot_data <- plt_dat %>% 
    dplyr::select(-one_of(c("siteID", "taxonID",  "band_site","band_species", "flightpath"))) 
  plot_data <- plot_data %>% dplyr::select(-one_of("individualID")) %>%
    t %>%
    data.frame
  colnames(plot_data) = unlist(plt_dat[1]) # the first row will be the header
  plot_data <- data.frame(bnd = 1:dim(plot_data)[1], plot_data)
  ggdat <- tidyr::gather(plot_data, treeID, Reflectance,-bnd)
  
  return(ggplot(ggdat, aes(x = bnd, y = Reflectance)) + 
           geom_line(aes(color = factor(treeID), alpha= 1), size = 0.2) +
           theme_bw()+
           theme(legend.position="none"))
  
}

convert_stei <- function(dat){
  library(rgdal)
  what_to_keep <- which(dat$easting > 270000)
  tmp_keep <- dat[what_to_keep,]
  dat$utmZone <- "15N"
  dat$siteID <- "CHEQ"
  coordinates(dat) <- c("easting", "northing")
  proj4string(dat) <- CRS("+init=epsg:32616") # WGS 84
  dat <- spTransform(dat, CRS("+init=epsg:32615"))
  new_dat <- cbind(dat@data, dat@coords)
  new_dat[what_to_keep,] <- tmp_keep
  return(new_dat)
}

get_cophenetic_distance <- function(){
  library(ape)
  library(data.table)
  library(tidyverse)
  library(treeio)
  library(ggplot2)
  library(ggtree)
  library(rotl)
  options(stringsAsFactors=FALSE)
  
  
  tree <- ape::read.tree(file="~/Documents/Data/OpenTree/opentree10.4_tree/labelled_supertree/labelled_supertree_ottnames.tre") # just using a small one for ease
  
  #tree <- ape::read.tree(file="~/Documents/Chapter2/")
  
  NEON_species <- read_csv("./TOS_retriever/out/field_data.csv") %>%
    dplyr::filter(height > 3) %>%
    dplyr::select(scientificName) %>% unique
  FIA_species <- read_csv("~/Documents/Data/FIA/Chapter2_FIA_nm.csv")
  full_list <- read_csv("./indir/Full_list2.csv")
  
  resolved_names_1 <- tnrs_match_names(full_list$species_sc[1:250]) %>%
    select(search_string, unique_name, ott_id)
  resolved_names_2 <- tnrs_match_names(full_list$species_sc[251:nrow(full_list)]) %>%
    select(search_string, unique_name, ott_id)
  
  resolved_names_2 <- resolved_names_2[complete.cases(resolved_names_2), ]
  resolved_names_1 <- resolved_names_1[complete.cases(resolved_names_1), ]
  
  otts_names <- rbind.data.frame(resolved_names_1, resolved_names_2)
  my_tree <- tol_induced_subtree(ott_ids = otts_names$ott_id)
  
  #not sure if this needed any more
  NEON_species <- strsplit(unlist(NEON_species), split = " ")
  NEON_species  <- do.call(rbind.data.frame, NEON_species)
  NEON_species <- lapply(1:nrow(NEON_species), function(x) paste(NEON_species[x,1:2], collapse = "_"))
  NEON_species <- do.call(rbind.data.frame, NEON_species)
  colnames(NEON_species) <- "speciesName"
  
  FIA_extra <- FIA_species$sp_name[!FIA_species$sp_name %in% unlist(NEON_species)]
  all_species <- c(unlist(FIA_extra), unlist(NEON_species)) %>% unique
  all_species[all_species %in% c("Aquifoliaceae_spp.", "Bourreria_spp.", "Quercus_spp.",  "Ulmus_spp.")] = NA
  all_species[all_species %in% c("Tree_broadleaf", "Tree_evergreen", "Tree_unknown", "Unknown_hardwood", "Unknown_plant")] = NA
  all_species[all_species == "Betula_glandulosa/nana"] = "Betula_glandulosa"
  all_species[all_species == "Carya_carolinae-septentrionalis"] = "Carya_carolinae"
  all_species <- all_species[complete.cases(all_species)]
  
  resolved_names <- tnrs_match_names(NEON_species)
  
  node = lapply(1:length(all_species), function(x) tree$tip.label[tree$tip.label %like% all_species[x]])
  
  nnodes = NULL
  for(ii in 1:length(node)){
    nnodes = c(nnodes, node[[ii]][1])
  }
  write_csv(data.frame(nnodes), "./indir/all_clades.csv")
  
  foo <- strsplit(nnodes, split = "_")
  foo <- do.call(rbind.data.frame, foo)[1:2]
  colnames(foo) <- c("genus", "species")
  foo$nnodes <- nnodes
  foo <- foo[complete.cases(foo),]
  foo <- foo %>% group_by(genus, species) %>%  filter(row_number()==1)
  
  sub_tree <- keep.tip(tree, foo$nnodes)
  
  
  #get species names
  tip_names <- sub_tree$tip.label %>%
    strsplit(split = "_") 
  tip_names <- do.call(rbind.data.frame, tip_names)
  colnames(tip_names) <- 1:ncol(tip_names)
  
  new_labels <- lapply(1:nrow(tip_names), function(x) paste(tip_names[x,1:2], collapse = "_"))
  new_labels <- do.call(rbind.data.frame, new_labels)
  colnames(new_labels) <- "speciesName"
  # pairwise distance matrix
  
  tree <- compute.brlen(sub_tree, 1)
  PatristicDistMatrix<-cophenetic.phylo(tree)
  PatristicDistMatrix
  #PatristicDistMatrix2<- distTips(tree)
  covphilo <- (max(PatristicDistMatrix) - PatristicDistMatrix ) / max(PatristicDistMatrix)
  
  #covphilo <-cophenetic(sub_tree)
  nms <- colnames(covphilo) %>%
    strsplit(split = "_")
  nms <- do.call(rbind.data.frame, nms)[1:2]
  colnames(nms) <- c("genus", "species")
  #nms$species[nms$species == "sp."] <- ""
  new_labels <- lapply(1:nrow(nms), function(x) paste(nms[x,1:2], collapse = " "))
  new_labels <- do.call(rbind.data.frame, new_labels)
  colnames(new_labels) <- "species_sc"
  
  newnames <- left_join(new_labels, full_list, by = "species_sc")
  colnames(covphilo) <- newnames$taxonID
  rownames(covphilo)<- newnames$taxonID
  saveRDS(covphilo, "./indir/phylogenic_covariance.rds")
}



merge_phylogenies <- function(my_tree, tree_2, root){
  library(ape)
  library(stringi)
  
  #get phylogeny for spermatophyta
  lemay <- ape::read.tree(file="~/Documents/Chapter2/trees/lemay_2012tre.txt")
  parfrey <- ape::read.tree(file="~/Documents/Chapter2/trees/parfrey.txt")
  Qiu <- ape::read.tree(file="~/Documents/Chapter2/trees/Qiu.tre.txt")
  ruhfel <- ape::read.tree(file="~/Documents/Chapter2/trees/ruhfeltre.txt")
  #contains most/all gynosperms
  leslie <- ape::read.tree(file="~/Documents/Chapter2/trees/leslietre.txt")
  
  
  #broadleaves
  foster <- ape::read.tree(file="~/Documents/Chapter2/trees/foster.tre.txt")
  
  n = nchar(root)
  
  
  #get_list_study
  
  lspcs <- gsub(" ", "", full_list$species_sc, fixed = TRUE)
  lspcs <- gsub("sp.", "",lspcs, fixed = TRUE)
  
  #tree_1 <- ape::read.tree(file="~/Documents/Chapter2/trees/pg_412.tre.txt")
  cut_1 <- ruhfel[[2]]
  cut_1 <- keep.tip(cut_1, cut_1$tip.label[startsWith(cut_1$tip.label, root)])
  stri_sub(cut_1$tip.label, n+1, n) <- " " 
  lspcs[!lspcs %in% leslie$tip.label]
  
  cut_2 <-  merge_trees <- ape::consensus(ruhfel[[1]], ruhfel[[2]])
  cut_2 <- keep.tip(cut_2, cut_2$tip.label[startsWith(cut_2$tip.label, root)])
  #rename tip the way tou want it
  stri_sub(cut_2$tip.label, n+1, n) <- " "
  cut_2$tip.label
  
  merge_trees <- ape::consensus(ruhfel[[1]], ruhfel[[2]])
  if(root == "Pinus"){
    main_root <- "Pinus subgenus Pinus"
  }
  tree_1 <- drop.tip(my_tree, 
        tip =  my_tree$tip_label[startsWith(my_tree$tip.label, cut_2$tip.label)])
  plot(cut_2)
  merged_tree <- bind.tree(tree_1, cut_2, where = main_root)
}



full_list$species_sc %>% sort





