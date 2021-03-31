get_species_effect_sigma <- function(data = data, nex_tree = "./indir/phylogeny/chapter2_nexus.txt", dat_pt = "./TOS_retriever/out/field_traits_dataset.csv",
                                     dict_pt = "./indir/phylogeny/Full_list.csv", unkn = "./indir/phylogeny/unknowns.csv"){
  library(tidyverse)
  library(ape)
  library(castor)
  library(adephylo)
  unknown <- readr::read_csv(unkn)
  dictionary <- readr::read_csv(dict_pt)
  field_Data <- data %>%
   select(scientificName, taxonID)  
  field_Data$scientificName[field_Data$scientificName == "Carya ovalis"] = "Carya ovata" 
  field_Data$scientificName[field_Data$scientificName == "Acer barbatum"] = "Acer saccharum" 
  field_Data$scientificName[field_Data$scientificName == "Acer saccharinum"] = "Acer saccharum" 
  field_Data$scientificName[field_Data$scientificName == "Carya ovalis"] = "Carya ovata" 
  
  #print(paste(field_Data$scientificName, collapse = ","))
  field_links <- strsplit(field_Data$scientificName, split = "\\s+")
  field_links  <- do.call(rbind, field_links)[,1:2]
  field_links[field_links[,2]=="sp.",2] <- ""
  field_links <- apply(field_links,1, function(x)paste(x, collapse = " " ))  
  field_links <- trimws(field_links)
  
  field_links <- field_links %>% data.frame(stringsAsFactors = F)
  colnames(field_links) <- "species_sc"
  #field_data <- cbind(field_links, field_Data$taxonID) %>% unique
  #field_data <- field_data[!duplicated(field_data[,1]),]
  
  #get taxonID to species link from USDA US tree species

  new_lbls <- left_join(field_links, dictionary)
  lbls_na <- new_lbls[is.na(new_lbls$taxonID), "species_sc"]
  new_lbls[new_lbls$species_sc=="Carya ovalis", ] <- c("Carya ovata", 407, "CAOV2")
  # load phylogeny tree of data
  tree <- ape::read.nexus(nex_tree)
  #tree <- ape::read.nexus("~/Desktop/phyloT_2_nexus.txt")
  tree$node.label <- paste(1:length(tree$node.label))
  #tree$node.label<-paste(1:171)
  #distances <- get_all_node_depths(tree, as_edge_count=FALSE)#/trr$max_distance
  #d2 <-     get_all_distances_to_root(tree, as_edge_count=FALSE)
  tree <- compute.brlen(tree, 1)
  PatristicDistMatrix<-cophenetic.phylo(tree)
  #PatristicDistMatrix2<- distTips(tree)
  pat_corr_mat <- (max(PatristicDistMatrix) - PatristicDistMatrix ) / max(PatristicDistMatrix)
  rename_cols <- (strsplit(colnames(pat_corr_mat), split = "_"))

  rename_cols <- do.call(rbind, rename_cols)
  rename_cols <- apply(rename_cols,1, function(x) paste(x, collapse = " ")) %>%
    data.frame(stringsAsFactors = F)
  colnames(rename_cols) <- "species_sc"
  
  new_names <- left_join(rename_cols, dictionary, by = "species_sc") 
  new_names[is.na(new_names$taxonID),]
  colnames(pat_corr_mat) <- rownames(pat_corr_mat) <- new_names$taxonID
  
  return(list(cov_mat = pat_corr_mat, dict = dictionary, data_lbls = new_lbls, in_lbls = field_links))
}
