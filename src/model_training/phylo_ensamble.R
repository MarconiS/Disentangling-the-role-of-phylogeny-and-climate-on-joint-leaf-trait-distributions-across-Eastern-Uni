## from http://tr.im/hH5A
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

get_closest_species <- function(raw_dat, cov_ranef, list_train, nclose = 3){
  # raw_dat$taxonID[raw_dat$taxonID == "DIVIS"] = "DIVI5"
  # raw_dat$taxonID[raw_dat$taxonID == "QUERCUS"] = "QULA2"
  list_to_predict <- raw_dat$taxonID %>% unique
  unknown_species <- list_to_predict[!(list_to_predict %in% list_train)]
  sample_cov_phylo <- cov_ranef$taxonID[rownames(cov_ranef$taxonID) %in% list_train, 
                                        colnames(cov_ranef$taxonID) %in% unknown_species] #%>%
  #data.frame
  take_five <- list()
  for(sp in 1:length(raw_dat$taxonID)){
    if(raw_dat$taxonID[sp] %in% colnames(sample_cov_phylo)){
      weights <- tail(sort(sample_cov_phylo[,(raw_dat$taxonID[sp])]),nclose)
      taxID <- rownames(sample_cov_phylo)[which(sample_cov_phylo[,(raw_dat$taxonID[sp])] %in% weights)]
      foo <- data.frame(mefa:::rep.data.frame(raw_dat[sp,], length(taxID)), 
                        weight = sample_cov_phylo[taxID, raw_dat$taxonID[sp]])
      foo$taxonID <- taxID
      foo$weight <- softmax(foo$weight)
      take_five[[sp]] <- foo
    }else{
      take_five[[sp]] <- data.frame(raw_dat[sp,], weight = 1)
    }
  }
  res <- do.call(rbind.data.frame, take_five) 
  return(res)
}

