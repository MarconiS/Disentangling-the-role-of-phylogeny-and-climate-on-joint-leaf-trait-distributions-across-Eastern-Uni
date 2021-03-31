library(tidyverse)
library(brms)
full = readRDS("./models/ch2perc_student_full_noml.rds")
predictions = predict(full$mod, newdata = test_data)
predictions_dims <- list()
for(tr in 1:dim(predictions)[3]){
  tr_lb <- dimnames(predictions)[[3]][tr]
  rescaled_tr <- predictions[, 1, tr] * full$scale$scale[[tr]] + full$scale$mean[[tr]]
  rescaled_tr <- exp(rescaled_tr)
  predictions_dims[[tr_lb]] <- rescaled_tr
}
predictions_dims <- do.call(cbind.data.frame, predictions_dims)

LES = cor(predictions_dims)
colnames(LES)= rownames(LES) = c("N%", "C%", "lignin%", "cell%", "LMA", "ChlA", "ChlB", "Carot")
corrplot::corrplot(LES, method = "circl", type="lower", order="hclust", tl.srt=45, addCoef.col = "black",
                   sig.level = 0.01, diag = F)
#get residuals correlation
res_preds = predictions_dims
for (tr in c(1:4,6:8)){
  mlm1=lm(as.formula(paste(colnames(predictions_dims)[tr], "~",
                     "leafMassPerArea",sep = "")),
    data=predictions_dims)
  test_res = residuals(mlm1)
  res_preds[tr] = test_res
}


spellman.cor = cor(LES,use="pairwise.complete.obs")
spellman.dist <- as.dist(1 - LES)
spellman.tree <- hclust(spellman.dist, method="complete")
plot(spellman.tree)

resLES = cor(res_preds)
colnames(resLES)= rownames(resLES) = c("N%", "C%",  "lignin%", "cell%","LMA", "ChlA", "ChlB", "Carot")
resLES=resLES[spellman.tree$order,spellman.tree$order]
corrplot::corrplot(resLES, method = "circle", type="lower", order="original", tl.srt=45, addCoef.col = "black",
                   sig.level = 0.01, diag = F)




spellman.cor = cor(resLES,use="pairwise.complete.obs")
spellman.dist <- as.dist(1 - resLES)
spellman.tree <- hclust(spellman.dist, method="complete")
plot(spellman.tree)
