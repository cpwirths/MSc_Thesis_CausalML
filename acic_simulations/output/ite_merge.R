# combines the ITE estimators from all methods and computes PEHE
rm(list=ls())
library(Metrics)
library(MLmetrics)

setwd("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/output")

# Initialize results output -----------------------------------------------

analyze <- seq(1:77)

methods <- c("DR MOM RF", "Causal Forest", "DR MOM Lasso", "DR MOM Boost", "DR MOM Trees", "DR MOM Nnet", "DR MOM BART",  "BART MChains") 
pehe_matrix <- matrix(NA, length(analyze), length(methods))
colnames(pehe_matrix) <- methods
knob <- rep(0, length(analyze))

rep <- 100

# Compute PEHE ------------------------------------------------------------

for(k in analyze){
  print(k)
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/dr_mom_ites/",k,"_dr_mom_ites_acic.RData", sep="",collapse = NULL))
  r <- r[1:rep,"iates_mat"]
  names <- colnames(r[[1]])[1:length(methods)]
  mse_sd <- matrix(NA, (length(names)), length(r))
  rownames(mse_sd) <- names
  for(i in 1:length(r)){
    for(j in 1:length(names)){
      all_methods <- r[[i]]
      mse_sd[j,i] <-  rmse(all_methods[,"tau_val_mu"], all_methods[,names[j]])/sd(all_methods[,"y_val"])
    }
  }
  index <- which(analyze==k)
  knob[index] <-  k
  pehe_matrix[index,] <- rowMeans(mse_sd)
  
}

pehe_matrix <- as.data.frame(round(pehe_matrix,2))

save(pehe_matrix, file="ite_final_output.RData")
