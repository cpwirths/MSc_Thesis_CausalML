#Function returning the ATE and GATE estimate for the simulation
bart_estimates <- function(data,y,d){
  source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/bart_group_effects.R")
  response = data[,y]
  treatment = data[,d]
  confounders = data[,!colnames(data)%in%c(y,d)]
  
  bart_sd <- bartc(response, treatment, confounders, estimand=c("ate"),commonSup.rule=c("sd")) #Overlap is addressed by the standard deviation rule
  bart_sd_ate <- summary(bart_sd)$estimates["estimate"]
  bart_sd_ate_lower <- summary(bart_sd)$estimates["ci.lower"]
  bart_sd_ate_upper <- summary(bart_sd)$estimates["ci.upper"]
  idx_sd <- bart_sd$commonSup.sub
  
  if(!is.null(bart_sd$samples.indiv.diff)){
    samples.indiv.diff <- bart_sd$samples.indiv.diff  
  }else{
    samples.indiv.diff <- bart_sd$fit.rsp$yhat.train - bart_sd$fit.rsp$yhat.test
    for(i in 1:dim(samples.indiv.diff)[1]){
      for(j in 1:dim(samples.indiv.diff)[2])
        samples.indiv.diff[i,j, which(treatment==0)] <- bart_sd$fit.rsp$yhat.test[i,j, which(treatment==0)] - bart_sd$fit.rsp$yhat.train[i,j, which(treatment==0)]
    }
    
    samples.indiv.diff <- aperm(samples.indiv.diff,perm=c(3,1,2))
  }
  
  bart_gate_sd <- gates_bart(samples.indiv.diff[idx_sd,,])
  
  ate_sd <- list("bart_sd_ate"=bart_sd_ate, "bart_sd_ate_lower"=bart_sd_ate_lower, "bart_sd_ate_upper"=bart_sd_ate_upper, "bart_gate_sd"=bart_gate_sd, "idx_sd"=idx_sd)

  return(list("bart_Y"=response, "bart_D"=treatment,  "ate_sd"=ate_sd))
}
