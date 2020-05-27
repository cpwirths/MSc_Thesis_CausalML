#functions compute the DML output for the simulations: point estimates and confidence intervals
source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/output_dml.R")

summary_rmse_sd <- function(input, sd=T){
  sd_rmse <- rep(0, (ncol(input)))
  if(sd==T){
    for(i in 1:ncol(input)){
      f <- !is.na(input[,i])
      sd_rmse[i] <- sqrt(mean(((input[f,1] - input[f,i])/input[f,2])^2))
    }
  }else{
    for(i in 1:ncol(input)){
      f <- !is.na(input[,i])
      sd_rmse[i] <- sqrt(mean(((input[f,1] - input[f,i]))^2))
    }
  }
 
  names(sd_rmse) <- names(input)
  return(round(sd_rmse[3:ncol(input)], digits=3))
}

output_dml_splits <- function(r,n=100, error, critic_v=1.96){
  library(aciccomp2016)
  output_raw <- r[,colnames(r) %in% c("table")] 
  output <- list()
  for(i in 1:n){
    output[[i]] <- output_p(output_raw[,i], adj=T)
  }
  
  DML_Lasso <- rep(0, n)
  DML_Lasso_lb  <- rep(0, n)
  DML_Lasso_ub  <- rep(0, n)
  
  DML_Trees <- rep(0, n)
  DML_Trees_lb <- rep(0, n)
  DML_Trees_ub <- rep(0, n)
  
  DML_Boosting <- rep(0, n)
  DML_Boosting_lb <- rep(0, n)
  DML_Boosting_ub <- rep(0, n)
  
  DML_Forest <- rep(0, n)
  DML_Forest_lb <- rep(0, n)
  DML_Forest_ub <- rep(0, n)
  
  DML_Nnet <- rep(0, n)
  DML_Nnet_lb <- rep(0, n)
  DML_Nnet_ub <- rep(0, n)
  
  DML_Bart<- rep(0, n)
  DML_Bart_lb  <- rep(0, n)
  DML_Bart_ub  <- rep(0, n)
  
  for(i in 1:n){
    
    DML_Lasso[i] <-  output[[i]]$output[[1]]["Median ATE","Lasso"] 
    DML_Lasso_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Lasso"]-output[[i]]$output[[1]][error,"Lasso"]*critic_v
    DML_Lasso_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Lasso"]+output[[i]]$output[[1]][error,"Lasso"]*critic_v
    
    DML_Trees[i] <-  output[[i]]$output[[1]]["Median ATE","Trees"] 
    DML_Trees_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Trees"]-output[[i]]$output[[1]][error,"Trees"]*critic_v
    DML_Trees_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Trees"]+output[[i]]$output[[1]][error,"Trees"]*critic_v
    
    DML_Boosting[i] <-  output[[i]]$output[[1]]["Median ATE","Boosting"] 
    DML_Boosting_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Boosting"]-output[[i]]$output[[1]][error,"Boosting"]*critic_v
    DML_Boosting_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Boosting"]+output[[i]]$output[[1]][error,"Boosting"]*critic_v
    
    DML_Forest[i] <-  output[[i]]$output[[1]]["Median ATE","Forest"] 
    DML_Forest_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Forest"]-output[[i]]$output[[1]][error,"Forest"]*critic_v
    DML_Forest_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Forest"]+output[[i]]$output[[1]][error,"Forest"]*critic_v
    
    DML_Nnet[i] <-  output[[i]]$output[[1]]["Median ATE","Nnet"]  
    DML_Nnet_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Nnet"]-output[[i]]$output[[1]][error,"Nnet"]*critic_v
    DML_Nnet_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Nnet"]+output[[i]]$output[[1]][error,"Nnet"]*critic_v
    
    DML_Bart[i] <-  output[[i]]$output[[1]]["Median ATE","Bart"] 
    DML_Bart_lb[i] <-  output[[i]]$output[[1]]["Median ATE","Bart"]-output[[i]]$output[[1]][error,"Bart"]*critic_v
    DML_Bart_ub[i] <-  output[[i]]$output[[1]]["Median ATE","Bart"]+output[[i]]$output[[1]][error,"Bart"]*critic_v
    
  
    }
  estimate <- as.data.frame(cbind(DML_Lasso, DML_Trees, DML_Boosting, DML_Forest, DML_Nnet, DML_Bart))
  estimate_lb <- as.data.frame(cbind(DML_Lasso_lb, DML_Trees_lb, DML_Boosting_lb, DML_Forest_lb, DML_Nnet_lb, DML_Bart_lb)) 
  estimate_ub <- as.data.frame(cbind(DML_Lasso_ub, DML_Trees_ub, DML_Boosting_ub, DML_Forest_ub, DML_Nnet_ub, DML_Bart_ub))
  
  return(list("estiamte"=estimate, "estimate_lb"=estimate_lb, "estimate_ub"=estimate_ub ))
}

coverage_m <- function(lb_m, ub_m, true){
  
  coverage <- rep(0, (ncol(lb_m)))
  for(i in 1:ncol(lb_m)){
    coverage[i] <- table(lb_m[,i] < true & true < ub_m[,i])["TRUE"]/nrow(lb_m)
    if(is.na(coverage[i])){coverage[i] <- 0}
  }
  names(coverage) <- names(lb_m)
  return(coverage)
}

coverage_per_sim <- function(l,u,t, sim, r){
  n <- nrow(l)
  coverage <- rep(0, n)
  for(i in 1:n){
    coverage[i] <- table(l[i,] < t[i,] & t[i,] < u[i,])[2]/length(u[i,])
    if(is.na(coverage[i])){coverage[i] <- 0}
  }
  return(coverage)
}


coverage <- function(l,u,t, sim, r){
  return(table(l < t & t < u)[2]/(sim*r))
}


rmse_s <- function(true, est, sd=NULL){
  if(is.null(sd)){
    return(sqrt(mean((true-est)^2)))
  }else{
    return(sqrt(mean(((true-est)/sd)^2)))
  }
}

coverage_s <- function(lb, ub, true){
  coverage <- 0
  coverage <- table(lb < true & true < ub)["TRUE"]/length(lb)
  if(is.na(coverage)){coverage <- 0}
  return(coverage)
}
