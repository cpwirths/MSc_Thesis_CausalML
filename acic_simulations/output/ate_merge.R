# combines the ATE estimators from all methods and computes RMSE, coverage, bias and interval length

rm(list=ls())
library("MLmetrics")
library("hydroGOF")
library("aciccomp2016")
source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/output_hte_ge_simulation.R")
source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/output_dml_simulation.R")

setwd("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/output")

# Initialize results output -----------------------------------------------

analyze <- seq(77) #specifiy knobs 
knobs_all <- read.csv("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/knobs_all.csv", sep=";")
knobs <- knobs_all[analyze,] 

names_methods <- c("OLS","DML Lasso", "DML Trees", "DML Boost", "DML RF","DML Nnet", "DML BART", "Causal Forest", "GE Lasso" ,"GE Trees" ,"GE Boost","GE RF","GE Nnet", "GE BART", "BART MC") #methods that are used to estimate the ATE

rmse_matrix <- data.frame(matrix(NA, length(analyze), length(names_methods))) #each row represents the rmse for all methods per knob
coverage_matrix <- data.frame(matrix(NA, length(analyze), length(names_methods))) #each row represents thecoverage for all methods per knob
bias_matrix <- data.frame(matrix(NA, length(analyze), length(names_methods))) #each row represents the bias for all methods per knob
interval_matrix <- data.frame(matrix(NA, length(analyze), length(names_methods))) #each row represents the interval length of the confidence interval per knob for all methods per knob

colnames(rmse_matrix) <-names_methods 
colnames(coverage_matrix) <- names_methods
colnames(bias_matrix) <- names_methods
colnames(interval_matrix) <- names_methods

rep <- 100 #number of simulations per knob
alpha <- 0.05
critic_v <- qnorm(1-alpha/2) #critical value for confidence intervals

# Compute bias, rmse, coverage and interval length per knob  -------------------------------
for(k in analyze){
  print(k)
  
  #load true ate, ols estimate and standard deviation of the response variable
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/acic_",k,"_data.RData", sep="",collapse = NULL))
  ate_matrix <- data.frame(true_ate_mu[1:rep])
  ate_matrix$sd_response <- sd_response[1:rep]
  ate_matrix$ols <- ate_ols_coeff[1:rep]
  ols_lb <- ate_ols_sum[1:rep,1]-ate_ols_sum[1:rep,2]*critic_v
  ols_ub <- ate_ols_sum[1:rep,1]+ate_ols_sum[1:rep,2]*critic_v
  lb_matrix <- ols_lb
  ub_matrix <- ols_ub
  ols_interval <- ub_matrix - lb_matrix
  
  #load DML output
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/dml/",k,"_dml_acic.RData", sep="",collapse = NULL))
  dml <- output_dml_splits(r, n=rep, error="se(adj)", critic_v=critic_v)
  dml_interval <- dml$estimate_ub - dml$estimate_lb
  
  #load ATE estimate of Causal Forest
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/cf/",k,"_cf_acic_results.RData", sep="",collapse = NULL))
  ate_cf_o <- rep(0,rep)
  ate_cf_lb_o <- rep(0,rep)
  ate_cf_ub_o <- rep(0,rep)
  
  for (i in 1:rep){
    ate_cf_o[i] <- r[, "cf_ate_o_raw"][[i]]["estimate"]
    ate_cf_lb_o[i] <- r[, "cf_ate_o_raw"][[i]]["estimate"] - r[, "cf_ate_o_raw"][[i]]["std.err"]*critic_v
    ate_cf_ub_o[i] <- r[, "cf_ate_o_raw"][[i]]["estimate"] + r[, "cf_ate_o_raw"][[i]]["std.err"]*critic_v
  }
  cf_interval_o <- ate_cf_ub_o - ate_cf_lb_o
  
  #load ATE estimate of the Generic method
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/generic/",k,"_generic_acic.RData", sep="",collapse = NULL))
  gates_generic_output <- output_generic(r)
  generic_ate <- gates_generic_output$ATE$ATE_est[1:rep,]
  generic_ate_lb <- gates_generic_output$ATE$ATE_lb[1:rep,]
  generic_ate_ub <- gates_generic_output$ATE$ATE_ub[1:rep,]
  generic_ate_interval <- generic_ate_ub - generic_ate_lb

  #load ATE estimate of BART MChains
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/bart_mc/",k,"_bart_acic_results.RData", sep="",collapse = NULL))
  sd_response_b <- rep(NA,rep)
  true_sd <- rep(NA,rep)
  bart_ate_sd <- rep(NA,rep)
  bart_ate_sd_lb <- rep(NA,rep)
  bart_ate_sd_ub <- rep(NA,rep)

  for(i in 1:rep){  
    true_sd[i] <- mean(true_cate_mu[,i][which(r[, "ate_sd"][[i]]$idx_sd==T)])
    sd_response_b[i] <- sd(r[, "bart_Y"][[i]][which(r[, "ate_sd"][[i]]$idx_sd==T)])  
    bart_ate_sd[i] <- as.numeric(r[, "ate_sd"][[i]]$bart_sd_ate)
    bart_ate_sd_lb[i] <- as.numeric(r[, "ate_sd"][[i]]$bart_sd_ate_lower)
    bart_ate_sd_ub[i] <- as.numeric(r[, "ate_sd"][[i]]$bart_sd_ate_upper)
  }  
  
  rmse_bart_sd_mc <-  rmse_s(true_sd, bart_ate_sd, sd_response_b) 
  coverage_bart_sd_mc <- coverage_s(bart_ate_sd_lb, bart_ate_sd_ub, true_sd)
  bias_bart_sd_mc <-  mean((true_sd-bart_ate_sd)/sd_response_b)
  interval_bart_sd_mc <- mean((bart_ate_sd_ub-bart_ate_sd_lb)/sd_response_b)

  #combine all estimates and calculate bias, rmse, coverage and interval length 
  ate_matrix <- data.frame(cbind(ate_matrix,dml$estiamte, ate_cf_o, generic_ate))
  lb_matrix <- data.frame(cbind(ols_lb,dml$estimate_lb,ate_cf_lb_o,  generic_ate_lb))
  ub_matrix <- data.frame(cbind(ols_ub,dml$estimate_ub, ate_cf_ub_o, generic_ate_ub))
  bias <-  (ate_matrix[,1] - ate_matrix[3:length(ate_matrix)])/sd_response[1:rep] 
  intervals <- data.frame(cbind(ols_interval,dml_interval, cf_interval_o, generic_ate_interval))
  intervals <-  intervals/sd_response[1:rep]
  
  rmse_sd_all <- c(summary_rmse_sd(ate_matrix, sd=T), round(rmse_bart_sd_mc, digits=3)) 
  coverage_all <- c(coverage_m(lb_matrix, ub_matrix, ate_matrix$true_ate_mu), coverage_bart_sd_mc)
  interval_all <- c(round(colMeans(intervals),3), round(interval_bart_sd_mc,3))
  bias_all <- c(round(colMeans(bias),4), round(bias_bart_sd_mc,4))
  
  index <- which(analyze==k)
  rmse_matrix[index,] <- rmse_sd_all
  coverage_matrix[index,] <-  coverage_all
  bias_matrix[index,] <- bias_all
  interval_matrix[index,] <-interval_all 
}

#save as finaloutput rmse, covergae, bias and interval matrix
save(rmse_matrix, coverage_matrix, bias_matrix, interval_matrix, file="ate_final_output.RData")


