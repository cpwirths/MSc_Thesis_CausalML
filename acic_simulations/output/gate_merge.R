# combines the GATE estimators from all methods and computes RMSE, coverage, bias and interval length

rm(list=ls())
library("MLmetrics")
library("hydroGOF")
library("aciccomp2016")
library("gridExtra")
library("grid")
library("ggplot2")
source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/output_hte_ge_simulation.R")

setwd("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/output")

# Initialize results output -----------------------------------------------

all_groups <- c("G1", "G2", "G3", "G4", "G5")
low <- "G1"
high <- "G5"
groups <- c(low, high)
analyze <- seq(1:77)

methods_gates <- c("GE Lasso" ,"GE Trees" ,"GE Boost","GE RF","GE Nnet", "GE BART", "Causal Forest", "BART MC")
rep <- 100

low_matrix_rmse <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
low_matrix_coverage <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
low_matrix_interval <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
low_matrix_bias <- data.frame(matrix(NA, length(analyze), length(methods_gates)))

high_matrix_rmse <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
high_matrix_coverage <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
high_matrix_interval <- data.frame(matrix(NA, length(analyze), length(methods_gates)))
high_matrix_bias <- data.frame(matrix(NA, length(analyze), length(methods_gates)))

colnames(low_matrix_rmse) <- methods_gates
colnames(low_matrix_coverage) <- methods_gates
colnames(low_matrix_bias) <- methods_gates
colnames(low_matrix_interval) <- methods_gates

colnames(high_matrix_rmse) <- methods_gates
colnames(high_matrix_coverage) <- methods_gates
colnames(high_matrix_bias) <- methods_gates
colnames(high_matrix_interval) <- methods_gates


# Compute bias, rmse, coverage and interval length per knob  -------------------------------
for(k in analyze){
  print(k)
  #true
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/acic_",k,"_data.RData", sep="",collapse = NULL))
  gates_true <- gates_persim(true_cate_mu)
  
  #generic
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/generic/",k,"_generic_acic.RData", sep="",collapse = NULL))
  gates_generic_output <- output_generic(r) 
  gates_generic_est <- gates_generic_output$GATE$GATE_hl
  gates_generic_lb <- gates_generic_output$GATE$GATE_hl_lb
  gates_generic_ub <- gates_generic_output$GATE$GATE_hl_ub
  
  #Causal Forest
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/cf/",k,"_cf_acic_results.RData", sep="",collapse = NULL))
  gates_cf <- matrix(NA, rep, length(all_groups))
  gates_cf_lb <- matrix(NA, rep, length(all_groups))
  gates_cf_ub <- matrix(NA, rep, length(all_groups))
  colnames(gates_cf) <- all_groups
  colnames(gates_cf_ub) <- all_groups
  colnames(gates_cf_lb) <- all_groups
  
  for (i in 1:rep){
    gates_cf[i,] <- r[, "cf_gates_o_raw"][[i]]["estimate",]
    gates_cf_lb[i,] <-  r[, "cf_gates_o_raw"][[i]]["estimate",] - r[, "cf_gates_o_raw"][[i]]["std.err",]*qnorm(1-0.025)
    gates_cf_ub[i,] <- r[, "cf_gates_o_raw"][[i]]["estimate",] + r[, "cf_gates_o_raw"][[i]]["std.err",]*qnorm(1-0.025)
  }
  
  #BART
  load(paste("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/bart_mc/",k,"_bart_acic_results.RData", sep="",collapse = NULL))
  sd_response_b <- rep(NA, rep)
  gates_bart_sd <- matrix(NA,rep,length(all_groups))
  gates_bart_lb_sd <- matrix(NA,rep,length(all_groups))
  gates_bart_ub_sd <- matrix(NA,rep,length(all_groups))
  colnames(gates_bart_sd) <- all_groups
  colnames(gates_bart_lb_sd) <- all_groups
  colnames(gates_bart_ub_sd) <- all_groups
  true_cates_sd <- true_cate_mu
  
  for(i in 1:rep){
    true_cates_sd[which(r[,"ate_sd"][[i]]$idx_sd==F),i] <- NA
    sd_response_b[i] <- sd(r[, "bart_Y"][[i]][which(r[, "ate_sd"][[i]]$idx_sd==T)])  
    gates_bart_sd[i,]  <- r[,"ate_sd"][[i]]$bart_gate_sd[,"est"]
    gates_bart_lb_sd[i,]  <- r[,"ate_sd"][[i]]$bart_gate_sd[,"ci.lower"]
    gates_bart_ub_sd[i,]  <- r[,"ate_sd"][[i]]$bart_gate_sd[,"ci.upper"]
  }
  gates_true_sd <- gates_persim(true_cates_sd) 
  
  array_all_gates <- array(matrix(NA,rep,length(groups)), dim=c(rep,length(groups),length(methods_gates)))
  dimnames(array_all_gates)[[2]] <- groups
  dimnames(array_all_gates)[[3]] <- c(dimnames(gates_generic_est)[[3]], "CF GATE", "BART MC GATE") 
  array_all_gates[,,1:6] <- gates_generic_est[1:rep,,]
  array_all_gates[,,7] <- gates_cf[1:rep,groups]
  array_all_gates[,,8] <- gates_bart_sd[1:rep,groups]
  
  array_all_gates_lb <- array(matrix(NA,rep,length(groups)), dim=c(rep,length(groups),length(methods_gates)))
  dimnames(array_all_gates_lb)[[2]] <- groups
  dimnames(array_all_gates_lb)[[3]] <- c(dimnames(gates_generic_lb)[[3]], "CF GATE", "BART MC GATE") 
  array_all_gates_lb[,,1:6] <- gates_generic_lb[1:rep,,]
  array_all_gates_lb[,,7] <- gates_cf_lb[1:rep,groups]
  array_all_gates_lb[,,8] <- gates_bart_lb_sd[1:rep,groups]
  
  array_all_gates_ub <- array(matrix(NA,rep,length(groups)), dim=c(rep,length(groups),length(methods_gates)))
  dimnames(array_all_gates_ub)[[2]] <- groups
  dimnames(array_all_gates_ub)[[3]] <- c(dimnames(gates_generic_ub)[[3]], "CF GATE", "BART MC GATE") 
  array_all_gates_ub[,,1:6] <- gates_generic_ub[1:rep,,]
  array_all_gates_ub[,,7] <- gates_cf_ub[1:rep,groups]
  array_all_gates_ub[,,8] <- gates_bart_ub_sd[1:rep,groups]
  
  index <- which(analyze==k)
  
  low_matrix_rmse[index,]  <- c(rmse_gates(array_all_gates[,,1:7], gates_true[1:rep,groups], sd_response[1:rep], sd=T)$rmse_gates[, low],rmse_gates_single(array_all_gates[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep], sd=T)$rmse_gates[low])  
  low_matrix_coverage[index,] <- c(coverage_gates(array_all_gates_lb[,,1:7], array_all_gates_ub[,,1:7], gates_true[1:rep,groups], sd_response[1:rep])$coverage_gates[, low],coverage_gates_single(array_all_gates_lb[,,8], array_all_gates_ub[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep])$coverage_gates[low])
  low_matrix_bias[index,] <- c(rmse_gates(array_all_gates[,,1:7], gates_true[1:rep,groups], sd_response[1:rep], sd=T)$bias_gates[, low],rmse_gates_single(array_all_gates[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep], sd=T)$bias_gates[low])  
  low_matrix_interval[index,] <- c(coverage_gates(array_all_gates_lb[,,1:7], array_all_gates_ub[,,1:7], gates_true[1:rep,groups], sd_response[1:rep])$interval_gates[, low],coverage_gates_single(array_all_gates_lb[,,8], array_all_gates_ub[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep])$interval_gates[low])
  
  
  high_matrix_rmse[index,]  <- c(rmse_gates(array_all_gates[,,1:7], gates_true[1:rep,groups], sd_response[1:rep], sd=T)$rmse_gates[, high],rmse_gates_single(array_all_gates[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep], sd=T)$rmse_gates[high])  
  high_matrix_coverage[index,] <- c(coverage_gates(array_all_gates_lb[,,1:7], array_all_gates_ub[,,1:7], gates_true[1:rep,groups], sd_response[1:rep])$coverage_gates[, high],coverage_gates_single(array_all_gates_lb[,,8], array_all_gates_ub[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep])$coverage_gates[high])
  high_matrix_bias[index,] <- c(rmse_gates(array_all_gates[,,1:7], gates_true[1:rep,groups], sd_response[1:rep], sd=T)$bias_gates[, high],rmse_gates_single(array_all_gates[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep], sd=T)$bias_gates[high])  
  high_matrix_interval[index,] <- c(coverage_gates(array_all_gates_lb[,,1:7], array_all_gates_ub[,,1:7], gates_true[1:rep,groups], sd_response[1:rep])$interval_gates[, high],coverage_gates_single(array_all_gates_lb[,,8], array_all_gates_ub[,,8], gates_true_sd[1:rep,groups], sd_response_b[1:rep])$interval_gates[high])
  

}

save(low_matrix_rmse, low_matrix_coverage, low_matrix_bias, low_matrix_interval, high_matrix_rmse, high_matrix_coverage, high_matrix_bias, high_matrix_interval, file="gate_final_ouput.RData")