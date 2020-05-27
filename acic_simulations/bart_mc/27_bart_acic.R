#this script creates the ATE and GATE bart estimates for the ACiC competition per knob

# Loading Packages --------------------------------------------------------
rm(list=ls())
library(bartCause)
library(iterators)
library(parallel)
library(doParallel)

# Load functions, data and specify input parameters -----------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/bart_mc" 
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 
  
k <- 27  #specify knob

setwd(dir_path)
load(paste(data_path,"acic_",k,"_data.RData", sep = "", collapse = NULL)) #load simulation replications for specified knob
source(paste(functions_path,"bart_estimate_sim.R", sep = "", collapse = NULL))

options(warn = -1)
set.seed(1211)
cl <- makeCluster(14) 
registerDoParallel(cl)
data <- datasets[[1]]
rep <- 100 #simulation replications

y  <- "y" #outcome variable
d  <- "z" #treatment variable

#BARTMC across simulation replications ----------------------------------------------
start_time <- Sys.time()

r <- foreach(i = 1:rep, .combine = 'rbind', .inorder=TRUE,.packages = c("bartCause")) %dopar% {
  b <- bart_estimates(data=datasets[[i]], y=y,d=d)
}

end_time <- Sys.time()
time <- end_time - start_time

#saving output data across all replications
save(r, time, file = paste(k,"_bart_acic_results.RData", sep = "", collapse = NULL))

