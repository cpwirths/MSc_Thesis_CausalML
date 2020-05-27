#this script create the causal forest estimates for the acic competition per datasets
#load dataset of interest
# Loading Packages --------------------------------------------------------
rm(list=ls())
library(grf)
library(iterators)
library(parallel)
library(doParallel)

# Load functions, data and specify input parameters -----------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/cf" 
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 

k <- 27 #specify knob

setwd(dir_path)
load(paste(data_path,"acic_",k,"_data.RData", sep = "", collapse = NULL)) #load simulation replications for specified knob
source(paste(functions_path,"cf_final_sim.R", sep = "", collapse = NULL))
  
options(warn = -1)
set.seed(1211)
cl <- makeCluster(14) 
registerDoParallel(cl)

start_time <- Sys.time()

y  <- "y" #outcome variable
d  <- "z" #treatment variable
rep <- 100 #simulation replications

#Causal Forest across simulation replications ----------------------------------------------
start_time <- Sys.time()

r <- foreach(i = 1:rep, .combine = 'rbind', .inorder=TRUE,.packages = c("grf")) %dopar% {
 cf_output <- cf(data=datasets[[i]], y=y, d=d)
}

end_time <- Sys.time()
time <- end_time - start_time 

#saving output data across all replications
save(r, k, file = paste(k,"_cf_acic_results.RData", sep = "", collapse = NULL))


