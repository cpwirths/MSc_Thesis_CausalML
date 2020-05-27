rm(list=ls())
# Loading Packages --------------------------------------------------------
library(caret)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(grf)
library(glmnet)
library(nnet)
library(rpart)
library(gbm)
library(MLmetrics)
library(BART)
library(bartCause)

# Load functions, data and specify input parameters -----------------------

dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/dr_mom_ites/"
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 
ML_parameters <- "~/Documents/GitHub/MSc_Thesis_CausalML/globals_parameters.R"

k <- 27 #specify knob

setwd(dir_path)
load(paste(data_path,"acic_",k,"_data.RData", sep = "", collapse = NULL)) #load simulation replications for specified knob
source(paste(functions_path,"ites_all_arguments_final.R", sep = "", collapse = NULL))
source(ML_parameters)

y = "y"  #outcome variable
d = "z"  #treatment variable
trim <- c(0.01, 0.99) #trim based on propensity score
rep <- 100 #simulation replications

options(warn=-1)
set.seed(1211);
cl <- makeCluster(14) 
registerDoParallel(cl)

data <- datasets[[1]] #take one dataset for intitialization
x <- paste(colnames(data[,!colnames(data)%in% c(y,d)]), collapse=" + ") #controls 
xl <- paste("poly(",colnames(data[,!colnames(data)%in% c(y,d)]),",3,raw=TRUE)",  collapse=" + ", sep = "") #controls used for Lasso based methods

#MOM DR across simulation replications ----------------------------------------------
start_time <- Sys.time()
r <- foreach(i = 1:rep, .combine = 'rbind', .inorder=TRUE,.packages = c("grf", "gbm", "glmnet", "nnet","rpart", "glmnet", "MLmetrics", "BART", "bartCause")) %dopar% {
  iates <- iates_all(data=datasets[[i]],y=y,d=d, xx=x, xl=xl, arguments=arguments, trim=trim)
}

end_time <- Sys.time()
time <- end_time - start_time

#saving output data across all replications
save(r,time, file = paste(k,"_dr_mom_ites_acic.RData", sep = "", collapse = NULL))
