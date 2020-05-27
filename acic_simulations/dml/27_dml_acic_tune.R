# References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     
#             "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016 
# Copyright 2017 <Victor Chernozhukov, Denis Chetverikov, Mert 
# Demirer, Esther Duflo, Whitney Newey, James Robins>
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#Modifications: The following code provides additionally a DML BART estimate and is tailored to the interactive moment equation. Additionally tuning of ML input parameters is enabled.

# Loading Packages --------------------------------------------------------
rm(list = ls()) 
library(foreign);
library(quantreg);
library(mnormt);
library(gbm);
library(glmnet);
library(MASS);
library(rpart);
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(sandwich);
library(hdm);
library(randomForest);
library(nnet)
library(neuralnet)
library(matrixStats)
library(quadprog)
library(xtable)
library(ivmodel)
library(Hmisc)
library(dummy)
library(BART) 

# Load functions, data and specify input parameters -----------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/dml"
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 
tune_parameters <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/tune/"

k <- 27 #specify knob 

setwd(dir_path)
load(paste(data_path,"acic_",k,"_data.RData", sep = "", collapse = NULL)) #load simulation replications for specified knob
source(paste(functions_path,"Moment_Functions_int_tune.R", sep = "", collapse = NULL))
source(paste(functions_path,"ML_Functions_b_tune.R", sep = "", collapse = NULL)) 
load(paste(tune_parameters,k,"_tune_p.RData", sep = "", collapse = NULL)) #load tuned parameters per replication

options(warn=-1)
set.seed(1211);
cl <- makeCluster(14) 
registerDoParallel(cl)

#fixed arguments that are not tuned
Bart <- list()
Boosting     <- list(bag.fraction = .5, train.fraction = 1.0, n.cores=1, cv.folds=2, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=500, na.action=na.omit, replace=TRUE)
Lasso       <- list(s = "lambda.min",intercept = TRUE)
RLasso       <- list(penalty = list(homoscedastic = FALSE, X.dependent.lambda =FALSE, lambda.start = NULL, c = 1.1), intercept = TRUE)
Nnet         <- list(trace = FALSE, MaxNWts=10000, maxit=1000)
Trees        <- list(reg_method="anova", clas_method="class")

arguments_general    <- list(Bart=Bart,Boosting=Boosting, Forest=Forest, RLasso=RLasso, Lasso=Lasso, Nnet=Nnet, Trees=Trees)
methods      <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet", "Bart") #ML methods used
split        <- 10
nfold 	     <- 2
rep <- 100 #simulation replications
trim <- c(0.01,0.99) #trim based on propensity score

start_time <- Sys.time()

data <- datasets[[1]] #dataset for initialization of the names
y  <- "y" #outcome variable
d  <- "z" #treatment variable

x <- paste(colnames(data[,!colnames(data)%in% c(y,d)]), collapse=" + ")  #controls 
xl <- paste("poly(",colnames(data[,!colnames(data)%in% c(y,d)]),",3,raw=TRUE)",  collapse=" + ", sep = "") #controls used for Lasso based methods 

#DML across simulation replications ----------------------------------------------
start_time <- Sys.time()

r <- foreach(i=1:rep, .combine='cbind', .inorder=TRUE) %:%
  foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','neuralnet','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet', 'BART')) %dopar% { 
    
    dml <- DoubleML_int_tune(data=datasets[[i]], y=y, d=d, xx=x, xL=xl, methods=methods, DML="DML2", nfold=nfold, arguments_general=arguments_general, tune_p=tune_p[i,], silent=FALSE, trim=trim) 
    
  }

end_time <- Sys.time()
time <- end_time - start_time

data_info <- list("n"=nrow(datasets[[1]]), "p_col" = ncol(datasets[[1]]), "y"= y, "d"=d, "x"=x, "xl"=xl, "methods"=methods, "time"= time, "split"=split, "nfold"=nfold, "arguments" = arguments_general)

#saving output data across all replications and splits
save(r, data_info, time, file = paste(k, "_dml_acic_tune.RData",sep="", collapse = NULL))
