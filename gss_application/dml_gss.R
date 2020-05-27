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


# Modifications: The following code provides additionally a DML BART estimate and is tailored to the interactive moment equation. 

# Returns the DML output for the GSS application
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

dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/gss_application/"
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 
ML_parameters <- "~/Documents/GitHub/MSc_Thesis_CausalML/globals_parameters.R"

setwd(dir_path)
load(paste(dir_path,"GSS_input.RData", sep = "", collapse = NULL))
source(paste(functions_path,"Moment_Functions_int.R", sep = "", collapse = NULL))
source(paste(functions_path,"ML_Functions_b.R", sep = "", collapse = NULL))
source(paste(functions_path,"output_dml.R", sep = "", collapse = NULL))
source(ML_parameters)

data <- model.matrix(~.-1,data)
data <- data.frame(data)
data <- na.omit(data)
ols <- lm(response~., data )

options(warn=-1)
set.seed(1211);
cl <- makeCluster(14)
registerDoParallel(cl)

methods      <- c("Bart","Lasso", "Trees", "Boosting", "Forest", "Nnet")  
split        <- 100
nfold 	     <- 2

start_time <- Sys.time()

y  <- "response" # Outcome Variable
d  <- "treatment" # Outcome Variable

x <- paste(colnames(data[,!colnames(data)%in% c(y,d)]), collapse=" + ")  #controls 
xl <- paste("(", x, ")^2", sep="") #controls used for Lasso based methods

############## Arguments for DoubleML function:

start_time <- Sys.time()

r <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','neuralnet','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet','BART')) %dopar% { 
	  
  dml <- DoubleML_int(data=data, y=y, d=d, xx=x, xL=xl, methods=methods, DML="DML2", nfold=nfold, arguments=arguments, silent=FALSE, trim=c(0.01,0.99)) 
  
}

end_time <- Sys.time()
time <- end_time - start_time

data_info <- list("n"=nrow(data), "p_col" = ncol(data), "y"= y, "d"=d, "x"=x, "xl"=xl, "methods"=methods, "time"= time, "split"=split, "nfold"=nfold, "ols" = ols$coefficients, "arguments" = arguments)

#calulculate DML output by taking the median across splits
output = output_p(r[,"table"])
save(r, output, data_info, file = "dml_gss_final.RData")










