#Tunes the ML inputs parameter over every grid for one knob over all simulation replications

# Loading Packages --------------------------------------------------------
rm(list=ls())
vec.pac= c("foreign", "quantreg", "mnormt", "gbm", "glmnet","MASS", "rpart", "doParallel", "sandwich", "randomForest","nnet", "matrixStats", "xtable", "readstata13", "car", "lfe", "doParallel","caret", "foreach","iterators", "parallel", "doParallel",  "multcomp","cowplot", "MASS", "sandwich", "hdm", "neuralnet", "quadprog", "Hmisc")
lapply(vec.pac, require, character.only = TRUE)

# Load functions, data and specify input parameters -----------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/tune"
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 

k <- 27 #specify knob

setwd(dir_path)
load(paste(data_path,"acic_",k,"_data.RData", sep = "", collapse = NULL)) #load simulation replications for specified knob
rep <- 100 #simulation replications

y <- "y"  #outcome variable
d <- "z"  #treatment variable

data <- datasets[[1]] #for initialization 

x <- paste(colnames(data[,!colnames(data)%in% c(y,d)]), collapse=" + ") #controls 
xl <-  paste("poly(",colnames(data[,!colnames(data)%in% c(y,d)]),",3,raw=TRUE)",  collapse=" + ", sep = "") #controls used for Lasso based methods


methods      <- c("glmnet", "rpart", "gbm", "rf", "nnet")   
method_names <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet")

options(warn=-1)
set.seed(1211);
cl <- makeCluster(14) 
registerDoParallel(cl)


# Fixed Parameters and input for tuning grid ------------------------------

args         <- list(glmnet=list(lambda.min=0.02), gbm=list(verbose=FALSE, bag.fraction = .5, train.fraction = 1.0, distribution="gaussian"), rf=list(nodesize=5, ntree=500, replace=TRUE), nnet=list(linout = TRUE, trace = FALSE, MaxNWts=10000, maxit=1000))

methodML   <- c("repeatedcv", "repeatedcv", "repeatedcv", "repeatedcv", "repeatedcv") 
tune       <- c(10, 10, NA, NA, NA) 
proces     <- c(NA, NA, NA, NA, NA)                  # pre-processing method
proces <- list("glmnet"= NULL, "rpart"=NULL, "gbm"=NULL, "rf"=NULL, "nnet"="range")
select     <- c("best", "best", "best", "best", "best")   
rep        <- c(2, 2, 2, 2, 2)
cv         <- c(2, 2, 2, 2, 2)

tune_param       <- list(0)
tune_param[[1]]  <- 0
tune_param[[2]]  <- 0
tune_param[[3]]  <- expand.grid(interaction.depth = c(2,3,4),
				                                n.trees = c(600,1000),
								                                shrinkage = c(0.1),
								                                n.minobsinnode = c(1,5))
tune_param[[4]]  <- data.frame(mtry=c(round((ncol(data)-2)/3)-5,round((ncol(data)-2)/3),round((ncol(data)-2)/3)+5))
tune_param[[5]]  <- expand.grid(size=c(2,4,8),decay=c(0.01,0.02))

start_time <- Sys.time()

#Tune ML parameters across simulation replications ----------------------------------------------
tune_p <- foreach(j = 1:rep, .combine='rbind', .inorder=T, .packages=vec.pac) %dopar%{ 
  set.seed(j)
  tune_param_yz1 <- list()
  tune_param_yz0 <- list()
  tune_param_d <- list()
  
  data <- datasets[[j]]
  for(l in 1:length(methods)){
   print(methods[l])
  
  if(tune_param[[l]]==0){f=NULL}else{f = tune_param[[l]]}
  
  if(methods[l]=="glmnet"){xx <- xl}else{xx<-x}
  form_y           <- as.formula(paste(y,"~",xx,sep=""));
  form_d           <- as.formula(paste(d,"~",xx,sep=""));
  
  fitControl   <- trainControl(method = methodML[l], number = cv[l], repeats = rep[l], allowParallel = FALSE, verboseIter=FALSE, search="grid", selectionFunction=select[l])
  arg          <- c(list(form=form_y, data = data[ind_u,],  method = methods[l],  tuneGrid = f, trControl = fitControl, preProcess=proces[[l]], tuneLength=tune[l]), args[[methods[l]]])
  fit.yz1      <- suppressWarnings(do.call(caret::train, arg))
  tune_param_yz1[[l]]  <- fit.yz1$bestTune

  fitControl   <- trainControl(method = methodML[l], number = cv[l], repeats = rep[l], allowParallel = FALSE, verboseIter=FALSE, search="grid", selectionFunction=select[l])
  arg          <- c(list(form=form_y, data = data[-ind_u,],  method = methods[l], tuneGrid = f, trControl = fitControl, preProcess=proces[[l]], tuneLength=tune[l]), args[[methods[l]]])
  fit.yz0      <- suppressWarnings(do.call(caret::train, arg))
  tune_param_yz0[[l]]  <- fit.yz0$bestTune
  
  fitControl   <- trainControl(method = methodML[l], number = cv[l], repeats = rep[l], allowParallel = FALSE, verboseIter=FALSE, search="grid", selectionFunction=select[l])
  arg          <- c(list(form=form_d, data = data,  method = methods[l], tuneGrid = f, trControl = fitControl, preProcess=proces[[l]], tuneLength=tune[l]), args[[methods[l]]])
  fit.d      <- suppressWarnings(do.call(caret::train, arg))
  tune_param_d[[l]]  <- fit.d$bestTune

  }
  names(tune_param_yz1) <- method_names
  names(tune_param_yz0) <- method_names
  names(tune_param_d)   <- method_names
  tune_p <- list("tune_param_yz1"=tune_param_yz1, "tune_param_yz0"=tune_param_yz0, "tune_param_d"=tune_param_d)
}

end_time <- Sys.time()
time <- end_time - start_time

#saving output data across all replications and splits
save(tune_p, time, file=paste(k, "_tune_p.RData",sep="", collapse = NULL))
