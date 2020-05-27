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


#Modifications: The following code provides additionally a DML BART estimate and is tailored to the interactive moment equation. Additionally ML input parameters are enabled.

source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/ML_Functions_b_tune.R")    
source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/cond_comp_int_tune.R")

DoubleML_int_tune <- function(data, y, d, xx, xL, methods, DML, nfold, arguments_general, tune_p, silent=FALSE, trim){

  K         <- nfold
  TE        <- matrix(0,1,length(methods))
  STE       <- matrix(0,1,length(methods))
  result    <- matrix(0,2,length(methods))
  result2   <- matrix(0,2,length(methods))
  MSE1      <- matrix(0,length(methods),K)
  MSE2      <- matrix(0,length(methods),K)
  MSE3      <- matrix(0,length(methods),K)
  cond.comp <- matrix(list(),length(methods),K)
  
  dpool     <- vector("list", length(methods))
  ypool     <- vector("list", length(methods))
  zpool     <- vector("list", length(methods))
  z1pool    <- vector("list", length(methods))
  z0pool    <- vector("list", length(methods))
  
  binary    <- as.numeric(checkBinary(data[,d]))
  
  if(!(binary==1)){
    print("treatment is not binary")
    stop()
  } 
  
  
  flag      <- 0
  split     <- runif(nrow(data))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
  
  for(k in 1:length(methods)){   
    
    if(silent==FALSE){
      cat(methods[k],'\n')
    }
    
    if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==methods[k])){
      x=xL
    } else {
      x=xx
    }
    
    for(j in 1:K){   
      
      if(silent==FALSE){
        cat('  fold',j,'\n')
      }
      
      ii  <- cvgroup == j
      nii <- cvgroup != j
      
      if(K==1){
        
        ii  <- cvgroup == j
        nii <- cvgroup == j
        
      }
      
      datause <- as.data.frame(data[nii,])
      dataout <- as.data.frame(data[ii,]) 
      
      if(length(methods)>0){
        
        cond.comp[[k,j]] <- cond_comp_int_tune(datause=datause, dataout=dataout, y=y, d=d, x=x, method=methods[k], xL=xL, arguments_general=arguments_general, tune=tune_p)
        
        MSE1[k,j]               <- cond.comp[[k,j]]$err.yz0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yz1
        MSE3[k,j]               <- cond.comp[[k,j]]$err.z 
        
        drop                   <- which(cond.comp[[k,j]]$mz_x>trim[1] & cond.comp[[k,j]]$mz_x<trim[2])      
        mz_x                   <- cond.comp[[k,j]]$mz_x[drop]
        my_z1x                 <- cond.comp[[k,j]]$my_z1x[drop]
        my_z0x                 <- cond.comp[[k,j]]$my_z0x[drop]
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]
        
        TE[1,k]                <- ATE(yout, dout, my_z1x, my_z0x, mz_x)/K + TE[1,k];
        STE[1,k]               <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,k];
        
        ypool[[k]]             <- c(ypool[[k]], yout)
        dpool[[k]]             <- c(dpool[[k]], dout)
        zpool[[k]]             <- c(zpool[[k]], mz_x)
        z1pool[[k]]            <- c(z1pool[[k]], my_z1x)
        z0pool[[k]]            <- c(z0pool[[k]], my_z0x)
        
      }
    }
  }
  
  rownames(MSE1) <- c(methods)
  rownames(MSE2) <- c(methods)
  rownames(MSE3) <- c(methods)
  
  MSE <- list("MSE_Y_X_D0"=MSE1, "MSE_Y_X_D1"=MSE2, "MSE_D_X"=MSE3)
  
  TE_pool        <- matrix(0,1,length(methods))
  STE_pool       <- matrix(0,1,length(methods))

  for(k in 1:length(methods)){ 
  TE_pool[1,(k)]         <- ATE(ypool[[k]], dpool[[k]], z1pool[[k]], z0pool[[k]], zpool[[k]])
  STE_pool[1,(k)]        <- ((SE.ATE(ypool[[k]], dpool[[k]], z1pool[[k]], z0pool[[k]], zpool[[k]]))^2)
  }

  if(length(methods)==1){
    
    TE_pool[1,(length(methods)+1)]  <- TE_pool[1,(length(methods))]
    STE_pool[1,(length(methods)+1)] <- STE_pool[1,(length(methods))]
    
  }
  colnames(result)   <- c(methods) 
  colnames(result2)  <- c(methods) 
  rownames(MSE1)     <- c(methods) 
  rownames(MSE2)     <- c(methods) 
  rownames(MSE3)     <- c(methods) 
  rownames(result)   <- c("ATE", "se")
  rownames(result2)  <- c("ATE", "se")
  
  if(DML=="DML1"){
    result[1,]         <- colMeans(TE)
    result[2,]         <- sqrt((STE))
  }
  
  if(DML=="DML2"){
    result[1,]         <- colMeans(TE_pool)
    result[2,]         <- sqrt((STE_pool))
  }
  

  table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3))   
  rownames(table)[3:5]   <- c("MSE[Y|X, D=0]", "MSE[Y|X, D=1]", "MSE[D|X]")
  
  table <- data.frame(t(table[1,]), t(table[2,]), t(table[3,]), t(table[4,]))
  return(list("table"=table, "MSE"=MSE))
  
}  

