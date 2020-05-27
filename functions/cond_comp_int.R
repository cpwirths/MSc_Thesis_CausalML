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


#Modifications: The following code provides additionally a DML BART estimate and is tailored to the interactive moment equations. 
  
cond_comp_int <- function(datause, dataout, y, d, x, method,xL, arguments){
  
  form_y   <- y
  form_d   <- d
  form_x   <- x
  form_xL  <- xL
  ind_u    <- which(datause[,d]==1)
  ind_o    <- which(dataout[,d]==1)
  err.yz1  <- NULL
  err.yz0  <- NULL
  my_z1x   <- NULL
  my_z0x   <- NULL
  fit.yz1  <- NULL
  fit.yz0  <- NULL
  fit.z    <- NULL
  mis.z    <- NULL
  
  ########################## Boosted  Trees ###################################################;
  
  if(method=="Bart"){
    
    option <- arguments[[method]]
    arg    <- option
    binary_y    <- as.numeric(checkBinary(datause[,form_y]))
      
    if(binary_y==0){
      fit.yz1        <- bartF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, type="cont",  option=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err 
      my_z1x         <- predict(fit.yz1$model, dataout[,colnames(dataout)%in%colnames(fit.yz1$mode$varcount)]) 
      my_z1x <- colMeans(my_z1x)
      
      fit.yz0        <- bartF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x=form_x,  form_y=form_y, type="cont",  option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout[,colnames(dataout)%in%colnames(fit.yz0$mode$varcount)]) 
      my_z0x <- colMeans(my_z0x)
    }else{
      fit.yz1        <- bartF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, type="binary",  option=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout[,colnames(dataout)%in%colnames(fit.yz1$mode$varcount)]) 
      my_z1x <- my_z1x$prob.test.mean
      
      fit.yz0        <- bartF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x=form_x,  form_y=form_y, type="binary",  option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout[,colnames(dataout)%in%colnames(fit.yz0$mode$varcount)]) 
      my_z0x <- my_z0x$prob.test.mean
    }
  
      fit.z          <- bartF(datause=datause, dataout=dataout, form_x=form_x,  form_y=form_d, type="binary",  option=arg)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
  }
  
  if(method=="Boosting"){
    
    option <- arguments[[method]]
    arg    <- option
    arg[which(names(arg) %in% c("clas_dist","reg_dist"))] <-  NULL
    
      fit.yz1        <- boost(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, distribution=option[['reg_dist']], option=arg)
      
          
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, n.trees=fit.yz1$best, dataout, type="response") 
      
      fit.yz0        <- boost(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, n.trees=fit.yz0$best, dataout, type="response") 
      
      fit.z          <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, distribution=option[['clas_dist']], option=arg)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
  }  
  
  ########################## Neural Network(Nnet Package) ###################################################;   
  
  if(method=="Nnet"){
    
    option <- arguments[[method]]
    arg    <- option
    
      fit.yz1        <- nnetF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz1$f] <- as.data.frame(scale(dataouts[,!fit.yz1$f], center = fit.yz1$min, scale = fit.yz1$max - fit.yz1$min))
      my_z1x         <- predict(fit.yz1$model, dataouts)*(fit.yz1$max[fit.yz1$k]-fit.yz1$min[fit.yz1$k])+fit.yz1$min[fit.yz1$k] 
      fit.yz0        <- nnetF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz0$f] <- as.data.frame(scale(dataouts[,!fit.yz0$f], center = fit.yz0$min, scale = fit.yz0$max - fit.yz0$min))
      my_z0x         <- predict(fit.yz0$model, dataouts)*(fit.yz0$max[fit.yz0$k]-fit.yz0$min[fit.yz0$k])+fit.yz0$min[fit.yz0$k] 
      
      fit.z          <- nnetF(datause=datause, dataout=dataout, form_x=form_x, form_y=form_d, clas=TRUE, arg=arg, scaling=FALSE)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
  } 
  
 
  ########################## Lasso and Post Lasso(Glmnet) Package) ###################################################;    
  
  if(method=="Ridge" || method=="Lasso" || method=="Elnet"){
    
    if(method=="Ridge"){ alp=0 }
    if(method=="Lasso"){ alp=1 }
    if(method=="Elnet"){ alp=0.5 }
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("s"))] <-  NULL
    s         <- option[['s']]
      
      fit.yz1        <- lassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z1x         <- predict(fit.yz1$model, newx=fit.p$x[,-1] ) 
      
      fit.yz0        <- lassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z0x         <- predict(fit.yz0$model,  newx=fit.p$x[,-1])   
    
      fit.z          <- lassoF(datause=datause, dataout=dataout,  form_x, form_d, logit=TRUE, alp=alp, arg=arg, s=s)
      mis.z          <- error(fit.z$yhatout, dataout[,d])$mis
  }    
  
  ############# Random Forest ###################################################;
  
  if(method=="Forest" | method=="TForest"){
    
    tune = FALSE
    if(method=="TForest"){tune=TRUE}
    
    option    <- arguments[[method]]
    
    arg       <- option
    arg[which(names(arg) %in% c("clas_nodesize","reg_nodesize"))] <-  NULL
   
      
      fit.yz1        <- RF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout, type="response") 
      
      fit.yz0        <- RF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout, type="response")
      
  
      fit.z          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, nodesize=option[["clas_nodesize"]], arg=arg, reg=TRUE, tune=tune)
      mis.z          <- error(as.numeric(fit.z$yhatout), dataout[,y])$mis
    
  }
  
  ########################## Regression Trees ###################################################;     
  
  if(method=="Trees"){
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("reg_method","clas_method"))] <-  NULL
    
      
      fit.yz1        <- tree(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout) 
      
      fit.yz0        <- tree(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model,dataout)   
    
      
      fit.z          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d, method=option[["clas_method"]], arg=arg)
      mis.z          <- error(as.numeric(fit.z$yhatout), dataout[,y])$mis
    
  }
  
  err.z          <- error(fit.z$yhatout, dataout[,d])$err
  mz_x           <- fit.z$yhatout       
  rz             <- fit.z$resout 
  err.z          <- error(fit.z$yhatout, dataout[,d])$err 
  
  return(list(method = method, my_z1x=my_z1x, my_z0x=my_z0x, mz_x= mz_x, err.yz1=err.yz1, err.yz0=err.yz0,  err.z = err.z, mis.z=mis.z, rz=rz, fit.z= fit.z,  fit.yz1out= fit.yz1$yhatout,  fit.yz0out= fit.yz0$yhatout));
}  

