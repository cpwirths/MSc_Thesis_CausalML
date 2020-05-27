#References: Knaus, Lechner, Strittmatter (2018). Machine Learning Estimation of Heterogeneous Causal Effects: Empirical Monte Carlo Evidence
#Copyright <2018> <Knaus, Lechner, Strittmatter> R-packages: Package: CATEs, MIT license

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Modifications: Addition of BART MC, Tree MOM DR, Boosting MOM DR, Neural Net MOM DR, BART MOM DR, DML 2 estimation procedure and propensity score trimming


# BART functions ----------------------------------------------------------

nuisance_cf_bart <- function(y,d,x,index, args_method,
                               args_p=NULL,
                               args_y=NULL,
                               args_y0=NULL,
                               args_y1=NULL) {
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {
    # P-score
    fit_p = do.call(pbart,c(list(x.train=x[-index[[i]],,drop=F],
                                     y.train=d[-index[[i]]]), args_method,
                                args_p))

    
    predict_p <- predict(fit_p, x[index[[i]],colnames(x)%in%colnames(fit_p$varcount),drop=F])
    predict_p <- predict_p$prob.test.mean
    
    np[index[[i]],1] = predict_p
    
    # Outcome
    fit_y = do.call(gbart,c(list(x.train=x[-index[[i]],,drop=F],
                                     y.train=y[-index[[i]]]), args_method,
                                args_y))
    predict_y <- predict(fit_y, x[index[[i]],colnames(x)%in%colnames(fit_y$varcount),drop=F])
    predict_y <- colMeans(predict_y)
    np[index[[i]],2] = predict_y 
    
    # Outcome of non-treated
    fit_y0 = do.call(gbart,c(list(x.train=x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F],
                                      y.train=y[-index[[i]]][d[-index[[i]]] == 0]), args_method,
                                 args_y0))
    predict_y0 <- predict(fit_y0, x[index[[i]],colnames(x)%in%colnames(fit_y0$varcount),drop=F])
    predict_y0 <- colMeans(predict_y0) 
    np[index[[i]],3] = predict_y0
    
    # Outcome of non-treated
    fit_y1 = do.call(gbart,c(list(x.train=x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F],
                                  y.train=y[-index[[i]]][d[-index[[i]]] == 1]), args_method,
                             args_y1))
    predict_y1 <- predict(fit_y1, x[index[[i]],colnames(x)%in%colnames(fit_y1$varcount),drop=F])
    predict_y1 <- colMeans(predict_y1) 
    np[index[[i]],4] = predict_y1
    
  }
  return(np)
}

mom_dr_bart = function(y,d,x,np,xnew,form_x=NULL, args_tau=NULL, trim) {
  
  
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  fit_tau = do.call(gbart,c(list(x.train=x[drop,],y.train=mo),args_tau))
  iate = predict(fit_tau,xnew[,colnames(xnew)%in%colnames(fit_tau$varcount)])
  iate <- colMeans(iate)
  return(iate)
}


# Lasso fucntions ---------------------------------------------------------

# Estimation of nusaince parameters

nuisance_cf_glmnet <- function(y,d,x,index, args_method, s,
                               args_p=NULL,
                               args_y=NULL,
                               args_y0=NULL,
                               args_y1=NULL) {
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {
    # P-score
    fit_p = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F],
                                     y=d[-index[[i]]],
                                     family="binomial"), args_method,
                                args_p))
    np[index[[i]],1] = predict(fit_p,x[index[[i]],,drop=F], s = s, type = "response")
    
    # Outcome
    fit_y = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F],
                                     y=y[-index[[i]]]), args_method,
                                args_y))
    np[index[[i]],2] = predict(fit_y,x[index[[i]],,drop=F], s = s)
    
    # Outcome of non-treated
    fit_y0 = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F],
                                      y=y[-index[[i]]][d[-index[[i]]] == 0]), args_method,
                                 args_y0))
    np[index[[i]],3] = predict(fit_y0,x[index[[i]],,drop=F], s = s)
    
    # Outcome of non-treated
    fit_y1 = do.call(cv.glmnet,c(list(x=x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F],
                                      y=y[-index[[i]]][d[-index[[i]]] == 1]), args_method,
                                 args_y1))
    np[index[[i]],4] = predict(fit_y1,x[index[[i]],,drop=F], s = s)
  }
  return(np)
}

# DR MOM Outcome transformation Lasso 

mom_dr_glmnet = function(y,d,x,np,xnew,form_x=NULL, args_tau=NULL, trim) {
  
  s         <- args_tau[['s']]
  args_tau[which(names(args_tau) %in% c("s"))] <-  NULL
  
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  fit_tau = do.call(cv.glmnet,c(list(x=x[drop,],y=mo),args_tau))
  iate = predict(fit_tau,xnew, s = s)
  return(iate)
}


# Tree functions ----------------------------------------------------------

# Estimation of nusaince parameters
nuisance_cf_trees <- function(dep_form,treat_form,df, d,index, args_method,
                              args_p=NULL,
                              args_y=NULL,
                              args_y0=NULL,
                              args_y1=NULL) {
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {

    fit_p <- do.call(rpart,c(list(formula=treat_form, data=df[-index[[i]],,drop=F]), method=args_method[["clas_method"]], args_p))
    
    bestcp         <- fit_p$cptable[which.min(fit_p$cptable[,"xerror"]),"CP"]
    pfit_p         <- prune(fit_p,cp=bestcp)
    np[index[[i]],1] <- predict(pfit_p, newdata=df[index[[i]],,drop=F])[,2]
    
    
    fit_y <- do.call(rpart,c(list(formula=dep_form, data=df[-index[[i]],,drop=F]), method=args_method[["reg_method"]], args_y))
    
    bestcp         <- fit_y$cptable[which.min(fit_y$cptable[,"xerror"]),"CP"]
    pfit_y         <- prune(fit_y,cp=bestcp)
    np[index[[i]],2] <- predict(pfit_y, newdata=df[index[[i]],,drop=F])
    
    
    fit_y0 <- do.call(rpart,c(list(formula=dep_form, data=df[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F]), method=args_method[["reg_method"]], args_y))
    
    bestcp         <- fit_y0$cptable[which.min(fit_y0$cptable[,"xerror"]),"CP"]
    pfit_y0         <- prune(fit_y0,cp=bestcp)
    np[index[[i]],3] <- predict(pfit_y0, newdata=df[index[[i]],,drop=F])
    
    fit_y1 <- do.call(rpart,c(list(formula=dep_form, data=df[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F]), method=args_method[["reg_method"]], args_y))
    
    bestcp         <- fit_y1$cptable[which.min(fit_y1$cptable[,"xerror"]),"CP"]
    pfit_y1         <- prune(fit_y1,cp=bestcp)
    np[index[[i]],4] <- predict(pfit_y1, newdata=df[index[[i]],,drop=F])
    
    
  }
  
  return(np)
}

mom_dr_trees = function(y,d,x,np,xnew,form_x,args_tau=NULL, trim) {
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  df_train <- as.data.frame(cbind(x[drop,], mo))
  df_new <- as.data.frame(xnew)
  form         <- as.formula(paste("mo ~", form_x))
  
  fit_tau <- do.call(rpart,c(list(formula=form, data=df_train), method="anova"))
  bestcp         <- fit_tau$cptable[which.min(fit_tau$cptable[,"xerror"]),"CP"]
  pfit_tau         <- prune(fit_tau,cp=bestcp)
  
  iate <- predict(pfit_tau, newdata=df_new)
  return(iate)
}


# Boosting functions ------------------------------------------------------


nuisance_cf_boost <- function(dep_form,treat_form,df, d,index, args_method, 
                              args_p=NULL,
                              args_y=NULL,
                              args_y0=NULL,
                              args_y1=NULL) {
  

  args    <- args_method
  args[which(names(args) %in% c("clas_dist","reg_dist"))] <-  NULL
  
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {
    fit_p <- do.call(gbm, c(list(formula=treat_form, data=df[-index[[i]],,drop=F]), args, distribution=args_method[['clas_dist']]))
    
    if( args[['cv.folds']]>0)  {best <- gbm.perf(fit_p,plot.it=FALSE,method="cv")}else{best <- gbm.perf(fit_p,plot.it=FALSE,method="OOB")}
    
    np[index[[i]],1] <- predict(fit_p, n.trees=best, newdata=df[index[[i]],,drop=F],  type="response")
    
    
    fit_y <- do.call(gbm, c(list(formula=dep_form, data=df[-index[[i]],,drop=F]), args, distribution=args_method[['reg_dist']], args_p))
    
    if( args[['cv.folds']]>0)  {best           <- gbm.perf(fit_y,plot.it=FALSE,method="cv")}else{best <- gbm.perf(fit_y,plot.it=FALSE,method="OOB")}
    
    np[index[[i]],2] <- predict(fit_y, n.trees=best, newdata=df[index[[i]],,drop=F],  type="response")
    
    
    fit_y0 <- do.call(gbm, c(list(formula=dep_form, data=df[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F]), args, distribution=args_method[['reg_dist']], args_p))
    
    if( args[['cv.folds']]>0)  {best<- gbm.perf(fit_y0,plot.it=FALSE,method="cv")}else{best<- gbm.perf(fit_y0,plot.it=FALSE,method="OOB")}
    
    np[index[[i]],3] <- predict(fit_y0, n.trees=best, newdata=df[index[[i]],,drop=F],  type="response")
    
    
    fit_y1 <- do.call(gbm, c(list(formula=dep_form, data=df[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F]), args, distribution=args_method[['reg_dist']], args_p))
    
    if( args[['cv.folds']]>0)  {best<- gbm.perf(fit_y1,plot.it=FALSE,method="cv")}else{best <- gbm.perf(fit_y1,plot.it=FALSE,method="OOB")}
    
    np[index[[i]],4] <- predict(fit_y1, n.trees=best, newdata=df[index[[i]],,drop=F],  type="response")
    
    
  }
  
  return(np)
}

mom_dr_boost = function(y,d,x,np,xnew,form_x,args_tau=NULL, trim) {
  
  args    <- args_tau
  args[which(names(args) %in% c("clas_dist","reg_dist"))] <-  NULL
  
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  df_train <- as.data.frame(cbind(x[drop,], mo))
  df_new <- as.data.frame(xnew)
  form         <- as.formula(paste("mo ~", form_x))
  fit_tau = do.call(gbm, c(list(formula=form, data=df_train), args, distribution=args_tau[['reg_dist']]))
  
  if(args_tau[['cv.folds']]>0)  {best           <- gbm.perf(fit_tau,plot.it=FALSE,method="cv")}else{best <- gbm.perf(fit_tau,plot.it=FALSE,method="OOB")}
  iate <- predict(fit_tau, n.trees=best, newdata=df_new,  type="response")
  return(iate)
}

# Random Forest functions -------------------------------------------------



nuisance_cf_grf <- function(y,d,x,index,
                            args_p=NULL,
                            args_y=NULL,
                            args_y0=NULL,
                            args_y1=NULL) {
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {
	  fit_p = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F],
                                             Y=d[-index[[i]]]),
                                        tune.parameters = "all", args_p))
  
    np[index[[i]],1] = predict(fit_p,x[index[[i]],,drop=F])$prediction
    
    fit_y = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F],
                                             Y=y[-index[[i]]]),
                                        tune.parameters = "all",
                                        args_y))
    np[index[[i]],2] = predict(fit_y,x[index[[i]],,drop=F])$prediction
    
    fit_y0 = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F],
                                              Y=y[-index[[i]]][d[-index[[i]]] == 0]),
                                         tune.parameters = "all",
                                         args_y0))
    np[index[[i]],3] = predict(fit_y0,x[index[[i]],,drop=F])$prediction
    
    fit_y1 = do.call(regression_forest,c(list(X=x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F],
                                              Y=y[-index[[i]]][d[-index[[i]]] == 1]),
                                         tune.parameters = "all",
                                         args_y1))
    np[index[[i]],4] = predict(fit_y1,x[index[[i]],,drop=F])$prediction
  }
  
  return(np)
}

mom_dr_grf = function(y,d,x,np,xnew, form_x = NULL, args_tau=NULL, trim) {
  
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  fit_tau = do.call(regression_forest,c(list(X=x[drop,],Y=mo),tune.parameters = "all",args_tau))
  iate = predict(fit_tau,xnew)$prediction

  return(iate)
}


# Neural Net functions ----------------------------------------------------


nuisance_cf_nnet <- function(y, d,x ,index, args_method,
                             args_p=NULL,
                             args_y=NULL,
                             args_y0=NULL,
                             args_y1=NULL){
  
  np = matrix(NA,length(d),4)
  colnames(np) = c("p_hat", "y_hat", "y0_hat", "y1_hat")
  for(i in 1:length(index)){
    
    np[index[[i]],1] <-  scaled_nnet(x_train= x[-index[[i]],,drop=F], y_train = d[-index[[i]]], x_test=x[index[[i]],,drop=F], arg_nnet=args_method, args_p, scale = F)
    
    np[index[[i]],2] <-  scaled_nnet(x_train= x[-index[[i]],,drop=F], y_train = y[-index[[i]]], x_test=x[index[[i]],,drop=F], arg_nnet=args_method, args_y)
    
    np[index[[i]],3] <-  scaled_nnet(x_train= x[-index[[i]],,drop=F][d[-index[[i]]] == 0,,drop=F], y_train = y[-index[[i]]][d[-index[[i]]] == 0], x_test=x[index[[i]],,drop=F], arg_nnet=args_method, arg_add = args_y0)
    
    np[index[[i]],4] <-  scaled_nnet(x_train= x[-index[[i]],,drop=F][d[-index[[i]]] == 1,,drop=F], y_train = y[-index[[i]]][d[-index[[i]]] == 1], x_test=x[index[[i]],,drop=F], arg_nnet=args_method, arg_add = args_y1)
    
  }
  return(np)
}


mom_dr_nnet = function(y,d,x,np,xnew,form_x=NULL, args_tau=NULL, trim) {
  drop <- which(np[,"p_hat"] > trim[1] & np[,"p_hat"] < trim[2])
  mo = np[drop,"y1_hat"] - np[drop,"y0_hat"] + d[drop] * (y[drop]-np[drop,"y1_hat"]) / np[drop,"p_hat"] - (1-d)[drop] * (y[drop]-np[drop,"y0_hat"]) / (1-np[drop,"p_hat"])
  iate = scaled_nnet(x_train=x[drop,], y_train=mo, x_test=xnew, arg_nnet=args_tau,arg_add=NULL, scale=TRUE)
  return(iate)
}

#max-min scaling of neural net
scaled_nnet <- function(x_train, y_train, x_test, arg_nnet, arg_add=NULL, scale=TRUE){
  if(scale==TRUE){
    x_train <- as.data.frame(x_train)
    x_test <- as.data.frame(x_test)
    f <- sapply(x_train,is.factor)
    f[which(lapply(sapply(x_train, unique), length)==1)] <- TRUE
    maxs <- apply(x_train[,!f], 2, max)
    mins <- apply(x_train[,!f], 2, min)
    x_train <- as.data.frame(scale(x_train[,!f], center = mins, scale = maxs - mins))
    x_test <- as.data.frame(scale(x_test[,!f], center = mins, scale = maxs - mins))
    max_dep <- max(y_train)
    min_dep <- min(y_train)
    y_train <- scale(y_train, center = min_dep, scale= max_dep-min_dep)
    fit = do.call(nnet, c(list(x=x_train, y=y_train), arg_nnet, arg_add))
    return(predict(fit,x_test)*(max_dep-min_dep)+min_dep)
  }else{
    fit = do.call(nnet , c(list(x=x_train, y=y_train), arg_nnet, arg_add))
    return(predict(fit,x_test))
  }
}

# Double robust estimation of IATE based on DML ---------------------------


cf_dml1 = function(est,y,d,x,np,xnew,index,form_x = NULL, args_tau=NULL, trim=NULL) {
  
  iate = matrix(0,length(d),1)
  
  for (i in 1:length(index)) {
    iate = iate + 1/length(index) *
      do.call(est,list(y[index[[i]]],
                       d[index[[i]]],
                       x[index[[i]],,drop=F],
                       np[index[[i]],],
                       xnew,
                       form_x, 
                       args_tau, trim))
  }
  return(iate)
}

cf_dml2 = function(est,y,d,x,np,xnew,index,form_x = NULL, args_tau=NULL, trim=NULL) {
  
  iate = matrix(0,length(d),1)
  iate = do.call(est,list(y,
                          d,
                          x,
                          np,
                          xnew,
                          form_x,
                          args_tau, trim))
  
  return(iate)
}







