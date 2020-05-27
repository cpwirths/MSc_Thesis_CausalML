#References: Knaus, Lechner, Strittmatter (2018). Machine Learning Estimation of Heterogeneous Causal Effects: Empirical Monte Carlo Evidence
#Copyright <2018> <Knaus, Lechner, Strittmatter> R-packages: Package: CATEs, MIT license

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#Modifications: Addition of BART MC, Tree MOM DR, Boosting MOM DR, Neural Net MOM DR, BART MOM DR, DML 2 estimation procedure and propensity score trimming

IATEs <- function(tr, val,  y, d, xx, xl, arguments, trim, total_index, tau_t, tau_val_mu, n,p) {
  
  d_t = tr[,d]
  y_t = tr[,y]
  y_val = val[,y]
  
  dep_form_x           <- as.formula(paste(y, "~", xx));
  treat_form_x         <- as.formula(paste(d, "~", xx))
  form_x <- xx
  form_xl           <- as.formula(paste("~", xl, "-1"));
  
  #moel matrix for all tree based methods
  x_t <- model.matrix(~.-1, tr[,!colnames(tr)%in%c(y, d)])
  x_v <- model.matrix(~.-1, val[,!colnames(val)%in%c(y, d)])
  
  methods <- c("RF MOM DR DML1","CF LC", "Lasso MOM DR DML1", "Boosting MOM DR DML1", "Trees MOM DR DML1", "Nnet MOM DR DML1","Bart MOM DR DML1", "BART MChains")
  ML_methods <- c("RF", "Lasso", "Boosting", "Trees", "Nnet","Bart")
  
  # Initialide matrix to store predicted IATEs of each method plus additional information of for evaluation pruposes
   add_info <- c("index_validation_sample", "tau_val_mu", "y_val")
  names <- c(methods, add_info)
  iates_mat = matrix(NA,n,length(names))
  colnames(iates_mat) =  names
  
  rmse <- matrix(NA, 3, length(ML_methods))
  colnames(rmse) <- ML_methods
  rownames(rmse) <- c("RMSE_Y1", "RMSE_Y0", "RMSE_D")
  
  mse <- matrix(NA, 3, length(ML_methods))
  colnames(mse) <- ML_methods
  rownames(mse) <- c("RMSE_Y1", "RMSE_Y0", "RMSE_D")
  
  #create index for sample splitting for nusiane parameters
  index = caret::createFolds(y_t, k = 2)
  
  # RF DR DML
  np = nuisance_cf_grf(y=y_t,d=d_t,x=x_t,index=index)
  rmse[, "RF"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "RF"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  iates_mat[,"RF MOM DR DML1"] <- do.call(cf_dml1,list(mom_dr_grf,y_t,d_t,x_t,np,x_v,index, form_x=NULL, args_tau=NULL, trim))
  
  # Causal Forest with local centering
  cf_lc = causal_forest(x_t, y_t, d_t, Y.hat = np[,2], W.hat = np[,1], tune.parameters= "all")
  iates_mat[,"CF LC"] = predict(cf_lc,x_v)$predictions
  
  # Estimate Lasso nuisance parameters with same index
  x_t_l <- model.matrix(form_xl,tr)
  x_v_l <- model.matrix(form_xl,val)
  
  args_glmnet <- arguments[["Lasso"]]
  args_glmnet$alpha <- 1
  s         <- args_glmnet[['s']]
  args <- args_glmnet
  args[which(names(args) %in% c("s"))] <-  NULL
  
  np = nuisance_cf_glmnet(y_t,d_t,x_t_l,index, args, s)
  rmse[, "Lasso"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "Lasso"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  # Lasso DR DML
  iates_mat[,"Lasso MOM DR DML1"] = do.call(cf_dml1,list(mom_dr_glmnet,y_t,d_t,x_t_l,np,x_v_l,index, form_x=NULL, args_tau=args_glmnet,trim))
  # Boosting based methods
  
  args_boost <- arguments[["Boosting"]]
  args    <- args_boost
  args[which(names(args) %in% c("clas_dist","reg_dist"))] <-  NULL
  
  np <-  nuisance_cf_boost(dep_form=dep_form_x,treat_form=treat_form_x,df=tr, d=d_t,index=index, args_method=args_boost)
  rmse[, "Boosting"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "Boosting"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  iates_mat[,"Boosting MOM DR DML1"] = do.call(cf_dml1,list(mom_dr_boost,y_t,d_t,x_t,np,x_v,index, form_x, args_tau = args_boost, trim))


  # Trees
  args_tree <- arguments[["Trees"]]
  
  np <-  nuisance_cf_trees(dep_form=dep_form_x,treat_form=treat_form_x,df=tr, d=d_t,index=index, args_method=args_tree)
  rmse[, "Trees"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "Trees"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  iates_mat[,"Trees MOM DR DML1"] = do.call(cf_dml1,list(mom_dr_trees,y_t,d_t,x_t,np,x_v,index, form_x, args_tau=args_tree, trim))
  
  #Neural Net
  
  args_nnet <- arguments[["Nnet"]]
  
  np = nuisance_cf_nnet(y_t,d_t,x_t,index, args_method=args_nnet)
  rmse[, "Nnet"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "Nnet"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  #DML DR Nnet
  iates_mat[,"Nnet MOM DR DML1"] = do.call(cf_dml1,list(mom_dr_nnet,y_t,d_t,x_t,np,x_v,index,form_x, args_tau = args_nnet, trim))
  
  
  args_bart <- arguments[["Bart"]]
  np = nuisance_cf_bart(y_t,d_t,x_t,index, args_method=args_bart)
  rmse[, "Bart"] <- c(RMSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), RMSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), RMSE(np[d_t==1, "p_hat"], d_t))
  mse[, "Bart"] <- c(MSE(np[d_t==1, "y1_hat"], y_t[d_t==1]), MSE(np[d_t==0, "y0_hat"], y_t[d_t==0]), MSE(np[d_t==1, "p_hat"], d_t))
  
  iates_mat[,"Bart MOM DR DML1"] <- do.call(cf_dml1,list(mom_dr_bart,y_t,d_t,x_t,np,x_v,index,form_x, args_tau = args_bart, trim))
  
  
  response = y_t
  treatment = d_t
  confounders = x_t
  
  bartmc <- bartc(response=response , treatment=treatment,confounders=confounders, keepTrees = TRUE, commonSup.rule="sd")
  library(bartCause)
  predict_ite <- predict(bartmc, newdata=x_v,type="ite") 
  predict_ite <- colMeans(predict_ite) 
  iates_mat[,"BART MChains"] <- predict_ite

  #additional information
  iates_mat[,"index_validation_sample"] = total_index
  iates_mat[,"tau_val_mu"] = tau_val_mu
  iates_mat[,"y_val"] = y_val

  return(list("iates"=iates_mat, "rmse"=rmse, "mse"=mse))
}

