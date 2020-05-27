#function which returns the ATE and the GATE point estimates and standard errors; GATEs are based on the five quantiles of the  heterogeneous treatment function
cf <- function(data,y,d){
	X <- data[,!colnames(data)%in%c(y,d)]
  X <- model.matrix(~.-1,X)
  Y <- data[,y]
  W <- data[,d]
  
  forest.W <- regression_forest(X, W, tune.parameters = "all")
  W.hat <- predict(forest.W)$predictions
  
  forest.Y <- regression_forest(X, Y, tune.parameters = "all")
  Y.hat <- predict(forest.Y)$predictions
  
  cf_raw = causal_forest(X, Y, W,
                         Y.hat = Y.hat, W.hat = W.hat,
                         tune.parameters = "all")
  
  cf_cate_raw <- predict(cf_raw, estimate.variance = TRUE)
  
  breaks_raw    <- quantile(cf_cate_raw$predictions, seq(0,1, 0.2),  include.lowest =TRUE, include.highest=TRUE)
  cf_ate_o_raw <- average_treatment_effect(cf_raw, target.sample = "overlap")
  
  cf_test_raw <- test_calibration(cf_raw)
  
  if(TRUE%in%duplicated(breaks_raw)){
    G1_raw <- mean(cf_cate_raw$predictions)
    G2_raw <- mean(cf_cate_raw$predictions)
    G3_raw <- mean(cf_cate_raw$predictions)
    G4_raw <- mean(cf_cate_raw$predictions)
    G5_raw <- mean(cf_cate_raw$predictions)
    cf_gates_raw <- cbind(G1_raw, G2_raw, G3_raw, G4_raw, G5_raw)
    cf_gates_o_raw <- cbind(G1_raw, G2_raw, G3_raw, G4_raw, G5_raw)
  }else{
    breaks_raw[1] <- breaks_raw[1]-0.01 
    breaks_raw[6] <- breaks_raw[6] +0.01
    SG_raw        <- cut(cf_cate_raw$predictions, breaks = breaks_raw)
    SGX_raw       <- model.matrix(~-1+SG_raw)
    SGX_raw       <- as.data.frame(SGX_raw)
    colnames(SGX_raw) <- c("G1", "G2", "G3", "G4", "G5")
  
    cf_gate_1_o_raw <- average_treatment_effect(cf_raw,subset=SGX_raw$G1==1, target.sample = "overlap")
    cf_gate_2_o_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G2==1, target.sample = "overlap")
    cf_gate_3_o_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G3==1, target.sample = "overlap")
    cf_gate_4_o_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G4==1, target.sample = "overlap")
    cf_gate_5_o_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G5==1, target.sample = "overlap")
    cf_gates_o_raw <- cbind(cf_gate_1_o_raw, cf_gate_2_o_raw, cf_gate_3_o_raw, cf_gate_4_o_raw, cf_gate_5_o_raw)
  }
  
	  return(list("cf_cate_raw"=cf_cate_raw, "cf_ate_o_raw" = cf_ate_o_raw,"cf_gates_o_raw"=cf_gates_o_raw, "hete_test_raw"=cf_test_raw))
}
