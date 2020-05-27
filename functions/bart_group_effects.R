# Function to compute group effects based on the quantiles of the heterogeneous treatment effect distribution

group_effects <- function(test_intervals, thres=0.2){
  breaks    <- quantile(test_intervals, seq(0,1, thres),  include.lowest =T, include.highest=T. ,na.rm=T)
  if(TRUE%in%duplicated(breaks)){
    G1 <- mean(test_intervals)
    G2 <- mean(test_intervals)
    G3 <- mean(test_intervals)
    G4 <- mean(test_intervals)
    G5 <- mean(test_intervals)
  }else{
    breaks[1] <- breaks[1] - 0.001
    breaks[6] <- breaks[6] + 0.001
    SG        <- cut(test_intervals, breaks = breaks)
    SGX       <- model.matrix(~-1+SG)
    SGX <- as.data.frame(SGX)
    colnames(SGX) <- c("G1", "G2", "G3", "G4", "G5")
    SGX$test_intervals <- na.omit(test_intervals)
    
    #Group Effects
    G1 = mean(SGX$test_intervals[SGX$G1==1])
    G2 = mean(SGX$test_intervals[SGX$G2==1])
    G3 = mean(SGX$test_intervals[SGX$G3==1])
    G4 = mean(SGX$test_intervals[SGX$G4==1])
    G5 = mean(SGX$test_intervals[SGX$G5==1])
  }
  return(c(G1, G2, G3, G4, G5))
}

gates_bart <- function(gates_array,groups=5, alpha=0.05){
  gates_est <- matrix(NA,groups,3)
  rownames(gates_est) <- c("G1", "G2", "G3", "G4", "G5")
  colnames(gates_est) <- c("est", "ci.lower", "ci.upper")
  gates_posterior <- array(NA, dim=c(nrow(gates_est),dim(gates_array)[2],dim(gates_array)[3]))
  for(i in 1:dim(gates_array)[2]){
    gates_posterior[,i,] <- apply(gates_array[,i,],2,group_effects)
  }
  for(i in 1:groups){
    gates_est[i,"est"] <- mean(as.vector(gates_posterior[i,,]))
    gates_est[i,"ci.lower"] <- mean(as.vector(gates_posterior[i,,])) - qnorm(1-alpha)*sd(as.vector(gates_posterior[i,,]))
    gates_est[i,"ci.upper"] <- mean(as.vector(gates_posterior[i,,])) +qnorm(1-alpha)*sd(as.vector(gates_posterior[i,,]))
  }
  return(gates_est)
}
