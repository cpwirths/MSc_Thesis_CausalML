# Reference: Dorie, Vincent, et al. "Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition." Statistical Science 34.1 (2019): 43-68.
# this file creates the datasets of the competition
# in total there are 77 simulation knobs
# each simulation knob conatins 100 replications
# additionally an OLS estimate is saved as benchmark
rm(list = ls())
library(aciccomp2016)
data_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/data_simulations/" 
setwd(data_path)

#create simulation datasets
for(k in seq(77)){  
print(k)
datasets <- list()
#stores ols estimates inlcuding std.error, t value and p value
ate_ols_coeff <- rep(NA, 100)
ate_ols_sum <- matrix(NA,100,4)
true_ate_mu <- rep(NA, 100)
true_cate_mu <-  matrix(NA,nrow(input_2016),100)
response <- matrix(NA,nrow(input_2016),100)
treatment <- matrix(NA,nrow(input_2016),100)
sd_response <-rep(NA,100)
 
for(i in 1:100){
  print(i)
  sim <- dgp_2016(input_2016, k, i)
  X <- model.matrix(~.-1,input_2016)
  X <- data.frame(X)
  df <- X
  df$y <- sim$y
  df$z <- sim$z
  response[,i] <- sim$y
  treatment[,i] <- sim$z
  sd_response[i] <- sd(sim$y)
  datasets[[i]] <- df
  
  true_ate_mu[i] <- mean(sim$mu.1 - sim$mu.0)
  true_cate_mu[,i] <- sim$mu.1 - sim$mu.0
  
  #save benchmark ols value
  fit <- lm(y~.,df)
  summary_fit <- summary(fit)
  ate_ols_coeff[i] <- coef(fit)["z"]
  ate_ols_sum[i,] <- summary_fit$coefficients["z",]
}
save.image(paste ("acic_",k,"_data.RData", sep = "", collapse = NULL))
}
      
