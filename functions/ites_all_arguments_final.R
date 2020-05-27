#References: Knaus, Lechner, Strittmatter (2018). Machine Learning Estimation of Heterogeneous Causal Effects: Empirical Monte Carlo Evidence
#Copyright <2018> <Knaus, Lechner, Strittmatter> R-packages: Package: CATEs, MIT license

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Modification Function which computes the individual treatment effects for the simulation; the ITEs are validated on a 50% hold out set
iates_all <- function(data,y,d, xx, xl, arguments, trim){
  
 source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/ites_functions_final.R")
 source("~/Documents/GitHub/MSc_Thesis_CausalML/functions/ites_specific_final.R")
 
  J <- 2
  split     <- runif(nrow(data))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/J)),include.lowest = TRUE))  
  tr = data[cvgroup==1,]
  val = data[cvgroup!=1,]
  total_index <- as.numeric(row.names(val))
  n <- nrow(val)
  p <- ncol(val)
  tau_tr_mu = true_cate_mu[cvgroup==1,i] 
  tau_val_mu = true_cate_mu[cvgroup!=1,i]
  trim <- trim
  iates_mat = IATEs(tr, val, y, d, xx, xl, arguments, trim, total_index, tau_tr_mu, tau_val_mu, n,p)
  return(list("iates_mat"=iates_mat$iates, "rmse"=iates_mat$rmse, "mse"=iates_mat$mse))
}
