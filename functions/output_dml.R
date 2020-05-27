# References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     
#             "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016 
# Copyright 2017 <Victor Chernozhukov, Denis Chetverikov, Mert 
# Demirer, Esther Duflo, Whitney Newey, James Robins>
#   
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Modifications: Returns p-values. 

#this function returns the median estimates of the DML method
output_p <- function(r, x = data_info$x, n = data_info$n , adj = FALSE, one_p = TRUE){  
  library(stringr)
  library(xtable)
  library(matrixStats)
  library(gtools)
  r <- do.call("rbind", r)
  
  methods <- na.omit(str_extract(colnames(r), "^[^.]+$"))
  r <- r[,1:(2*length(methods))] 
  result           <- matrix(0,5, length(methods))
  colnames(result) <- cbind(t(methods))
  rownames(result) <- cbind("Median ATE", "se(adj)",  "se(median)", "T-stat", "p-value")
  
  r <- data.matrix(r)
  
  result[1,]        <- colQuantiles(r[,1:length(methods)], probs=0.5)
  result[2,]        <- colQuantiles(sqrt(r[,(length(methods)+1):ncol(r)]^2+(r[,1:length(methods)] - colQuantiles(r[,1:length(methods)], probs=0.5))^2), probs=0.5)
  result[3,]        <- colQuantiles(r[,(length(methods)+1):ncol(r)], probs=0.5)
  
  
  if(adj == FALSE){
    result[4,] <- result[1,]/result[3,]
  }else{
    result[4,] <- result[1,]/result[2,]
  }
  
  no_controls = str_count(x, pattern = "\\+")+ 1
  
  p = no_controls + 1 
  
  if(one_p==FALSE){
    df <- n - p   
  }else{
    df <- n - 1 
  }
  
  result[5,] <- 2*pt(-abs(result[4,]),df=df) 
  
  result <- round(result, digits = 3)
  result_table <- result
  result_table[1,] <- paste(result_table[1,], stars.pval(result_table[5,]), sep="")
  result_table_n <- result_table
  
  for(i in 1:ncol(result_table)){
    for(j in seq(2,nrow(result_table),3)){
      
      result_table[j,i] <- paste("(", result_table[j,i], ")", sep="")
      
    }
    for(j in seq(3,nrow(result_table),3)){
      
      result_table[j,i] <- paste("(", result_table[j,i], ")", sep="")
      
    }
  }
  
  
  output <- list(result, xtable(result_table[1:3,]))
  return(list("table" = r, "output"=output))
}