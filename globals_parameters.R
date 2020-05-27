# Global ML parameters used across al Causal ML methods, when general Machine Learning algorithms are used as metalearner
Boosting     <- list(n.minobsinnode = 1, bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=2, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=250, na.action=na.omit, replace=TRUE)
Lasso       <- list(s = "lambda.min",intercept = TRUE, lambda.min=0.02)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")
Bart        <- list()
arguments    <- list(Boosting=Boosting, Forest=Forest, Lasso=Lasso, Nnet=Nnet, Trees=Trees, Bart=Bart)


