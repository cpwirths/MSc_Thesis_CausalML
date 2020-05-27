rm(list=ls())  
# Returns BART MChains output for the GSS application
library(xtable)
library(ggplot2)
library(bartCause)
# Load data and functions -------------------------------------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/gss_application/" 
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 

setwd(dir_path)
load(paste(dir_path,"GSS_input.RData", sep = "", collapse = NULL))
source(paste(functions_path,"bart_cates.R", sep = "", collapse = NULL))
source(paste(functions_path,"bart_group_effects.R", sep = "", collapse = NULL))

response <- data$response
treatment <- data$treatment
confounders <- data[, !colnames(data)%in%c("treatment", "response")]

# Estimation --------------------------------------------------------------
bart <- bartc(response, treatment, confounders, estimand = "ate", method.rsp=c("bart"), method.trt=c("none"))
bart_ate_output <- summary(bart) #ATE

ite <- fitted(bart,value= c("indiv.diff"))
hist(ite, xlab="ITE", main="Histogram BART MChains")

#Retrieve estimation of individual treatment effects
if(!is.null(bart$samples.indiv.diff)){
  samples.indiv.diff <- bart$samples.indiv.diff  
}else{
  samples.indiv.diff <- bart$mu.hat.obs - bart$mu.hat.cf
  for(i in 1:dim(samples.indiv.diff)[1]){
    for(j in 1:dim(samples.indiv.diff)[2])
      samples.indiv.diff[i,j, which(treatment==0)] <- bart$mu.hat.cf[i,j, which(treatment==0)] - bart$mu.hat.obs[i,j, which(treatment==0)]
  }
 samples.indiv.diff <- aperm(samples.indiv.diff,perm=c(3,1,2))
}

bart_gates_output <- gates_bart(samples.indiv.diff) #GATE

# Heterogeneity analysis of control variables ---------------------------

bart_rac_1_yes <- bart_cate(samples.indiv.diff, data$RACDIF1=="yes")
bart_rac_2_yes <- bart_cate(samples.indiv.diff, data$RACDIF2=="yes")
bart_rac_3_yes <- bart_cate(samples.indiv.diff, data$RACDIF3=="yes")
bart_rac_4_yes <- bart_cate(samples.indiv.diff, data$RACDIF4=="yes")

bart_rac_1_no <- bart_cate(samples.indiv.diff, data$RACDIF1=="no")
bart_rac_2_no <- bart_cate(samples.indiv.diff, data$RACDIF2=="no")
bart_rac_3_no <- bart_cate(samples.indiv.diff, data$RACDIF3=="no")
bart_rac_4_no <- bart_cate(samples.indiv.diff, data$RACDIF4=="no")

rac_bart <- cbind(bart_rac_1_yes, bart_rac_1_no, bart_rac_2_yes, bart_rac_2_no, bart_rac_3_yes, bart_rac_3_no, bart_rac_4_yes, bart_rac_4_no) 
rac_bart <- round(rac_bart, 2)

plot_cate_bart(array=samples.indiv.diff, column=data$year, numeric=T, title="BART MC", y_name="CATE", x_name="year")
plot_cate_bart(array=samples.indiv.diff,column=data$partyid, numeric=F, title="BART MC", y_name="CATE", x_name="PARTY ID", exclude ="OTHER PARTY", angle=90)
plot_cate_bart(array=samples.indiv.diff,column=data$polviews, numeric=F, title="BART MC", y_name="CATE", x_name="Political Views", angle=90)

save(bart_ate_output, bart_gates_output, file="bart_gss_ate_gates.RData")




