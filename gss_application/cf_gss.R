#Returns the Causal Forest output for the GSS application
rm(list=ls())  
library(grf)

# Load data and functions -------------------------------------------------
dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/gss_application/" 
functions_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/functions/" 

setwd(dir_path)
load(paste(dir_path,"GSS_input.RData", sep = "", collapse = NULL))
source(paste(functions_path,"cf_cates.R", sep = "", collapse = NULL))

X <-  data[, !colnames(data)%in%c("treatment", "response")]
X <- model.matrix(~.-1,X)
Y <- as.numeric(data$response)
W <- as.numeric(data$treatment)

# Estimation --------------------------------------------------------------
forest.W <- regression_forest(X, W, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions
forest.Y <- regression_forest(X, Y, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

cf_raw = causal_forest(X, Y, W,
                       tune.parameters = "all")

cf_cate_raw <- predict(cf_raw, estimate.variance = TRUE)
hist(cf_cate_raw$predictions, xlab= "ITE", main="Histogram Causal Forest")
test_calibration(cf_raw)

cf_ate_output <- average_treatment_effect(cf_raw) #ATE

breaks_raw    <- quantile(cf_cate_raw$predictions, seq(0,1, 0.2),  include.lowest =TRUE, include.highest=TRUE)
breaks_raw[1] <- breaks_raw[1]-0.01 
breaks_raw[6] <- breaks_raw[6] +0.01
SG_raw        <- cut(cf_cate_raw$predictions, breaks = breaks_raw)
SGX_raw       <- model.matrix(~-1+SG_raw)
SGX_raw       <- as.data.frame(SGX_raw)
colnames(SGX_raw) <- c("G1", "G2", "G3", "G4", "G5")

cf_gate_1_raw <- average_treatment_effect(cf_raw,subset=SGX_raw$G1==1)
cf_gate_2_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G2==1)
cf_gate_3_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G3==1)
cf_gate_4_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G4==1)
cf_gate_5_raw <- average_treatment_effect(cf_raw, subset=SGX_raw$G5==1)
cf_gates_raw <- cbind(cf_gate_1_raw, cf_gate_2_raw, cf_gate_3_raw, cf_gate_4_raw, cf_gate_5_raw)
lb <- cf_gates_raw["estimate",] - cf_gates_raw["std.err",]*1.96
ub <- cf_gates_raw["estimate",] + cf_gates_raw["std.err",]*1.96

cf_gates_output <- round(rbind(cf_gates_raw, lb, ub),2) #GATE
colnames(cf_gates_output) <- c("G1","G2","G3", "G4", "G5")


# Heterogeneity analysis of control variables ---------------------------

varimp = variable_importance(cf_raw)
selected.idx = which(varimp > mean(varimp))

rac_1_no <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF1=="no")
rac_1_no <- c(rac_1_no["estimate"],rac_1_no["estimate"] - 1.96*rac_1_no["std.err"],rac_1_no["estimate"] + 1.96*rac_1_no["std.err"]) 
names(rac_1_no) <- c("est", "lb", "ub")

rac_1_yes <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF1=="yes")
rac_1_yes <- c(rac_1_yes["estimate"],rac_1_yes["estimate"] - 1.96*rac_1_yes["std.err"],rac_1_yes["estimate"] + 1.96*rac_1_yes["std.err"]) 
names(rac_1_yes) <- c("est", "lb", "ub")

rac_2_no <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF2=="no")
rac_2_no <- c(rac_2_no["estimate"],rac_2_no["estimate"] - 1.96*rac_2_no["std.err"],rac_2_no["estimate"] + 1.96*rac_2_no["std.err"]) 
names(rac_2_no) <- c("est", "lb", "ub")

rac_2_yes <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF2=="yes")
rac_2_yes <- c(rac_2_yes["estimate"],rac_2_yes["estimate"] - 1.96*rac_2_yes["std.err"],rac_2_yes["estimate"] + 1.96*rac_2_yes["std.err"]) 
names(rac_2_yes) <- c("est", "lb", "ub")

rac_3_no <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF3=="no")
rac_3_no <- c(rac_3_no["estimate"],rac_3_no["estimate"] - 1.96*rac_3_no["std.err"],rac_3_no["estimate"] + 1.96*rac_3_no["std.err"]) 
names(rac_3_no) <- c("est", "lb", "ub")

rac_3_yes <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF3=="yes")
rac_3_yes <- c(rac_3_yes["estimate"],rac_3_yes["estimate"] - 1.96*rac_3_yes["std.err"],rac_3_yes["estimate"] + 1.96*rac_3_yes["std.err"]) 
names(rac_3_yes) <- c("est", "lb", "ub")

rac_4_no <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF4=="no")
rac_4_no <- c(rac_4_no["estimate"],rac_4_no["estimate"] - 1.96*rac_4_no["std.err"],rac_4_no["estimate"] + 1.96*rac_4_no["std.err"]) 
names(rac_4_no) <- c("est", "lb", "ub")

rac_4_yes <- average_treatment_effect(cf_raw, target.sample = "overlap", subset = data$RACDIF4=="yes")
rac_4_yes <- c(rac_4_yes["estimate"],rac_4_yes["estimate"] - 1.96*rac_4_yes["std.err"],rac_4_yes["estimate"] + 1.96*rac_4_yes["std.err"]) 
names(rac_4_yes) <- c("est", "lb", "ub")
rac_cf <- cbind(rac_1_yes, rac_1_no,rac_2_yes, rac_2_no, rac_3_yes, rac_3_no, rac_4_yes, rac_4_no)
rac_cf <- round(rac_cf,2) 
cis <- function(input){apply(input,2,function(x){
  paste0("(",paste(x, collapse=", "),")")
})}

rac_cf <- t(rbind(rac_cf[1,], cis(rac_cf[2:3,])))

plot_cate_cf(cf=cf_raw, column=data$year, title="CF", y_name="CATE", x_name="year", angle=90)
plot_cate_cf(cf=cf_raw, column=data$partyid, numeric=F, title="CF", y_name="CATE", x_name="PARTY ID", exclude = "OTHER PARTY", angle=90)
plot_cate_cf(cf=cf_raw, column=data$polviews, numeric=F, title="CF", y_name="CATE", x_name="Political Views", angle=90)

save(cf_ate_output, cf_gates_output, file="cf_gss_ate_gates.RData")





