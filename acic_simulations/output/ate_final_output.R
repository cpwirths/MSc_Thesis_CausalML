# Plots the rmse, coverage, bias and interval length graphs of the ATE results

rm(list=ls())
library("gridExtra")
library("grid")
library("ggplot2")

load("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/output/ate_final_output.RData")
knobs <- read.csv("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/knobs_all.csv", sep=";")
b_index <-c(27,55,40,21,31,29,28) #Benchmark analysis

# Tables ------------------------------------------------------------------

#summary overview
order <- order(colnames(rmse_matrix))
xtable(rbind(colMeans(bias_matrix[b_index,order(colnames(bias_matrix))]), colMeans(rmse_matrix[,order(colnames(rmse_matrix))]), colMeans(coverage_matrix[,order(colnames(coverage_matrix))]), colMeans(interval_matrix[,order(colnames(interval_matrix))])))

#Benchmark analysis
xtable(rmse_matrix[b_index,order(colnames(rmse_matrix))])
xtable(coverage_matrix[b_index,order(colnames(coverage_matrix))])
xtable(round(bias_matrix[b_index,order(colnames(bias_matrix))], digits=3))
xtable(interval_matrix[b_index,order(colnames(interval_matrix))])


# plotting results --------------------------------------------------------
plot <-function(data, y_name, methods, title, min=0, lim=0.13, shape=18,fill="blue", r=100, pos=0.12, line=0, blank=F, marker_size=3, textsize=17){
  means <- data.frame(round(colMeans(data,na.rm = T),3))
  colnames(means) <- y_name
  means$methods <- as.factor(methods)
  means$methods <- factor(means$methods, levels = means$methods)
  means$methods <- methods
  library(ggplot2)
  if(blank==F){
    ggplot(means, aes(x=methods, y=means[,1])) + theme_bw() +
      geom_point(shape=shape, color=fill, size=marker_size)  +ggtitle(title)+ylab(y_name)+ylim(min,lim)+ theme(axis.text.x = element_text(angle = 90, size=textsize))+  theme(axis.text.y = element_text(size=textsize))+ theme(axis.title.x=element_blank())+annotate("text", label = paste("k = ", nrow(data),", r = ",r, sep="") , x = "DML Lasso", y = pos, size=5)+ theme(plot.title = element_text(size = textsize, face = "bold")) + theme(axis.title = element_text(size = textsize, face = "bold"))+geom_hline(yintercept=line)
  }else{
    ggplot(means, aes(x=methods, y=means[,1])) + theme_bw() +
      geom_point(shape=shape, fill=fill, size=marker_size) +ylab(y_name)+ylim(min,lim)+theme(axis.title.x=element_blank()) + theme(axis.text.y = element_text(size=textsize)) + annotate("text", label = paste("k = ", nrow(data),", r = ",r, sep=""), x = "DML Lasso", y = pos, size=5)+ theme(plot.title = element_text(size = textsize, face = "bold")) + theme(axis.title = element_text(size = textsize, face = "bold"))+geom_hline(yintercept=line)+ theme(axis.text.x = element_blank())+scale_x_discrete(position = "top")
  }
}

#overview
all_rmse <- plot(rmse_matrix, "rmse", colnames(rmse_matrix), "RMSE ATE Overview", min=0, lim=0.11, pos=0.1, blank = F, marker_size = 4)
all_coverage <-  plot(coverage_matrix, "coverage", colnames(coverage_matrix), "Coverage ATE Overview",min=0.3, lim=1, pos=0.31, shape=17, fill="red", line=0.95)
all_bias <- plot(bias_matrix, "bias", colnames(bias_matrix), "Bias ATE Overview", min=-0.04, lim=0.02, pos=0.01, shape=15, fill="orange")
all_interval <- plot(interval_matrix, "int length", colnames(interval_matrix), "Interval ATE Overview",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
grid.arrange(all_bias, all_rmse, ncol=2)
grid.arrange(all_coverage, all_interval, ncol=2)


#rmse
plot1_r <- plot(rmse_matrix[which(knobs$Response.Model=="linear" & knobs$Treatment.Model=="linear") ,], "rmse ATE", colnames(rmse_matrix), "Rsp. & Trt. model linear",marker_size = 4)
plot2_r <- plot(rmse_matrix[which(knobs$Response.Model!="linear" & knobs$Treatment.Model!="linear"),], "rmse ATE", colnames(rmse_matrix), "Rsp. & Trt. model non-linear", marker_size = 4)
plot3_r <- plot(rmse_matrix[which(knobs$Percent.treated=="low"),], "rmse ATE", colnames(rmse_matrix), "% treated low", marker_size = 4)
plot4_r <- plot(rmse_matrix[which(knobs$Percent.treated=="high"),], "rmse ATE", colnames(rmse_matrix), "% treated high", marker_size = 4)
plot5_r <- plot(rmse_matrix[which(knobs$Overlap=="full"),], "rmse ATE", colnames(rmse_matrix), "Full overlap", marker_size = 4)
plot6_r <- plot(rmse_matrix[which(knobs$Overlap=="penalize"),], "rmse ATE", colnames(rmse_matrix), "Penalized Overlap", marker_size = 4)
plot7_r <- plot(rmse_matrix[which(knobs$Trt.Rsp.Alignment=="low"),], "rmse ATE", colnames(rmse_matrix), "Low Alignment", marker_size = 4)
plot8_r <- plot(rmse_matrix[which(knobs$Trt.Rsp.Alignment=="high"),], "rmse ATE", colnames(rmse_matrix), "High Alignment", marker_size = 4)
plot9_r <- plot(rmse_matrix[which(knobs$heterogneity=="low"),], "rmse ATE", colnames(rmse_matrix), "Low Heterogeneity", marker_size = 4)
plot10_r <- plot(rmse_matrix[which(knobs$heterogneity=="high"),], "rmse ATE", colnames(rmse_matrix), "High Heterogeneity", marker_size = 4)
grid.arrange(plot1_r, plot2_r, ncol=2)
grid.arrange(plot3_r, plot4_r, ncol=2)
grid.arrange(plot5_r, plot6_r, ncol=2)
grid.arrange(plot7_r, plot8_r, ncol=2)
grid.arrange(plot9_r, plot10_r, ncol=2)

#coverage
plot1_c <- plot(coverage_matrix[which(knobs$Response.Model=="linear" & knobs$Treatment.Model=="linear") ,], "coverage ATE", colnames(coverage_matrix), "Rsp. & Trt. model linear", min=0.1, lim=1, shape=17, fill="red", line=0.95, blank=F, pos=0.2)
plot2_c <- plot(coverage_matrix[which(knobs$Response.Model!="linear" & knobs$Treatment.Model!="linear"),], "coverage ATE", colnames(coverage_matrix), "Rsp. & Trt. model non-linear",min=0.1,lim=1, shape=17, fill="red", line=0.95, blank=F, pos=0.2)
plot3_c <- plot(coverage_matrix[which(knobs$Percent.treated=="low"),], "coverage ATE", colnames(coverage_matrix), "% treated low", min=0.25, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank=F)
plot4_c <- plot(coverage_matrix[which(knobs$Percent.treated=="high"),], "coverage ATE", colnames(coverage_matrix), "% treated high", min=0.25, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank=F)
plot5_c <- plot(coverage_matrix[which(knobs$Overlap=="full"),], "coverage ATE", colnames(coverage_matrix), "Full overlap", min=0.2, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank=F)
plot6_c <- plot(coverage_matrix[which(knobs$Overlap=="penalize"),], "coverage ATE", colnames(coverage_matrix), "Penalized Overlap", min=0.2, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank=F)
plot7_c <- plot(coverage_matrix[which(knobs$Trt.Rsp.Alignment=="low"),], "coverage ATE", colnames(coverage_matrix), "Low Alignment", min=0.2, pos=0.3, lim=1, shape=17, fill="red",line=0.95, blank=F)
plot8_c <- plot(coverage_matrix[which(knobs$Trt.Rsp.Alignment=="high"),], "coverage ATE", colnames(coverage_matrix), "High Alignment", min=0.2,pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank = F)
plot9_c <- plot(coverage_matrix[which(knobs$heterogneity=="low"),], "coverage ATE", colnames(coverage_matrix), "Low Heterogeneity", min=0.2, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank = F)
plot10_c <- plot(coverage_matrix[which(knobs$heterogneity=="high"),], "coverage ATE", colnames(coverage_matrix), "High Heterogeneity", min=0.2, pos=0.3, lim=1, shape=17, fill="red", line=0.95, blank = F)
grid.arrange(plot1_c, plot2_c, ncol=2)
grid.arrange(plot3_c, plot4_c, ncol=2)
grid.arrange(plot5_c, plot6_c, ncol=2)
grid.arrange(plot7_c, plot8_c, ncol=2)
grid.arrange(plot9_c, plot10_c, ncol=2)


#bias
plot1_m <- plot(bias_matrix[which(knobs$Response.Model=="linear" & knobs$Treatment.Model=="linear") ,], "bias ATE", colnames(bias_matrix), "Rsp. & Trt. model linear", min=-0.04, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot2_m <- plot(bias_matrix[which(knobs$Response.Model!="linear" & knobs$Treatment.Model!="linear"),], "bias ATE", colnames(bias_matrix), "Rsp. & Trt. model non-linear",min=-0.04, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot3_m <- plot(bias_matrix[which(knobs$Percent.treated=="low"),], "bias ATE", colnames(bias_matrix), "% treated low",min=-0.04, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot4_m <- plot(bias_matrix[which(knobs$Percent.treated=="high"),], "bias ATE", colnames(bias_matrix), "% treated high",min=-0.04, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot5_m <- plot(bias_matrix[which(knobs$Overlap=="full"),], "bias ATE", colnames(bias_matrix), "Full overlap",min=-0.05, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot6_m <- plot(bias_matrix[which(knobs$Overlap=="penalize"),], "bias ATE", colnames(bias_matrix), "Penalized Overlap",min=-0.05, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot7_m <- plot(bias_matrix[which(knobs$Trt.Rsp.Alignment=="low"),], "bias ATE", colnames(bias_matrix), "Low Alignment",min=-0.06, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot8_m <- plot(bias_matrix[which(knobs$Trt.Rsp.Alignment=="high"),], "bias ATE", colnames(bias_matrix), "High Alignment",min=-0.06, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot9_m <- plot(bias_matrix[which(knobs$heterogneity=="low"),], "bias ATE", colnames(bias_matrix), "Low Heterogeneity",min=-0.05, lim=0.03, pos=0.01, shape=15, fill="darkorange")
plot10_m <- plot(bias_matrix[which(knobs$heterogneity=="high"),], "bias ATE", colnames(bias_matrix), "High Heterogeneity",min=-0.05, lim=0.03, pos=0.01, shape=15, fill="darkorange")

grid.arrange(plot1_m, plot2_m, ncol=2)
grid.arrange(plot3_m, plot4_m, ncol=2)
grid.arrange(plot5_m, plot6_m, ncol=2)
grid.arrange(plot7_m, plot8_m, ncol=2)
grid.arrange(plot9_m, plot10_m, ncol=2)

#interval
plot1_i <- plot(interval_matrix[which(knobs$Response.Model=="linear" & knobs$Treatment.Model=="linear") ,], "int length ATE", colnames(interval_matrix), "Rsp. & Trt. model linear", min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot2_i <- plot(interval_matrix[which(knobs$Response.Model!="linear" & knobs$Treatment.Model!="linear"),], "int length ATE", colnames(interval_matrix), "Rsp. & Trt. model non-linear",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot3_i <- plot(interval_matrix[which(knobs$Percent.treated=="low"),], "int length ATE", colnames(interval_matrix), "% treated low",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot4_i <- plot(interval_matrix[which(knobs$Percent.treated=="high"),], "int length ATE", colnames(interval_matrix), "% treated high",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue",line=NULL)
plot5_i <- plot(interval_matrix[which(knobs$Overlap=="full"),], "int length ATE", colnames(interval_matrix), "Full overlap",min=0.03, lim=0.11, pos=0.1, shape=19, fill="steelblue", line=NULL)
plot6_i <- plot(interval_matrix[which(knobs$Overlap=="penalize"),], "int length ATE", colnames(interval_matrix), "Penalized Overlap",min=0.03, lim=0.13, pos=0.1, shape=19, fill="steelblue", line=NULL)
plot7_i <- plot(interval_matrix[which(knobs$Trt.Rsp.Alignment=="low"),], "int length ATE", colnames(interval_matrix), "Low Alignment",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot8_i <- plot(interval_matrix[which(knobs$Trt.Rsp.Alignment=="high"),], "int length ATE", colnames(interval_matrix), "High Alignment",min=0.04, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot9_i <- plot(interval_matrix[which(knobs$heterogneity=="low"),], "int length ATE", colnames(interval_matrix), "Low Heterogeneity",min=0.03, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)
plot10_i <- plot(interval_matrix[which(knobs$heterogneity=="high"),], "int length ATE", colnames(interval_matrix), "High Heterogeneity",min=0.03, lim=0.12, pos=0.095, shape=19, fill="steelblue", line=NULL)

grid.arrange(plot1_i, plot2_i, ncol=2)
grid.arrange(plot3_i, plot4_i, ncol=2)
grid.arrange(plot5_i, plot6_i, ncol=2)
grid.arrange(plot7_i, plot8_i, ncol=2)
grid.arrange(plot9_i, plot10_i, ncol=2)


