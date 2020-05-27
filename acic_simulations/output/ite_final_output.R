# Plots the  PEHE of the ITE results
rm(list=ls())
library("gridExtra")
library("grid")
library("ggplot2")

knobs <- read.csv("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/knobs_all.csv", sep=";")
load("~/Documents/GitHub/MSc_Thesis_CausalML/acic_simulations/output/ite_final_output.RData")

plot <-function(data, y_name, methods, title, r=100, pos=0.75, lim=0.8, text_size=17){
  k <- nrow(data)
  data <- data.frame(round(colMeans(data),2))
  colnames(data) <- y_name
  ggplot(data, aes(x=methods, y=data[,1])) +
    theme_bw() +geom_point(shape=18, color="blue", size=4)+ggtitle(title)+ylab(y_name)+ylim(0,lim)+theme(axis.text.x = element_text(angle = 90, size=text_size))+theme(axis.title.x=element_blank())+annotate("text", label = paste("k =", k,", r = ",r, sep=""), x = "DR MOM RF", y = pos, size=5)+ theme(plot.title = element_text(size = text_size, face = "bold")) + theme(axis.title = element_text(size = text_size, face = "bold"))+ theme(axis.text.y = element_text(size=text_size))+theme(axis.text.x = element_text(size=text_size))
}

plot(pehe_matrix, "pehe", colnames(pehe_matrix), "Overview", r=100)

plot1  <- plot(pehe_matrix[which(knobs$Treatment.Model=="polynomial") ,], "pehe", colnames(pehe_matrix), "Trt. model polynomial")
plot2 <- plot(pehe_matrix[which(knobs$Treatment.Model=="step") ,], "PEHE", colnames(pehe_matrix), "Trt. model step")
plot3 <- plot(pehe_matrix[which(knobs$Percent.treated=="low"),], "PEHE", colnames(pehe_matrix), "% treated low")
plot4 <- plot(pehe_matrix[which(knobs$Percent.treated=="high"),], "PEHE", colnames(pehe_matrix), "% treated high")
plot5 <- plot(pehe_matrix[which(knobs$Overlap=="full"),], "PEHE", colnames(pehe_matrix), "Full overlap")
plot6 <- plot(pehe_matrix[which(knobs$Overlap=="penalize"),], "PEHE", colnames(pehe_matrix), "Penalized Overlap")
plot7 <- plot(pehe_matrix[which(knobs$Trt.Rsp.Alignment=="low"),], "PEHE", colnames(pehe_matrix), "Low Alignment")
plot8 <- plot(pehe_matrix[which(knobs$Trt.Rsp.Alignment=="high"),], "PEHE", colnames(pehe_matrix), "High Alignment")
plot9 <- plot(pehe_matrix[which(knobs$heterogneity=="low"),], "PEHE", colnames(pehe_matrix), "Low Heterogeneity")
plot10 <- plot(pehe_matrix[which(knobs$heterogneity=="high"),], "PEHE", colnames(pehe_matrix), "High Heterogeneity")

grid.arrange(plot1, plot2, nrow=1, ncol=2)
grid.arrange(plot3, plot4, nrow=1, ncol=2)
grid.arrange(plot5, plot6, nrow=1, ncol=2)
grid.arrange(plot7, plot8, nrow=1, ncol=2)
grid.arrange(plot9, plot10, nrow=1, ncol=2)


