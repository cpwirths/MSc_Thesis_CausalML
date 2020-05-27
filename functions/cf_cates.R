#Function which plots the CATEs point estimates and confidence intervals of the Causal Forest 
#CATEs are based on every distinct value in a column; e.g if the column represents year 2010-2018, the function plots the CATE for each year and the corresponding confidence interval

library(ggplot2)

plot_cate_cf <-function(cf,column, numeric=T, x_name="x_name", y_name="y_name", title="title", textsize=12, exclude=NULL,angle=0){
  if(numeric==T){
    if(is.factor(column)){
      vec <- as.numeric(levels(column))    
    }else{
      vec <- unique(column) 
    }
  }else{
    vec <- levels(column)
  }
  if(!is.null(exclude)){
    which(vec%in%exclude)
    vec <- vec[-which(vec%in%exclude)]  
  }
  matrix <- matrix(NA,length(vec), 4)
  for(i in vec){
    index <- which(i==vec)
    print(i)
    est <- average_treatment_effect(cf, target.sample = "overlap", subset = column==i)
    matrix[index,1:2] <- est
    matrix[index,3] <- est[1] - 1.96 * est[2]
    matrix[index,4] <- est[1] + 1.96 * est[2]
  }
  matrix <- data.frame(matrix)
  colnames(matrix) <- c("est", "std.error", "lb", "ub") 
  if(numeric==T){
    x_axis <- vec  
  }else{
    x_axis <- factor(vec, levels=vec) 
  }
  if(numeric){
    p <- ggplot(matrix, aes(x=x_axis, y=est))+
      geom_line(data=matrix,group=1)+geom_point()+
      geom_ribbon(data=matrix,aes(ymin=lb,ymax=ub), alpha=0.3,group=1)+theme_bw()+ggtitle(title)+  theme(axis.text.x = element_text(angle = angle, size=textsize, vjust = 0.5)) + scale_x_continuous(x_name, labels = as.character(x_axis), breaks = x_axis)+ylab(y_name)+ xlab(x_name)+  theme(axis.text.y = element_text(size=textsize))+theme(plot.title = element_text(size = textsize, face = "bold"))+theme(axis.title = element_text(size = textsize, face = "bold"))+xlab("yes")
  }else{p <- ggplot(matrix, aes(x=x_axis, y=est))+
    geom_line(data=matrix,group=1)+geom_point()+
    geom_ribbon(data=matrix,aes(ymin=lb,ymax=ub), alpha=0.3,group=1)+theme_bw()+ggtitle(title)+  theme(axis.text.x = element_text(angle = angle, size=textsize, vjust = 0.5)) +ylab(y_name)+ xlab(x_name) + theme(axis.text.y = element_text(size=textsize))+theme(plot.title = element_text(size = textsize, face = "bold"))+theme(axis.title = element_text(size = textsize, face = "bold"))
  }
  return(p)
}