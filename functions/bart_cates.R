#functions to compute BART MChains CATE estimates; requires an array of the individual treatment effects as input with the dimension number of draws, number of chains, number of observations; 
#CATEs are based on every distinct value in a column; e.g if the column represents year 2010-2018, the function plots the CATE for each year and the corresponding confidence interval

library(ggplot2)
#
# Function returning point estimates and confidence interval --------------
bart_cate <- function(array, index, alpha=0.05){
  effect <- matrix(NA, dim(array)[2],dim(array)[3])
  for(i in 1:dim(array)[2]){
    for(j in 1:dim(array)[3]){
      effect[i,j] <- mean(array[index,i,j])  
    }
  }
  est <- mean(as.vector(effect))
  sd_bart <- sd(as.vector(effect))
  lower <- quantile(as.vector(effect), alpha)
  upper <- quantile(as.vector(effect), 1-alpha)
  output <- c(est, lower,upper)
  names(output) <- c("est", "lower", "upper")
  return(output)  
}


# Function plotting point estimates and confidence interval ---------------
plot_cate_bart <-function(array, column, numeric=T, x_name="x_name", y_name="y_name", title="title", textsize=12, angle=90, exclude=NULL){
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
  matrix <- matrix(NA,length(vec), 3)
  for(i in vec){
    index <- which(i==vec)
    matrix[index,1:3] <-  bart_cate(array, column==i)
  }
  matrix <- data.frame(matrix)
  colnames(matrix) <- c("est", "lb", "ub") 
  if(numeric==T){
    x_axis <- vec  
  }else{
    x_axis <- factor(vec, levels=vec) 
  }
  if(numeric){
    p <- ggplot(matrix, aes(x=x_axis, y=est))+
      geom_line(data=matrix,group=1)+geom_point()+
      geom_ribbon(data=matrix,aes(ymin=lb,ymax=ub), alpha=0.3,group=1)+theme_bw()+ggtitle(title)+  theme(axis.text.x = element_text(angle = angle, size=textsize, vjust = 0.5)) + scale_x_continuous(x_name, labels = as.character(x_axis), breaks = x_axis)+ylab(y_name)+ xlab(x_name)+  theme(axis.text.y = element_text(size=textsize))+theme(plot.title = element_text(size = textsize, face = "bold"))+theme(axis.title = element_text(size = textsize, face = "bold"))+xlab("xlab")
  }else{p <- ggplot(matrix, aes(x=x_axis, y=est))+
    geom_line(data=matrix,group=1)+geom_point()+
    geom_ribbon(data=matrix,aes(ymin=lb,ymax=ub), alpha=0.3,group=1)+theme_bw()+ggtitle(title)+  theme(axis.text.x = element_text(angle = angle, size=textsize, vjust = 0.5)) +ylab(y_name)+ xlab(x_name) + theme(axis.text.y = element_text(size=textsize))+theme(plot.title = element_text(size = textsize, face = "bold"))+theme(axis.title = element_text(size = textsize, face = "bold"))
  }
  return(p)
}

