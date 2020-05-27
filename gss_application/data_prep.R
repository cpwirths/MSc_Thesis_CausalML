rm(list=ls())
library(tidyverse)
library(foreign)

#prepares the GSS raw data for the estimation, download rawwdata from https://gss.norc.org/

dir_path <- "~/Documents/GitHub/MSc_Thesis_CausalML/gss_application" 
setwd(dir_path)
data <- read.dta("GSS7218_R1.DTA")
data$welfare <- data$natfare
data$assistance <- data$natfarey 

data <- data[, c("welfare","assistance", "partyid", "polviews", "age", "educ", "year", "RACDIF1", "RACDIF2", "RACDIF3", "RACDIF4")]
data <- data[rowSums(is.na(data[, c("welfare", "assistance")])) != ncol(data[, c("welfare", "assistance")]), ]
data$treatment <- ifelse(is.na(data$assistance), 1, 0)
data$treatment <- as.factor(ifelse(is.na(data$assistance), "welfare", "assistance"))
data$response <- rep(0, nrow(data))
data$response[which(data$welfare== "TOO MUCH" | data$assistance== "TOO MUCH")] <- 1
data <- data[,!colnames(data)%in%c("welfare", "assistance")]
data$treatment <- ifelse(data$treatment == "welfare", 1, 0)
data <- data[data$year!=1977,] 
data <- na.omit(data)
data$partyid <- droplevels(data$partyid)
data$polviews <- droplevels(data$polviews)
data$RACDIF1 <- droplevels(data$RACDIF1)
data$RACDIF2 <- droplevels(data$RACDIF2)
data$RACDIF3 <- droplevels(data$RACDIF3)
data$RACDIF4 <- droplevels(data$RACDIF4)
summary(data)

save(data, file="GSS_input.RData")

