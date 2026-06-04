
### This is the Nonlinear Shift Example #####


require(energy)
require(mice)
print("mice loaded")
require(data.table)
print("data.table loaded")
require(ranger)
print("ranger loaded")
require(parallel)
print("parallel loaded")
require(missMDA)
print("missMDA loaded")
library(tidyverse)
print("tidyverse loaded")
library(softImpute)
print("softImpute loaded")
library(Matrix)
print("Matrix loaded")
library(readxl)
print("readxl loaded")
library(plotrix)
print("plotrix loaded")
library(missForest)
print("missForest loaded")
#library(naniar)
#print("naniar loaded")
require(kernlab)
print("kernlab loaded")
library(MASS)
print("MASS loaded")
library(AER)
print("AER loaded")
library(missForest)
print("missForest loaded")
library(scoringRules)
library(miceDRF)


source("helpers.R")
#source("Iscores_new.R")



methods <- c("cart") # "rf" "cart", "GAIN", "MIWAE"


if (any(c("GAIN", "MIWAE") %in% methods)){
  
  #install.packages("reticulate")
  library(reticulate)
  

##Laptop
Sys.setenv("gain_env" =  "C:/Users/jeffr/anaconda3/envs/gain_env")
#use_python("C:/Users/jeffr/anaconda3/envs/gain_env/bin/python")
use_condaenv("C:/Users/jeffr/anaconda3/envs/gain_env")

## Write 
# conda activate gain_env
## in the terminal!

py_config()

###MIWAE method
torch <- import("torch") 
torchvision <- import("torchvision")
numpy <- import("numpy")
scipy <- import("scipy")
pandas <- import("pandas")
sklearn<- import("sklearn")

###GAIN method
#Required Python Packages
tensorflow <- import("tensorflow")
numpy <- import("numpy")
tqdm <- import("tqdm")
keras <- import("keras")
argparse <- import("argparse")  #pip install argparse
sys<- import("sys")



reticulate::source_python("gain.py") #there will be  warning but don't worry 
reticulate::source_python("MIWAE_Pytorch.py") #there will be  warning but don't worry




py_config()

#Required Python Packages
tensorflow <- import("tensorflow")
numpy <- import("numpy")
tqdm <- import("tqdm")
keras <- import("keras")
argparse <- import("argparse")  #pip install argparse
sys<- import("sys")
torch <- import("torch") 
numpy <- import("numpy")
scipy <- import("scipy")
pandas <- import("pandas")
sklearn<- import("sklearn")
torchvision <- import("torchvision")

reticulate::source_python("gain.py") #there will be  warning but don't worry 
reticulate::source_python("MIWAE_Pytorch.py") #there will be  warning but don't worry
}


nrep.total<-10
m<-1


set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)




Results<-list()


for (s in 1:10){
  set.seed(seeds[s])
  
  
  
  X0 <- readRDS(paste0("datasets/", "scm1d", ".RDS"))
  N <- nrow(X0)
  
  
  ###Only choose numerical variables ####
  #X0 <- X0[, apply(X0, 2, function(x)
  #  length(unique(x))) / N > 0.1 , drop = F]
  
  #if (ncol(X0) <= 5) {
  # next
  #}
  
  n <- round(N / 2) - 1
  indextrain <- sample(1:N, size = n , replace = F)
  
  Xtrain <- X0[indextrain, ]
  
  
  d <- ncol(X0)
  npattern <- sample(3:max(round(n / 100), 3), size = 1)
  
  patterns <- matrix(
    sample(c(0, 1), size = npattern * (d - 1), replace = T),
    nrow = npattern,
    ncol = d - 1,
    byrow = T
  )
  patterns <- cbind(patterns, rep(1, nrow(patterns)))
  
  ##add fully observed pattern
  if (all(rowSums(patterns) < d)) {
    patterns <- rbind(patterns, rep(1, d))
    
  }
  
  
  # Separate the all-ones (complete) pattern
  all_ones <- apply(patterns, 1, function(r)
    all(r == 1))
  incomplete_patterns <- patterns[!all_ones, , drop = FALSE]
  n_pats <- nrow(incomplete_patterns)
  
  n_total <- nrow(Xtrain)
  n_complete <- round(n_total / (n_pats + 1))  # ~1/21 of rows stay complete
  
  # Randomly assign rows to the "complete" pattern
  complete_idx <- sample(n_total, n_complete)
  Xtrain_to_ampute <- Xtrain[-complete_idx, , drop = FALSE]
  
  # Ampute the remaining rows with equal frequency across incomplete patterns
  tmp <- ampute(
    Xtrain_to_ampute,
    patterns  = incomplete_patterns,
    freq      = rep(1 / n_pats, n_pats),
    prop      = 0.99,
    # ampute all but 1/(n_pats+1) fraction
    mech      = "MAR",
    bycases   = TRUE
  )
  
  # Reassemble, preserving original row order
  X.NA <- rbind(Xtrain[complete_idx, , drop = FALSE], tmp$amp)
  X.NA <- X.NA[order(c(complete_idx, seq_len(nrow(Xtrain))[-complete_idx])), , drop = FALSE]
  
  M <- is.na(X.NA) * 1
  colnames(X.NA) <- paste0("X", 1:ncol(Xtrain))
  
  
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  ## Deactivate standardization for MIWAE here!!!!
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=1)
  #methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  escore<-rep(0, length(methods))
  RMSE<-rep(0, length(methods))
  names(escore)<-methods
  names(RMSE)<-methods
  

  
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]
      
      colnames(Ximp)<-paste0("X",1:ncol(Xtrain))
      escore[method]<-escore[method]+eqdist.e( rbind(Xtrain,Ximp), c(nrow(Xtrain), nrow(Ximp))  )*(2*n)/(n^2)
      #escore[method]<-
      #  escore[method]+ 0.5*scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))- owndistance(X,Ximp)
      
      
      RMSE[method] <-
        RMSE[method] -  round(mean(apply(Xtrain - Ximp,1,function(x) norm(as.matrix(x), type="F"  ) )),2)
      
    }
    escore[method] <- -1/m*escore[method]
    RMSE[method] <- 1/m*RMSE[method]
  }
  
  print("e-score")
  print(sort( round(escore,3) , decreasing=T))
  #print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(energy.score=escore, RMSE=RMSE)
  
  saveRDS(Results, file = paste0("results_", paste0(methods, collapse="_"), ".rds"))
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  
}


## Analysis



energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
energydata<-energydata[,!(colnames(energydata) %in% "sample")]


## Standardize
## Analysis
energydata<-energydata[,!(colnames(energydata) %in% "sample")]
energydata<-(energydata - max(energydata))/abs(min(energydata)- max(energydata))
meanvalsenergy<- colMeans(energydata)

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$RMSE)))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))


png(filename = "Application_5_EnergyDistance_RMSE.png", 
    width = 1700,    # Width in pixels
    height = 800,    # Height in pixels
    res = 120)       # Resolution in dpi


par(mfrow=c(1,1))
## Setup
boxplot(energydata[,order(meanvalsenergy)], boxfill = NA, border = NA, ylim=c(-1,0),cex.axis=1.5,cex.lab=1.5) #invisible boxes - only axes and plot area
##

boxplot(energydata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="white",
        boxwex=0.25, at = 1:ncol(energydata) + 0.15) #shift these left by -0.15


meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="gray",
        boxwex=0.25, at = 1:ncol(scoredata) - 0.15) #shift to the right by +0.15

abline(v=1:(ncol(scoredata)-1)+0.50,lty = 2)

# Close the PNG device
dev.off()





filename =paste0("Application_3_", paste0(methods, collapse="_"))

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))












