
### This is the Uniform Shift Example #####


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
require(kernlab)
print("kernlab loaded")
library(MASS)
print("MASS loaded")
library(AER)
print("AER loaded")
library(missForest)
print("missForest loaded")
library(Iscores)
print("IScores loaded")
library(scoringRules)
library(miceDRF)

source("helpers.R")
#source("Iscores_new.R")



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

## Add MIWAE here:
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "GAIN", "MIWAE")
nrep.total<-10
m<-1


set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)


n<-5000
d<-5




#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

for (s in 1:1){
  set.seed(seeds[s])
  
  # independent uniform
  X<-matrix(runif(n=d*n), nrow=n, ncol=d)
  # uniform with Gaussian copula
  #X <- gaussian_copula_uniform_sim(n = n, d = d)$uniform_data
  
  vectors <- matrix(c(
    rep(0, d),
    0, 1, rep(0,d-2),
    1, rep(0,d-1)
  ), nrow = 3, byrow = TRUE)
  
  
  # Generate random draws
  # sample() will generate indices, which we use to select rows from the matrix
  M <- vectors[apply(X,1, function(x) sample(1:3, size = 1, prob=c(x[1]/3, 2/3-x[1]/3, 1/3), replace = TRUE)), ]
  
  X.NA<-X
  X.NA[M==1]<-NA
  
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:d)
  
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  
  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  escore<-rep(0, length(methods))
  RMSE<-rep(0, length(methods))
  names(escore)<-methods
  names(RMSE)<-methods
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]

      colnames(Ximp)<-paste0("X",1:ncol(X))
      escore[method]<-escore[method]+eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
      #escore[method]<-
      #  escore[method]+ 0.5*scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))- owndistance(X,Ximp)
      
      
      RMSE[method] <-
        RMSE[method] -  round(mean(apply(X - Ximp,1,function(x) norm(as.matrix(x), type="F"  ) )),2)
      
    }
    escore[method] <- -1/m*escore[method]
    RMSE[method] <- 1/m*RMSE[method]
  }
  
  #print("drf-score2:")
  #print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(energy.score=escore, RMSE=RMSE)
  

    #Plotting:
    # Open a PNG device for saving the plot
    png(filename = paste0("Application_2_X1Imputation_", "s=",s, ".png"), 
        width = 1200,    # Width in pixels
        height = 400,    # Height in pixels
        res = 120)       # Resolution in dpi
    
    # Set up the plotting layout
    par(mfrow=c(1,3))
    
    # Create the three histograms
    hist(X[is.na(X.NA[,1]),1], 
         freq=FALSE, 
         main="Truth", 
         xlab="")
    
    hist(imputations[["GAIN"]][[1]][is.na(X.NA[,1]),1], 
         freq=FALSE, 
         main="GAIN", 
         xlab="")
    
    hist(imputations[["cart"]][[1]][is.na(X.NA[,1]),1], 
         freq=FALSE, 
         main="mice-cart", 
         xlab="")
    
    # Close the PNG device
    dev.off()
    
    
    

  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  
}





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


png(filename = "Application_2_EnergyDistance_RMSE.png", 
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

filename =paste0("Application_2_", paste0(methods, collapse="_"))

#filename ="Application_1_withGAINMIWAE"

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))












