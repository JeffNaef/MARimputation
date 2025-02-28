
### This is the Gaussian Shift Example #####


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
library(Iscores)
print("IScores loaded")
library(scoringRules)
library(miceDRF)

source("helpers.R")
#source("Iscores_new.R")


#install.packages("reticulate")
library(reticulate)

## Better code
#Sys.setenv("gain_env" =  path.expand("~/anaconda3/envs/gain_env"))
#use_python(path.expand("~/opt/anaconda3/envs/gain_env/bin/python"))
#use_condaenv(path.expand("~/opt/anaconda3/envs/gain_env"))

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


methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "GAIN", "MIWAE")


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



nrep.total<-10



dataset <- "multivariateGaussian"

## TO DO
## 1.Create Nice multivariate Gaussian dataset with 3-4 patterns with changing distribution!! (say 5 fully observed values and
## that change their distribution in each pattern + always the same conditional distribution)
## => Show how hard imputation can be + score proper + energy distance useful
## 2. Use real dataset with ampute MAR
## => Show score proper + energy distance useful
## 3. Use spam dataset and impute with DRF block + GAIN
m <- 1


set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)



d<-3
Beta<-diag(d)#matrix(1, nrow=d,ncol=d)


N<-500

C<-matrix(0, nrow=d, ncol=d)
for (i in 1:3){
  for (j in 1:3){
    
    
    C[i,j] <- 0.5^(abs(i-j))
    
  }
  
  
}

## Build the observed data
Xobs1<- MASS::mvrnorm(n = 10*N, mu = rep(5, d), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs2<- MASS::mvrnorm(n = 10*N, mu = rep(-5, d), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs3<- MASS::mvrnorm(n = 10*N, mu = rep(0, d), Sigma = C)#genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train

# matrix(Xobs1, nrow=nrow(Xobs1), )
X.NA1 <- Xobs1%*%Beta+ matrix( rnorm(n=3*10*N,sd=2), nrow=10*N, ncol= 3   )
X.NA2 <- Xobs2%*%Beta+ matrix( rnorm(n=3*10*N,sd=2), nrow=10*N, ncol= 3   )
X.NA3 <- Xobs3%*%Beta + matrix( rnorm(n=3*10*N,sd=2), nrow=10*N, ncol= 3   )

#Fulldata<-cbind( rbind( X.NA1, X.NA2, X.NA3 ), rbind(Xobs1, Xobs2, Xobs3)    )

M1<-matrix(c(1,0,0), nrow=10*N, ncol=3, byrow=T)
M2<-matrix(c(0,1,0), nrow=10*N, ncol=3, byrow=T)
M3<-matrix(c(0,0,1), nrow=10*N, ncol=3, byrow=T)

#FullM<-cbind( rbind( M1, M2, M3 ), matrix(0, nrow=3*10000, ncol=3)    )

#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

for (s in 1:10){
  set.seed(seeds[s])
  
  
  
  
  X<-cbind( rbind( X.NA1[ (N*(s-1)+1):(N*s),], 
                   X.NA2[ (N*(s-1)+1):(N*s),], X.NA3[ (N*(s-1)+1):(N*s),]), 
            rbind(Xobs1[ (N*(s-1)+1):(N*s),], Xobs2[ (N*(s-1)+1):(N*s),], Xobs3[ (N*(s-1)+1):(N*s),])    )
  M<-cbind( rbind( M1[ (N*(s-1)+1):(N*s),], M2[ (N*(s-1)+1):(N*s),], M3[ (N*(s-1)+1):(N*s),] ), 
            matrix(0, nrow=3*N, ncol=3)    )
  X.NA<-X
  X.NA[M==1] <- NA
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:6)
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  ## Deactivate standardization for MIWAE here!!!!
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
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
  
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(energy.score=escore, RMSE=RMSE)
  
  
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


png(filename = "Application_3_EnergyDistance_RMSE.png", 
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












