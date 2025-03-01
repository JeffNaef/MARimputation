
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

## Add MIWAE here:
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "sample", "GAIN", "MIWAE", "mipca", "mean")

nrep.total<-10




## Maybe use spam??

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




# Given the observed data build the data containing missing values
Beta<-matrix(runif(n=9, min=-2, max=2), nrow=3,ncol=3)

N<-1500

# ## Build the observed data
# Xobs1<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
# Xobs2<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
# Xobs3<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
# 
# X.NA1 <- Xobs1%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )
# X.NA2 <- Xobs2%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )
# X.NA3 <- Xobs3%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )
# 


Xstar<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=6)$train


## Add: 
pM6<-1
pM5<-1 #1/(1+exp(-( rowMeans(Xstar[,-5]) )))
pM4<-1 #1/(1+exp(-( rowMeans(Xstar[,-4]) )))
pM3<-1/(1+exp(-( rowMeans(Xstar[,-3]) )))
pM2<-1/(1+exp(-( rowMeans(Xstar[,-2]) )))
pM1<-1/(1+exp(-( rowMeans(Xstar[,-1]) )))

U6<-runif(n=length(pM3))
U5<-runif(n=length(pM2))
U4<-runif(n=length(pM1))
U3<-runif(n=length(pM3))
U2<-runif(n=length(pM2))
U1<-runif(n=length(pM1))

Mfull<-matrix(0, nrow=nrow(Xstar), ncol=ncol(Xstar))


case<-"MNAR"

if (case=="MNAR"){
  
  ## MNAR
  Mfull[,6]<-U6 < (1-pM6)
  Mfull[,5]<-U5 < (1-pM5)
  Mfull[,4]<-U4 < (1-pM4)
  Mfull[,3]<-U3 < (1-pM3)
  Mfull[,2]<-U2 < (1-pM2)
  Mfull[,1]<-U1 < (1-pM1)
}else{
  ## MCAR
  
  Mfull[,6]<-U6 < (1-pM6) #rep(mean(1-pM6), length(U6))
  Mfull[,5]<-U5 < (1-pM5) #rep(mean(1-pM5), length(U6))
  Mfull[,4]<-U4 < (1-pM4) #rep(mean(1-pM4), length(U6))
  Mfull[,3]<-U3 < rep(mean(1-pM3), length(U6))
  Mfull[,2]<-U2 < rep(mean(1-pM2), length(U6))
  Mfull[,1]<-U1 < rep(mean(1-pM1), length(U6))
  
  #Fulldata<-cbind( rbind( X.NA1, X.NA2, X.NA3 ), rbind(Xobs1, Xobs2, Xobs3)    )
}




#FullM<-cbind( rbind( M1, M2, M3 ), matrix(0, nrow=3*10000, ncol=3)    )



#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

for (s in 1:10){
  set.seed(seeds[s])
  
  
  
  X<-Xstar[(N*(s-1)+1):(N*s),]
  M<-Mfull[(N*(s-1)+1):(N*s),]
  X.NA<-X
  X.NA[M==1] <- NA
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:6)
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  
  imputationfuncs<-list()
  for (method in methods){
    ##Probably not the smartes solution
    imputationfuncs[[method]][[1]] <- method
    imputationfuncs[[method]][[2]] <- function(X,m, method){ 
      doimputation(X.NA=X, methods=method, m=m,print=F, visitSequence="arabic", maxit = 1)}
  }
  ################################## evaluations #########################################
  ########################################################################################
  
  #Step 1: Without access to true underlying data, check Iscore
  
  
  
  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  escore<-rep(0, length(methods))
  names(escore)<-methods
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      escore[method]<-escore[method]+eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
      #escore[method]<-
      #  escore[method]+ 0.5*scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))- owndistance(X,Ximp)
      
      
    }
    escore[method] <- -1/m*escore[method]
  }
  

  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list( energy.score=escore)
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  
}


## Analysis



par(mfrow=c(1,1))
energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
energydata<-energydata[,!(colnames(energydata) %in% "sample")]
meanvalsenergy<- colMeans(energydata)
boxplot(energydata[,order(meanvalsenergy)], cex.axis=1.5)







filename =paste0("Application_5_", paste0(methods, collapse="_"))

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))












