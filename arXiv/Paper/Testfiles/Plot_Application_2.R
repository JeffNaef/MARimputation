
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
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "GAIN", "MIWAE", "mipca")
#methods <- c(  "GAIN", "MIWAE")

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




# Given the observed data build the data containing missing values
Beta<-matrix(runif(n=9, min=-2, max=2), nrow=3,ncol=3)

N<-1000

## Build the observed data

## As done in the paper ##
 Xobs1<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
 Xobs2<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
 Xobs3<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train


# ##With fixed covariance matrix ##
# ############################
# d<-3
# # params<-list(...)
# # mu<- params$mu
# # lambda <- params$lambda
# # C<-params$C
# # L<-params$L
# mu1=runif(d*(d-1)/2, min=-5, max=5)#2*runif(d)-1
# mu2=runif(d*(d-1)/2, min=-5, max=5)
# mu3=runif(d*(d-1)/2, min=-5, max=5)
# 
# L<-diag(d)
# 
# Xobs1<- MASS::mvrnorm(n = 10*N, mu = mu1, Sigma = L)
# Xobs2<- MASS::mvrnorm(n = 10*N, mu = mu2, Sigma = L)
# Xobs3<- MASS::mvrnorm(n = 10*N, mu = mu3, Sigma = L)
# 
# ############################



X.NA1 <- Xobs1%*%Beta + matrix( rnorm(n=3*10*N, sd=1), nrow=10*N, ncol= 3   )
X.NA2 <- Xobs2%*%Beta + matrix( rnorm(n=3*10*N, sd=1), nrow=10*N, ncol= 3   )
X.NA3 <- Xobs3%*%Beta + matrix( rnorm(n=3*10*N, sd=1), nrow=10*N, ncol= 3   )

#Fulldata<-cbind( rbind( X.NA1, X.NA2, X.NA3 ), rbind(Xobs1, Xobs2, Xobs3)    )

M1<-matrix(c(1,0,0), nrow=10*N, ncol=3, byrow=T)
M2<-matrix(c(0,1,0), nrow=10*N, ncol=3, byrow=T)
M3<-matrix(c(0,0,1), nrow=10*N, ncol=3, byrow=T)

#FullM<-cbind( rbind( M1, M2, M3 ), matrix(0, nrow=3*10000, ncol=3)    )



#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

s<-1
  
  set.seed(seeds[s])
  
  
  ###something is not right here!!! ####
  
  
  X<-cbind( rbind( X.NA1[ (N*(s-1)+1):(N*s),], 
                   X.NA2[ (N*(s-1)+1):(N*s),], X.NA3[ (N*(s-1)+1):(N*s),]), 
            rbind(Xobs1[ (N*(s-1)+1):(N*s),], Xobs2[ (N*(s-1)+1):(N*s),], Xobs3[ (N*(s-1)+1):(N*s),])    )
  M<-cbind( rbind( M1[ (N*(s-1)+1):(N*s),], M2[ (N*(s-1)+1):(N*s),], M3[ (N*(s-1)+1):(N*s),] ), 
            matrix(0, nrow=3*N, ncol=3)    )
  X.NA<-X
  X.NA[M==1] <- NA
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:6)
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  ## Deactivate standardization for MIWAE here!!!!
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  
  # imputationfuncs<-list()
  # for (method in methods){
  #   ##Probably not the smartes solution
  #   imputationfuncs[[method]][[1]] <- method
  #   imputationfuncs[[method]][[2]] <- function(X,m, method){ 
  #     doimputation(X.NA=X, methods=method, m=m,print=F, visitSequence="arabic", maxit = 1)}
  # }
  ################################## evaluations #########################################
  ########################################################################################
  
  #Step 1: Without access to true underlying data, check Iscore
  
  
  # start_time <- Sys.time()
  # new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf", imputationfuncs=imputationfuncs)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  # 
  # 
  # 
  # new.score.drf <- unlist(lapply(new.score.list.drf, function(x) x$score))
  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  
  
  
  escore<-rep(0, length(methods))
  names(escore)<-methods
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]
      
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      escore[method]<-escore[method]+eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )
      
    }
    escore[method] <- -1/m*escore[method]
  }
  
  #print("drf-score2:")
  #print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(energy.score=escore)
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  



## Analysis
  
sort(escore)  

pattern <- c(rep(1,N), rep(2,N), rep(3,N))  
par(mfrow=c(3,length(methods) + 2 ))

#M1<-matrix(c(1:6), nrow=3,ncol=2, byrow = TRUE)
#M2<-matrix(c(7:(3*length(methods))), nrow=3,ncol=length(methods)-2, byrow = FALSE)
#layout(cbind(M1,M2))

M<-matrix(c(1:(3*(2+length(methods)))), nrow=3,ncol=length(methods)+2, byrow = FALSE)
layout(M)






for (method in c("truth", "observed", methods[order(escore, decreasing=T)])){

  if (method=="truth"){
    plot(X[,1:2], col=pattern, main="truth")
    plot(X[,2:3], col=pattern, main="truth")
    plot(X[,c(1,3)], col=pattern,main="truth")
  }else if (method=="observed"){
    
    
    plot(X.NA[,1:2], col=pattern, main="observed")
    plot(X.NA[,2:3], col=pattern, main="observed")
    plot(X.NA[,c(1,3)], col=pattern, main="observed")
    
    
  }else{

plot(imputations[[method]][[1]][,1:2], col=pattern, main=method, xlab="", ylab="")
plot(imputations[[method]][[1]][,2:3], col=pattern, main=method, xlab="", ylab="")
plot(imputations[[method]][[1]][,c(1,3)], col=pattern, main=method, xlab="", ylab="")
  }
}














