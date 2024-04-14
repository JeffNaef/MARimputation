
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
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "sample", "GAIN", "MIWAE")


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

## Add MIWAE here:
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "sample", "GAIN", "MIWAE")


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

N<-500

## Build the observed data
Xobs1<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs2<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs3<- genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train

Xobs1tranformed<-t(apply(Xobs1,1,function(x){ 
  c(x[3]*sin(x[1]*x[2]), x[2]*(x[2] > 0), atan(x[1])*atan(x[2]) ) }  ))
Xobs2tranformed<-t(apply(Xobs2,1,function(x){ 
  c(x[3]*sin(x[1]*x[2]), x[2]*(x[2] > 0), atan(x[1])*atan(x[2]) ) }  ))
Xobs3tranformed<-t(apply(Xobs3,1,function(x){ 
  c(x[3]*sin(x[1]*x[2]),x[2]*(x[2] > 0), atan(x[1])*atan(x[2]) ) }  ))

# matrix(Xobs1, nrow=nrow(Xobs1), )
X.NA1 <- Xobs1tranformed%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )
X.NA2 <- Xobs2tranformed%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )
X.NA3 <- Xobs3tranformed%*%Beta + matrix( rnorm(n=3*10*N), nrow=10*N, ncol= 3   )

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
  # ################################## evaluations #########################################
  # ########################################################################################
  # 
  # #Step 1: Without access to true underlying data, check Iscore
  # 
  # 
  # start_time <- Sys.time()
  # new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf", imputationfuncs=imputationfuncs)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  #new.score.drf <- unlist(lapply(new.score.list.drf, function(x) x$score))
  
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
  
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(energy.score=escore)
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  
}


## Analysis


par(mfrow=c(1,1))
energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
energydata<-energydata[,!(colnames(energydata) %in% "sample")]
meanvalsenergy<- colMeans(energydata)
boxplot(energydata[,order(meanvalsenergy)], cex.axis=1.5)







filename ="Application_3_withGAINMIWAE"

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))












