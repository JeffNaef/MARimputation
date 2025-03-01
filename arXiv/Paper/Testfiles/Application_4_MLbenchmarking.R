### Application Ready to Run #####


##Need a different dataset here.

# Potential datets:
#yeast

#######################################
###One real dataset with ampute MAR ###
#######################################

###Need to find a solution for binary variables!!!

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



# data sets we use in the paper
# real data sets: dataset = ...
# "airfoil","Boston","concrete.compression","connectionist.bench.vowel",
# "yacht","climate.model.crashes", "CASchools", "ecoli","wine","yeast", "ionosphere"
# "iris","concrete.slump","seeds","planning.relax",


parallelize <- TRUE
num.trees.per.proj <- 5 # need to be at least 2
min.node.size <- 10
m <- 1
n<-2000
real.data <- TRUE
n.cores <- 10
frac.subsample<-0.9


set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)

######################################################################################################
########################################### mask, imputation and evaluation ##########################
######################################################################################################

Results<-list()

#length(ids.jack)
#Results <- lapply(1:10, function(s){

for (s in 9:10){
  
  
  ### For s=9, the DRF imputation breaks down !!####
  ### Need to find the reason###
  
  set.seed(seeds[s])
  
  
  #Resample data
  X <- genData_real(dataset = "real_spam", n = n)
  X <- as.matrix(X)
  colnames(X)<-NULL
  d <- ncol(X)
  n.train <- nrow(X)
  
  colnames(X) <- paste0("X",1:ncol(X))
  
  #X.NA <- genMask(X, mech = missing.mech, pmiss=pmiss)
  ##Focus on this!
  
  #patterns<-ampute.default.patterns(d )[1:5,]
  
  
  
  
  #patterns[1,d-2]<-0
  #patterns[2,d-1]<-0
  #patterns[3,d]<-0
  
  M<-apply(X,2, function(x) sample(c(0,1), size=length(x), replace=T, prob = c(1-0.2,0.2) )  )
  X.NA<-X
  X.NA[M==1] <- NA
  
  X <- X[rowSums(is.na(X.NA)) < d,]
  X.NA <- X.NA[rowSums(is.na(X.NA)) < d,]
  
  #colnames(X)<-NULL
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  
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
  colnames(X)<-paste0("X",1:ncol(X))
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
  
  
  Results[[s]] <- list(energy.score=escore)
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
}





## Analysis


par(mfrow=c(1,1))
energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
energydata<-energydata[,!(colnames(energydata) %in% "sample")]
meanvalsenergy<- colMeans(energydata)
boxplot(energydata[,order(meanvalsenergy)], cex.axis=1.5)




filename =paste0("Application_4_", paste0(methods, collapse="_"))

#filename = "Application_4_withGAINMIWAE"
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))





