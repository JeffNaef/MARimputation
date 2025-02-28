### Application Ready to Run #####


#### Application of Section 4.2. Air Quality Data ##

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
library(transport)

source("helpers.R")
#source("Iscores_new.R")





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


methods <- c("DRF", "cart","norm.predict", "missForest", "norm.nob", "sample", "GAIN", "MIWAE")


pmiss <- 0.3
nrep.total<-10




dataset <- "real_airquality"
data.type <- "real"

m <- 1
n<-2000



set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)

######################################################################################################
########################################### mask, imputation and evaluation ##########################
######################################################################################################

Results<-list()

#length(ids.jack)
#Results <- lapply(1:10, function(s){

for (s in 1:10){
  
  set.seed(seeds[s])

  
  #Resample data
  X <- genData_real(dataset = dataset, n = n)
  #X<-X[,1:6]
  X <- as.matrix(X)
  
  X2 <- genData_real(dataset = dataset, n = n)
# X2<-X2[,1:6]
  X2 <- as.matrix(X2)
  
  rownames(X)<-NULL
  rownames(X2)<-NULL
  
  d <- ncol(X)
  colnames(X)<-paste0("X", 1:d)
  n.train <- nrow(X)
  colnames(X2)<-paste0("X", 1:d)
  
  X.NA<-X
  
  X.NA[,1:6]<-apply(X[,1:6],2, function(x) {
    U<-runif(length(x))
    ifelse(U <= pmiss, rep(NA, length(x)), x)
  })
  
  sum(is.na(X.NA))/(n*d)
  
  # Throw away observations which are only NA
  X <- X[rowSums(is.na(X.NA)) < d,]
  X.NA <- X.NA[rowSums(is.na(X.NA)) < d,]
  
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:ncol(X))
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  
  
  escore<-rep(0, length(methods))
  ediff<-rep(0, length(methods))
  RMSE <- rep(0, length(methods))
  W2res<-rep(0, length(methods))
  names(escore)<-methods
  names(ediff)<-methods
  names(RMSE) <- methods
  names(W2res)<-methods
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-as.matrix(imputations[[method]][[j]])
      
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      
      Y1<-pp(Ximp)
      Y2<-pp(X2[1:nrow(Ximp),])
      
      W2res[method]<-wasserstein(Y1,Y2, p=2)
      
      ## Energy-score
      #escore[method]<-
      #  escore[method]+ 0.5*scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))- owndistance(X,Ximp)
      
      ediff[method]<-ediff[method]-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
      
      
      RMSE[method] <-
        RMSE[method] -  round(mean(apply(X - Ximp,1,function(x) norm(as.matrix(x), type="F"  ) )),2)
      
      
      
    }
    escore[method] <- 1/m*escore[method]
    ediff[method] <- 1/m*ediff[method]
    RMSE[method] <- 1/m*RMSE[method]
  }
  
  Y1full<-pp(X)
  Y2full<-pp(X2[1:nrow(X),])
  W2res["truth"]<-wasserstein(Y1full,Y2full, p=2)
  
  
  #print("drf-score2:")
  #print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(W2res=W2res,energy.score=escore, ediff.score=ediff, RMSE=RMSE)
  
  
  
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
}


Wassersteinmeasure<-lapply(Results, function(r) {
 # The higher the better
  diff<- - abs(r$W2res-r$W2res["truth"])
 
 diff<-diff[!(names(diff) %in% "truth")] 
 
} )


png(filename = "Application_2_EnergyDistance_RMSE_Wasserstein.png", 
    width = 1500,    # Width in pixels
    height = 800,    # Height in pixels
    res = 120)       # Resolution in dpi


## Analysis
energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$ediff.score))
energydata<-energydata[,!(colnames(energydata) %in% "sample")]
energydata<-(energydata - max(energydata))/abs(min(energydata)- max(energydata))
meanvalsenergy<- colMeans(energydata)
## Setup
boxplot(energydata[,order(meanvalsenergy)], boxfill = NA, border = NA, ylim=c(-1,0),cex.axis=1.5,cex.lab=1.5) #invisible boxes - only axes and plot area
##

boxplot(energydata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="white",
        boxwex=0.25, at = 1:ncol(energydata) + 0.30) #shift these left by -0.15

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$RMSE)))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="gray",
        boxwex=0.25, at = 1:ncol(scoredata) - 0) #shift to the right by +0.15


scoredata<-t(sapply(1:length(Results), function(j)  unlist(Wassersteinmeasure[[j]])))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="gray",
        boxwex=0.25, at = 1:ncol(scoredata) - 0.30) #shift to the right by +0.15

abline(v=1:(ncol(scoredata)-1)+0.50,lty = 2)


dev.off()


filename =paste0("Application_1_", paste0(methods, collapse="_"))
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))









