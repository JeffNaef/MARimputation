
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


## Add MIWAE here:
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "GAIN", "MIWAE")
#methods <- c("DRF", "cart", "missForest")


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


n<-5000
d<-5




#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

for (s in 1:nrep.total){
  set.seed(seeds[s])
  
  # independent uniform
  #X<-matrix(runif(n=d*n), nrow=n, ncol=d)
  # uniform with Gaussian copula
  # X <- gaussian_copula_uniform_sim(n = n, d = 2)$uniform_data
  X<-simulate_fgm(n=n, alpha=1)
  X<-cbind(X,matrix(runif( (d-2)*n ), nrow=n, ncol=d-2 ))
  
  vectors <- matrix(c(
    rep(0, d),
    0, 1, rep(0,d-2),
    1, rep(0,d-1)
  ), nrow = 3, byrow = TRUE)
  
  
  # Generate random draws
  # sample() will generate indices, which we use to select rows from the matrix
  M <- vectors[apply(X,1, function(x) sample(1:3, size = 1, prob=c((x[1]+x[2])/3, (2-x[1])/3, (1-x[2])/3), replace = TRUE)), ]
  
  X.NA<-X
  X.NA[M==1]<-NA
  
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:d)
  colnames(X.NA)<-paste0("X",1:d)
  
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  
  imputations[["truth"]][[1]]<-X.NA
  imputations[["truth"]][[1]][M[,2]==1, 2]<- sapply(imputations[["truth"]][[1]][M[,2]==1, 1], function(u1) sample_u2_given_u1(u1, alpha=1))
  imputations[["truth"]][[1]][M[,1]==1, 1]<- sapply(imputations[["truth"]][[1]][M[,1]==1, 2], function(u2) sample_u2_given_u1(u2, alpha=1))

  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  escore<-rep(0, length(methods)+1)
  RMSE<-rep(0, length(methods)+1)
  names(escore)<-c(methods,"truth")
  names(RMSE)<-c(methods,"truth")
  for (method in c(methods,"truth")){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      escore[method]<-escore[method]+eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)

      
      ### Do an actual test with independently generated X2
      #eqdist.etest(  rbind(X2, imputations[["DRF"]][[1]]), sizes=c(n,n) , R=50)
      
      
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


png(filename = "Application_Revision_EnergyDistance_RMSE.png", 
    width = 1800,    # Width in pixels
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






#######
#For plotting
#######

d<-2
n<-20000

#X <- gaussian_copula_uniform_sim(n = n, d = 2)$uniform_data
X<-simulate_fgm(n=n, alpha=1)


X<-cbind(X,matrix(runif( (d-2)*n ), nrow=n, ncol=d-2 ))

vectors <- matrix(c(
  rep(0, d),
  0, 1, rep(0,d-2),
  1, rep(0,d-1)
), nrow = 3, byrow = TRUE)


# Generate random draws
# sample() will generate indices, which we use to select rows from the matrix
M <- vectors[apply(X,1, function(x) sample(1:3, size = 1, prob=c((x[1]+x[2])/3, (2-x[1])/3, (1-x[2])/3), replace = TRUE)), ]

## both variables 
plot(X[,1:2], pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
     xlab = "U1", ylab = "U2",
     #main = paste0("FGM Copula (alpha = ", alpha, ")")
)

## both variabels in the fully observed pattern:
#plot(X[M[,1]==0& M[,2]==0,1:2])

plot(X[M[,1]==0& M[,2]==0,1:2], pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3),
           xlab = "U1", ylab = "U2",
           #main = paste0("FGM Copula (alpha = ", alpha, ")")
     )

# X_1 is observed
plot(X[M[,1]==0,1:2])
# X_2 is observed
plot(X[M[,2]==0,1:2])



### Quantile Estimation


methods <- c( "DRF", "cart", "missForest")

set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)


n<-5000
d<-5
alpha<-0.1
nrep.total<-20

Resultsquantile<-list()

for (s in 1:nrep.total){
  set.seed(seeds[s])
  
  # independent uniform
  #X<-matrix(runif(n=d*n), nrow=n, ncol=d)
  # uniform with Gaussian copula
  # X <- gaussian_copula_uniform_sim(n = n, d = 2)$uniform_data
  X<-simulate_fgm(n=n, alpha=1)
  X<-cbind(X,matrix(runif( (d-2)*n ), nrow=n, ncol=d-2 ))
  
  vectors <- matrix(c(
    rep(0, d),
    0, 1, rep(0,d-2),
    1, rep(0,d-1)
  ), nrow = 3, byrow = TRUE)
  
  
  # Generate random draws
  # sample() will generate indices, which we use to select rows from the matrix
  M <- vectors[apply(X,1, function(x) sample(1:3, size = 1, prob=c((x[1]+x[2])/3, (2-x[1])/3, (1-x[2])/3), replace = TRUE)), ]
  
  X.NA<-X
  X.NA[M==1]<-NA
  
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:d)
  colnames(X.NA)<-paste0("X",1:d)
  
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations

  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  quantile<-rep(0, length(methods)+1)
  names(quantile)<-c(methods,"observedonly")
  for (method in c(methods)){
    
      
      Ximp<-imputations[[method]][[1]]
      
      colnames(Ximp)<-paste0("X",1:ncol(X))
      quantile[method]<-quantile(Ximp[,1], probs=alpha)
    
  }
  
  quantile["observedonly"]<-quantile(X.NA[!is.na(X.NA[,1]),1], probs=alpha)
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Resultsquantile[[s]] <- quantile
  
  
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
  
}


##Quantile of X_1 \mid M_1=0
-7 + sqrt(49+15*alpha)




png(filename = "Application_Revision_quantile.png", 
    width = 1700,    # Width in pixels
    height = 800,    # Height in pixels
    res = 120)       # Resolution in dpi


par(mfrow=c(1,1))
## Setup
quantiledata<-t(sapply(1:length(Resultsquantile), function(j) Resultsquantile[[j]]))
quantiledatamtruth<-abs(quantiledata-alpha)

meanvalsquantiles<-colMeans(quantiledatamtruth)

boxplot(quantiledata[,order(meanvalsquantiles, decreasing = T)],,cex.axis=1.5,cex.lab=1.5)
abline(h=alpha)

# Close the PNG device
dev.off()




filename =paste0("Application_Revision_", paste0(methods, collapse="_"))

#filename ="Application_1_withGAINMIWAE"

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))





