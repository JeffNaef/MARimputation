
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

source("helpers.R")
#source("Iscores_new.R")


##Add drf!!
#methods <- c("DRF", "cart","norm.predict", "missForest", "norm.nob")
## Add MIPCA:
#methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob", "mipca")
## Add mice-sample
#methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob")

methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob")


#methods <- c("norm.predict", "norm.nob")

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
#Beta<-matrix(runif(n=9, min=-2, max=2), nrow=3,ncol=3)
d<-3
Beta<-matrix(c(0.5,1,1.5), nrow=d,ncol=d, byrow=T)


N<-500

C<-matrix(0, nrow=d, ncol=d)
for (t in 1:d){
  for (s in 1:d){
    
    
    C[t,s] <- 0.5^(abs(t-s))
    
  }
  
  
}

## Build the observed data
Xobs1<- MASS::mvrnorm(n = 10*N, mu = rep(3, d), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs2<- MASS::mvrnorm(n = 10*N, mu = rep(-3, d), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train
Xobs3<- MASS::mvrnorm(n = 10*N, mu = rep(0, d), Sigma = C)#genDataNoNA_synthetic(dataset = dataset, n.train = 10*N, d=3)$train


X.NA1 <- Xobs1%*%Beta + matrix( rnorm(n=3*10*N,sd=2), nrow=10*N, ncol= 3   )
X.NA2 <- Xobs2%*%Beta + matrix( rnorm(n=3*10*N,sd=2), nrow=10*N, ncol= 3   )
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
  
  
  # start_time <- Sys.time()
  # new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf2", imputationfuncs=imputationfuncs)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  # 
  # start_time <- Sys.time()
  # new.score.list.imp <- Iscores_new(X.NA,imputations,score="mulitpleimp2", imputationfuncs=imputationfuncs)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  
  # start_time <- Sys.time()
  # new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf", imputationfuncs=imputationfuncs)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  
  
  # new.score.list.drf <- Iscores(imputations,
  #                      methods,
  #                      X.NA,
  #                      num.proj=1, rescale=F, projection.function = function(X){1:ncol(X)})
  
  
  new.score.list.drf <- Iscores(imputations,
                                methods,
                                X.NA,
                                num.proj=20, rescale=F)
  
  start_time <- Sys.time()
  new.score.list.imp <- Iscores_new(X.NA,imputations=imputations, imputationfuncs=imputationfuncs)
  end_time <- Sys.time()
  
  end_time-start_time
  
  
 DRIScore <- unlist(new.score.list.drf)
 names(DRIScore) <- colnames(new.score.list.drf)
 new.score.imp <- unlist(lapply(new.score.list.imp, function(x) x$score))

    #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  escore<-rep(0, length(methods))
  ediff<-rep(0, length(methods))
  RMSE <- rep(0, length(methods))
  names(escore)<-methods
  names(ediff)<-methods
  names(RMSE) <- methods
  for (method in methods){
    
    for (j in 1:m){
    
    Ximp<-as.matrix(imputations[[method]][[j]])
    
    
    colnames(Ximp)<-paste0("X",1:ncol(X))
   
    ## Energy-score
    escore[method]<-
      escore[method]+ 0.5*scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))- owndistance(X,Ximp)

    ediff[method]<-ediff[method]-eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )*(2*n)/(n^2)
    
    
    RMSE[method] <-
      RMSE[method] -  round(mean(apply(X - Ximp,1,function(x) norm(as.matrix(x), type="F"  ) )),2)

    
  
    }
    escore[method] <- 1/m*escore[method]
    ediff[method] <- 1/m*ediff[method]
    RMSE[method] <- 1/m*RMSE[method]
  }

  # print("drf-score2:")
  # print( sort( round( unlist(DRIScore)/sum(unlist(DRIScore)),3) , decreasing=T)   )
  # print("m-score2:")
  # print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
  # print("e-score")
  # print( sort( round(escore/sum(escore),3) , decreasing=T)   )

  print(paste0("nrep ",s, " out of ", nrep.total ))

  Results[[s]] <- list(new.score.imp = new.score.imp, DRIScore=DRIScore , energy.score=escore, ediff.score=ediff, RMSE=RMSE)
  
  
  #return(list(new.score.imp = new.score.imp,DRIScore=DRIScore , energy.score=escore))


}

##Standardize




### Version 1 ####

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
        boxwex=0.25, at = 1:ncol(scoredata)) #shift to the right by +0.15

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$new.score.imp)))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="darkgray",
        boxwex=0.25, at = 1:ncol(scoredata) - 0.30) #shift to the right by +0.15

abline(v=1:(ncol(scoredata)-1)+0.50,lty = 2)




### Version with DR-I-Score ####

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

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$DRIScore)))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="gray",
        boxwex=0.25, at = 1:ncol(scoredata)) #shift to the right by +0.15

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$new.score.imp)))
scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
scoredata<-(scoredata - max(scoredata))/abs(min(scoredata)- max(scoredata))
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="darkgray",
        boxwex=0.25, at = 1:ncol(scoredata) - 0.30) #shift to the right by +0.15

abline(v=1:(ncol(scoredata)-1)+0.50,lty = 2)

# 
# ### Version 2 ####
# ## Analysis
# energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
# energydata<-energydata[,!(colnames(energydata) %in% "sample")]
# energydata<--(energydata - max(energydata))/min(energydata)
# meanvalsenergy<- colMeans(energydata)
# ## Setup
# boxplot(energydata[,order(meanvalsenergy)], boxfill = NA, border = NA, ylim=c(-1,0),cex.axis=1.5,cex.lab=1.5) #invisible boxes - only axes and plot areas
# ##
# 
# scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$RMSE)))
# scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
# scoredata<--(scoredata - max(scoredata))/min(scoredata)
# meanvalsnewscore<- colMeans(scoredata)
# boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="white", 
#         boxwex=0.25, at = 1:ncol(scoredata)-0.15) #shift to the right by +0.15
# 
# scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$new.score.imp)))
# scoredata<-scoredata[,!(colnames(scoredata) %in% "sample")]
# scoredata<--(scoredata - max(scoredata))/min(scoredata)
# meanvalsnewscore<- colMeans(scoredata)
# boxplot(scoredata[,order(meanvalsenergy)],ylab="",yaxt="n", xaxt = "n", add = TRUE, boxfill="darkgray", 
#         boxwex=0.25, at = 1:ncol(scoredata) + 0.15) #shift to the right by +0.15
# 
# abline(v=1:(ncol(scoredata)-1)+0.50,lty = 2)
# 


filename ="Application_2_withsample_O"

assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))












