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
library(scoringRules)

source("helpers.R")
#source("Iscores_new.R")



#install.packages("reticulate")
library(reticulate)


##Add drf!!
#methods <- c("DRF", "cart", "missForest")

#methods <- c("DRF","cart","norm.predict", "missForest", "norm.nob", "pmm")
## Add sample:
#methods <- c("DRF","cart","norm.predict", "missForest", "norm.nob", "mipca")
## Add MIPCA:
methods <- c("DRF","cart","norm.predict", "missForest", "norm.nob", "sample")


#methods <- c("pmm", "midastouch","mipca", "cart", "sample", "norm.predict",
#             "mean","rf", "missForest")
missing.mech <- "MAR"
pmiss <- 0.6
nrep.total<-10




## Maybe use spam??

#dataset <- "multivariateGaussian"
#data.type <- "simulated"
dataset <- "real_airquality"
#dataset <- "yeast"
#dataset <- "real_birthdata"
data.type <- "real"
## TO DO
## 1.Create Nice multivariate Gaussian dataset with 3-4 patterns with changing distribution!! (say 5 fully observed values and
## that change their distribution in each pattern + always the same conditional distribution)
## => Show how hard imputation can be + score proper + energy distance useful
## 2. Use real dataset with ampute MAR
## => Show score proper + energy distance useful
## 3. Use spam dataset and impute with DRF block + GAIN


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
 
for (s in 1:10){
 
  set.seed(seeds[s])
  
  
  #Resample data
  X <- genData_real(dataset = dataset, n = n)
  X <- as.matrix(X)
  colnames(X)<-NULL
  d <- ncol(X)
  n.train <- nrow(X)
  #X.NA <- genMask(X, mech = missing.mech, pmiss=pmiss)
  ##Focus on this!
  
  #patterns<-ampute.default.patterns(d )[1:5,]
  
  
  patterns<-matrix(1, nrow=4, ncol=d)
  patterns[1,1]<-0
  patterns[2,2]<-0
  patterns[3,3]<-0
  patterns[4,4]<-0
  
  #patterns[1,d-2]<-0
  #patterns[2,d-1]<-0
  #patterns[3,d]<-0
  
  X.NA<-ampute(X,
               prop = pmiss,
               patterns = patterns,
               freq = NULL,
               mech = missing.mech)$amp
  
  # Throw away observations which are only NA
  X <- X[rowSums(is.na(X.NA)) < d,]
  X.NA <- X.NA[rowSums(is.na(X.NA)) < d,]
  
  
  colnames(X)<-NULL
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
      doimputation(X.NA=X, methods=method, m=m,print=T, visitSequence="arabic", maxit = 1)}
  }
  
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
  #                               methods,
  #                               X.NA,
  #                               num.proj=1, rescale=F, projection.function = function(X){1:ncol(X)})
  
  new.score.list.drf <- Iscores(imputations,
                                methods,
                                X.NA,
                                num.proj=20, rescale=F, projection.function = NULL)
  
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
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
}




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

##Standardization:

# score in [min(all scores),max(all scores)]
# (score - max(all scores))/|min(all scores)-max(all scores)| in [-1,0]

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


# ### Version 2 ####
# ## Analysis
# energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
# energydata<-energydata[,!(colnames(energydata) %in% "sample")]
# energydata<--(energydata - max(energydata))/min(energydata)
# meanvalsenergy<- colMeans(energydata)
# ## Setup
# boxplot(energydata[,order(meanvalsenergy)], boxfill = NA, border = NA, ylim=c(-0.8,0),cex.axis=1.5,cex.lab=1.5) #invisible boxes - only axes and plot areas
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



filename = "Application_1_withsample_O"
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))





