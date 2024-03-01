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


##Add drf!!
#methods <- c("DRF", "cart", "missForest")

#methods <- c("DRF","cart","norm.predict", "missForest", "norm.nob", "pmm")
## Add MIPCA:
methods <- c("DRF","cart","norm.predict", "missForest", "norm.nob", "mipca")

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
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m)
  methods<-imputations$methods
  
  imputations <-imputations$imputations
  
  # append the truth, that is X
  
  # for ( i in 1:m){
  #   
  #   if (data.type == "simulated"){
  #     imputations[["truth"]][[i]] <- as.data.frame(X)
  #   }else{ imputations[["truth"]][[i]] <- as.data.frame(X)}
  #   
  #   names(imputations[["truth"]])[[i]] <- paste0(i)
  # }
  # 
  # methods <- c(methods, "truth")
  # names(imputations) <- methods
  
  imputationfuncs<-list()
  for (method in methods){
    ##Probably not the smartes solution
    imputationfuncs[[method]][[1]] <- method
    imputationfuncs[[method]][[2]] <- function(X,m, method){ 
      doimputation(X.NA=X, methods=method, m=m,print=FALSE, visitSequence="arabic", maxit = 1)}
  }
  ################################## evaluations #########################################
  ########################################################################################
  
  #Step 1: Without access to true underlying data, check Iscore
  
  
  start_time <- Sys.time()
  new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf2", imputationfuncs=imputationfuncs)
  end_time <- Sys.time()
  
  end_time-start_time
  
  start_time <- Sys.time()
  new.score.list.imp <- Iscores_new(X.NA,imputations,score="mulitpleimp2", imputationfuncs=imputationfuncs)
  end_time <- Sys.time()
  
  end_time-start_time
  
  
  new.score.drf <- unlist(lapply(new.score.list.drf, function(x) x$score))
  new.score.imp <- unlist(lapply(new.score.list.imp, function(x) x$score))
  
  
  #Step 2: With access to the full data, check energy score:
  # So far only for m=1!!!
  escore<-rep(0, length(methods))
  names(escore)<-methods
  for (method in methods){
    
    for (j in 1:m){
      
      Ximp<-imputations[[method]][[j]]
      escore[method]<-escore[method]+eqdist.e( rbind(X,Ximp), c(nrow(X), nrow(Ximp))  )
      
    }
    escore[method] <- -1/m*escore[method]
  }
  
#  print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
#  print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
#  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
 
  print("drf-score2:")
  print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
  print("m-score2:")
  print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
   
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore)
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
}





## Analysis





##Aspect Ratio: 1000 x 1000
par(mfrow=c(2,1))

scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$new.score.drf)))
scoredata<-scoredata[,!(colnames(scoredata)%in%"pmm")]
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsnewscore)], cex.axis=1.5)


scoredata<-t(sapply(1:length(Results), function(j)  unlist(Results[[j]]$new.score.imp)))
scoredata<-scoredata[,!(colnames(scoredata)%in%"pmm")]
meanvalsnewscore<- colMeans(scoredata)
boxplot(scoredata[,order(meanvalsnewscore)], cex.axis=1.5)


par(mfrow=c(1,1))
energydata<-t(sapply(1:length(Results), function(j) Results[[j]]$energy.score))
energydata<-energydata[,!(colnames(energydata)%in%"pmm")]
meanvalsenergy<- colMeans(energydata)
boxplot(energydata[,order(meanvalsenergy)], cex.axis=1.5)






filename = "Application_2_withMIPCA"
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))





