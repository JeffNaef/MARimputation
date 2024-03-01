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

##TO DO:
##- Try to include GAIN into doimputation!
##- Fix DRF imputation



##Add drf!!
#methods <- c("DRF", "cart", "missForest")

#methods <- c("DRF", "cart")
#methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob")
## Add MIPCA:
methods <- c( "DRF", "cart","norm.predict", "missForest", "norm.nob")



m<-1




# data sets we use in the paper
# real data sets: dataset = ...
# "airfoil","Boston","concrete.compression","connectionist.bench.vowel",
# "yacht","climate.model.crashes", "CASchools", "ecoli","wine","yeast", "ionosphere"
# "iris","concrete.slump","seeds","planning.relax",




set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)

######################################################################################################
########################################### mask, imputation and evaluation ##########################
######################################################################################################

Results<-list()

  
  
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
  

  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add drf
  
  imputations <- doimputation(X.NA=X.NA, methods=methods, m=m) # , blocksize=2
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
      doimputation(X.NA=X, methods=method, m=m,print=F, visitSequence="arabic", maxit = 1)}
  }
  ################################## evaluations #########################################
  ########################################################################################
  
  # #Step 1: Without access to true underlying data, check Iscore
  # start_time <- Sys.time()
  # new.score.list.drf <- Iscores_new(X.NA,imputations,score="drf2", imputationfuncs=imputationfuncs, projections=T, projectionsize=10)
  # end_time <- Sys.time()
  # 
  # end_time-start_time
  
  start_time <- Sys.time()
  new.score.list.imp <- Iscores_new(X.NA,imputations,score="mulitpleimp2", imputationfuncs=imputationfuncs, projections=T, projectionsize=10, maxlength = 10)
  end_time <- Sys.time()
  
  end_time-start_time
  
  
  
  #new.score.drf <- unlist(lapply(new.score.list.drf, function(x) x$score))
  new.score.imp <- unlist(lapply(new.score.list.imp, function(x) x$score))
  
  
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
  
  
  #  print( sort( round( unlist(new.score.drf)/sum(unlist(new.score.drf)),3) , decreasing=T)   )
  #  print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
  #  print( sort( round(escore/sum(escore),3) , decreasing=T)   )

  print("m-score2:")
  print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
  print("e-score")
  print( sort( round(escore/sum(escore),3) , decreasing=T)   )

  
  Results <- list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore)
  
  #return(list(new.score.imp = new.score.imp,new.score.drf=new.score.drf , energy.score=escore))
  
#}






filename = "Application_3"
assign(filename, Results)
save(Results, file=paste(filename, ".Rda",sep=""))



RMSE<-list()
## RMSE 
for (method in methods){
  
  
  RMSE[[method]]<-norm(as.matrix(imputations[[method]][[1]] - X) , type="F")
  
}
sort(unlist(RMSE), decreasing=T)
res3<-unlist(RMSE)



res1<-sort(unlist(new.score.imp), decreasing=T )
res2<-escore[order(-res1)]

res2<-sort(escore, decreasing=T)
res1<-new.score.imp[order(-escore)]


Restable<-cbind(res1,res2)
colnames(Restable)<-c("m-I-score", "Full Data Score")


xtable(Restable)




