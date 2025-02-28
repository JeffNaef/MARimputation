
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

methods <- c("norm.nob")


#methods <- c("norm.predict", "norm.nob")

nrep.total<-10


###Notes: Whatever I try to think of is essentially just a large sample version 
### of the score

### Take the imputation H_n. For each variable j, learn the parameters of X_j | X_{-j} from the set 
### where X_j is not observed. This is H_{n, X_j \mid X_{-j}}. Now generate a large number of
### observations from the patterns where X_j is observed (here just a Gaussian) and 
### use H_{n, X_j \mid X_{-j}}for each observation of X_{-j} in the other group to score 


## Maybe use spam??

dataset <- "multivariateGaussian"


m <- 1


set.seed(2) #1
seeds <- sample(c(0:2000),100,replace = FALSE)




# Given the observed data build the data containing missing values
#Beta<-matrix(runif(n=9, min=-2, max=2), nrow=3,ncol=3)
d<-3
Beta<-matrix(c(0.5,1,1.5), nrow=d,ncol=d, byrow=T)


n<-1000

C<-matrix(0, nrow=d, ncol=d)
for (t in 1:d){
  for (s in 1:d){
    
    
    C[t,s] <- 0.5^(abs(t-s))
    
  }
  
  
}

### Data generating function

datageneration <- function(n, Ntest, Beta, C){

## Build the observed data
Xobs1<- MASS::mvrnorm(n = Ntest+n, mu = rep(3, 3), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*n, d=3)$train
Xobs2<- MASS::mvrnorm(n = Ntest+n, mu = rep(-3, 3), Sigma = C) #genDataNoNA_synthetic(dataset = dataset, n.train = 10*n, d=3)$train
Xobs3<- MASS::mvrnorm(n = Ntest+n, mu = rep(0, 3), Sigma = C)#genDataNoNA_synthetic(dataset = dataset, n.train = 10*n, d=3)$train


X.NA1 <- Xobs1%*%Beta + matrix( rnorm(n=3*(Ntest+n),sd=2), nrow=Ntest+n, ncol= 3   )
X.NA2 <- Xobs2%*%Beta + matrix( rnorm(n=3*(Ntest+n),sd=2), nrow=Ntest+n, ncol= 3   )
X.NA3 <- Xobs3%*%Beta + matrix( rnorm(n=3*(Ntest+n),sd=2), nrow=Ntest+n, ncol= 3   )

Fulldata<-cbind( rbind( X.NA1, X.NA2, X.NA3 ), rbind(Xobs1, Xobs2, Xobs3)    )

M1<-matrix(c(1,0,0), nrow=(Ntest+n), ncol=3, byrow=T)
M2<-matrix(c(0,1,0), nrow=(Ntest+n), ncol=3, byrow=T)
M3<-matrix(c(0,0,1), nrow=(Ntest+n), ncol=3, byrow=T)

FullM<-cbind( rbind( M1, M2, M3 ), matrix(0, nrow=3*(Ntest+n), ncol=3)    )




X<-cbind( rbind( X.NA1[ 1:n,], 
                 X.NA2[ 1:n,], X.NA3[ 1:n,]), 
          rbind(Xobs1[1:n,], Xobs2[ 1:n,], Xobs3[ 1:n,])    )
M<-cbind( rbind( M1[ 1:n,], M2[ 1:n,], M3[ 1:n,] ), 
          matrix(0, nrow=3*n, ncol=3)    )

Xtest<-cbind( rbind( X.NA1[ -(1:n),], 
                     X.NA2[-(1:n),], X.NA3[ -(1:n),]), 
              rbind(Xobs1[ -(1:n),], Xobs2[ -(1:n),], Xobs3[ -(1:n),])    )
Mtest<-cbind( rbind( M1[ -(1:n),], M2[ -(1:n),], M3[ -(1:n),] ), 
              matrix(0, nrow=dim(Xtest)[1], ncol=3)    )

return(list(Fulldata=Fulldata, FullM=FullM, X=X, M=M,Xtest=Xtest, Mtest=Mtest  ))

}


nvec<-c(200, 400,1000,2000,4000, 8000)
Ntest<-10000


#length(ids.jack)
#Results <- lapply(1:10, function(s){
Results<-list()

for (s in 1:length(nvec)){
  set.seed(seeds[s])
  
  n <-nvec[s]
  
  tmp<-datageneration(n=n, Ntest=Ntest,Beta=Beta, C=C)
  
  X = tmp$X
  M = tmp$M
  Xtest=tmp$Xtest
  Mtest=tmp$Mtest
  
  d<-ncol(X)

  
  X.NA<-X
  X.NA[M==1] <- NA
  

  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:6)
  
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
  

  start_time <- Sys.time()
  new.score.list.imp <- Iscores_new(X.NA,imputations=imputations, imputationfuncs=imputationfuncs, N=100)
  end_time <- Sys.time()
  
  end_time-start_time
  

  new.score.imp <- unlist(lapply(new.score.list.imp, function(x) x$score))
  Hn<-imputations[["norm.nob"]][[1]]
  colnames(Hn)<-paste0("X",1:ncol(X))
  scorej<-c()
  #### Check for j=3 ######
  #########################
  
  ###
  for (j in 1:3){
  #Fulldata[FullM[,j]!=1,]
  

  ### Need to change this to O_j!!
  Oj<-apply(X.NA[M[,j]==0,-j , drop=F],2,function(x) !any(is.na(x)) )
  indexOj<-(1:d)[-j][Oj]
  
  ## Learn H_{n, X_j \mid X_{O_j}}
  formula <- as.formula(paste(paste0("X", j), "~", paste(paste0("X", indexOj), collapse = " + "))) 
  #as.formula(paste(paste0("X", j), "~", paste0("X", (1:d)[-j]  )))
  regj<-lm(formula, data=Hn[M[,j]==1,] )

  
  
  ## Now evaluate H_{n, X_j \mid X_{O_j}} on a big independent test sample
  

  
  # only keep points where X_j is observed (those are the points we want to score)
  Xtestj<-Xtest[Mtest[,j]==0,]
  # Get the mean 
  meanX<-regj$coefficients[1] + Xtestj[,indexOj]%*%regj$coefficients[-1] 
  # Get the variance
  varhat<-var(regj$residuals)
  

  #Xmatrix<-do.call(cbind, lapply(1:nrow(Xtest), function(l)  rnorm(n, mean = 0, sd = 1)  ))
  ## Score with a larger test sample
  #scorej[j] <- -mean(sapply(1:2000, function(l)  {  crps_sample(y = Xtest[l,j], dat = Xmatrix[l,]) }))
  
  scorej[j] <- -mean(sapply(1:nrow(Xtestj), function(l)  { 
    Xmatrixl <-  rnorm(n=500, mean = meanX[l], sd = sqrt(varhat))
    crps_sample(y = Xtestj[l,j], dat = Xmatrixl) })) #Xmatrix[l,]

  }
  
  print(scorej)
  print(new.score.list.imp$norm.nob$scorelist)
  
  # nrow(Xtest)
  # print("drf-score2:")
  # print( sort( round( unlist(DRIScore)/sum(unlist(DRIScore)),3) , decreasing=T)   )
  # print("m-score2:")
  # print( sort( round( unlist(new.score.imp)/sum(unlist(new.score.imp)),3) , decreasing=T)   )
  # print("e-score")
  # print( sort( round(escore/sum(escore),3) , decreasing=T)   )
  
  print(paste0("nrep ",s, " out of ", nrep.total ))
  
  Results[[s]] <- list(scoretheory = scorej, scoreestimate = new.score.list.imp$norm.nob$scorelist )
  
  
  #return(list(new.score.imp = new.score.imp,DRIScore=DRIScore , energy.score=escore))
  
  
}



tmp<-do.call( rbind, lapply(Results, function(x) c(mean(x$scoretheory), mean(x$scoreestimate)) ))
colnames(tmp) <- c("Theory", "Estimate")
tmp


tmpreal<-datageneration(n=n, Ntest=Ntest, C=C, Beta=Beta)
Fulldata=tmpreal$Fulldata
FullM=tmpreal$FullM

##Problem: H_n changes with the data!

scorejtrue<-c()
###
for (j in 1:3){
  #Fulldata[FullM[,j]!=1,]
  
  
  ### Need to change this to O_j!!
  indexOj<-c(4,5,6)
  
  ## Learn H_{n, X_j \mid X_{O_j}}
  formula <- as.formula(paste(paste0("X", j), "~", paste(paste0("X", indexOj), collapse = " + "))) 
  #as.formula(paste(paste0("X", j), "~", paste0("X", (1:d)[-j]  )))
  regj<-lm(formula, data=data.frame(Fulldata[FullM[,j]==1,]) )
  ### Here we could also derive the actual true values! But the estimates should be
  ### so close it should not be an issue
  
  # only keep points where X_j is observed (those are the points we want to score)
  Xtestj<-Fulldata[FullM[,j]==0,]
  # Get the mean 
  meanX<-regj$coefficients[1] + Xtestj[,indexOj]%*%regj$coefficients[-1] 
  # Get the variance
  varhat<-var(regj$residuals)
  
  
  #Xmatrix<-do.call(cbind, lapply(1:nrow(Xtest), function(l)  rnorm(n, mean = 0, sd = 1)  ))
  ## Score with a larger test sample
  #scorej[j] <- -mean(sapply(1:2000, function(l)  {  crps_sample(y = Xtest[l,j], dat = Xmatrix[l,]) }))
  
  scorejtrue[j] <- -mean(sapply(1:nrow(Xtestj), function(l)  { 
    Xmatrixl <-  rnorm(n=500, mean = meanX[l], sd = sqrt(varhat))
    crps_sample(y = Xtestj[l,j], dat = Xmatrixl) })) #Xmatrix[l,]
  
}



res<-cbind(nvec,tmp, mean(scorejtrue))

res

filename ="Application_2_consistencycheck"

assign(filename, res)
save(res, file=paste(filename, ".Rda",sep=""))



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













