



n<-2000
X<-mvrnorm(n=2000, mu=c(0,0), Sigma=matrix( c(1,0.7,0.7,1), nrow=2, byrow =T   ))



XnoNA<-X[1:round(n/2), ]
XNA<-X[ (round(n/2) +1 ):n, ]




####Test
# linear<-lm(XnoNA[,1]~XnoNA[,2])
# mean <-  cbind(1,XNA[,2])%*%linear$coefficients
# varreg <- var(linear$residuals)
# cdftrue <- function(t) {pnorm(t, mean = mean, sd = sqrt(varreg), lower.tail = TRUE, log.p = FALSE)}
# tmp2<-lapply( seq(min(Ytrain), max(Ytrain), by=0.1), function(t)  (cdftrue(t)-1*(mean < t) )^2) 
# res2<-do.call(cbind,tmp2)
# mean(colMeans(res2))
# 
# tmp2<-lapply( seq(min(Ytrain), max(Ytrain), by=0.1), function(t)  (cdftrue(t)-1*(XNA[j,1] < t) )^2) 
# res2<-do.call(cbind,tmp2)
# mean(colMeans(res2))
# 
# 
# 
# tmp2<-lapply( seq(min(Ytrain), max(Ytrain), by=0.1), function(t)  (cdftrue(t)-1*(rep(100,nrow(XNA)) < t) )^2) 
# res2<-do.call(cbind,tmp2)
# mean(colMeans(res2))
# 
# 
# #scorej[i] <- -mean(do.call(cbind,tmp))
# #scorej[i] <- -mean(do.call(cbind,tmp2))
# ####
# 
# #scorej[i] <- mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = c(Ytrain), w=Fhat[j,]) }))
# 
# # Truth
# sctruth <-mean(sapply(1:nrow(XNA), function(j)  { crps_sample(y = XNA[j,1], dat = rnorm(n=1000, mean=mean[j], sd=sqrt(varreg)) ) }))
# scimp <-mean(sapply(1:nrow(XNA), function(j)  { crps_sample(y = mean[j], dat = rnorm(n=1000, mean=mean[j], sd=sqrt(varreg)) ) }))
# scimp <-mean(sapply(1:nrow(XNA), function(j)  { crps_sample(y = mean[j], dat = rnorm(n=1000, mean=100, sd=0.05 ) )}))

# (bad) imputation
linear<-lm(XnoNA[,1]~XnoNA[,2])

#Learn model on imputed data...
Ytraintruue<-XNA[,1]
Ytrainnormpredict<-cbind(1,XNA[,2])%*%linear$coefficients
Xtrain<-XNA[,2]
#... and evaluate on true data
Ytest<- XnoNA[,1]
Xtest<- XnoNA[,2]

###correct way:
linear<-lm(Ytrainnormpredict~Xtrain)
mean <-  cbind(1,Xtest)%*%linear$coefficients
varreg <- var(linear$residuals)
#cdftrue <- function(t) {pnorm(t, mean = mean, sd = sqrt(varreg), lower.tail = TRUE, log.p = FALSE)}
scmean <-mean(sapply(1:length(Ytest), function(j)  { crps_sample(y = Ytest[j], dat = rnorm(n=1000, mean=mean[j], sd=sqrt(varreg)) ) }))

linear<-lm(Ytraintruue~Xtrain)
mean <-  cbind(1,Xtest)%*%linear$coefficients
varreg <- var(linear$residuals)
#cdftrue <- function(t) {pnorm(t, mean = mean, sd = sqrt(varreg), lower.tail = TRUE, log.p = FALSE)}
sctruth <-mean(sapply(1:length(Ytest), function(j)  { crps_sample(y = Ytest[j], dat = rnorm(n=1000, mean=mean[j], sd=sqrt(varreg)) ) }))



scnonsensical <-mean(sapply(1:nrow(XNA), function(j)  { crps_sample(y = Ytest[j], dat = rnorm(n=1000, mean=100, sd=sqrt(varreg)) ) }))







X<-X.NA

Ximp<-imputations$truth[[1]]
Ximp<-imputations$missForest[[1]]
Ximp<-imputations$norm.predict[[1]]



colnames(X) <- colnames(Ximp) <- paste0("X", 1:ncol(X))

X<-as.matrix(X)
Ximp<-as.matrix(Ximp)

n<-nrow(X)
p<-ncol(X)

##Step 1: Reoder the data according to the number of missing values
## (least missing first)
numberofmissingbyj<-sapply(1:p, function(j)  sum(is.na(X[,j]))  )
X0<-X
Ximp0<-Ximp
X<-X[,order(numberofmissingbyj, decreasing=T)]
Ximp<-Ximp[,order(numberofmissingbyj, decreasing=T)]


## Done in the function
M<-1*is.na(X)
colnames(M) <- colnames(X)

# Order first according to most missing values

# Get dimensions with missing values (all other are not important)
dimwithNA<-(colSums(M) > 0)
dimwithNA<-dimwithNA[dimwithNA==TRUE]
index<-1:ncol(X)
scorej<-matrix(NA, nrow= sum(dimwithNA), ncol=1)
weight<-matrix(NA, nrow= sum(dimwithNA), ncol=1)
i<-0

for (j in names(dimwithNA)[1]){
  

  ## TO DO: 
  ##- Need to do something about class-imbalance!
  ##- Need to weigh the resulting score according to number of obs!
  
  print( paste0("Dimension ", i, " out of ", sum(dimwithNA) )   ) 
  
  
  require(drf)
  
  Ximp1<-Ximp[M[,j]==1, ]
  Ximp0<-Ximp[M[,j]==0, ]
  
  n1<-nrow(Ximp1)
  n0<-nrow(Ximp0)
 
  nmin<-min(n0,n1)
  
  if (nmin==n1){
    # Train DRF on observed data
    Xtrain<-Ximp0[,!colnames(Ximp0) %in% j, drop=F]
    Ytrain<-Ximp0[,j, drop=F]
    
    #Evaluate on imputed data
    Xtest <- Ximp1[,!colnames(Ximp1) %in% j, drop=F]
    Ytest <-Ximp1[,j, drop=F]
    
  }else{
    # Train DRF on imputed data
    Xtrain<-Ximp1[,!colnames(Ximp1) %in% j, drop=F]
    Ytrain<-Ximp1[,j, drop=F]
    
    # Evaluate on observed data
    Xtest <- Ximp0[,!colnames(Ximp0) %in% j, drop=F]
    Ytest <-Ximp0[,j, drop=F]
    
  }
  
  fit <- drf(X=Xtrain, 
             Y=Ytrain,  splitting.rule = "FourierMMD", num.trees = 2000, num.features=20, min.node.size = 15)
  
  
  # Predict on the test data
  Fhat <- predict(fit, newdata=Xtest  )$weights
  
  
  
  linear<-lm(Ytrain[,1]~Xtrain)
  mean <-  cbind(1,Xtest)%*%linear$coefficients
  varreg <- var(linear$residuals)
  cdftrue <- function(t) {pnorm(t, mean = mean, sd = sqrt(varreg), lower.tail = TRUE, log.p = FALSE)}
  
  
  # Calculate Energy score
  nmax<-max(n0,n1)
  
  
  ##Energy score
  # ##########################
  # K<-abs(outer( c(Ytrain), c(Ytrain), '-'))
  # Ktest<-abs(outer(c(Ytrain),c(Ytest), '-'))
  # normalizingpart<-1/2*apply(Fhat,1, function(x) x%*%K%*%x )
  # firstpart <- diag(Fhat%*%Ktest   )
  # scorej[i]<- mean(normalizingpart - firstpart)
  # ##########################
  
  ##Energy score
  # ##########################
  cdf<-function(t){

    Fhat%*%(1*(Ytrain < t))

  }
  tmp<-lapply( seq(min(Ytrain), max(Ytrain), by=0.1), function(t)  (cdf(t)-1*(Ytest < t) )^2)
  
  res1<-do.call(cbind,tmp)
  
  mean(res1)
  
  
  tmp2<-lapply( seq(min(Ytrain), max(Ytrain), by=0.1), function(t)  (cdftrue(t)-1*(Ytest < t) )^2)
  
  res2<-do.call(cbind,tmp2)
  
  #scorej[i] <- -mean(do.call(cbind,tmp))
  scorej[i] <- -mean(do.call(cbind,tmp2))
  
  # ##########################

  ##
  mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = c(Ytrain), w=Fhat[j,]) }))
  mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = Ytrain[sample(1:nrow(Ytrain), size=1000, replace=T, prob=  Fhat[j,])] ) }))
  mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = rnorm(n=1000, mean=mean, sd=sqrt(varreg)) ) }))
  
  ##same same!
  
  
  #es_sample(y = Ytest[1,, drop=F], dat = Ytrain, w=Fhat[1,, drop=F])
  
  
  es_sample(y = obs[1], dat = fc_sample[1, , drop = FALSE])
  crps_sample(y = obs[1], dat = fc_sample[1, ])
  
  
  # # Kernel Score
  # # ##########################
  # Ytrainb<-scale(Ytrain)
  # Ytestb <- scale(Ytest)
  # bandwidth_Y <- drf:::medianHeuristic(Ytrainb)
  # k_Y <- rbfdot(sigma = 1/(2*bandwidth_Y^2))
  # K <- kernelMatrix(k_Y, Ytrainb, y = Ytrainb)
  # 
  # Ktest<-kernelMatrix(k_Y, Ytrainb, y = Ytestb)
  # normalizingpart<-1/2*apply(Fhat,1, function(x) x%*%K%*%x )
  # firstpart <- diag(Fhat%*%Ktest   )
  # scorej[i]<- mean(normalizingpart - firstpart)
  # # ##########################
  
  
  
  
  # Ytrainb<-Ytrain
  # bandwidth_Y <- drf:::medianHeuristic(Ytrainb)
  # k_Y <- rbfdot(sigma = 1/(2*bandwidth_Y^2))
  # K <- kernelMatrix(k_Y, Ytrainb, y = Ytrainb)
  # 
  # 
  # 
  # scorej[i]<- mean( sapply(1:nrow(Fhat) ,  function(l){
  #   
  #   Z1<-Ytrain[sample(1:nrow(Ytrain), size=max(nmin, nmax/2), replace = T, Fhat[l,]),]
  #   Z2<-Ytrain[sample(1:nrow(Ytrain), size=max(nmin, nmax/2), replace = T, Fhat[l,]),]
  #   
  #   K <- kernelMatrix(k_Y, Z1, y = Z2)
  #   k <- kernelMatrix(k_Y, Z1, y=Ytest[l,])
  #   
  #   return(0.5* mean(K)  - mean(k) )
  #   
  # } )
  # )
  
  
  # OOB: Ask herb!
  #Fhat <- get_causal_sample_weights(fit, newdata = NULL, newtreatment = matrix(rep(0, n0+n1), ncol=1), g = rep(0.5, n0+n1))
  #Hhat <- get_causal_sample_weights(fit, newdata = NULL, newtreatment = matrix(rep(1, n0+n1), ncol=1), g = rep(0.5, n0+n1))
  
   
}
  
