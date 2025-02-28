library(mice)
library(miceDRF)
library(drf)
library(MASS)
library(miceDRF)
source("helpers.R")

###Example 1: Showing that all patterns need to be used for imputation##
set.seed(10)
n<-10000

Xstar<-matrix( runif(n=3*n), nrow=n, ncol=3    )
Mindex<-sapply(1:n, function(i) 
  sample(1:3, size=1, replace=F, prob=c(Xstar[i,1]/3, 2/3-Xstar[i,1]/3, 1/3    )))
M<-t(sapply(1:n, function(i) if (Mindex[i]==1){c(0,0,0)}else if (Mindex[i]==2) {c(0,1,0)} else if  (Mindex[i]==3) {c(1,0,0)}) )


X<-Xstar
X[Mindex==2,2]<-NA
X[Mindex==3,1]<-NA

head(cbind(X,M,Mindex))



hist(Xstar[,2])




# What we want to impute = unconditional distribution
hist(Xstar[Mindex==3,1], xlab="", main="", probability=T)

# Fully observed pattern:
hist(Xstar[Mindex==1,1], xlab="", main="", probability=T)


## Learning from both patterns
hist(Xstar[Mindex==1 | Mindex==2,1], xlab="",main="", probability=T)

#left: Distribution we want to impute X_2 \mid M=m_3
#middle: Distribution of X_2 in the fully observed pattern (X_2 \mid M=m_1)
#right: Distribution in all observed patterns (Mixture of X_2 \mid M=m_1 and X_2 \mid M=m_2)
# Close the PNG device
dev.off()



###Example 2: Showing MAR allows for changes in observed distribution
library(MASS)
source("helpers.R")
require(drf)
n<-10000
U<-runif(n)

## In both cases: X_1|X_2 ~ N(X_2, 1)

X2m1<-mvrnorm(n= sum(U <= 1/2),  mu=0, Sigma=1  )
X1m1<-1*X2m1 + rnorm(n= sum(U <= 1/2))
X2m2<-mvrnorm(n= sum(U > 1/2),  mu=5, Sigma=1  )
X1m2<-1*X2m2 + rnorm(n= sum(U > 1/2))

Xstar<-rbind( cbind(X1m1, X2m1), cbind(X1m2, X2m2)   )

par(mfrow=c(1,2))
hist( X2m1, xlab="", main="")
hist( X2m2, xlab="", main="")

#Left: Distribution of the observed X_2 in pattern 1 (X_2 \mid M=m_1)
#Right: Distribution of the observed X_2 in pattern 2 (X_2 \mid M=m_2)





##Try to impute this Example ##
##This is not ideal however, since RF tends to be bad with only one X!
set.seed(10)
n<-2000
U<-runif(n)

## In both cases: X_1|X_2 ~ N(X_2, 1)
##Careful: Compared to the paper the two patterns are reversed:
## m_1=(1,0), m_2=(0,1)
## (X_1, X_2) | M=m_1 ~ N((5,5), Sigma  )
## (X_1, X_2) | M=m_2 ~ N((0,0), Sigma  )
X2m1<-mvrnorm(n= sum(U <= 1/2),  mu=5, Sigma=1  )
X1m1<-1*X2m1 + rnorm(n= sum(U <= 1/2))
X2m2<-mvrnorm(n= sum(U > 1/2),  mu=0, Sigma=1  )
X1m2<-1*X2m2 + rnorm(n= sum(U > 1/2))


Xstar<-rbind( cbind(X1m1, X2m1), cbind(X1m2, X2m2)   )

X.NA<-Xstar
X.NA[1:nrow(X1m1),1]<-NA

imputations<- doimputation(X.NA=X.NA, methods="cart", parallelize=F, m=1, min.node.size=15)
impcart<-imputations$imputations$cart$'1'

# impdrf (reuse code!)
imputations<- doimputation(X.NA=X.NA, methods="DRF", parallelize=F, m=1, min.node.size=15)
impDRF<-imputations$imputations$DRF$'1'


reg<- lm(X1m2~X2m2)
Ysample<-reg$coefficients[1] + reg$coefficients[2]*X2m1+ rnorm(nrow(X2m1),0, sd=var(reg$residuals))
impreg<-X.NA
impreg[1:nrow(X1m1),1]<-Ysample



# DRF<-drf(Y=X1m2, X=X2m2, num.trees=10, min.node.size=15, num.features = 100,
#          compute.oob.predictions = F)
# HhatDRF<-predict(DRF, newdata=X2m1)$weights
# 
# Ygen<-sapply(1:nrow(X1m1), function(i) X1m2[sample(1:nrow(X1m2), size=1, replace = T, HhatDRF[i,])])
# impDRF<-X.NA
# impDRF[1:nrow(X1m1),1]<-Ygen

# ## Robust DRF trial = > Doesnt really work ###
# X2m2rob<- (X2m2 - mean(X2m2))/abs(X2m2-mean(X2m2))
# DRFrob<-drf(Y=X1m2, X=X2m2rob, num.trees=2000, min.node.size=15)
# X2m1rob<- (X2m1 - mean(X2m1))/abs(X2m1-mean(X2m1))
# HhatDRFrob<-predict(DRFrob, newdata=X2m1rob)$weights
# 
# Ygen<-sapply(1:nrow(X1m1), function(i) X1m2[sample(1:nrow(X1m2), size=1, replace = T, HhatDRFrob[i,])])
# impDRFrob<-X.NA
# impDRFrob[1:nrow(X1m1),1]<-Ygen






# engr = engression(X=X2m2,Y=X1m2,num_epochs = 100)
# Ysample = predict(engr,X2m1,type="sample",nsample=1)
# impeng<-X.NA
# impeng[1:nrow(X1m1),1]<-Ysample

#linreg<-lm(X1m2~ X2m2 )
#X1m2rob<- X1m2-cbind(rep(1,nrow(X2m2)), X2m2)%*%linreg$coefficients
#DRFlin<-drf(Y=X1m2rob, X=X2m2, num.trees=2000, min.node.size=15)
#HhatDRFlin<-predict(DRFlin, newdata=X2m1)$weights

#Ygen0<-sapply(1:nrow(X1m1), function(i) X1m2rob[sample(1:nrow(X1m2rob), size=1, replace = T, HhatDRFlin[i,])])
#Ygen <- Ygen0 + cbind(rep(1,nrow(X2m1)), X2m1)%*%linreg$coefficients
#HhatDRFlin<-X.NA
#HhatDRFlin[1:nrow(X1m1),1]<-Ygen

# Open a PNG device for saving the plot
png(filename = "Example2.png", 
    width = 1200,    # Width in pixels
    height = 700,    # Height in pixels
    res = 120)       # Resolution in dpi

par(mfrow=c(2,2))
hist( X1m1, xlab="", ylab="", main="Truth", probability = T, ylim=c(0,1.2), xlim=c(-5,10))
# Almost perfect imputation!
hist( impreg[1:nrow(X1m1),1], xlab="", ylab="", main="mice-norm.nob", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
hist( impcart[1:nrow(X1m1),1], xlab="", ylab="", main="mice-cart", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
hist( impDRF[1:nrow(X1m1),1], xlab="", ylab="", main="mice-DRF", probability = T, xlim=c(-5,10), ylim=c(0,1.2))
#hist( impDRFrob[1:nrow(X1m1),1], xlab="", ylab="", main="DRF Robust", probability = T, xlim=c(-5,10), ylim=c(0,0.3))
#hist( impeng[1:nrow(X1m1),1], xlab="", ylab="", main="Engression", probability = T, xlim=c(-5,10), ylim=c(0,0.3))
# hist( impDRF[1:nrow(X1m1),1],probability = T, add=T, col=rgb(0,1,0,0.3))
# hist( impcart[1:nrow(X1m1),1], probability = T, add=T, col=rgb(1,0,1,0.3))
# hist( impeng[1:nrow(X1m1),1], probability = T, add=T, col=rgb(1,0.5,1,0.3))


dev.off()