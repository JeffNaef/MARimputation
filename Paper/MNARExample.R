#| warning: false
###Code adapted from Stef van Buuren's book
###https://stefvanbuuren.name/fimd/sec-true.html
library(mice)
library(energy)

create.X <- function(beta = 1, sigma2 = 1, n = 50,
                     run = 1) {
  set.seed(seed = run)
  x <- rnorm(n)
  y <- beta * x + rnorm(n, sd = sqrt(sigma2))
  z <- beta * x + beta*y + rnorm(n, sd = sqrt(sigma2))
  cbind(X1 = x, X2 = y, X3=z)
}


set.seed(10)
n<-8000

Xstar <- create.X(run = 1, n=n, beta=0)

plot(Xstar[,1:2], cex=0.8, col="darkblue")

## MCAR Example
#M<-cbind(sample(c(0,1), size=nrow(Xstar), replace=T, prob = c(1-0.2,0.2) ), 0  )
#X<-Xstar
#X[M==1] <- NA
## MAR Example
#M<-matrix(0, nrow=nrow(Xstar), ncol=ncol(Xstar))
#M[Xstar[,2] > 0, 1] <- sample(c(0,1), size=sum(Xstar[,2] > 0), replace=T, prob = c(1-0.8,0.8) )
#colnames(M) <- cbind("M1", "M2")
## MNAR Example (actually nicely a combination of the two!)
p<-0.5
M<-matrix(0, ncol=ncol(Xstar), nrow=nrow(Xstar))
M[,2]<-sample(c(0,1), size=nrow(Xstar), replace=T, prob = c(1-p,p))
M[Xstar[,2] > 0, 1] <- sample(c(0,1), size=sum(Xstar[,2] > 0), replace=T, prob = c(1-0.8,0.8) )


X<-Xstar
X[M==1] <- NA



## (0) Mean Imputation ##

## (2) Gaussian Imputation ##

blub <- mice(X, method = "norm.nob", m = 1) # visitSequence="arabic"
impnorm<-mice::complete(blub, action="all")[[1]]

# 1. Estimate Regression
#lmodelX1X2<-lm(X1~X2, X=as.X.frame(X[!is.na(X[,1]),])   )
# (same as before)
##RMSE calculation

RMSEcalc<-function(impX){
  
  round(mean(apply(Xstar - impX,1,function(x) norm(as.matrix(x), type="F"  ) )),2)
  
}

energycalc <- function(impX){
  
  round(eqdist.e( rbind(Xstar,impX), c(nrow(Xstar), nrow(impX))  ),2)
  
}



###Test 1: There should be no bias when learning p(x_j|x_{-j}) for all patterns where x_j is observed
lmodelX2<-lm(X2~X1+X3, data=data.frame(X1=Xstar[!is.na(X[,2]),1], X2=X[!is.na(X[,2]),2], X3=Xstar[!is.na(X[,2]),3])   )

summary(lmodelX2)


lmodelX2X3full<-lm(X2~X1+X3, data=as.data.frame(Xstar)   )

summary(lmodelX2X3full)

##same same!!!! This confirms that h^*(x_2 | x_1, x_3 )=p^*(x_2 | x_1, x_3)
## and we learn the right thing! However, the difference comes in imputation for X_2: we need the correction factors 
# P(M=m_1 \mid x)/P(M=m_1 \mid x_1,x_3) = 2\mathbf{1}\{ x_2 > 0\}
# P(M=m_2 \mid x)/P(M=m_2 \mid x_1,x_3) = 2\mathbf{1}\{ x_2 <= 0\}
# P(M=m_3 \mid x)/P(M=m_3 \mid x_1,x_3) = 2\mathbf{1}\{ x_2 <= 0\}

## This is why X_2 is off here! It has negative values when it shouldn't have!
impnormtheory<-Xstar
meanx<-predict(lmodelX2, newdata= as.data.frame(impnormtheory[is.na(X[,2]),])  )
var <- var(lmodelX2$residuals)
impnormtheory[is.na(X[,2]),2] <-rnorm(n=length(meanx), mean = meanx, sd=sqrt(var) )

par(mfrow=c(2,1))
plot(Xstar[!is.na(X[,2]),c("X2","X1")], main="Truth", col="darkblue", cex.main=1.5)
points(Xstar[is.na(X[,2]),c("X2","X1")], col="darkgreen", cex=0.8 )



plot(impnormtheory[!is.na(X[,2]),c("X2","X1")], main=paste("Gaussian Imputation","\nRMSE", RMSEcalc(impnormtheory), "\nEnergy", energycalc(impnormtheory)), col="darkblue", cex.main=1.5)
points(impnormtheory[is.na(X[,2]),c("X2","X1")], col="darkred", cex=0.8 )



par(mfrow=c(2,1))
hist(Xstar[,2])
hist(impnormtheory[,2])



##Instead, lets impute by the corrected distribution: (only focussing on X_2)
impnormtheory2<-cbind(X1=Xstar[,1], X2=X[,2], X3=Xstar[,3] )
meanx<-predict(lmodelX2, newdata= as.data.frame(impnormtheory2[is.na(X[,2]),])  )
var <- var(lmodelX2$residuals)

# # M=(0,0,0) => correction(x)=2\mathbf{1}\{ x_2 <= 0\}
# n1<-nrow(impnormtheory2[M[,1]==0 & M[,2]==0 &  M[,3]==0,])
# impnormtheory2[M[,1]==0 & M[,2]==0 &  M[,3]==0,2] <- -rhalfnorm(n1, theta=sqrt( (pi/2)/var ))
# # M=(1,0,0) => correction(x)=2\mathbf{1}\{ x_2 <= 0\}
# n2<-nrow(impnormtheory2[M[,1]==1 & M[,2]==0 &  M[,3]==0,])
# impnormtheory2[M[,1]==1 & M[,2]==0 &  M[,3]==0,2] <- -rhalfnorm(n2, theta=sqrt( (pi/2)/var ))
# # M=(1,1,0) => correction(x)=2\mathbf{1}\{ x_2 > 0\}
# n3<-nrow(impnormtheory2[M[,1]==1 & M[,2]==1 &  M[,3]==0,])
# impnormtheory2[M[,1]==1 & M[,2]==1 &  M[,3]==0,2] <- rhalfnorm(n3, theta=sqrt( (pi/2)/var ))
# # M=(0,1,0) => correction(x)=2\mathbf{1}\{ x_2 > 0\}
# n4<-nrow(impnormtheory2[M[,1]==0 & M[,2]==1 &  M[,3]==0,])
# impnormtheory2[M[,1]==0 & M[,2]==1 &  M[,3]==0,2] <- rhalfnorm(n4, theta=sqrt( (pi/2)/var ))



# plot(impnormtheory2[!is.na(X[,1]),c("X2","X1")], main=paste("Gaussian Imputation","\nRMSE", RMSEcalc(impnormtheory2), "\nEnergy", energycalc(impnormtheory2)), col="darkblue", cex.main=1.5)
# points(impnormtheory2[is.na(X[,1]),c("X2","X1")], col="darkred", cex=0.8 )


# plot(impnormtheory2[!is.na(X[,1]),c("X2","X1")], main=paste("Gaussian Imputation","\nRMSE", RMSEcalc(impnormtheory2), "\nEnergy", energycalc(impnormtheory2)), col="darkblue", cex.main=1.5)
# points(impnormtheory2[is.na(X[,1]),c("X2","X1")], col="darkred", cex=0.8 )

par(mfrow=c(1,2))
## If beta=0 this is an actual halfnormal distribution!!!
hist(Xstar[is.na(X[,1])&is.na(X[,2]),2], main="", xlab="", freq=F)
hist(Xstar[is.na(X[,2]),2], main="", xlab="", freq=F)

## But this is not acturately captured by the imputation: (which is again gaussian)
#hist(impnormtheory[is.na(X[,1])&is.na(X[,2]),2])


par(mfrow=c(2,2))
plot(Xstar[,c("X2","X1")], main="Truth", col="darkblue", cex.main=2, xlim=c(-4,4))
# plot(Xstar[!is.na(X[,1])& !is.na(X[,2]),c("X2","X1")], main="Truth", col="darkblue", cex.main=1.5)
# points(Xstar[is.na(X[,1])& !is.na(X[,2]),c("X2","X1")], col="darkred", cex=0.8 )
# points(Xstar[!is.na(X[,1])& is.na(X[,2]),c("X2","X1")], col="black", cex=0.8 )
plot(impnorm[,c("X2","X1")], main="Gaussian Imputation", col="darkblue", cex.main=2, xlim=c(-4,4))
hist(Xstar[,2], main="", xlab="", freq=F,xlim=c(-4,4))
hist(impnorm[,2], main="", xlab="", freq=F,xlim=c(-4,4))
