###Code adapted from Stef van Buuren's book
###https://stefvanbuuren.name/fimd/sec-true.html
library(mice)
source("helpers.R")


create.data <- function(beta = 1, sigma2 = 2, n = 50,
                        run = 1) {
  set.seed(seed = run)
  x <- rnorm(n)
  y <- beta * x + rnorm(n, sd = sqrt(sigma2))
  cbind(x = x, y = y)
}


make.missing <- function(data, p = 0.5){
  rx <- rbinom(nrow(data), 1, p)
  data[rx == 0, "x"] <- NA
  data
}

######Figure out why the Causal-DRF example works, but the IntroExample is NaN #####



set.seed(10)

truedata <- create.data(run = 1, n=2000)

#data <- genMask(truedata, mech = "MAR", pmiss=0.2)
data <- make.missing(truedata)

imputations<- doimputation(X.NA=data, methods="norm.predict", parallelize=F, m=1, visitSequence="arabic")
impnormpredict<-imputations$imputations$norm.predict$'1'

par(mfrow=c(1,2))
plot(impnormpredict[!is.na(data[,1]),], main="Regression Imputation", cex=0.8, col="darkblue", cex.main=1.5)
points(impnormpredict[is.na(data[,1]),], col="darkred", cex=0.8 )

# ###use correct norm prediction instead!!
# imputations<- doimputation(X.NA=data, methods="cart", parallelize=F, m=1)
# impcart<-imputations$imputations$cart$'1'
# 
# plot(impcart[!is.na(data[,1]),], main="CART Imputation", col="darkblue", cex.main=1.5)
# points(impcart[is.na(data[,1]),], col="darkred", cex=0.8 )
# ###
imputations<- doimputation(X.NA=data, methods="norm.nob", parallelize=F, m=1)
impcart<-imputations$imputations$norm.nob$'1'

plot(impcart[!is.na(data[,1]),], main="Gaussian Imputation", col="darkblue", cex.main=1.5)
points(impcart[is.na(data[,1]),], col="darkred", cex=0.8 )


### For Presentation
par(mfrow=c(1,3))
plot(impnormpredict[!is.na(data[,1]),], main="Regression Imputation", cex=0.8, col="darkblue", cex.main=1.5)
points(impnormpredict[is.na(data[,1]),], col="darkred", cex=0.8 )
plot(impcart[!is.na(data[,1]),], main="Gaussian Imputation", col="darkblue", cex.main=1.5)
points(impcart[is.na(data[,1]),], col="darkred", cex=0.8 )

#plot(truedata[!is.na(data[,1]),], main="Truth", col="darkblue", cex.main=1.5)
#points(truedata[is.na(data[,1]),], col="green", cex=0.8 )
plot(truedata, main="Truth", col="darkblue", cex.main=1.5)




##Biased beta estimates

lm(y~x, data=impnormpredict)
lm(y~x, data=impcart)


sum(sqrt(colMeans((truedata - impnormpredict)^2)))
sum(sqrt(colMeans((truedata - impcart)^2)))


lm(y~x, data=impnormpredict)

lm(y~x, data=impcart)

##Complete.cases
lm(y~x, data=as.data.frame(data[complete.cases(data),])   )



