library(mice)
source("helpers.R")
library(MASS)
library(energy)


set.seed(10)

p <- 10  
A <- matrix(runif(p^2)*2-1, ncol=p) 
Sigma <- t(A) %*% A

mu<-runif(p)

truedata <- mvrnorm(n=2000, mu=mu, Sigma=Sigma)

data <- genMask(truedata, mech = "MAR", pmiss=0.1)


imputations<- doimputation(X.NA=data, methods=c("cart","missForest"), parallelize=F, m=1)



imputations <-imputations$imputations
imputations$truth[[1]] <- truedata



## Score function
# unconditional score is not valid!
#Iscores_new(X,imputations, score="unconditional")
Iscores_new(data,imputations, score="conditional")



## Energy Distance

