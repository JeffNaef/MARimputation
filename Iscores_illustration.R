source("test_helpers2.R")


d<-3
patterns <- matrix(rep(NA,d),nrow=1)

set.seed(1)

X <- genDataNoNA_synthetic(n.train = 10,
                           n.test = 10,
                           d=d,
                           dataset = "spiral",
                           patterns = patterns)$train
X <- as.matrix(X)
colnames(X) <- NULL


missing.mech <- "MCAR"

X.NA <- genMask(X, mech = missing.mech, pmiss=0.3)
# Throw away observations which are only NA
X <- X[rowSums(is.na(X.NA)) < d,]
X.NA <- X.NA[rowSums(is.na(X.NA)) < d,]
View( round( X.NA*100 ))
View( round( X*100 ))

X.full.observed <- X.NA[complete.cases(X.NA[,2:3]), 2:3]
X.miss <- X[!complete.cases(X.NA[,2:3]), 2:3]
View( round( X.full.observed*100 ))
View( round( X.miss*100 ))
