n <- 1000
X <- cbind(rnorm(n),rnorm(n),rnorm(n))
X.NA <- X
X.NA[,1] <- ifelse(stats::runif(n)<=0.2, NA, X[,1])
X.NA[,2] <- ifelse(stats::runif(n)<=0.2, NA, X[,2])
X.NA[,3] <- ifelse(stats::runif(n)<=0.5, NA, X[,3])

# Remove all observations with all NA

ind<-as.vector(apply(X.NA, 1 , function(x) sum(is.na(x)) == 3))
X.NA<-X.NA[!ind,]


imputations <- list()

imputations[[1]] <- lapply(1:5, function(i) {
 X.loc <- X.NA
 X.loc[is.na(X.NA[,1]),1] <- mean(X.NA[,1],na.rm=TRUE)
 X.loc[is.na(X.NA[,2]),2] <- mean(X.NA[,2],na.rm=TRUE)
 X.loc[is.na(X.NA[,3]),3] <- mean(X.NA[,3],na.rm=TRUE)

 return(X.loc)
})

imputations[[2]] <- lapply(1:5, function(i) {
 X.loc <- X.NA
 X.loc[is.na(X.NA[,1]),1] <- sample(X.NA[!is.na(X.NA[,1]),1],
 size = sum(is.na(X.NA[,1])), replace = TRUE)
 X.loc[is.na(X.NA[,2]),2] <- sample(X.NA[!is.na(X.NA[,2]),2],
                                    size = sum(is.na(X.NA[,2])), replace = TRUE)
 X.loc[is.na(X.NA[,3]),3] <- sample(X.NA[!is.na(X.NA[,3]),3],
                                    size = sum(is.na(X.NA[,3])), replace = TRUE)
 return(X.loc)
})

methods <- c("mean","sample")


### Old
install.packages("Iscores")

library(Iscores)

Iscores(imputations,
        methods,
        X.NA,
        num.proj=5
)


library(parallel)
library(kernlab)

source("R/helpers.R")
source("R/Iscores_new.R")

Iscores_new(imputations,
methods,
X.NA,
num.proj=5
)





