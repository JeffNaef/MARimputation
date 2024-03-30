
### 2d imputation example ##

n<-1000
library(data.table)
library(mlbench)
library(mice)
#dat<-data.table(l=sample(x=c(0,1), size=n, replace=T),Y1=rnorm(n), Y2=rnorm(n), X1=rnorm(n), X2=rnorm(n))
dat<- mlbench::mlbench.spirals(n, sd = 0.05)$x


tmp<-ampute(dat, prop=0.3)
X.NA<-tmp$amp

imputations<-list()

#blocklist<-list()
#blocklist[[1]]<-colnames(datNA)[1]
#blocklist[[2]]<-colnames(datNA)[2]

imputedDRF <- blockmiceDRF(X.NA, blocksize=1) #mice(X.NA, method="myfunc", num.features=10, m=1, blocks = blocklist ) ## lets be honest, "cart" will do the job
imputations[["DRF"]][[1]]<-complete(imputedDRF)

imputedcart = mice(X.NA, method="cart", m=1) ## lets be honest, "cart" will do the job
imputations[["cart"]][[1]]<-complete(imputedcart)

imputations[["true.data"]][[1]]<-dat
X<-dat



layout(mat = matrix(c(1, 2, 3), 
                    nrow = 1, 
                    ncol = 3),
       heights = c(1, 1),    # Heights of the two rows
       widths = c(1, 1))     # Widths of the two columns

# Plot 1: Scatterplot
#par(mar = c(0.1, 0.4, 0.5, 0.4))

ind.candidates <- which(!complete.cases(X.NA))
ind.candidates <- sort(ind.candidates)

par(mfrow=c(1,3), mai=c(0.4, 0.2, 0.4, 0.01))

for (method in c("true.data", "DRF", "cart")){
  
  if (method == "true.data"){
    imp <- X[ind.candidates,]
  }else{
    imp<-imputations[[method]][[1]][ind.candidates,]
  }
  # Just change the naming
  if (method=="true.data"){
    nameofmethod <- "truth"
  }
  if (method=="cart"){
    nameofmethod <- "mice-cart"
  }
  if (method=="DRF"){
    nameofmethod <- "mice-DRF"
  }
  if (method=="loess"){
    nameofmethod <- "loess"
  }
  
  # score in the title
  plot(X.NA[-ind.candidates,],pch=19,col="gray87",cex=1.2,yaxt="n",
       xlab = "", ylab="" ,font.lab=1, 
       main=nameofmethod,cex.main=2, xaxt="n")
  Axis(side=1, at=c(-1, 0, 1), labels=c(-1, 0, 1), cex.axis=2,padj=0.6
  )# paste0(method, ", ", round(unname(CI[2,method]),2)," (", round(unname(CI[1,method]),2), ",",round(unname(CI[3,method]),2) , ")")
  points(imp,pch=0,col="gray40",cex=1.2,cex.lab=1,font.lab=1,
  )
  #paste0(method, ", ", round(unname(CI[2,method]),2)," (", round(unname(CI[1,method]),2), ",",round(unname(CI[3,method]),2) , ")")
  #paste0(method, ", ", round(unname(scoredr[method,]),2))
  # no score in the title
  #plot(X.NA[-ind.candidates,],pch=19,col="darkblue",cex=0.5, xlab = expression(x[1]), ylab = expression(x[2]),cex.lab=1.75,font.lab=1,main=method )
  #points(imp,pch=19,col="grey",cex=0.5, xlab = expression(x[1]), ylab = expression(x[2]),cex.lab=1.75,font.lab=1)
  #print(paste("Method=", method, ", RMSE for mean: ", sum( ( colMeans( imputations$imputations[[method]][[1]]) - colMeans(X) )^2)   ))
  #print(paste("Method=", method, ", RMSE for var: ", sum( ( colVars( imputations$imputations[[method]][[1]]) - colVars(X) )^2)   ))
  #print(paste("Method=", method, ", RMSE for var: ", sum( ( colVars( imputations$imputations[[method]][[1]]) - colVars(X) )^2)   ))
  
}


## Large-Dimensional Example

library(mvtnorm)
library(missMDA)
library(Iscores)
library(mice)

n<-1000
d<-12
set.seed(1)



## Find better datasets!!!
Sigma<-matrix(0, nrow=d, ncol=d)

for (i in 1:d){
  
  for (j in 1:d){
    
    Sigma[i,j]=0.7^(abs(i-j))
  }
}


# library(datamicroarray)
# data(christensen)
# X <-christensen$x[,1:d]
# colnames(X)<-NULL


X<-rmvnorm(n=n, sigma=Sigma )#matrix(rnorm(n*d), nrow=n, ncol=d)





X.NA<-ampute(X, prop=0.01)$amp


imputations<-list()

imputedDRF<-blockmiceDRF(X.NA, blocksize=round(d/4))
imputations[["DRF"]][[1]]<-complete(imputedDRF)

imputedmean = mice(X.NA, method="mean", m=1) 
imputations[["mean"]][[1]]<-complete(imputedmean)

imputedmean = mice(X.NA, method="cart", m=1) 
imputations[["cart"]][[1]]<-complete(imputedmean)


imputations[["true.data"]][[1]]<-X

methods<-names(imputations)


## Scoring the imputations
scores<-Iscores(imputations = imputations,
        methods = methods,
        m = 1,
        X.NA = X.NA,
        num.proj=200
)


sort(scores)






#### Old ###
############





### Small DRF code example ###
library(drf)
set.seed(1)

n<-1000
d<-5
p<-3


# generate X and a test point x
X<-matrix(rnorm(n*d), nrow=n, ncol=d)
x<-matrix(rnorm(d), nrow=1, ncol=d)

# generate some multivariate Y
Y<-matrix(0, nrow=n, ncol=d)
Y[,1]<-X[,1]^2 + 2*exp(X[,2]) + rnorm(n)
Y[,2] <- Y[,1]^2 + rnorm(n)

# Fit drf
DRF<-drf(X,Y)

# Sample from the estimated conditional distribution of Y|X=x
w<-predict(DRF, newdata=x)$weights[1,]
Yxest<-Y[sample(1:n, size=n, replace=T, prob=w), ]


# Simulate from true conditional distribution
Yx<-matrix(0, nrow=n, ncol=d)
Yx[,1]<-x[,1]^2 + 2*exp(x[,2]) + rnorm(n)
Yx[,2]<-Yx[,1]^2 + rnorm(n)

plot(Yx, col="darkblue")
points(Yxest, col="red")







