#############################
### 2d imputation example ###
#############################


n<-1000
library(data.table)
library(mlbench)
library(mice)

path<-getwd()
source(paste0(path,"/mice_DRF.R"))


#dat<-data.table(l=sample(x=c(0,1), size=n, replace=T),Y1=rnorm(n), Y2=rnorm(n), X1=rnorm(n), X2=rnorm(n))
dat<- mlbench::mlbench.spirals(n, sd = 0.05)$x


tmp<-ampute(dat, prop=0.3)
X.NA<-tmp$amp

imputations<-list()

# New imputation method
imputedDRF <- blockmiceDRF(X.NA, blocksize=1) #mice(X.NA, method="myfunc", num.features=10, m=1, blocks = blocklist ) ## lets be honest, "cart" will do the job
imputations[["DRF"]][[1]]<-complete(imputedDRF)

imputedcart = mice(X.NA, method="cart", m=1) ## lets be honest, "cart" will do the job
imputations[["cart"]][[1]]<-complete(imputedcart)

imputations[["true.data"]][[1]]<-dat
X<-dat


# Just plotting the results
###############################################
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
  Axis(side=1, at=c(-1, 0, 1), labels=c(-1, 0, 1), cex.axis=2,padj=0.6)
  points(imp,pch=0,col="gray40",cex=1.2,cex.lab=1,font.lab=1)
}
###############################################


#################################
## Higher Dimensional Example ### 
#################################



library(mvtnorm)
library(missMDA)
library(Iscores)
library(mice)

path<-getwd()
source(paste0(path,"/mice_DRF.R"))


n<-1000
d<-200
set.seed(1)



## Some Sigma matrix
Sigma<-matrix(0, nrow=d, ncol=d)

for (i in 1:d){
  
  for (j in 1:d){
    
    Sigma[i,j]=0.7^(abs(i-j))
  }
}


X<-rmvnorm(n=n, sigma=Sigma )#matrix(rnorm(n*d), nrow=n, ncol=d)





X.NA<-ampute(X, prop=0.2)$amp


imputations<-list()

imputedDRF<-blockmiceDRF(X.NA, blocksize=round(d/4))
Ximputed<-complete(imputedDRF)



# Just plotting the results
###############################################
## Plot some random two components
ind.candidates <- which(!complete.cases(X.NA))
ind.candidates <- sort(ind.candidates)

# The farther way i,j the smaller the correlation.
par(mfrow=c(1,3))

i<-1
j<-2


plot(X.NA[-ind.candidates,c(i,j)],pch=19,col="gray87",cex=1.2,yaxt="n",
     xlab = "", ylab="" ,font.lab=1,cex.main=2, xaxt="n")
points(Ximputed[ind.candidates,c(i,j)],pch=0,col="gray40",cex=1.2,cex.lab=1,font.lab=1,
)


i<-1
j<-5


plot(X.NA[-ind.candidates,c(i,j)],pch=19,col="gray87",cex=1.2,yaxt="n",
     xlab = "", ylab="" ,font.lab=1,cex.main=2, xaxt="n")
points(Ximputed[ind.candidates,c(i,j)],pch=0,col="gray40",cex=1.2,cex.lab=1,font.lab=1,
)

i<-1
j<-10


plot(X.NA[-ind.candidates,c(i,j)],pch=19,col="gray87",cex=1.2,yaxt="n",
     xlab = "", ylab="" ,font.lab=1,cex.main=2, xaxt="n")
points(Ximputed[ind.candidates,c(i,j)],pch=0,col="gray40",cex=1.2,cex.lab=1,font.lab=1,
)
###############################################
