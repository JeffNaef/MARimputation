d <- 10  # number of dimensions
m <- 50  # number of samples from multivariate forecast distribution

# parameters for multivariate normal example
mu0 <- rep(0, d)
mu <- rep(1, d)
S0 <- S <- diag(d)
S0[S0==0] <- 0.2
S[S==0] <- 0.1

# generate samples from multivariate normal distributions
obs <- drop(mu0 + rnorm(d) %*% chol(S0))
fc_sample <- t(replicate(m, drop(mu + rnorm(d) %*% chol(S))))

# compute Energy Score
es_sample(y = obs, dat = t(fc_sample))


# First part
# sqrt(colSums((fc_sample - matrix(obs, nrow=nrow(fc_sample), ncol=ncol(fc_sample), byrow = T))^2))

(firstpart<-mean( sapply(1:m, function(j)  norm(fc_sample[j,, drop=F]- obs,type="F")    )  ))

##same as 

scoringRules:::esC_xy(obs, t(fc_sample), rep(1/m,m))

# Second part
#secondpart<- mean( sapply(1:(m-1), function(j)  norm(fc_sample[,j+1, drop=F]- fc_sample[,j, drop=F],type="F")    ))

D<-matrix(NaN, nrow=m, ncol=m)
for (i in 1:m){
  
  for (j in 1:m){
    
    D[i,j]<- norm(fc_sample[,j, drop=F] - fc_sample[,i, drop=F], type="F")
    
  }
  
}



w<-rep(1/m,m)

secondpart<-w%*%D%*%w

scoringRules:::esC_xx(fc_sample, rep(1/m,m))

