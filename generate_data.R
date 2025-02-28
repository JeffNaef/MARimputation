genData <- function(dataset = "synthetic1", n = 5000, p = 10, meanShift = 0.8, sdShift = 1, std=1) {
  
  if (dataset == "meanshift") {
    x <- runif(n,-1,1)
    y <- rnorm(n, mean = meanShift*(x > 0))
    X <- matrix(runif(n * (p-1)), ncol = p-1)
    X <- cbind(x,X)
    return(list(y=y,X=X))
  } else if (dataset == "sdshift") {
    x <- runif(n,-1,1)
    y <- rnorm(n, sd = 1 + sdShift*(x > 0))
    X <- matrix(runif(n * (p-1)), ncol = p-1)
    X <- cbind(x,X)
    return(list(y=y,X=X))
  }else if (dataset == "motivatingexample") {
    x1 <- runif(n,-1,1)
    x2 <- runif(n,-1,1)
    y <- matrix(rnorm(n,mean = meanShift*(x1 > 0), sd = 1 + sdShift*(x2 > 0)), nrow=n, ncol=1)
    X <- cbind(x1,x2)
    return(list(y=y,X=X)) 
  } else if (dataset == "distshift") {
    x <- runif(n,-1,1)
    y <- ifelse(x >= 0, rexp(n = n, 1), rnorm(n, 1, 1))
    X <- matrix(runif(n * (p-1)), ncol = p-1)
    X <- cbind(x,X)
    return(list(y=y,X=X))
  } else if (dataset == "bivariatesynthetic"){
    x1 <- runif(n,0,1)
    x2 <- runif(n,0,1)
    y1 <- runif(n,x1,x1+1)
    y2 <- runif(n,0,x2)
    X <- matrix(runif(n * (p-2)), ncol = p-2)
    X <- cbind(x1,x2,X)
    y<- cbind(y1,y2)
    return(list(y=y,X=X))
  } else if (dataset == "copulasynthetic" ){
    X  <- matrix(runif(n*p, min = -1, max = 1), ncol = p)
    d<-2
    gen = function(xx) {
      # copula
      normCop <- normalCopula(param=c(xx[1]), dim = 2)
      # margins
      paramMargins = list(list(mean = 0, sd = 0.5), list(mean = 0, sd = 0.5))
      mdvNorm <- mvdc(copula=normCop, margins=c("norm", "norm"), paramMargins=paramMargins)
      # gen
      c(rMvdc(n = 1, mvdc = mdvNorm), rnorm(d-2))
    }
    y <- t(apply(X, 1, gen))
    names = c()
    for(i in 1:d){
      names = c(names, paste0('Y', i))
    }
    colnames(y) <- names
    
    return(list(y=y,X=X))
  }else if (dataset=="GP"){
    
    X  <- matrix(runif(n*p, min = 0, max = 1), ncol = p)
    d<-10
    
    # Create a grid of points
    t <- seq(-5, 5, length.out = d)
    y<-sapply(1:n, function(i){
      # Set the parameters
      lengthscale <- X[i,2]  # lengthscale parameter for the RBF kernel depending on X_2
      
      # Compute the covariance matrix using the RBF kernel
      K <-   kernelMatrix(rbfdot(sigma =lengthscale), t, y = t)
      
      # Simulate from the Gaussian process
      rmvnorm(1, mean = rep(X[i,1],d), sigma = K)
    })
    ## Plot the simulation
    #plot(x, y, type = "l", main = "Gaussian Process Simulation with Gaussian Kernel")
    
    return(list(y=t(y),X=X))
    
  }else if (dataset=="real_wagedata") {
    
    load("~/GitHub/DRFvarimporance/applications/wage_data/data/datasets/wage_benchmark.Rdata")
    
    index<-sample(1:nrow(X), size = n,replace = F)
    
    
    X<-cbind(Y[index,],X[index,])
    # X<-X[index,]
    # Y<-Y[index,]
    # 
    # X<-cbind(X,Y[,"male"])
    # colnames(X)[ncol(X)]<-"male"
    # Y<-Y[,1, drop=F]
    
    
    return(X)
    
    
  } else if(dataset=="real_birthdata") {
    
    load("~/GitHub/DRFvarimporance/applications/births_data/data/datasets/births_benchmark.Rdata")
    
    #load("~/GitHub/DRFvarimporance/applications/births_data/data/datasets/births_benchmark2.Rdata")
    
    index<-sample(1:nrow(X), size = n,replace = F)
    
    X<-cbind(Y[index,],X[index,])
    #Y<-Y[index,]
    
    return(X)
    
  } else if(dataset=="real_airquality") {
    
    load(paste0(getwd(), "/datasets/air.Rdata"))
    
    
    index<-sample(1:nrow(X), size = n,replace = F)
    
    X<-cbind(Y[index,],X[index,-12])
    #X<-X[index,]
    #Y<-Y[index,]
    
    return(X)
    
  }else if (dataset == "synthetic4") {
    x <- runif(n)
    y <- sin(4*pi*x) + ifelse(x>=.5, rnorm(n), rnorm(n, sd=2))
    X <- matrix(runif(n * p), ncol = p)
    X <- cbind(x,X)
    return(list(y=y,X=X))
  } else if (dataset == "friedman1") {
    d <- mlbench::mlbench.friedman1(n,std)
    return(list(y=d$y,X=d$x))
  } else if (dataset == "friedman2") {
    d <- mlbench::mlbench.friedman2(n,std)
    return(list(y=d$y,X=d$x))
  } else if (dataset == "friedman3") {
    d <- mlbench::mlbench.friedman3(n,std)
    return(list(y=d$y,X=d$x))
  } else if (dataset == "Abalone") {
    data("abalone", package = "AppliedPredictiveModeling")
    response <- "Rings"
    abalone[[response]] <- as.numeric(abalone[[response]])
    return(list(y=abalone$Rings, X=model.matrix(~.-Rings-1, data = abalone)))
  } else if (dataset == "Boston") {
    data("BostonHousing2", package = "mlbench")
    response <- "cmedv"
    BostonHousing2[[response]] <- as.numeric(BostonHousing2[[response]])
    return(list(y=BostonHousing2$cmedv,X=model.matrix(~crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + b + lstat-1, data = BostonHousing2)))
  } else if (dataset == "BigMac") {
    data("BigMac2003", package = "alr3")
    response <- "BigMac"
    BigMac2003[[response]] <- as.numeric(BigMac2003[[response]])
    return(list(y=BigMac2003$BigMac,X=model.matrix(~Bread + Rice + FoodIndex + Bus + Apt + TeachGI + 
                                                     TeachNI + TaxRate + TeachHours-1, BigMac2003)))
  } else if (dataset == "Ozone") {
    data("Ozone", package = "mlbench")
    Ozone <- subset(Ozone, complete.cases(Ozone))
    Ozone <- as.data.frame(lapply(Ozone, function(x) {
      x <- x[, drop = TRUE]
      if (is.factor(x)) return(as.ordered(x))
      x
    }))
    response <- "V4"
    Ozone[[response]] <- as.numeric(Ozone[[response]])
    return(list(y=Ozone$V4 ,X=model.matrix(~V1 + V2 + V3 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13-1, Ozone)))
  }
}

