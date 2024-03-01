


###Notes to be deleted
###- Implemented a new sampling for DRF, with sometimes weird results, need to check this
###- The args argument when blocksize is specified does not work at all, need to check this!


Iscores_new<-function(X, imputations, num.trees = 1000, min.node.size=5, score="drf", imputationfuncs=NULL, maxlength=NULL, projections=FALSE,...){

  require(drf)
  require(grf)
  require(kernlab)
  require(Matrix)
  require(scoringRules)
  
  
  numberofmissingbyj<-sapply(1:ncol(X), function(j)  sum(is.na(X[,j]))  )
  print("Number of missing values per dimension:")
  print(paste0(numberofmissingbyj, collapse=",")  )
  
  # imputations: List with
  # imputations$nameofmethod being the imputation method
  # imputations$nameofmethod[[j]] being the jth imputation out of m, of imputation nameofmethod
  
  methods<-names(imputations)
  
  m<-length(imputations[[1]])
  score_all<-list()
  
  for (method in methods) {
    print(paste0("Evaluating method ", method))
   
    
    if (score=="drf"){
    ##Reinstate!! Or make even better
    # tmp<-rep(NaN, m)
    # ##Calculate score for method j
    # for (j in 1:m){
    # tmp[j]<-Iscores_new_perimp(X, imputations[[method]][[j]], num.trees = num.trees, min.node.size=min.node.size,score=score)
    #  
    # }
    
    
    tmp<-Iscores_new_perimp(X, imputations[[method]][[1]], num.trees = num.trees, min.node.size=min.node.size,score=score, maxlength=maxlength, projections=projections,...)
    score_all[[method]] <- tmp  
    }else{
      
      
      tmp<-Iscores_new_perimp(X, imputations[[method]][[1]], num.trees = num.trees, min.node.size=min.node.size,score=score, imputationfuncs[[method]], maxlength=maxlength, projections=projections, ...)
      score_all[[method]] <- tmp  
      
    }
    
  }
  
  return(score_all)
  
}


Iscores_new_perimp <- function(X, Ximp, num.trees = 1000, min.node.size=5, score="drf", imputationfunc=NULL, maxlength=NULL, projections=FALSE,...){

colnames(X) <- colnames(Ximp) <- paste0("X", 1:ncol(X))

args<-list(...)

X<-as.matrix(X)
Ximp<-as.matrix(Ximp)

n<-nrow(X)
p<-ncol(X)

##Step 1: Reoder the data according to the number of missing values
## (least missing first)
numberofmissingbyj<-sapply(1:p, function(j)  sum(is.na(X[,j]))  )
 #X0<-X
 #Ximp0<-Ximp
 #X<-X[,order(numberofmissingbyj, decreasing=T)]
 #Ximp<-Ximp[,order(numberofmissingbyj, decreasing=T)]


## Done in the function
M<-1*is.na(X)
colnames(M) <- colnames(X)

# if (any(colSums(M)==0)){
# 
#   indexfull<-colnames(X)[colSums(M)==0]
# 
# }else{
# 
#   indexfull<-colnames(X)
# 
#   warning("No fully observed variable, using ad-hoc solution")
# 
# }

indexfull<-colnames(X)


# Order first according to most missing values

# Get dimensions with missing values (all other are not important)
dimwithNA<-(colSums(M) > 0)
dimwithNA <- dimwithNA[order(numberofmissingbyj, decreasing=T)]
dimwithNA<-dimwithNA[dimwithNA==TRUE]

if (is.null(maxlength)){maxlength<-sum(dimwithNA) }

if (sum(dimwithNA) < maxlength){
  warning("maxlength was set smaller than sum(dimwithNA)")
  maxlength<-sum(dimwithNA)
}


index<-1:ncol(X)
scorej<-matrix(NA, nrow= min(sum(dimwithNA), maxlength), ncol=1)
weight<-matrix(NA, nrow= min(sum(dimwithNA), maxlength), ncol=1)
i<-0

for (j in names(dimwithNA)[1:maxlength]){
  
  i<-i+1
  
  ## TO DO: 
  ##- Need to do something about class-imbalance!
  ##- Need to weigh the resulting score according to number of obs!
  
  print( paste0("Dimension ", i, " out of ", maxlength )   ) 

      
    require(drf)
    
    Ximp1<-Ximp[M[,j]==1, ]
    Ximp0<-Ximp[M[,j]==0, ]
      
    n1<-nrow(Ximp1)
    n0<-nrow(Ximp0)
      
    
    if (n1 < 200){
      scorej[i]<-NA
      
      warning('Sample size of missing and nonmissing too small for nonparametric distributional regression, setting to NA')
      
    }else{
      
      if (score=="drf"){
        
        
        
        # Problem with Categorical variables!!
        
        # Train DRF on imputed data
        Xtrain<-Ximp1[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]
        Ytrain<-Ximp1[,j, drop=F]
        
        # Evaluate on observed data
        Xtest <- Ximp0[,!(colnames(Ximp0) %in% j) &  (colnames(Ximp0) %in% indexfull), drop=F]
        Ytest <-Ximp0[,j, drop=F]
        
        
        
        if (all(Ytest %in% 0:1)){
          
          if (!all(Ytrain %in% 0:1)){
            
            warning("Imputation Method should impute with binary variables, but has numerical values => Imputation Method should not be used at all. Using quickfix for now")
            
            Ytrain<-(Ytrain > 0.5 )
            
          }
          
          
          fit <- probability_forest(X=Xtrain,
                                    Y=as.factor(Ytrain), num.trees = 2000, min.node.size = 5)
          
          Fhat <- predict(fit, newdata=Xtest  )$predictions
          
          
          scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = sample(c(0,1), size=2000, replace=T, prob=c(Fhat[j,1],Fhat[j,2]) ) ) }))
          
        }else{
          
          #####Fit DRF##########
          fit <- drf(X=Xtrain,
                     Y=Ytrain,  splitting.rule = "FourierMMD", num.trees = 2000, num.features=20, min.node.size = 15, honesty=F)
          
          
          # Predict on the test data
          Fhat <- predict(fit, newdata=Xtest  )$weights
          scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = c(Ytrain), w=Fhat[j,]) }))
          #####Fit DRF##########
        }
        

      }else if (score=="drf2"){
        
        # generate h(x_{-j})
        HXmj<-imputationfunc[[2]](X=X[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]  , m=1, method= imputationfunc[[1]])$imputations[[1]][[1]]
        
        HXmj1<-HXmj[M[,j]==1, ]
        HXmj0<-HXmj[M[,j]==0, ]
        
        
        
        # Train DRF on imputed data
        Xtrain<-HXmj1
        Ytrain<-Ximp1[,j, drop=F]
        
        # Evaluate on observed data
        Xtest <- HXmj0
        Ytest <-Ximp0[,j, drop=F]
        
        
        if (all(Ytest %in% 0:1)){
          
          if (!all(Ytrain %in% 0:1)){
            
            warning("Imputation Method should impute with binary variables, but has numerical values => Imputation Method should not be used at all. Using quickfix for now")
            
            Ytrain<-(Ytrain > 0.5 )
            
          }
          
          
          fit <- probability_forest(X=Xtrain,
                                    Y=as.factor(Ytrain), num.trees = 2000, min.node.size = 5)
          
          Fhat <- predict(fit, newdata=Xtest  )$predictions
          
          
          scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = sample(c(0,1), size=2000, replace=T, prob=c(Fhat[j,1],Fhat[j,2]) ) ) }))
          
        }else{
          
          #####Fit DRF##########
          fit <- drf(X=Xtrain,
                     Y=Ytrain,  splitting.rule = "FourierMMD", num.trees = 2000, num.features=20, min.node.size = 15, honesty=F)
          
          
          # Predict on the test data
          Fhat <- predict(fit, newdata=Xtest  )$weights
          scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = c(Ytrain), w=Fhat[j,]) }))
          #####Fit DRF##########
        }
        
      }else if (score=="mulitpleimp"){

        
        # Problem with Categorical variables!!
        
        # Train DRF on imputed data
        Xtrain<-Ximp1[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]
        Ytrain<-Ximp1[,j, drop=F]
        
        # Evaluate on observed data
        Xtest <- Ximp0[,!(colnames(Ximp0) %in% j) &  (colnames(Ximp0) %in% indexfull), drop=F]
        Ytest <-Ximp0[,j, drop=F]
        
        Xartificial<-cbind(c(rep(NA,nrow(Ytest)),c(Ytrain)),rbind(Xtest, Xtrain)   )
        colnames(Xartificial)<-c(colnames(Ytrain), colnames(Xtrain))
        
        m<-200
        ##reverse the ordering of the mice imputation here!!! Will likely make it much faster, as 
        ##we always have all columns fully observed except the ones we want to impute by construction
        ## This is weird it shouldn't take that long!
        Imputationlist<-imputationfunc[[2]](X=Xartificial  , m=m, method= imputationfunc[[1]])
        
        Ymatrix<-do.call(cbind, lapply(Imputationlist$imputations[[1]], function(x)  x[1:nrow(Ytest),1]  ))
        
        scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = Ymatrix[j,]) }))
        
      
      }else if (score=="mulitpleimp2"){
        
        # Step 1: Impute using all except variable j!
        
        # Problem with Categorical variables!!
        
        
        if (projections==T){
          
          if (!("projectionsize" %in% names(args) )){
            projectionsize <- sample(c(1:ncol(X.NA)), 1)
          }else{
          projectionsize <- args$projectionsize
          }
          
          A<-sample( indexfull[!(colnames(Ximp1) %in% j)], projectionsize  )
          
        }else{
          
          A <- indexfull 
        }
        
        
        
        
        
        # generate h(x_{-j})
        HXmj<-imputationfunc[[2]](X=X[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% A), drop=F]  , m=1, method= imputationfunc[[1]])$imputations[[1]][[1]]
        
        HXmj1<-HXmj[M[,j]==1, ]
        HXmj0<-HXmj[M[,j]==0, ]
        
        
        
        # Train DRF on imputed data
        Xtrain<-HXmj1
        Ytrain<-Ximp1[,j, drop=F]
        
        # Evaluate on observed data
        Xtest <- HXmj0
        Ytest <-Ximp0[,j, drop=F]
        
        Xartificial<-cbind(c(rep(NA,nrow(Ytest)),c(Ytrain)),rbind(Xtest, Xtrain)   )
        colnames(Xartificial)<-c(colnames(Ytrain), colnames(Xtrain))
        
        m<-100
        ##reverse the ordering of the mice imputation here!!! Will likely make it much faster, as 
        ##we always have all columns fully observed except the ones we want to impute by construction
        ## This is weird it shouldn't take that long!
        Imputationlist<-imputationfunc[[2]](X=Xartificial  , m=m, method= imputationfunc[[1]])
        
        Ymatrix<-do.call(cbind, lapply(Imputationlist$imputations[[1]], function(x)  x[1:nrow(Ytest),1]  ))
        
        scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = Ymatrix[j,]) }))
        
      }
    }
  
    

  #weighted
  weight[i]<-(n1/n)*(n0/n)
  
}

scorelist<-c(scorej)
names(scorelist) <- names(dimwithNA)[1:maxlength]
weightlist<-c(weight)
names(weightlist) <- names(dimwithNA)[1:maxlength]

weightedscore<-scorej*weight/(sum(weight, na.rm=T))

## Weight the score according to n0/n * n1/n!!
return( list(score= sum(weightedscore, na.rm=T), scorelist=scorelist, weightlist=weightlist)  )
}






blockmiceDRF<-function(data, blocksize=1, blocklist=NULL,num.trees=10, num.features=10, honesty=T, min.node.size=5, m=1, robust=F,... ){
  
  n<-nrow(data)
  d<-ncol(data)
  
  if (is.null(colnames(data))){
    colnames(data)<-paste0("X.",c(1:d))
  }
  
  if (is.null(blocklist) | length(blocklist)==0){
    
    if (blocksize > d){
      warning("blocksize too large, taking blocksize = 1")
      blocksize<-1
    }else{
      
      # build arbitrary groups of data
      sequence<-seq(1,d,by=blocksize)
      # Remove the last point
      sequence<-sequence[-length(sequence)]
      
      blocklist<-lapply(sequence, function(j) {colnames(data)[j:(j+blocksize-1)]})
      
      # Use whathever is left over => Maybe better to have an overlap!
      blocklist[[length(blocklist)+1]]<-colnames(data)[((sequence[length(sequence)])+blocksize):d]
      
    }
    
  }
  
  if (robust==T){
    warning("robust is set to True")
  }
  
  ## Make sure the blocks actually contain NA (important to not crash the function!)
  listindex<-lapply(blocklist, function(x) any(is.na(data[,x]))  )
  blocklist <- blocklist[ do.call(rbind,listindex)]
  
  return(mymice(data, method="drffuncmult", num.trees=num.trees, num.features=num.features, min.node.size=min.node.size, m=m, blocks = blocklist, robust=robust, ... ))
  
# ### This needs to be improved!
#   if (blocksize > 1){
#   return(mymice(data, method="drffuncmult", num.trees=num.trees, num.features=num.features, min.node.size=min.node.size, m=m, blocks = blocklist, robust=robust, args ))
#   }else{
#     return(mymice(data, method="drffuncmult", num.trees=num.trees, num.features=num.features, min.node.size=min.node.size, m=m, blocks = blocklist, robust=robust, ... ))
#     #return(mice(data, method="drffuncuniv", num.trees=num.trees, num.features=num.features, min.node.size=min.node.size, m=m, blocks = blocklist, robust=robust, ... ))
#   }

}


mice.impute.drffuncmult<- function(data, type=NULL, min.node.size=5, robust=F, num.trees=10,num.features=20, ...) #function (y, ry, x, wy = NULL, minbucket = 5, cp = 1e-04,min.node.size=5, robust=F , ...) 
{
  require(drf)
  
  type=NULL
  
  imputeindex<-names( which( (colSums(apply(data,2, function(x) is.na(x))) > 0)) )
  
  #Careful, we might not have enough yobs!!
  yobs<- data[complete.cases(data[,imputeindex]), imputeindex, drop=F]
  xobs<- data[complete.cases(data[,imputeindex]), !(colnames(data)  %in% imputeindex), drop=F]
  xmis<- data[!complete.cases(data[,imputeindex]), !(colnames(data)  %in% imputeindex), drop=F]
  
  if (robust==T){
    #xobs<-t(apply(xobs,1, function(x){x/sqrt(sum(x^2))}))
    #xmis<-t(apply(xmis,1, function(x){x/sqrt(sum(x^2))}))
    
    ## Linear Version
    # linreg<-lm(yobs~ xobs )
    # yobs<- yobs-cbind(rep(1,nrow(xobs)), xobs)%*%linreg$coefficients
    # fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    # DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    # impute0 <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(1))  # sample one observation per xmis
    # impute <- impute0 + cbind(rep(1,nrow(xmis)), xmis)%*%linreg$coefficients 
    
    
    ## RF Version
    RF<-drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=5, splitting.rule="CART")
    meanpred0<-predict(RF, newdata=xobs, functional="mean")
    yobs<- yobs-meanpred0
    fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    impute0 <- t(vapply(1:nrow(xmis), function(s) unlist(yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ]), numeric(ncol(yobs))))  # sample one observation per xmis
    meanpred<-predict(RF, newdata=xmis, functional="mean")
    
    if (nrow(meanpred)==ncol(impute0)){
    impute0<-t(impute0)  
    }
    
    impute <- impute0 + meanpred
    
    
    
  }else{
    

    
    fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    #impute <- t(vapply(1:nrow(xmis), function(s) unlist(yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ]), numeric(ncol(yobs))))  # sample one observation per xmis
    
    ## Much faster, but is it correct??? Seem to get ambigious resutls
    indexes <- apply(DRFw, 1, function(probs) sample(1:nrow(yobs), size=1, replace=T, prob=probs))
    impute <- yobs[indexes, , drop=F]
    
    
  }
  
  dataimp<-data
  
  # in case we need to transpose
  if (nrow(dataimp[!complete.cases(data[,imputeindex]), imputeindex, drop=F])==ncol(impute)){
    impute <- t(impute)
    
  }
  
  
  dataimp[!complete.cases(data[,imputeindex]), imputeindex] <- impute 
  ## Note: With these we overwrite true values, but its not a problem, this is just an auxiliary filtering step
  
  imputelist<-list()
  for (j in imputeindex){
    
    # Here we take out the actual imputed values
    imputelist[[j]] <-  dataimp[!complete.cases(data[,j]),j]
    
  }
  
  return(imputelist)
  
}




mice.impute.drffuncuniv<- function (y, ry, x, wy = NULL, minbucket = 5, cp = 1e-04,min.node.size=5, robust=F , ...)
{
  #install.on.demand("drf", ...)
  require(drf)
  if (is.null(wy)) {
    wy <- !ry
  }
  minbucket <- max(1, minbucket)
  if (dim(x)[2] == 0) {
    x <- cbind(x, 1)
    dimnames(x) <- list(NULL, "int")
  }
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- as.matrix(y[ry])

  args<-list(...)

  if ("num.trees" %in% names(args)){
    num.trees<-args$num.trees
  }else{
    num.trees=10
  }

  if ("num.features" %in% names(args)){
    num.features<-args$num.features
  }else{
    num.features=10
  }


  if (robust==T){
    #xobs<-t(apply(xobs,1, function(x){x/sqrt(sum(x^2))}))
    #xmis<-t(apply(xmis,1, function(x){x/sqrt(sum(x^2))}))

    ## Linear Version
    # linreg<-lm(yobs~ xobs )
    # yobs<- yobs-cbind(rep(1,nrow(xobs)), xobs)%*%linreg$coefficients
    # fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    # DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    # impute0 <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(1))  # sample one observation per xmis
    # impute <- impute0 + cbind(rep(1,nrow(xmis)), xmis)%*%linreg$coefficients


    ## RF Version
    #RF<-drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=5, splitting.rule="CART")
    RF<-drf(Y=yobs, X=xobs,num.trees=num.trees, compute.oob.predictions = F, min.node.size=1, splitting.rule="CART")
    meanpred<-predict(RF, newdata=NULL, functional="mean")$mean
    yobs<- yobs-meanpred
    fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    impute0 <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(1))  # sample one observation per xmis
    impute <- impute0 + c(predict(RF, newdata=xmis, functional="mean")$mean)



  }else{


    fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
    DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
    impute <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(1))  # sample one observation per xmis



  }


  # if (!is.factor(yobs)) {
  #   fit <- rpart::rpart(yobs ~ ., data = cbind(yobs, xobs),
  #                       method = "anova", control = rpart::rpart.control(minbucket = minbucket,
  #                                                                        cp = cp, ...))
  #   leafnr <- floor(as.numeric(row.names(fit$frame[fit$where,
  #   ])))
  #   fit$frame$yval <- as.numeric(row.names(fit$frame))
  #   nodes <- predict(object = fit, newdata = xmis)
  #   donor <- lapply(nodes, function(s) yobs[leafnr == s])
  #   impute <- vapply(seq_along(donor), function(s) sample(donor[[s]],
  #                                                         1), numeric(1))
  # }
  #   cat.has.all.obs <- table(yobs) == sum(ry)
  #   if (any(cat.has.all.obs)) {
  #     return(rep(levels(yobs)[cat.has.all.obs], sum(wy)))
  #   }
  #   xy <- cbind(yobs, xobs)
  #   xy <- droplevels(xy)
  #   fit <- rpart::rpart(yobs ~ ., data = xy, method = "class",
  #                       control = rpart::rpart.control(minbucket = minbucket,
  #                                                      cp = cp, ...))
  #   nodes <- predict(object = fit, newdata = xmis)
  #   impute <- apply(nodes, MARGIN = 1, FUN = function(s) {
  #     sample(colnames(nodes), size = 1, prob = s)
  #   })
  # }
  impute
}











generateJackKifeIndices <- function(n = 100, n.rep = 6, n.groups = 2) {
  l <- c(list(1:n),
         unlist(lapply(1:n.rep,
                       function(i)
                       {ids <- cut(sample(1:n),
                                   breaks = n.groups,
                                   labels = FALSE);
                       lapply(1:n.groups, function(j) which(ids!=j))}),
                recursive = FALSE))
  return(l)
}

genDataNoNA_synthetic <- function(n.train = 500,
                                  n.test = 10000,
                                  dataset = "bivariateGaussian",
                                  patterns = matrix(c(NA,NA),nrow=1),
                                  pattern_type = pattern_type,
                                  d=2,
                                  ...) {

  if (nrow(patterns)<1) {
    stop("patterns should contains at least one element.")
  }

  if (dataset == "bivariateGaussian") {

    # description: bivariate Gaussian with a specific covariance matrix
    d <- 2

    # train set
    train <- MASS::mvrnorm(n = n.train,
                           mu = rep(0, d),
                           Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example

    # tests list
    tests <- list()

    # iterate over patterns
    for (i in 1:nrow(patterns)) {

      # unconditional
      if (sum(is.na(patterns[i,]))==2) {

        tests[[i]]  <- MASS::mvrnorm(n = n.test,
                                     mu = rep(0, d),
                                     Sigma = diag(0.3,d,d)+matrix(0.7,d,d)) # naive matrix example

        # second missing
      } else if (is.na(patterns[i,2])) {

        tests[[i]] <- cbind(patterns[i,1], rnorm(n = n.test,
                                                 mean = 0.7*patterns[i,1],
                                                 sd = 1-0.7^2))

        # first missing
      } else if (is.na(patterns[i,1])) {

        tests[[i]] <- cbind(rnorm(n = n.test,
                                  mean = 0.7*patterns[i,2],
                                  sd = 1-0.7^2), patterns[i,2])
      } else {
        stop("one at least should be NA.")
      }
    }

    return(list(train = train,tests = tests, patterns=patterns))

  } else if (dataset == "InverseSmiley") {
    # description: inverse parabola in
    d.alt <- d
    d <- 2
    sigma <- 0.5

    # train
    x <- rnorm(n.train)
    train <- cbind(x, -5*x^2+ rnorm(n.train)*0.5)

    # test list
    tests <- list()

    for (i in 1:nrow(patterns)) {

      if (sum(is.na(patterns[i,]))==d.alt) {
        x <- rnorm(n.test)
        tests[[i]] <- cbind(x, -x^2+1 + rnorm(n.test)*sigma)

      } else if (sum(is.na(patterns[i,]))==0) {
        stop("should contain at least one missing value.")

      } else  {
        condition.on.ind <- which(!is.na(patterns[i,]))
        condition.on.val <- patterns[i,condition.on.ind]
        dimensions.to.fill <- setdiff(1:d,condition.on.ind)

        n <- max(100000, n.test)
        epsilon <- 0.05
        x <- rnorm(n)
        data.proposed <- cbind(x, -x^2+1 + rnorm(n)*sigma)

        ids <- which(data.proposed[,condition.on.ind]<= condition.on.val + epsilon
                     & data.proposed[,condition.on.ind] >= condition.on.val - epsilon)
        if(length(ids)>=n.test){
          ids <- ids[1:n.test]
        }
        data.proposed <- data.proposed[ids,]
        data.proposed[,condition.on.ind] <- matrix(rep(condition.on.val,length(ids)), byrow=T,ncol=length(condition.on.val))
        tests[[i]] <- data.proposed

      }
    }

    if(d.alt >2){
      train <- cbind(train, MASS::mvrnorm(n = n.train, mu = rep(0, d.alt-2), Sigma = diag(d.alt-2)))
    }
    return(list(train = train, tests = tests, patterns=patterns))

  } else if (substr(dataset, start = 1, stop = 13)== "multivariateT" || dataset == "HGH" || dataset == "multivariateGaussian"){
    # description: Jeff's HGH model, can be used in any dimensions


    # params<-list(...)
    # mu<- params$mu
    # lambda <- params$lambda
    # C<-params$C
    # L<-params$L
    mu=runif(d*(d-1)/2, min=-5, max=5)#2*runif(d)-1
    lambda=rep(10,d) # For the Gaussian case, lambda is just ignored
    C=diag(d)
    ##Generate Lower triangular matrix
    a <- 2*runif(d*(d-1)/2)-1
    L <- diag(d)
    i.upr <- which(lower.tri(L, diag = FALSE), arr.ind=TRUE)
    L[i.upr] <- a


    train <- matrix(NA, nrow=n.train, ncol=d)

    for (j in 1:n.train){
      if (dataset == "multivariateGaussian"){
        G<-rep(1,d)
        D <- diag(sqrt(G))
      } else if (dataset == "HGH"){
        G <- sapply(lambda,function(x) {rgamma(1,x, rate = 1)})
        D <- diag(sqrt(G))
      } else if (substr(dataset, start = 1, stop = 13)== "multivariateT"){
        G<-rchisq(1, df=lambda)/lambda
        D <- diag(sqrt(G)^(-1))
      }
      Z <- MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = C)
      train[j,] <- mu + L%*%D%*%Z
    }


    tests <- list()

    # iterate over patterns
    for (i in 1:nrow(patterns)) {

      # unconditional
      if (sum(is.na(patterns[i,]))==d) {
        test <- matrix(NA, nrow=n.test, ncol=d)
        for (j in 1:n.test){
          if (dataset == "multivariateGaussian"){
            G<-rep(1,d)
            D <- diag(sqrt(G))
          } else if (dataset == "HGH"){
            G <- sapply(lambda,function(x) {rgamma(1,x, rate = 1)})
            D <- diag(sqrt(G))
          } else if (substr(dataset, start = 1, stop = 13)== "multivariateT"){
            G<- rchisq(1, df=lambda)/lambda
            D <- diag(sqrt(G)^(-1))
          }
          Z <- MASS::mvrnorm(n = 1, mu = rep(0, d), Sigma = C)
          test[j,] <- mu + L%*%D%*%Z
        }
        tests[[i]] <- test
      }
    }

    return(list(train =train, tests= tests, patterns=patterns))

  }else if (dataset == "spiral"){

    if (!require(mlbench)) {
      stop("mlbench should be installed.")
    }
    # description: Here we take examples directly from the mlbench package

    train <- mlbench::mlbench.spirals(n.train, sd = 0.05)$x
    if(d >2){
      train <- cbind(train, MASS::mvrnorm(n = n.train, mu = rep(0, d-2), Sigma = diag(d-2)))
    }

    tests <- list()

    # iterate over patterns
    for (i in 1:nrow(patterns)) {

      # unconditional
      if (sum(is.na(patterns[i,]))==d) {
        X <- mlbench.spirals(n.test,...)$x
        if( d==2){
          tests[[i]] <- X
        } else{
          tests[[i]] <- cbind(X, MASS::mvrnorm(n = n.test, mu = rep(0, d-2), Sigma = diag(d-2)))
        }
        # second missing
      } else{

        dimensions.to.fill <- which(is.na(patterns[i,]))
        condition.on.ind <- which(!is.na(patterns[i,]))
        condition.on.val <- patterns[i,condition.on.ind]

        n <- 100000
        epsilon <- rep(0.05, length(condition.on.ind))

        X <- mlbench.spirals(n,...)$x
        if( d==2){
          data.proposed <- X
        } else{
          data.proposed <- cbind(X, MASS::mvrnorm(n = n, mu = rep(0, d-2), Sigma = diag(d-2)))
        }

        ids <- which(data.proposed[,condition.on.ind]<= condition.on.val + epsilon
                     & data.proposed[,condition.on.ind] >= condition.on.val - epsilon)
        if(length(ids)>=n.test) ids <- ids[1:n.test]

        data.proposed <- data.proposed[ids,]
        data.proposed[,condition.on.ind] <- matrix(rep(condition.on.val,length(ids)), byrow=T,ncol=length(condition.on.val))
        tests[[i]] <- data.proposed
      }
    }

  } else if(dataset =="swissroll"){

    train <- swissroll(n.train)
    if (d > 2){
      train <- cbind(train, MASS::mvrnorm(n = n.train, mu = rep(0, d-2), Sigma = diag(d-2)))
    }
    tests <- c()
    patterns <- patterns

  }else if(dataset =="four.clouds"){

    sigma <- 0.01
    Sigma <- matrix(c(sigma, 0, 0, sigma), ncol=2)
    mode1 <- mvrnorm(n.train/4, mu=c(1,0),Sigma = Sigma )
    mode2 <- mvrnorm(n.train/4, mu=c(0,1),Sigma = Sigma )
    mode3 <- mvrnorm(n.train/4, mu=c(-1,0),Sigma = Sigma )
    mode4 <- mvrnorm(n.train/4, mu=c(0,-1),Sigma = Sigma )

    train <- rbind(mode1,mode2,mode3,mode4)

    if (d > 2){
      train <- cbind(train, MASS::mvrnorm(n = n.train, mu = rep(0, d-2), Sigma = diag(d-2)))
    }
    tests <- c()
    patterns <- patterns

  }else{
    train <-  eval(parse(text=paste0("mlbench.", dataset,"(", "n=", "n.train)$x")))

    tests<-c()
    patterns<-c()


  }


  return(list(train=train,tests=tests, patterns=patterns))


}


genMask_introexample_MAR <- function(X, pmiss=0.3){

  X.NA <- matrix(NA, ncol=2, nrow=nrow(X))

  X.NA[,1] <- sapply(1:nrow(X),function(ind){

    ifelse(X[ind,2] > -0.3 & X[ind,2] < 0.3, rbinom(1, 1, prob=0.3), 0)

  })

  X.NA[,2] <- sapply(1:nrow(X),function(ind){

    ifelse(X[ind,1] < -0.3 | X[ind,1] > 0.3, rbinom(1, 1, prob=0.3), 0)

  })

  X.NA.final <- X
  X.NA.final[which(X.NA==1)]<-NA

  return(X.NA.final)


}

genMask <- function(X, mech = "MCAR", pmiss = 0.5) {

  create.patterns <- function(pmiss, d){
    patterns <- c()
    nar <- 0
    while(nar != ceiling(d/2)){
      nr.of.0 <- rbinom(prob=pmiss,n = 1, size=d)
      if(nr.of.0 > ceiling(d/2)){
        nr.of.0 <- d-nr.of.0
      }

      pat <- rep(1, d)
      pat[sample(1:d, nr.of.0)]<-0
      patterns <- rbind(patterns, pat)
      patterns <- patterns[!duplicated(patterns), , drop=F]
      nar <- ifelse(is.null(nrow(patterns)),0,nrow(patterns))

    }

    return(patterns)
  }

  n <- nrow(X)
  d <- ncol(X)
  # method = simple deletes any entry in the data matrix with probability pNA
  if(mech == "MCAR"){
    na <- apply(X,2,function(x) sample(c(0,1), size = length(x), prob = c(1-pmiss, pmiss), replace = TRUE))
    # inds.total.na <- which(apply(na, 1, sum)==d)
    X.NA <- X
    X.NA[na==1] <- NA
  }else if (mech %in% c("MAR", "MNAR")){

    res <- NA
    while(length(res)==1) {
      #num.patterns <- 0
      #while(num.patterns != ceiling(pmiss*d)){
      patterns <- create.patterns(pmiss,d)
      # num.patterns <- nrow(patterns)
      #}

      #patterns <- apply(X,2,function(x) sample(c(1,0), size = length(x), prob = c(1-pmiss, pmiss), replace = TRUE))
      #patterns <- patterns[!duplicated(patterns), , drop=F]

      freq <- rep(1/nrow(patterns), nrow(patterns))
      cont <- TRUE
      bycases <- FALSE
      run <- TRUE

      res <- tryCatch(ampute(X,
                             prop = pmiss,
                             patterns = patterns,
                             freq = freq,
                             mech = mech,
                             cont = cont,
                             bycases = bycases,
                             run = run)$amp, error = function(x) NA)
    }
    X.NA <- res
  }


  return(X.NA)
}

doimputation <- function(X.NA, methods, m=m,...){

  require(mice)
  
  counter<-0
  imputations <- list()
  methods0<-methods
  
  args<-list(...)

  indexbinary<-apply(X.NA,2, function(x) all(x[!is.na(x)] %in% 0:1) )


  for (method in methods0){
    
    
    start_time <- Sys.time()
    
    counter<-which(methods==method)

    ## MIPCA imputation
    if (method=="mipca"){

      imp <-tryCatch({
        nbdim <- estim_ncpPCA(X.NA, ncp.min=1) # estimate the number of dimensions to impute
        #imputations[[counter]] <-MIPCA(X.NA, ncp = nbdim$ncp, nboot = m)$res.MI
        #names(imputations)[counter]<-method
        MIPCA(X.NA, ncp = nbdim$ncp, nboot = m)$res.MI

      },error = function(e){
        warning(paste("Method", method, "could not be used for imputation"))
        print(paste("Method", method, "could not be used for imputation"))
        # Delete the method from the list
        return(NA)
      })
      # ,
      # warning=function(e) {
      #   message(paste("Method", method, "could not be used for imputation (warnings)"))
      #   message("Here's the original warning message:")
      #   message(e)
      #   return(NA)
      # })

      if(any(is.na(imp))){
        methods <- methods[-which(methods==method)]
      }else{
        imputations[[method]] <- imp
      }



    } else if (method=="knn"){
      
      # Load DMwR package
      require(DMwR2)
      
      # KNN imputation
      imputations[[method]] <- knnImputation(X.NA)
      
    } else if (method=="amelia"){

      imp <- tryCatch({
        amelia(X.NA, m = m)$imputations

      },error = function(e){
        warning(paste("Method", method, "could not be used for imputation"))
        print(paste("Method", method, "could not be used for imputation"))
        # Delete the method from the list
        return(NA)
      })
      # ,
      # warning=function(e) {
      #   message(paste("Method", method, "could not be used for imputation (warnings)"))
      #   message("Here's the original warning message:")
      #   message(e)
      #   return(NA)
      # })
      #

      if(any(is.na(imp))){
        methods <- methods[-which(methods==method)]
      }else{
        imputations[[method]] <- imp
      }



    } else if (method=="missForest"){
      
      require(missForest)

      imp <- tryCatch({

        lapply(1:m,function(i) missForest(X.NA)$ximp)
        #lapply(1:m,function(i) mymissForest(X.NA)$ximp)

      },error = function(e){
        warning(paste("Method", method, "could not be used for imputation"))
        print(paste("Method", method, "could not be used for imputation"))
        # Delete the method from the list
        return(NA)
      })
      # ,
      # warning=function(e) {
      #   message(paste("Method", method, "could not be used for imputation (warnings)"))
      #   message("Here's the original warning message:")
      #   message(e)
      #   return(NA)
      # })


      if(any(is.na(imp))){
        methods <- methods[-which(methods==method)]
      }else{
        imputations[[method]] <- imp
      }

    } else if (method=="DRF"){
      require(drf)
      
      #blub <- mice(data, method="drffunc", num.trees=5, num.features=50, min.node.size=5, m=m, blocks = blocklist, robust=robust, ... )
      
      if (!("blocksize" %in% names(args) )){
        blocksize=1
        blub <- blockmiceDRF(data=X.NA, blocksize=blocksize, blocklist=NULL, robust=F, num.trees=10, honesty=F, num.features=50, min.node.size=1, m=m, ...)
      }else{
        #stop("This does not work yet!")
        ##Check out https://stackoverflow.com/questions/16321760/r-change-value-of-an-argument-in-ellipsis-and-pass-ellipsis-to-the-other-functi
        blocksize <- args[["blocksize"]]
        args <- args[!(names(args) %in% "blocksize" )] 
        blub <- blockmiceDRF(data=X.NA, blocksize=blocksize, blocklist=NULL, robust=F, num.trees=10, honesty=F, num.features=50, min.node.size=1, m=m)
        #do.call(blockmiceDRF(data=X.NA, blocksize=blocksize, blocklist=NULL, robust=F, num.trees=10, honesty=F, num.features=50, min.node.size=1, m=m), args)
      }
      
      
      #blub <- blockmiceDRF(X.NA, blocksize=blocksize, robust=F, num.trees=10, honesty=F, num.features=20, min.node.size=5, m=m, ...)

      imputations[[method]]<-mice::complete(blub, action="all")
      
    }else if (method=="DRF-lin"){
      require(drf)
      blub <- blockmiceDRF(data=X.NA, blocksize=blocksize, robust=T, num.trees=2, num.features=50, min.node.size=15, m=m,...)
      imputations[[method]]<-mice::complete(blub, action="all")
      
    }else {

      imp <- tryCatch({

        # mice imputation
        blub <- mice(X.NA, method = method, m = m, ...) # visitSequence="arabic"
        mice::complete(blub, action="all")

      },error = function(e){
        warning(paste("Method", method, "could not be used for imputation"))
        print(paste("Method", method, "could not be used for imputation"))
        # Delete the method from the list
        return(NA)
      })
      # ,
      # warning=function(e) {
      #   message(paste("Method", method, "could not be used for imputation (warnings)"))
      #   message("Here's the original warning message:")
      #   message(e)
      #   return(NA)
      # })


      if(any(is.na(imp[[1]]))){
        methods <- methods[-which(methods==method)]
      }else{
        imputations[[method]] <- imp
      }


    }

    end_time <- Sys.time()
    
   
    
    print(paste0("Method ", method, " needed ",  round(end_time-start_time,3), " to impute" ))
    
    
  }


  
  return(list(imputations=imputations, methods=methods))


}

genData_real <- function(dataset= "dengue", n=NULL){


  if (dataset == "dengue"){
    # description: dengue fever in africa.
    # https://www.picostat.com/dataset/r-dataset-package-daag-dengue
    dat <- read.csv(file="dataset_dengue.csv")
    X <- as.matrix(data.frame(na.omit(dat)))
    colnames(X) <-NULL

  }else if (dataset== "Boston"){
    X <- as.matrix(Boston) ## For the moment: take out categorical variables
    colnames(X) <- NULL

  }else if (dataset == "network"){
    nv <- vcount(fblog)
    X <- as_adjacency_matrix(fblog)
    colnames(X) <- NULL
    rownames(X) <- NULL

  }else if (dataset == "iris"){
    X <- iris
    X <- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X) <- NULL

    # X <- as.matrix(model.matrix(~., iris))
    # X <- X[,-1]
    # colnames(X) <- NULL

  }else if (dataset == "wine"){

    X <- read.csv(file="wine.data", header=F)
    X <- as.matrix(X)
    X <- X[,-1]
    colnames(X) <- NULL

    # data("wine")
    # X <- wine
    # X$Class <- NULL
    # X <- as.matrix(X)
    # colnames(X) <- NULL

    # data("wine")
    # wine$Class <- as.factor(wine$Class)
    # X <- as.matrix(model.matrix(~., wine))
    # X <- X[,-1]
    # colnames(X) <- NULL

  }else if(dataset == "wisconsin"){

    dat <- read.csv(file="wisconsin.csv")
    dat<-dat[,-c(1,2,33)]
    X <- as.matrix(dat)
    colnames(X) <-NULL

    # dat <- read.csv(file="wisconsin.csv")
    # dat$diagnosis<-as.factor(dat$diagnosis)
    # dat<-dat[,c(-1,-33)]
    # X <- as.matrix(model.matrix(~., dat))
    # X<-X[,-1]
    # colnames(X) <-NULL

  }else if(dataset == "wine.quality.red"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/
    dat <- read.csv(file="winequality-red.csv", sep=";")
    dat$quality<-NULL
    X <- as.matrix(dat)
    colnames(X) <-NULL

    # dat <- read.csv(file="winequality-red.csv", sep=";")
    # dat$quality<-as.factor(dat$quality)
    # X <- as.matrix(model.matrix(~., dat))
    # X<-X[,-1]
    # colnames(X) <-NULL

  }else if(dataset == "wine.quality.white"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/
    dat <- read.csv(file="winequality-white.csv", sep=";")
    dat$quality<-NULL
    X <- as.matrix(dat)
    colnames(X) <- NULL

    # dat <- read.csv(file="winequality-white.csv", sep=";")
    # dat$quality<-as.factor(dat$quality)
    # X <- as.matrix(model.matrix(~., dat))
    # X<-X[,-1]
    # colnames(X) <- NULL

  }else if(dataset == "yacht"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/00243/
    X <- read.csv(file="yacht_hydrodynamics.data", sep="", header=F)
    X <- as.matrix(X)
    colnames(X) <-NULL

  }else if(dataset == "yeast"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/yeast/
    X <- read.csv(file=paste0(getwd(), "/datasets/yeast.data"), sep="", header=F)
    X <- X[,-c(1, ncol(X))]
    X <- as.matrix(X)
    colnames(X) <-NULL

    # X <- read.csv(file="yeast.data", sep="", header=F)
    # X <- X[,-1]
    # X <- as.matrix(model.matrix(~., X))
    # X<-X[,-1]
    # colnames(X) <-NULL

  }else if(dataset == "DNA"){
    data("DNA")
    dat <- DNA
    X <- as.matrix(model.matrix(~., dat))
    X<-X[,-1]
    colnames(X) <-NULL

  }else if(dataset == "airfoil"){
    X <- read.table(file="airfoil_self_noise.dat")
    colnames(X) <-NULL
    X <-as.matrix(X)

  }else if (dataset =="blood.transfusion"){
    X <- read.csv(file="transfusion.data")
    X<- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="transfusion.data")
    # X[, ncol(X)]<- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset == "concrete.slump"){
    X <- read.csv(file="slump_test.data")
    X <- as.matrix(X)
    X <- X[,-1]
    colnames(X)<- NULL

  }else if (dataset == "connectionist.bench.sonar"){

    X <- read.csv(file="sonar.all-data", header=F)
    X <- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="sonar.all-data", header=F)
    # X.last <- as.factor(X[, ncol(X)])
    # X <- apply(X[,1:(ncol(X)-1)], 2, as.numeric)
    # X <- data.frame(X, X.last)
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL
    #
  }else if (dataset == "connectionist.bench.vowel"){
    X <- read.csv(file="vowel-context.data", sep="", header=F)
    X <- X[,-c(1,2,3, ncol(X))]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="vowel-context.data", sep="", header=F)
    # inds.to.convert <- c(1,2,3, ncol(X))
    # for ( i in inds.to.convert){
    #   X[,i] <- as.factor(X[,i])
    # }
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset == "ecoli"){
    X <- read.csv(file="ecoli.data", sep="", header=F)
    X <- X[,-c(1, ncol(X))]
    X <- as.matrix(X)
    X <- X[,-7]# hihgly cor with var8
    X <- X[,-4]
    colnames(X)<- NULL

    # X <- read.csv(file="ecoli.data", sep="", header=F)
    # X <- data.frame(X)
    # X <- X[,-1]
    # X[,ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset=="glass"){
    X <- read.csv(file="glass.data",header=F)
    X <- X[,-c(1, ncol(X))]
    colnames(X)<- NULL

    # X <- read.csv(file="glass.data",header=F)
    # X <- X[,-1]
    # X <- data.frame(X)
    # X[,ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset=="ionosphere"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/
    X <- read.csv(file="ionosphere.data",header=F)
    X <- X[,-c(1,2, ncol(X))]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="ionosphere.data",header=F)
    # X <- X[,-2]
    # X <- data.frame(X)
    # X[,1] <- as.factor(X[,1])
    # X[,ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset=="libras"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/libras/
    X <- read.csv(file="movement_libras.data", header=F)
    X <- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="movement_libras.data", header=F)
    # X <- data.frame(X)
    # X[,ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL
    #
  }else if (dataset=="parkinsons"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/
    X <- read.csv(file="parkinsons.data", header=T)
    X$status <- NULL
    X$name <- NULL
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv(file="parkinsons.data", header=T)
    # X <- data.frame(X)
    # X <- X[,-1]
    # X$status <- as.factor(X$status)
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset=="planning.relax"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/00230/
    X <- read.csv("plrx.txt", sep="", header=F)
    X <- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv("plrx.txt", sep="", header=F)
    # X <- data.frame(X)
    # X[, ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL
    #
  }else if (dataset=="seeds"){ #https://archive.ics.uci.edu/ml/machine-learning-databases/00236/
    X <- read.csv("seeds_dataset.txt", sep="", header=F)
    X <- X[,-ncol(X)]
    X <- as.matrix(X)
    colnames(X)<- NULL

    # X <- read.csv("seeds_dataset.txt", sep="", header=F)
    # X <- data.frame(X)
    # X[, ncol(X)] <- as.factor(X[,ncol(X)])
    # X <- model.matrix(~., X)
    # X <- X[,-1]
    # colnames(X)<- NULL

  }else if (dataset=="climate.model.crashes"){
    X <- read.table(file="pop_failures.dat", header = TRUE)
    X <- X[,-c(1,ncol(X))]
    X<-as.matrix(X)
    colnames(X)<- NULL

    # X <- read.table(file="pop_failures.dat", header = F)
    # X <- as.data.frame(X)
    # colnames(X) <- NULL
    # X <- X[-1,]
    # inds.to.convert <- setdiff(1:ncol(X),c(1, ncol(X)))
    # blub <- lapply(inds.to.convert, function(i) if(is.factor(X[,i])) {
    #   X[,i]<- droplevels(X[,i])
    #   return(as.numeric(levels(X[,i]))[X[,i]])
    #   }else X[,i])
    # blub <- do.call(cbind, blub)
    # blub <- data.frame(blub)
    # blub <- cbind(blub, X[,1])
    # blub <- cbind(blub, X[,ncol(X)])
    # X <- model.matrix(~., blub)
    # X <- X[,-1]
    # colnames(X)<- NULL
    #
  }else if (dataset == "concrete.compression"){
    X <- read_excel("Concrete_Data.xls")
    X <- as.matrix(X)
    colnames(X)<- NULL

  }else if (dataset == "CASchools"){
    data("CASchools")
    X <- CASchools
    X <- X[,-c(1,2,3,4)]
    X <- as.matrix(X)
    colnames(X) <- NULL


  }else if (dataset=="real_wagedata") {
    
    #load("~/GitHub/DRFvarimporance/applications/wage_data/data/datasets/wage_benchmark.Rdata")
    
    load(paste0(getwd(), "/datasets/wage.Rdata"))
    
    index<-sample(1:nrow(X), size = n,replace = F)
    
    
    X<-cbind(data$Y[index,],data$X[index,])
    # X<-X[index,]
    # Y<-Y[index,]
    # 
    # X<-cbind(X,Y[,"male"])
    # colnames(X)[ncol(X)]<-"male"
    # Y<-Y[,1, drop=F]
    
    
    return(X)
    
    
  } else if(dataset=="real_birthdata") {
    
    #load("~/GitHub/DRFvarimporance/applications/births_data/data/datasets/births_benchmark.Rdata")
    load(paste0(getwd(), "/datasets/births.Rdata"))
    #load("~/GitHub/DRFvarimporance/applications/births_data/data/datasets/births_benchmark2.Rdata")
    
    index<-sample(1:nrow(data$X), size = n,replace = F)
    
    X<-cbind(data$Y[index,],data$X[index,])
    #Y<-Y[index,]
    
    return(X)
    
  } else if(dataset=="real_airquality") {
    
    load(paste0(getwd(), "/datasets/air.Rdata"))
    
    
    index<-sample(1:nrow(data$X), size = n,replace = F)
    
    ##Need to remove one coumn for Land.use and one for Location.Setting
    #X<-cbind(data$Y[index,],data$X[index, !(colnames(data$X)%in% c("Land.Use_RESIDENTIAL", "Location.Setting_URBAN AND CENTER CITY")  ) ])
    
    X<-cbind(data$Y[index,],data$X[index,])
    
    #X<-X[index,]
    #Y<-Y[index,]
    
    return(X)
    
  }else if(dataset=="real_spam") {
    
    X<-as.matrix(read.csv(paste0(getwd(), "/datasets/spam.csv")))
    
    #"C:/Users/Jeff/OneDrive/Today/MAR_project/Code/datasets/spam.csv"
  
    #index<-sample(1:nrow(data$X), size = n,replace = F)
    
    #X<-as.matrix(spam)
    #X<-X[index,]
    #Y<-Y[index,]
    
    return(X)
    
  }
  
  
  
  return(X)

}


###Continue here adapting this and playing around!!!


mymice <- function (data, m = 5, method = NULL, predictorMatrix, ignore = NULL, 
                    where = NULL, blocks, visitSequence = NULL, formulas, blots = NULL, 
                    post = NULL, defaultMethod = c("pmm", "logreg", "polyreg", 
                                                   "polr"), maxit = 5, printFlag = TRUE, seed = NA, data.init = NULL, 
                    ...) 
{
  call <- match.call()
  if (!is.na(seed)) 
    set.seed(seed)
  data <- mice:::check.dataform(data)
  m <- mice:::check.m(m)
  mp <- missing(predictorMatrix)
  mb <- missing(blocks)
  mf <- missing(formulas)
  if (mp & mb & mf) {
    blocks <- mice:::make.blocks(colnames(data))
    predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
    formulas <- mice:::make.formulas(data, blocks)
  }
  if (!mp & mb & mf) {
    predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
                                                    data)
    blocks <- mice:::make.blocks(colnames(predictorMatrix), partition = "scatter")
    formulas <- mice:::make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (mp & !mb & mf) {
    blocks <- mice:::check.blocks(blocks, data)
    predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
    formulas <- mice:::make.formulas(data, blocks)
  }
  if (mp & mb & !mf) {
    formulas <- mice:::check.formulas(formulas, data)
    blocks <- mice:::construct.blocks(formulas)
    predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & mf) {
    blocks <- mice:::check.blocks(blocks, data)
    z <- mice:::check.predictorMatrix(predictorMatrix, data, blocks)
    predictorMatrix <- z$predictorMatrix
    blocks <- z$blocks
    formulas <- mice:::make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (!mp & mb & !mf) {
    formulas <- mice:::check.formulas(formulas, data)
    predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
                                                    data)
    blocks <- mice:::construct.blocks(formulas, predictorMatrix)
    predictorMatrix <- mice:::make.predictorMatrix(data, blocks, 
                                                   predictorMatrix)
  }
  if (mp & !mb & !mf) {
    blocks <- mice:::check.blocks(blocks, data, calltype = "formula")
    formulas <- mice:::check.formulas(formulas, blocks)
    predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & !mf) {
    blocks <- mice:::check.blocks(blocks, data)
    formulas <- mice:::check.formulas(formulas, data)
    predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
                                                    data, blocks)
  }
  chk <- mice:::check.cluster(data, predictorMatrix)
  where <- mice:::check.where(where, data, blocks)
  user.visitSequence <- visitSequence
  visitSequence <- mice:::check.visitSequence(visitSequence, data = data, 
                                              where = where, blocks = blocks)
  predictorMatrix <- mice:::edit.predictorMatrix(predictorMatrix = predictorMatrix, 
                                                 visitSequence = visitSequence, user.visitSequence = user.visitSequence, 
                                                 maxit = maxit)
  method <- mice:::check.method(method = method, data = data, where = where, 
                                blocks = blocks, defaultMethod = defaultMethod)
  post <- mice:::check.post(post, data)
  blots <- mice:::check.blots(blots, data, blocks)
  ignore <- mice:::check.ignore(ignore, data)
  state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "", 
                             out = "")
  setup <- list(method = method, predictorMatrix = predictorMatrix, 
                visitSequence = visitSequence, post = post)
  setup <- mice:::edit.setup(data, setup, ...)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  nmis <- apply(is.na(data), 2, sum)
  imp <- mice:::initialize.imp(data, m, ignore, where, blocks, visitSequence, 
                               method, nmis, data.init)
  from <- 1
  to <- from + maxit - 1
  q <- mysampler(data, m, ignore, where, imp, blocks, method, 
                 visitSequence, predictorMatrix, formulas, blots, post, 
                 c(from, to), printFlag, ...)
  if (!state$log) 
    loggedEvents <- NULL
  if (state$log) 
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  midsobj <- list(data = data, imp = q$imp, m = m, where = where, 
                  blocks = blocks, call = call, nmis = nmis, method = method, 
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence, 
                  formulas = formulas, post = post, blots = blots, ignore = ignore, 
                  seed = seed, iteration = q$iteration, lastSeedValue = get(".Random.seed", 
                                                                            envir = globalenv(), mode = "integer", inherits = FALSE), 
                  chainMean = q$chainMean, chainVar = q$chainVar, loggedEvents = loggedEvents, 
                  version = packageVersion("mice"), date = Sys.Date())
  oldClass(midsobj) <- "mids"
  if (!is.null(midsobj$loggedEvents)) {
    warning("Number of logged events: ", nrow(midsobj$loggedEvents), 
            call. = FALSE)
  }
  midsobj
}


mysampler <- function (data, m, ignore, where, imp, blocks, method, visitSequence, 
                       predictorMatrix, formulas, blots, post, fromto, printFlag, 
                       ...) 
{
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  r <- !is.na(data)
  chainMean <- chainVar <- mice:::initialize.chain(blocks, maxit, 
                                            m)
  if (maxit < 1) 
    iteration <- 0
  if (maxit >= 1) {
    if (printFlag) {
      cat("\n iter imp variable")
    }
    for (k in from:to) {
      iteration <- k
      for (i in seq_len(m)) {
        if (printFlag) {
          cat("\n ", iteration, " ", i)
        }
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            y <- data[, j]
            ry <- r[, j]
            wy <- where[, j]
            data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], 
                                            i]
          }
        }
        for (h in visitSequence) {
          ct <- attr(blocks, "calltype")
          calltype <- ifelse(length(ct) == 1, ct[1], 
                             ct[h])
          b <- blocks[[h]]
          if (calltype == "formula") 
            ff <- formulas[[h]]
          else ff <- NULL
          type <- predictorMatrix[h, ]
          user <- blots[[h]]
          theMethod <- method[h]
          empt <- theMethod == ""
          univ <- F
          mult <- T
          pass <- !empt && mice:::is.passive(theMethod) && length(blocks[[h]]) == 
            1
          if (printFlag & !empt) 
            cat(" ", b)
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(it = k, im = i, dep = h, meth = theMethod, 
                           log = oldstate$log)
          assign("state", newstate, pos = parent.frame(), 
                 inherits = TRUE)
          if (univ) {
            for (j in b) {
              imp[[j]][, i] <- sampler.univ(data = data, 
                                            r = r, where = where, type = type, formula = ff, 
                                            method = theMethod, yname = j, k = k, 
                                            calltype = calltype, user = user, ignore = ignore, 
                                            ...)
              data[(!r[, j]) & where[, j], j] <- imp[[j]][(!r[, 
                                                              j])[where[, j]], i]
              cmd <- post[j]
              if (cmd != "") {
                eval(parse(text = cmd))
                data[(!r[, j]) & where[, j], j] <- imp[[j]][(!r[, 
                                                                j])[where[, j]], i]
              }
            }
          }
          if (mult) {
            mis <- !r
            mis[, setdiff(colnames(data), b)] <- FALSE
            data[mis] <- NA
            fm <- paste("mice.impute", theMethod, sep = ".")
            if (calltype == "formula") {
              imputes <- do.call(fm, args = list(data = data, 
                                                 formula = ff, ...))
            }
            else if (calltype == "type") {
              ###Problem: observed parts for the batch can be empty!!
              imputes <- do.call(fm, args = list(data = data, 
                                                 type = type, ...))
              
              ###Continue checking from here!!
            }
            else {
              stop("Cannot call function of type ", calltype, 
                   call. = FALSE)
            }
            if (is.null(imputes)) {
              stop("No imputations from ", theMethod, 
                   h, call. = FALSE)
            }
            for (j in names(imputes)) {
              imp[[j]][, i] <- imputes[[j]]
              data[!r[, j], j] <- imp[[j]][, i]
            }
          }
          if (pass) {
            for (j in b) {
              wy <- where[, j]
              ry <- r[, j]
              imp[[j]][, i] <- model.frame(as.formula(theMethod), 
                                           data[wy, ], na.action = na.pass)
              data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], 
                                              i]
            }
          }
        }
      }
      k2 <- k - from + 1L
      if (length(visitSequence) > 0L) {
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            if (!is.factor(data[, j])) {
              chainVar[j, k2, ] <- apply(imp[[j]], 2L, 
                                         var, na.rm = TRUE)
              chainMean[j, k2, ] <- colMeans(as.matrix(imp[[j]]), 
                                             na.rm = TRUE)
            }
            if (is.factor(data[, j])) {
              for (mm in seq_len(m)) {
                nc <- as.integer(factor(imp[[j]][, mm], 
                                        levels = levels(data[, j])))
                chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
                chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
              }
            }
          }
        }
      }
    }
    if (printFlag) {
      r <- get("loggedEvents", parent.frame(1))
      ridge.used <- any(grepl("A ridge penalty", r$out))
      if (ridge.used) {
        cat("\n * Please inspect the loggedEvents \n")
      }
      else {
        cat("\n")
      }
    }
  }
  list(iteration = maxit, imp = imp, chainMean = chainMean, 
       chainVar = chainVar)
}


