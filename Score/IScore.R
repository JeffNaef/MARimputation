Iscores_new<-function(X, N=100, imputations=NULL, num.trees = 1000, min.node.size=5, score="drf", imputationfuncs=NULL, maxlength=NULL, projections=FALSE,...){
  
  ## X: Data with NAs
  ## N: Number of samples from imputation distribution H
  ## imputations: Either NULL or a list of imputations for the methods considered, each imputed X saved as 
  ##              imputations[[method]], whereby method is a string
  ## num.trees: num.trees used for estimating H_{X_j \mid X_{-j}} by DRF. Not relevant for the m-I-Score
  ## min.node.size: min.node.size used for estimating H_{X_j \mid X_{-j}} by DRF. Not relevant for the m-I-Score
  
  
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
  
  methods<-names(imputationfuncs)
  
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
      if (is.null(imputations)){
        # If there is no prior imputation
        tmp<-Iscores_new_perimp(X, Ximp=NULL, N=N, num.trees = num.trees, min.node.size=min.node.size,score=score, maxlength=maxlength, projections=projections,...)
        
      }else{
        
        tmp<-Iscores_new_perimp(X, Ximp=imputations[[method]], N=N, num.trees = num.trees, min.node.size=min.node.size,score=score, maxlength=maxlength, projections=projections,...)
        
        
      }
      
       score_all[[method]] <- tmp  
    }else{
      
      
      tmp<-Iscores_new_perimp(X, Ximp=imputations[[method]], N=N, num.trees = num.trees, min.node.size=min.node.size,score=score, imputationfunc=imputationfuncs[[method]], maxlength=maxlength, projections=projections, ...)
      score_all[[method]] <- tmp  
      
    }
    
  }
  
  return(score_all)
  
}


Iscores_new_perimp <- function(X, Ximp, N=100, num.trees = 1000, min.node.size=5, score="mIScore2", imputationfunc=NULL, maxlength=NULL, projections=FALSE,...){
  
  if (is.null(Ximp)){
    # Impute, maxit should not be 1 here!
    Ximp<-imputationfunc[[2]](X=X  , m=1, method= imputationfunc[[1]])[[1]]
  }
  
  
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
    
    
    if (n1 < 10){
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
        
        # generate h(x_{-j}), careful maxit should not be 1 here!
        HXmj<-imputationfunc[[2]](X=X[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]  , m=1, method= imputationfunc[[1]])[[1]]
        
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
        
      }else if (score=="mIScore"){
        
        
        # Problem with Categorical variables!!
        
        # Train DRF on imputed data
        Xtrain<-Ximp1[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]
        Ytrain<-Ximp1[,j, drop=F]
        
        # Evaluate on observed data
        Xtest <- Ximp0[,!(colnames(Ximp0) %in% j) &  (colnames(Ximp0) %in% indexfull), drop=F]
        Ytest <-Ximp0[,j, drop=F]
        
        Xartificial<-cbind(c(rep(NA,nrow(Ytest)),c(Ytrain)),rbind(Xtest, Xtrain)   )
        colnames(Xartificial)<-c(colnames(Ytrain), colnames(Xtrain))
        
        Imputationlist<-imputationfunc[[2]](X=Xartificial  , m=N, method= imputationfunc[[1]], maxit=1)
        
        Ymatrix<-do.call(cbind, lapply(Imputationlist, function(x)  x[1:nrow(Ytest),1]  ))
        
        scorej[i] <- -mean(sapply(1:nrow(Ytest), function(j)  { crps_sample(y = Ytest[j,], dat = Ymatrix[j,]) }))
        
        
      }else if (score=="mIScore2"){
        
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
        
        
        
        
        
        # generate h(x_{-j}), careful maxit should not be 1 here!
        HXmj<-imputationfunc[[2]](X=X[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% A), drop=F]  , m=1, method= imputationfunc[[1]])[[1]]
        
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
        

        ##reverse the ordering of the mice imputation here!!! Will likely make it much faster, as 
        ##we always have all columns fully observed except the ones we want to impute by construction
        ## This is weird it shouldn't take that long!
        Imputationlist<-imputationfunc[[2]](X=Xartificial  , m=N, method= imputationfunc[[1]],maxit = 1)
        
        Ymatrix<-do.call(cbind, lapply(Imputationlist, function(x)  x[1:nrow(Ytest),1]  ))
        
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

