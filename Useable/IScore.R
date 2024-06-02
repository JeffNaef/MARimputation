Iscores_new<-function(X, N=50,  imputationfuncs=NULL, imputations=NULL, maxlength=NULL,...){
  
  ## X: Data with NAs
  ## N: Number of samples from imputation distribution H
  ## imputationfuncs: A list of functions, whereby each imputationfuncs[[method]] is a function that takes the arguments
  ## X,m and imputes X m times using method: imputations= imputationfuncs[[method]](X,m).
  ## imputations: Either NULL or a list of imputations for the methods considered, each imputed X saved as 
  ##              imputations[[method]], whereby method is a string
  ## maxlength: Maximum number of variables X_j to consider, can speed up the code
  
  
  require(Matrix)
  require(scoringRules)
  
  
  numberofmissingbyj<-sapply(1:ncol(X), function(j)  sum(is.na(X[,j]))  )
  print("Number of missing values per dimension:")
  print(paste0(numberofmissingbyj, collapse=",")  )

  methods<-names(imputationfuncs)

  score_all<-list()
  
  for (method in methods) {
    print(paste0("Evaluating method ", method))
    
    
    # }
    if (is.null(imputations)){
      # If there is no prior imputation
      tmp<-Iscores_new_perimp(X, Ximp=NULL, N=N, imputationfunc=imputationfuncs[[method]], maxlength=maxlength,...)
      score_all[[method]] <- tmp  
      
      
    }else{
      
      tmp<-Iscores_new_perimp(X, Ximp=imputations[[method]][[1]], N=N, imputationfunc=imputationfuncs[[method]], maxlength=maxlength, ...)
      score_all[[method]] <- tmp  
      
    }
    
    
    
  }
  
  return(score_all)
  
}


Iscores_new_perimp <- function(X, Ximp, N=50, imputationfunc, maxlength=NULL,...){
  
  if (is.null(Ximp)){
    # Impute, maxit should not be 1 here!
    Ximp<-imputationfunc(X=X  , m=1)[[1]]
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

  ## Done in the function
  M<-1*is.na(X)
  colnames(M) <- colnames(X)
  
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

    
    print( paste0("Dimension ", i, " out of ", maxlength )   ) 
    
  
    
    # H for all missing values of X_j
    Ximp1<-Ximp[M[,j]==1, ]
    
    # H for all observed values of X_j
    Ximp0<-Ximp[M[,j]==0, ]
    
    X0 <-X[M[,j]==0, ]
    
    n1<-nrow(Ximp1)
    n0<-nrow(Ximp0)
    
    
    if (n1 < 10){
      scorej[i]<-NA
      
      warning('Sample size of missing and nonmissing too small for nonparametric distributional regression, setting to NA')
      
    }else{
      
      
      # Evaluate on observed data
      Xtest <- Ximp0[,!(colnames(Ximp0) %in% j) &  (colnames(Ximp0) %in% indexfull), drop=F]
      Oj<-apply(X0[,!(colnames(Ximp0) %in% j) &  (colnames(Ximp0) %in% indexfull), drop=F],2,function(x) !any(is.na(x)) )
      # Only take those that are fully observed
      Xtest<-Xtest[,Oj, drop=F]
      
      Ytest <-Ximp0[,j, drop=F]
      
      if (is.null(Xtest)){
        scorej[i]<-NA
        #weighted
        weight[i]<-(n1/n)*(n0/n)
        warning("Oj was empty")
        next
      }
      
      ###Test 1:
      # Train DRF on imputed data
      Xtrain<-Ximp1[,!(colnames(Ximp1) %in% j) & (colnames(Ximp1) %in% indexfull), drop=F]
      # Only take those that are fully observed
      Xtrain<-Xtrain[,Oj, drop=F]
      
      Ytrain<-Ximp1[,j, drop=F]

      
      Xartificial<-cbind(c(rep(NA,nrow(Ytest)),c(Ytrain)),rbind(Xtest, Xtrain)   )
      colnames(Xartificial)<-c(colnames(Ytrain), colnames(Xtrain))
      
      Imputationlist<-imputationfunc(X=Xartificial  , m=N)
      
      Ymatrix<-do.call(cbind, lapply(Imputationlist, function(x)  x[1:nrow(Ytest),1]  ))
      
      scorej[i] <- -mean(sapply(1:nrow(Ytest), function(l)  { crps_sample(y = Ytest[l,], dat = Ymatrix[l,]) }))
      
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
