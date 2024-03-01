
blockmiceDRF<-function(data, blocksize, blocklist=NULL,num.trees=10, num.features=10, m=1 ){
  
  n<-nrow(data)
  d<-ncol(data)
  
  if (is.null(colnames(data))){
    colnames(data)<-paste0("X.",c(1:d))
  }
  
  if (is.null(blocklist)){
  
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
  
  
  return(mice(data, method="myfunc", num.trees=num.trees, num.features=num.features, m=m, blocks = blocklist ))
}






mice.impute.myfunc<-function (y, ry, x, wy = NULL, minbucket = 5, cp = 1e-04, ...) 
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
  xobs <- data.frame(x[ry, , drop = FALSE])
  xmis <- data.frame(x[wy, , drop = FALSE])
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
  
    
    
  
  fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F)
  DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
  impute <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(1))  # sample one observation per xmis
  
  
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
  # else {
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