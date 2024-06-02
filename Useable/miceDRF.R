



mice.impute.DRF<- function (y, ry, x, wy = NULL,min.node.size=1, num.features=10,  num.trees=10 , ...)
{
    #install.on.demand("drf", ...)
  require(drf)
  if (is.null(wy)) {
    wy <- !ry
  }
  if (dim(x)[2] == 0) {
    x <- cbind(x, 1)
    dimnames(x) <- list(NULL, "int")
  }
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- as.matrix(y[ry])
  
  # args<-list(...)
  # if ("m" %in% names(args)){
  #   m<-args$m
  # }else{
  #   m=1
  # }
  m<-1
  

  fit <- drf(Y=yobs, X=xobs,num.trees=num.trees, num.features=num.features, compute.oob.predictions = F, min.node.size=min.node.size)
  DRFw <- predict(fit, newdata=xmis)$weights # These are the nodes now
  impute <- vapply(1:nrow(xmis), function(s) yobs[sample(1:nrow(yobs), size=1, replace=T, prob=DRFw[s,]), ], numeric(m))  # sample one observation per xmis
  
  
  

  impute
}
