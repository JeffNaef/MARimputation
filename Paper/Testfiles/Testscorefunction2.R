


firstpart<-
  mean(apply(X,1,
      function (x) {
        mean( sapply(1:nrow(X), function(j)  norm(Ximp[j,, drop=F]- x,type="F")    )  )
      }))


firstpart<-
  mean(apply(X,1,
             function (x) {
               scoringRules:::esC_xy(x, t(Ximp), rep(1/nrow(Ximp),nrow(Ximp)))
             }))



  
  
secondpart<-scoringRules:::esC_xx(t(Ximp), rep(1/nrow(Ximp),nrow(Ximp)))
  

firstpart - 0.5*secondpart


es_sample(y = X[l,], dat = t(Ximp))
scoringRules:::esC_xy(X[l,], t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp))) - 
  0.5 * scoringRules:::esC_xx(t(Ximp), w=rep(1/nrow(Ximp),nrow(Ximp)))
###Why does this not give the same????