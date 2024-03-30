# evaluation of results of main_test_newcode.R



datasets <- c("yeast","airfoil","concrete.slump","seeds","Boston","concrete.compression",
             "CASchools","connectionist.bench.vowel","ecoli","iris","yacht",
             "wine", "planning.relax","climate.model.crashes")


  # ,,
  # ,,  , "ionosphere"

num.with.errors.new <- rep(NA, length(datasets))
names(num.with.errors.new)<- datasets

frac.first.new <- rep(NA, length(datasets))
names(frac.first.new)<- datasets

num.with.errors.old <- rep(NA, length(datasets))
names(num.with.errors.old)<- datasets

frac.first.old <- rep(NA, length(datasets))
names(frac.first.old)<- datasets

for(use.score in c("old", "new")){

  i <- 1

  for(dataset in datasets){


    if(use.score=="new"){
      if(dataset =="climate.model.crashes"){
        load(paste0(dataset,"_MAR_pmiss0.2_m5_nrep30_upsamplingTRUE_perpatternTRUE_min.node.size10_num.proj200_num.trees.per.proj5.Rda"))

      }else{
        load(paste0(dataset,"_MAR_pmiss0.2_m5_nrep30_upsamplingTRUE_perpatternTRUE_min.node.size10_num.proj100_num.trees.per.proj5.Rda"))

      }
    }else{
      if(dataset =="climate.model.crashes"){
        load(paste0(dataset,"_old_MAR_pmiss0.2_m5_nrep30_upsamplingTRUE_perpatternTRUE_min.node.size10_num.proj200_num.trees.per.proj5.Rda"))

      }else{
        load(paste0(dataset,"_old_MAR_pmiss0.2_m5_nrep30_upsamplingTRUE_perpatternTRUE_min.node.size10_num.proj100_num.trees.per.proj5.Rda"))

      }
    }



    Results <- unlist(Results, recursive = F)
    Results <- lapply(Results,function(l){

      if(length(l)>1){
        return(l)
      }else{

      }

    })

    Results <- Results[lapply(Results,length)>0]

    if(use.score=="new"){
      num.with.errors.new[i] <- 30- length(Results)
    }else{
      num.with.errors.old[i] <- 30- length(Results)
    }


    # num.datasets <- unlist(lapply(1:length(Results), function(i){
    #   length(Results[[i]]$new.score)
    #
    # }))

    ranks <- lapply(1:length(Results), function(i){
      inds.ordered <- order(Results[[i]], decreasing=TRUE)
      res <- Results[[i]][,inds.ordered, drop=F]
      print(res)
      rank <- which(colnames(res)=="truth")
      rank
    })

    if(use.score=="new"){
      frac.first.new[i] <- sum( unlist(ranks)==1 )/length(Results)

    }else{
      frac.first.old[i] <- sum( unlist(ranks)==1 )/length(Results)

    }

    i <- i+1
    rm(Results)


  }
}



cbind(frac.first.new, frac.first.old, num.with.errors.new, num.with.errors.old)

