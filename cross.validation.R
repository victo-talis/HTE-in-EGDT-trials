

# oest.cv.fun takes arguments:
# 1. The cleaned list of datasets (if multiply imputed); list of length 1 if data are complete
# 2. The proportion to use for splitting; if proportion= -x, a random proportion between x and 1-x is used.
# 3. The name of the outcome variable
# 4. The name of the treatment variable
# 5. A character vector indicating the names of the covariates to model
# 6. The fixed value of the randomization probability to the treatment group (known propensity score)
# 6. A character vector of the function names to call, each producing a score 
# 7. The number of cores to pass ot mclapply (if cores=1, uses lapply)

# The output values are:
# 1. A list of datasets with length(mods)+2 variables: outcome, treatment, score1,...,scorex, 
#     and nrow equal to the number of patients in the CV-test set.
# 2. A list of variable importance scores (where applicable), of the same length as 1.


library(parallel)

oest.cv.fun<-function(train,
                      cv.tr.pr=0.5,
                      cv_k=50,
                      trids0=NULL,
                      y,
                      w,
                      x,
                      pi=0.5,
                      mods,
                      cores=1,
                      seed=15223,
                      ...){
 
  M<-length(train)
  J<-length(mods)

  train1<-train[[1]]
  #browser()
  iy<-which(colnames(train1)==y)
  iw<-which(colnames(train1)==w)
  ix<-which(colnames(train1) %in% x)
  
  # Create a matrix with the indices for each of the cv_k training sets
  iw0<-as.numeric(rownames(train1[train1[,iw]==0,]))
  iw1<-as.numeric(rownames(train1[train1[,iw]==1,]))
  
  
  # have option to do CV or not
  if(is.null(trids0) & cv_k>1){
    cv.tr.ids<-matrix(,floor(length(iw0)*cv.tr.pr)+floor(length(iw1)*cv.tr.pr), cv_k)
    set.seed(seed) #15213
    for (i in 1:cv_k){
      cv.w0<-sample(iw0,floor(length(iw0)*cv.tr.pr))
      cv.w1<-sample(iw1,floor(length(iw1)*cv.tr.pr))
      cv.tr.ids[,i]<-c(cv.w0,cv.w1)
    }
  } else if (!is.null(trids0)) {
    cv.tr.ids<-matrix(trids0,length(trids0), 1)
  } else {
    stop("cv_k=0 but no training set identified through trids0")
  }
  
  VImat<-matrix(0,ncol(as.matrix(train1[,-c(iy,iw)])),0)
  rownames(VImat)<-names(train1[,-c(iy,iw)])
  
  # Define a function that produces the output for a given cv split
  cvfun<-function(k,
                  M,
                  cv.tr.ids,
                  mthd,
                  train,
                  y,w,x,pii=pi,
                  ...){ 
    
    #browser()
    train1<-train[[1]]
    cvt<-cv.tr.ids[,k]
    cve<-c(1:nrow(train1))[-cvt]
    #browser()
    
    fun.exec<-match.fun(mthd)(train,cvt,cve,M,y,w,x,pii,...)
    
    #if( !(mthd %in% c("ModCovLogitLASSO.CATE","ModCovLogitLASSO.CATE.xinter"))){
    #  fun.exec<-match.fun(mthd)(train,cvt,cve,M,y,w,x,pii)
    #} else {
    #  #browser()
    #  fun.exec<-match.fun(mthd)(train,cvt,cve,M,y,w,x,pii,xinter=FALSE,rmat)
    #}
    
    score<-fun.exec[[1]]
    importance<-fun.exec[[2]]
    
    colnames(train1)[iy]<-"Y"
    colnames(train1)[iw]<-"W"
    tempeval<-cbind(train1[cve,c("Y","W")], score)
    list(tempeval, importance)
  }
  
  res.allmods<-list()
  for (j in 1:J){
    if (cores>1){
      res.allmods[[j]]<-mclapply(1:cv_k, cvfun, M, cv.tr.ids, mods[j], train, y, w, x, pi, mc.cores=cores,...)
    } else {
      res.allmods[[j]]<-lapply(1:cv_k, cvfun, M, cv.tr.ids, mods[j], train, y, w, x, pi,...)
    }
  }
  
  names(res.allmods)<-mods
  
  res.allmods # a list of length J, each item alist of length cv_k, each item of length 2
    
}



###############################



# Fast version of the above that does cv_k CV iterations, one on each of M=cv_k MI datasets
# So, length(train) must equal cv_k

oest.cv.fast<-function(train,
                      cv.tr.pr=0.5,
                      cv_k=50,
                      trids0=NULL,
                      y,
                      w,
                      x,
                      pi=0.5,
                      mods,
                      cores=1,
                      seed=15212,
                      ...){
  
  if(cv_k>1){
    if(length(train)!=cv_k) stop(paste0("train object does not contain ",cv_k," datasets."))
  }
  
  
  M<-length(train)
  J<-length(mods)
  
  train1<-train[[1]]
  #browser()
  iy<-which(colnames(train1)==y)
  iw<-which(colnames(train1)==w)
  ix<-which(colnames(train1) %in% x)
  
  # Create a matrix with the indices for each of the cv_k training sets
  iw0<-as.numeric(rownames(train1[train1[,iw]==0,]))
  iw1<-as.numeric(rownames(train1[train1[,iw]==1,]))
  
  # have option to do CV or not
  if(is.null(trids0) & cv_k>1){
    cv.tr.ids<-matrix(,floor(length(iw0)*cv.tr.pr)+floor(length(iw1)*cv.tr.pr), cv_k)
    set.seed(seed) #15213
    for (i in 1:cv_k){
      cv.w0<-sample(iw0,floor(length(iw0)*cv.tr.pr))
      cv.w1<-sample(iw1,floor(length(iw1)*cv.tr.pr))
      cv.tr.ids[,i]<-c(cv.w0,cv.w1)
    }
  } else if (!is.null(trids0)) {
    cv.tr.ids<-matrix(trids0,length(trids0), 1)
  } else {
    stop("cv_k=0 but no training set identified through trids0")
  }
  
  
  VImat<-matrix(0,ncol(as.matrix(train1[,-c(iy,iw)])),0)
  rownames(VImat)<-names(train1[,-c(iy,iw)])
  
  # Define a function that produces the output for a given cv split
  cvfun<-function(k,
                  M,
                  cv.tr.ids,
                  mthd,
                  train,
                  y,w,x,pii=pi,
                  ...){ 
    
    #browser()
    train1<-train[[1]]
    cvt<-cv.tr.ids[,k]
    cve<-c(1:nrow(train1))[-cvt]
    #browser()
    
    fun.exec<-match.fun(mthd)(list(train[[k]]), cvt, cve, 1, y, w, x, pii,...)
    
    #if( !(mthd %in% c("ModCovLogitLASSO.CATE","ModCovLogitLASSO.CATE.xinter"))){
    #  fun.exec<-match.fun(mthd)(train,cvt,cve,M,y,w,x,pii)
    #} else {
    #  #browser()
    #  fun.exec<-match.fun(mthd)(train,cvt,cve,M,y,w,x,pii,xinter=FALSE,rmat)
    #}
    
    score<-fun.exec[[1]]
    importance<-fun.exec[[2]]
    
    colnames(train1)[iy]<-"Y"
    colnames(train1)[iw]<-"W"
    tempeval<-cbind(train1[cve,c("Y","W")], score)
    list(tempeval, importance)
  }
  
  res.allmods<-list()
  for (j in 1:J){
    if (cores>1){
      res.allmods[[j]]<-mclapply(1:cv_k, cvfun, M, cv.tr.ids, mods[j], train, y, w, x, pi, mc.cores=cores,...)
    } else {
      res.allmods[[j]]<-lapply(1:cv_k, cvfun, M, cv.tr.ids, mods[j], train, y, w, x, pi,...)
    }
  }
  
  #browser()
  
  names(res.allmods)<-mods
  
  res.allmods # a list of length J, each item alist of length cv_k, each item of length 2
  
}


