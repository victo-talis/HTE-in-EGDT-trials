
### ~~~~~~~~~~~~~~~~~~~~~~~~ Section 7 ~~~~~~~~~~~~~~~~~~~~~~~###

# Optimal enrichment score threshold (OEST)
# Chooses the method associated with the threshold leading to the highest
# estimated Pr(p<0.05 | c) for possible score values c, where the p-value
# is testing H0: CATE=0 within the subgroup defined as score>c.

# First I define methods for CATE estimation: 
# 1. Random forest s-learner for binary outcomes
# 2. Causal forest
# 3. BART s-learner for binary outcomes
# 4. Ridge logistic, s-learner
# 5. LASSO logistic, s-learner
# 6. FindIt s-learner

# I also define methods for risk estimation in control group: 
# 1. Logistic regression
# 2. Penalized logistic, ridge
# 3. Penalized logistic, LASSO
# 4. Random forest for binary outcomes
# 5. BART for binary outcomes



library(randomForest)
library(BART)
library(grf)
library(glmnet)
devtools::load_all(paste0(machine,"_BoxMigration/HTE analysis programs/miselect-master/miselect-master"),recompile=TRUE)
devtools::load_all(paste0(machine,"_BoxMigration/HTE analysis programs/causalTree-master/causalTree-master"),recompile=TRUE)
devtools::load_all(paste0(machine,"_BoxMigration/HTE analysis programs/decor-main/decor-main"),recompile=TRUE)
library(miselect)
#library(XBART)
#source("C:/Users/vit13/Box Sync/Faraaz PROWESS/Programs/OEST Programs/causalMARS_mod1.R")
library(caret)
library(recipes)
library(causalTree)


rloss<-function(data,lev=NULL,model=NULL){
  #browser()
  # w<-data$w
  # y<-data$y
  # mx<-data$mx
  # ex<-data$ex
  # pred<-data$pred
  rm<-mean((data$y-data$mx-(data$w-data$ex)*data$pred)^2)
  c(rloss=rm)
}


nnet.CATE.rlearn<-function(train, cvt, cve, M, y, w, x, pi=0.5){
  
  iw<-which(colnames(train[[1]])==w)
  wname<-names(train[[1]])[iw]
  iy<-which(colnames(train[[1]])==y)
  
  y <- w <- list()
  xmat <- list()
  xpred<-list()
  mxmat<-matrix(,length(cvt),M)
  
  for(m in 1:M){
    
    ds<-train[[m]]
    ds0 <- ds[,-c(iy,iw)]
    
    if(identical(cvt,cve)){
      y[[m]]<-ds[,iy]
      w[[m]]<-ds[,iw]
      xmat[[m]]<-xpred[[m]]<-as.matrix(ds0)  #excludes x
    } else {
      dst<-ds0[cvt,]
      dse<-ds0[cve,]
      y[[m]]<-ds[cvt,iy]
      w[[m]]<-ds[cvt,iw]
      xmat[[m]]<-as.matrix(dst)
      xpred[[m]]<-as.matrix(dse)
    }
    
    mfit<-regression_forest(X=as.matrix(dst),
                          Y=y[[m]])
    mxmat[,m]<-predict(mfit)$predictions
    
  }
  
  p <- ncol(xmat[[1]]) 
  wt <- rep(1/M,length(cvt))*((w[[1]]-0.5)^2)
  
  xmat.all <- matrix(numeric(0),,p)
  psy.all<-wt.all<-y.all<-mx.all<-ex.all<-w.all<-c()
  for(m in 1:M){
    xmat.all<-rbind(xmat.all,xmat[[m]])
    y.all<-c(y.all,y[[m]])
    w.all<-c(w.all,w[[m]])
    mx.all<-c(mx.all,mxmat[,m])
    ex.all<-c(ex.all,0.5)
    psy.all<-c(psy.all,(y[[m]]-mxmat[,m])/(w[[m]]-0.5))
    wt.all<-c(wt.all,wt)
  }
  
  browser()
  trainDat<-as.data.frame(cbind(psy.all, xmat.all, mx=mx.all, ex=ex.all, y=y.all, w=w.all, wt=wt.all))
  
  nnr.recipe<-recipe(psy.all~., data=trainDat) %>%
    update_role(mx, new_role = "performance var") %>%
    update_role(ex, new_role = "performance var") %>%
    update_role(y, new_role = "performance var") %>%
    update_role(w, new_role = "performance var") %>%
    update_role(wt, new_role = "case weight")
  
  #browser()
  nnet.cv<-train(nnr.recipe,
                 data=trainDat,
                 method="nnet",
                 #weights=wt.all,
                 metric="rloss",
                 #tuneLength = 5,
                 trControl = trainControl(summaryFunction = rloss),
                 nnet.trace=FALSE,
                 linout=TRUE)
  
  browser()
  nnet.fin<-nnet.cv$finalModel
  
  predmat<-matrix(,length(cve),M)
  for(m in 1:M){
    predmat[,m]<- predict(nnet.fin, newdata=xpred[[m]])
  }
  
  #browser()
  cate<-apply(predmat,1,mean)

  #browser()
  list(cate)
  
}



nnet.CATE.ModCov<-function(train, cvt, cve, M, y, w, x, pi=0.5){
 
  iw<-which(colnames(train[[1]])==w)
  wname<-names(train[[1]])[iw]
  iy<-which(colnames(train[[1]])==y)
  
  y <- list()
  xmat <- list()
  xpred<-list()
  xcols<-which(names(train[[1]]) %in% x)
  xcols2<-c()
  
  for(m in 1:M){
    
    ds<-train[[m]]
    ds[,iw]<-ds[,iw]-0.5
    
    lastcol<-ncol(ds)
    newj<-lastcol+1
    
    if(length(x)>1){
      for(j1 in xcols){
        k<-which(xcols==j1)
        ds[,newj]<-ds[,j1]*ds[,iw]
        NewName<-paste0(names(ds)[j1],"_trans")
        names(ds)[newj]<-NewName
        if(m==1){
          #rmat<-cbind(rmat,pmax(rmat[,j1],rmat[,iw]))
          #colnames(rmat)[newj]<-NewName
          xcols2<-c(xcols2,newj)
        }
        newj=newj+1
      }
    } else {
      j1<-xcols[1]
      ds[,newj]<-ds[,j1]*ds[,iw]
      NewName<-paste0(names(ds)[j1],"_trans")
      names(ds)[newj]<-NewName
      if(m==1){
        #rmat<-cbind(rmat,pmax(rmat[,j1],rmat[,iw]))
        #colnames(rmat)[newj]<-NewName
        xcols2<-c(xcols,newj)
      }
    }

    ds0 <- ds[,c(iw,(lastcol+1):(newj-1))]
    
    if(identical(cvt,cve)){
      y[[m]]<-ds[,iy]
      xmat[[m]]<-xpred[[m]]<-as.matrix(ds0)
    } else {
      dst<-ds0[cvt,]
      dse<-ds0[cve,]
      y[[m]]<-ds[cvt,iy]
      xmat[[m]]<-as.matrix(dst)
      xpred[[m]]<-as.matrix(dse)
    }
    
  }
  
  
  p <- ncol(xmat[[1]]) 
  wt <- rep(1/M,length(cvt))

  xmat.all <- matrix(numeric(0),,p)
  y.all<-wt.all<-c()
  for(m in 1:M){
    xmat.all<-rbind(xmat.all,xmat[[m]])
    y.all<-c(y.all,y[[m]])
    wt.all<-c(wt.all,wt)
  }
  
  #browser()
  nnet.cv<-train(x=xmat.all,
            y=y.all,
            method="nnet",
            weights=wt.all,
            metric="RMSE",
            nnet.trace=FALSE)
  
  nnet.fin<-nnet.cv$finalModel
  
  predmat<-matrix(,length(cve),M)
  for(m in 1:M){
    xpred0 <-xpred[[m]]
    wpred0 <-xpred0[,1]
    for(j in 1:p){
      xpred0[,j]<-xpred0[,j]/wpred0
    }
    predmat[,m]<- predict(nnet.fin, newdata=xpred0)
  }

  pred.mean<-apply(predmat,1,mean)
  cate<-(exp(pred.mean/2)-1)/(exp(pred.mean/2)+1)
  #phat<-exp(pred.mean)/(1+exp(pred.mean))
  
  #browser()
  list(cate,
       pred.mean)
   
}


nnet.CATE.slearn<-function(train, cvt, cve, M, y, w, x, pi=0.5){
  
  iw<-which(colnames(train[[1]])==w)
  wname<-names(train[[1]])[iw]
  iy<-which(colnames(train[[1]])==y)
  
  y <- list()
  xmat <- list()
  xpred<-list()
  xcols<-which(names(train[[1]]) %in% x)
  xcols2<-c()
  
  for(m in 1:M){
    
    ds<-train[[m]]
    ds0 <- ds[,-iy]
    
    if(identical(cvt,cve)){
      y[[m]]<-ds[,iy]
      xmat[[m]]<-xpred[[m]]<-as.matrix(ds0)
    } else {
      dst<-ds0[cvt,]
      dse<-ds0[cve,]
      y[[m]]<-ds[cvt,iy]
      xmat[[m]]<-as.matrix(dst)
      xpred[[m]]<-as.matrix(dse)
    }
    
  }

  p <- ncol(xmat[[1]]) 
  wt <- rep(1/M,length(cvt))
  
  xmat.all <- matrix(numeric(0),,p)
  y.all<-wt.all<-c()
  for(m in 1:M){
    xmat.all<-rbind(xmat.all,xmat[[m]])
    y.all<-c(y.all,y[[m]])
    wt.all<-c(wt.all,wt)
  }
  
  #browser()
  nnet.cv<-train(x=xmat.all,
                 y=y.all,
                 method="nnet",
                 weights=wt.all,
                 metric="RMSE",
                 nnet.trace=FALSE)
  
  nnet.fin<-nnet.cv$finalModel

  predmat0<-predmat1<-predmat<-matrix(,length(cve),M)
  for(m in 1:M){
    xpred0<-xpred1<-xpred[[m]]
    xpred0[,1]<-0
    xpred1[,1]<-1
    predmat0[,m]<- predict(nnet.fin, newdata=xpred0)
    predmat1[,m]<- predict(nnet.fin, newdata=xpred1)
    predmat[,m]<- predict(nnet.fin, newdata=xpred[[m]])
  }
  
  #browser()
  cate<-apply(predmat1-predmat0,1,mean)
  phat<-apply(predmat,1,mean)
  
  #browser()
  list(cate,
       phat)
  
}




rf.CATE.xlearn<-function(train, cvt, cve, M, y, w, x, pi=0.5, ntrees=1000,hon=TRUE){
  
  iw<-which(colnames(train[[1]])==w)
  iy<-which(colnames(train[[1]])==y)
  rf.x.pred.all<-import.all<-c()
  
  for (m in 1:M) {
    ds<-train[[m]]
    dst<-ds[cvt,]

    if (identical(ds[cvt,],ds[cve,])){
      perm<-order(dst[,w]) #to pop back, use order(perm)
    }

    dst.0<-dst[dst[,w]==0,]
    dst.1<-dst[dst[,w]==1,]

    # Step A: Fit rf models to each arm separately, and make predictions on the opposite arm
    
    #browser()
    
    rf.0<-regression_forest(X=as.matrix(dst.0[,x]),
                            Y=dst.0[,y],
                            num.trees=ntrees,
                            mtry=floor(sqrt(length(x))),
                            honesty=hon)                    #Devault is sqrt(p) for survival tree

    # Is this oob?
    pred.0.tr<-predict(rf.0,
                       newdata = as.matrix(dst.1[,x]))[,1]
    
    rf.1<-regression_forest(X=as.matrix(dst.1[,x]),
                             Y=dst.1[,y],
                             num.trees=ntrees,
                             mtry=floor(sqrt(length(x))),
                            honesty=hon)
    pred.1.tr<-predict(rf.1,
                       newdata = as.matrix(dst.0[,x]))[,1]
    
    
    # Step B: Fit rf models to each pseudo-outcome ps.1=Y1-pred.0.tr and ps.0=pred.1.tr-Y0
    
    dst.0$ps<-pred.1.tr-dst.0[,y]
    dst.1$ps<-dst.1[,y]-pred.0.tr
    
    tau<-list()
    
    tau[[1]]<-regression_forest(X=as.matrix(dst.0[,x]),
                                Y=dst.0[,"ps"],
                                num.trees=ntrees,
                                 mtry=floor(sqrt(length(x))),
                                honesty=hon)
    tau[[2]]<-regression_forest(X=as.matrix(dst.1[,x]),
                                Y=dst.1[,"ps"],
                                num.trees=ntrees,
                                mtry=floor(sqrt(length(x))),
                                honesty=hon)
    
    dse<-ds[cve,]

    # Step C: Make predictions from both tau.0 and tau.1 for every x from the test set
    # This is what needs to be OOB if I'm using for LOO cross-fitting
    
    # If train==test, and W==1, then tau[[2]] should be run so that I have OOB predictions
    # for the W==1 set.
    # If train==test, and W==0, then tau[[1]] should be run so that I get OOB predictions 
    
    if (identical(cvt,cve)){
      pred.tau0.w0<-predict(tau[[1]])[,1]
      pred.tau0.w1<-predict(tau[[1]],
                             newdata = as.matrix(dst.1[,x]))[,1]
      pred.tau0<-c(pred.tau0.w0,pred.tau0.w1)[order(perm)]
      
      pred.tau1.w0<-predict(tau[[2]],
                            newdata=as.matrix(dst.0[,x]))[,1]
      pred.tau1.w1<-predict(tau[[2]])[,1]
      pred.tau1<-c(pred.tau1.w0,pred.tau1.w1)[order(perm)]
      
    } else {
      pred.tau0<-predict(tau[[1]],
                         newdata = as.matrix(dse[,x]))[,1]
      pred.tau1<-predict(tau[[2]],
                         newdata = as.matrix(dse[,x]))[,1]
      
    }

    #browser()
    score<-pi*pred.tau0 + (1-pi)*pred.tau1
    rf.x.pred.all<-cbind(rf.x.pred.all,score)
    
    import<-apply(cbind(variable_importance(tau[[1]]),variable_importance(tau[[2]])),1,mean)
    import.all<-cbind(import.all,import)
  }

  
  #browser()
  list("preds"=apply(rf.x.pred.all,1,mean),
       "import"=apply(import.all,1,mean))
}






cf.lin.CATE<-function(train, cvt, cve, M, y, w, x, pi=0.5,ntrees=5000,tune=FALSE,predY=FALSE,
                  forest.tune.params=c("sample.fraction","mtry","min.node.size","alpha"), 
                  lin.vars=NULL){
  
  if(is.null(lin.vars)) lin.vars<-x
  
  cf.pred.all<-cf.imp.all<-cf.predy.all<-c()
  for (m in 1:M){
    
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]
    
    #browser()
    if(tune){
      #browser()
      cf<-causal_forest(X=as.matrix(dst[,x]),
                        Y=dst[,y],
                        W=dst[,w],
                        Y.hat=mean(dst[,y]),
                        W.hat=pi,
                        tune.parameters=forest.tune.params,
                        num.trees=ntrees,
                        honesty=TRUE,   
                        honesty.fraction=0.5,
                        num.trees=ntrees)
    } else {
      hyp=list("params"=c(5,0.5,ceiling(2*length(x)/3),0.05,0))
      names(hyp[[1]])<-c("min.node.size","sample.fraction",
                         "mtry","alpha","imbalance.penalty")
      cf<-causal_forest(X=as.matrix(dst[,x]),
                        Y=dst[,y],
                        W=dst[,w],
                        Y.hat=NULL, #mean(dst[,y]),     #Don't need to fit a y.hat model
                        W.hat=pi,                #otherwise silently fits an RF model to estimate propensity
                        honesty=TRUE,   
                        honesty.fraction=0.5,
                        num.trees=ntrees,
                        sample.fraction=as.numeric(hyp$params["sample.fraction"]),
                        min.node.size=as.numeric(hyp$params["min.node.size"]),
                        alpha=as.numeric(hyp$params["alpha"]),
                        imbalance.penalty = as.numeric(hyp$params["imbalance.penalty"]),
                        mtry=as.numeric(hyp$params["mtry"]))
    }
    
    browser()
    tuned.lambda <- tune_ll_causal_forest(forest)
    pwc<-predict(cf, 
                 linear.correction.variables = lin.vars,
                 ll.lambda = tuned.lambda$lambda.min)
    
    
    #browser()
    if (identical(dst,dse)){
      #browser()
      cf.pred<-predict(cf, 
                       linear.correction.variables = lin.vars,
                       ll.lambda = tuned.lambda$lambda.min)$predictions  # oob predictions
    } else {
      cf.pred<-predict(cf, 
                       linear.correction.variables = lin.vars,
                       ll.lambda = tuned.lambda$lambda.min,
                       newdata=as.matrix(dse[,x]))$predictions
    }
    
    cf.pred.all<-cbind(cf.pred.all,cf.pred)
    
    cf.imp<-variable_importance(cf)
    cf.imp.all<-cbind(cf.imp.all,cf.imp)
    
    if(predY){
      #browser()
      a_fromCF<-get_forest_weights(cf,
                                   newdata=as.matrix(dse[,x]))
      
      sum_aY_mu0<-a_fromCF[dse$w==0,dst$w==0] %*% matrix(dst[dst$w==0,y],,1)
      sum_a_mu0<-apply(a_fromCF[dse$w==0,dst$w==0],1,sum)
      preds_mu0_fromCF<-as.numeric(sum_aY_mu0/sum_a_mu0)
      
      sum_aY_mu1<-a_fromCF[dse$w==1,dst$w==1] %*% matrix(dst[dst$w==1,y],,1)
      sum_a_mu1<-apply(a_fromCF[dse$w==1,dst$w==1],1,sum)
      preds_mu1_fromCF<-as.numeric(sum_aY_mu1/sum_a_mu1)
      
      pred.y<-c(preds_mu0_fromCF,preds_mu1_fromCF)[order(order(dse$w))]
      
      cf.predy.all<-cbind(cf.predy.all,pred.y)
      
    } else {
      cf.predy.all<-cbind(cf.predy.all,rep(0,length(cve)))
    }
    
    
  }
  
  
  #browser()
  list(apply(cf.pred.all,1,mean),
       apply(cf.imp.all,1,mean),
       apply(cf.predy.all,1,mean))
  
}


llrf.CATE.rlearn<-function(train, cvt, cve, M, y, w, x, pi=0.5, ntrees=5000, ...){

  args<-list(...)
  
  #browser()
  
  pred.all<-c()
  for (m in 1:M){
    
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]
    
    if("cluster" %in% names(args)){
      clust<-dst[,args$cluster]
    } else {
      clust<-NULL
    }
    
    if("mxvars" %in% names(args)){
      mxvars<-args$mxvars
    }
    
    mfit<-regression_forest(X=as.matrix(dst[,mxvars]),
                            Y=dst[,y],
                            cluster=clust)
    mx<-mfit$predictions[,1]
    yst<-(dst[,y]-mx)/(dst[,w]-0.5)
    rtest_ll<-ll_regression_forest(X=as.matrix(dst[,x]),
                                   Y=yst,
                                   enable.ll.split=TRUE,
                                   ll.split.weight.penalty=TRUE,
                                   ll.split.lambda=0.1,
                                   num.trees=ntrees,
                                   cluster=clust)
    
    if (identical(dst,dse)){
      #browser()
      pred<-predict(rtest_ll)$predictions  # oob predictions
    } else {
      pred<-predict(rtest_ll, 
                    newdata=as.matrix(dse[,x]))$predictions
    }
    
    pred.all<-cbind(pred.all,pred)
  }
    
  list(apply(pred.all,1,mean),
       NA,
       NA)
}



cf.CATE<-function(train, cvt, cve, M, y, w, x, pi=0.5,ntrees=5000,tune=FALSE,predY=FALSE){

  cf.pred.all<-cf.imp.all<-cf.predy.all<-c()
  for (m in 1:M){
  
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]

  #browser()
    if(tune){
      #browser()
      hyp<-tune_causal_forest(X=as.matrix(dst[,x]),
                              Y=dst[,y],
                              W=dst[,w],
                              Y.hat=mean(dst[,y]),
                              W.hat=pi,
                              tune.num.trees=1000,
                              tune.parameters=c("sample.fraction","mtry","min.node.size","alpha","imbalance.penalty"))
    } else {
      hyp=list("params"=c(5,0.5,ceiling(2*length(x)/3),0.05,0))
      names(hyp[[1]])<-c("min.node.size","sample.fraction",
               "mtry","alpha","imbalance.penalty")
    }

    cf<-causal_forest(X=as.matrix(dst[,x]),
                      Y=dst[,y],
                      W=dst[,w],
                      Y.hat=NULL, #mean(dst[,y]),     #Don't need to fit a y.hat model
                      W.hat=pi,                #otherwise silently fits an RF model to estimate propensity
                      honesty=TRUE,   
                      honesty.fraction=0.5,
                      num.trees=ntrees,
                      sample.fraction=as.numeric(hyp$params["sample.fraction"]),
                      min.node.size=as.numeric(hyp$params["min.node.size"]),
                      alpha=as.numeric(hyp$params["alpha"]),
                      imbalance.penalty = as.numeric(hyp$params["imbalance.penalty"]),
                      mtry=as.numeric(hyp$params["mtry"]))
    #browser()
    if (identical(dst,dse)){
      #browser()
      cf.pred<-predict(cf)$predictions  # oob predictions
    } else {
      cf.pred<-predict(cf, newdata=as.matrix(dse[,x]))$predictions
    }
    
    cf.pred.all<-cbind(cf.pred.all,cf.pred)
    
    cf.imp<-variable_importance(cf)
    cf.imp.all<-cbind(cf.imp.all,cf.imp)

    if(predY){
      #browser()
      a_fromCF<-get_forest_weights(cf,
                                   newdata=as.matrix(dse[,x]))

      sum_aY_mu0<-a_fromCF[dse$w==0,dst$w==0] %*% matrix(dst[dst$w==0,y],,1)
      sum_a_mu0<-apply(a_fromCF[dse$w==0,dst$w==0],1,sum)
      preds_mu0_fromCF<-as.numeric(sum_aY_mu0/sum_a_mu0)

      sum_aY_mu1<-a_fromCF[dse$w==1,dst$w==1] %*% matrix(dst[dst$w==1,y],,1)
      sum_a_mu1<-apply(a_fromCF[dse$w==1,dst$w==1],1,sum)
      preds_mu1_fromCF<-as.numeric(sum_aY_mu1/sum_a_mu1)
      
      pred.y<-c(preds_mu0_fromCF,preds_mu1_fromCF)[order(order(dse$w))]
      
      cf.predy.all<-cbind(cf.predy.all,pred.y)
      
    } else {
      cf.predy.all<-cbind(cf.predy.all,rep(0,length(cve)))
    }
    
    
  }
  
  
  #browser()
  list(apply(cf.pred.all,1,mean),
       apply(cf.imp.all,1,mean),
       apply(cf.predy.all,1,mean))
  
}

cf.CATE.2<-function(train, cvt, cve, M, y, w, x, pi=0.5){
  
  cf.pred.all<-cf.imp.all<-c()
  for (m in 1:M){
    
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]
    
    cf<-causal_forest(X=as.matrix(dst[,x]),
                      Y=dst[,y],
                      W=dst[,w],
                      Y.hat=NULL,     #Silently fits an RF model to estimate E[Y|X=x]
                      W.hat=pi,      #otherwise silently fits an RF model to estimate propensity
                      honesty=TRUE,   #No need to predict in the training sample
                      honesty.fraction=0.5,
                      sample.fraction=0.5,
                      num.trees=2000,
                      ci.group.size=2,
                      min.node.size=2,
                      alpha=0.2,
                      imbalance.penalty = 0.01,
                      mtry=ceiling(2*length(x)/3))
    
    if (identical(dst,dse)){
      cf.pred<-predict(cf)$predictions  # oob predictions
    } else {
      cf.pred<-predict(cf, newdata=as.matrix(dse[,x]))$predictions
    }
    
    cf.yhat<-cf$Y.hat
    
    cf.pred.all<-cbind(cf.pred.all,cf.pred)
    
    cf.yhat.all<-cbind(cf.imp.all,cf.yhat)
    
  }
  
  #browser()
  list(apply(cf.pred.all,1,mean),
       apply(cf.yhat.all,1,mean))
  
}



cf2018.CATE<-function(train, cvt, cve, M, y, w, x, pi=0.5){
  
  form<-as.formula(paste0(y,"~",paste(x,collapse="+")))
  
  cf.pred.all<-c()
  cf.imp.all<-data.frame(var=numeric())
  for (m in 1:M){
    
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]
    
    #browser()
    cf<-causalForest(formula=form,
                     data=dst,
                     treatment=dst[,w],
                     #double.sample=TRUE,
                     #split.honest=TRUE,
                     minsize=2,
                     mtry=ceiling(2*length(x)/3),
                     sample.size.total=floor(nrow(dst)/2),
                     sample.size.train.frac=0.5,
                     num.trees=2000,
                     ncov_sample=ceiling(2*length(x)/3),
                     ncolx=length(x),
                     oob.predict=TRUE)
    
    
    if (identical(dst,dse)){
      cf.pred<-predict(cf)  # oob predictions
    } else {
      cf.pred<-predict(cf,newdata=dse[,x])
    }
    
    #cf.yhat<-cf$Y.hat
    
    cf.pred.all<-cbind(cf.pred.all,cf.pred)
    
    cf.imp<-bind_rows(lapply(1:length(cf$trees),function(x){
      d1<-as.data.frame(cf$trees[[x]]$variable.importance)
      d1$var<-rownames(d1)
      colnames(d1)[1]<-"imp"
      rownames(d1)<-1:nrow(d1)
      d1
    }))
    
    cf.imp2<-as.data.frame(tapply(cf.imp$imp,cf.imp$var,mean))
    cf.imp2$var<-rownames(cf.imp2)
    colnames(cf.imp2)[1]<-paste0("imp",m)
    rownames(cf.imp2)<-1:nrow(cf.imp2)
    #browser()
    cf.imp.all<-merge(cf.imp.all,cf.imp2,by="var",all.x=TRUE,all.y=TRUE)

  }
  
  cfimpout<-data.frame(var=cf.imp.all[,1],
                       imp=apply(cf.imp.all[,-1],1,function(x){mean(x,na.rm=TRUE)}))
  
  #browser()
  list(apply(cf.pred.all,1,mean),
       cfimpout)
  
}



BART.CATE.slearn<-function(train, cvt, cve, M, y, w, x, pi){
  
  iw<-which(colnames(train[[1]])==w)
  b.pred.all<-c()
  
  for (m in 1:M){
    
    ds<-train[[m]]
    colnames(ds)[iw]<-"W"
    dst<-dst.cf<-ds[cvt,]
    dst.cf$W<- 1-dst$W
    dse<-dse.cf<-ds[cve,]
    dse.cf$W<- 1-dse$W
    dse.all<-rbind(dse, dse.cf)
    
    pb<-pbart(x.train=dst[,c("W",x)],
              y.train=as.numeric(dst[,y]),
              x.test=dse.all[,c("W",x)])
    
    browser()
    te.obs<-pb$prob.test[,1:nrow(dse)]
    te.cf<-pb$prob.test[,(nrow(dse)+1):nrow(dse.all)]
    
    te.t1<-te.obs
    te.t1[,which(dse$W==0)]<-te.cf[,which(dse$W==0)] # Take the counterfactuals for T=0 subjects
    te.t0<-te.obs
    te.t0[,which(dse$W==1)]<-te.cf[,which(dse$W==1)] # Take the counterfactuals for T=1 subjects
    
    b.pred<-te.t1-te.t0
    b.pred.all<-rbind(b.pred.all,b.pred)
    
  }
  
  list(apply(b.pred.all,1,mean),NA)
  
}



XBART.CATE.slearn<-function(train, cvt, cve, M, y, w, x, pi){
  
  iw<-which(colnames(train[[1]])==w)
  cate_all<-c()
  
  for (m in 1:M){
    
    ds<-train[[m]]
    colnames(ds)[iw]<-"W"
    dst<-dst.cf<-ds[cvt,]
    dst.cf$W<- 1-dst$W
    dse<-dse.cf<-ds[cve,]
    dse.cf$W<- 1-dse$W
    dse.all<-rbind(dse, dse.cf)
    
    pb<-XBART.Probit(X=as.matrix(dst[,c("W",x)]),
                     y=as.matrix(dst[,y]),
                     Xtest=as.matrix(dse.all[,c("W",x)]),
                     num_sweeps=40,
                     burnin=15,
                     num_trees=floor(0.25*log(nrow(dst))^log(log(nrow(dst))))
                     )
    
    pf_draws<-pnorm(pb$yhats_test[1:nrow(dse),])
    pcf_draws<-pnorm(pb$yhats_test[(nrow(dse)+1):nrow(dse.all),])

    cate_draws<-t((pf_draws-pcf_draws)*(dse$W*2 -1))
    
    #browser()
    cate_all<-rbind(cate_all,cate_draws)
    print(m)
  }
  
  list(colMeans(cate_all),NA)
  
}


cMARS.CATE<-function(train, cvt, cve, M, y, w, x, propensity=FALSE,binaryY=TRUE, pi){
  
  iw<-which(colnames(train[[1]])==w)
  pred.all<-pred.all.0<-pred.all.1<-c()
  
  for (m in 1:M){
    
    ds<-train[[m]]
    colnames(ds)[iw]<-"W"
    dst<-ds[cvt,]
    dse<-ds[cve,]
    
    #browser()
    cmars.test<-bagged.causalMARS(y=dst$Y,
                                  tx=dst$W,
                                  x=as.matrix(dst[,x]),
                                  nbag=20,
                                  propensity=propensity)
    #browser()
    
    if (binaryY){
      pred<-predict.bagged.causalMARS_binaryY(cmars.test,
                                              newx=as.matrix(dse[,x]),
                                              type="lp.mat.list")
      pred.all.0<-rbind(pred.all.0, pred[[1]])
      pred.all.1<-rbind(pred.all.1, pred[[2]])
      
    } else {
      pred.tau<-predict.bagged.causalMARS(cmars.test,
                                              newx=as.matrix(dse[,x]),
                                              type="average")
      pred.all<-rbind(pred.all,pred.tau)
    }

  }
  
  if (binaryY){
    pred.0<-colMeans(pred.all.0)
    pred.1<-colMeans(pred.all.1)
    
    #browser()
    print("done")
    list(pred.1-pred.0,NA)
  } else {
    print("done")
    list(pred.all,NA)
  }
  
}


# New version using stacked datasets and weighted observations in an elastic net
ModCovLogitLASSO.CATE<-function(train, cvt, cve, M, y, w, x, pi, xinter=FALSE){

  
  iw<-which(colnames(train[[1]])==w)
  wname<-names(train[[1]])[iw]
  iy<-which(colnames(train[[1]])==y)

  y <- list()
  xmat <- list()
  xpred<-list()
  xcols<-which(names(train[[1]]) %in% x)
  xcols2<-c()

  for(m in 1:M){
    ds<-train[[m]]
    ds[,iw]<-ds[,iw]-0.5
    
    lastcol<-ncol(ds)
    newj<-lastcol+1
    
    if(length(x)>1){
      for(j1 in xcols){
        k<-which(xcols==j1)
        ds[,newj]<-ds[,j1]*ds[,iw]
        NewName<-paste0(names(ds)[j1],"_trans")
        names(ds)[newj]<-NewName
        if(m==1){
          #rmat<-cbind(rmat,pmax(rmat[,j1],rmat[,iw]))
          #colnames(rmat)[newj]<-NewName
          xcols2<-c(xcols2,newj)
        }
        newj=newj+1
        if(xinter & length(xcols)>k){
          for(j2 in xcols[-c(1:k)]){
            ds[,newj]<-ds[,j1]*ds[,j2]*ds[,iw]
            NewName<-paste0(names(ds)[j1],"_trans"," x ",names(ds)[j2],"_trans")
            names(ds)[newj]<-NewName
            if(m==1){
              #rmat<-cbind(rmat,pmax(rmat[,j1],rmat[,j2],rmat[,iw]))
              #colnames(rmat)[newj]<-NewName
              xcols2<-c(xcols2,newj)
            }
            newj=newj+1
          }
        }
      }
    } else {
      j1<-xcols[1]
      ds[,newj]<-ds[,j1]*ds[,iw]
      NewName<-paste0(names(ds)[j1],"_trans")
      names(ds)[newj]<-NewName
      if(m==1){
        #rmat<-cbind(rmat,pmax(rmat[,j1],rmat[,iw]))
        #colnames(rmat)[newj]<-NewName
        xcols2<-c(xcols,newj)
      }
    }

    ds0 <- ds[,c(iw,(lastcol+1):(newj-1))]
    
    if(identical(cvt,cve)){
      y[[m]]<-ds[,iy]
      xmat[[m]]<-xpred[[m]]<-as.matrix(ds0)
    } else {
      dst<-ds0[cvt,]
      dse<-ds0[cve,]
      y[[m]]<-ds[cvt,iy]
      xmat[[m]]<-as.matrix(dst)
      xpred[[m]]<-as.matrix(dse)
    }
    
  }

  p <- ncol(xmat[[1]])                
  n <- nrow(xmat[[1]])
  #ruse<-rmat[cvt,xcols2]
  #wt  <- (1 - apply(ruse,1,sum)/p)/M         # Du et al don't see advantages for this weight
  wt <- rep(1/M,length(cvt))
  pf       <- rep(1, p)
  adWeight <- rep(1, p)
  alpha    <- c(0, .5 , 1)

  xmat.all <- matrix(numeric(0),,p)
  y.all<-wt.all<-c()
  for(m in 1:M){
    xmat.all<-rbind(xmat.all,xmat[[m]])
    y.all<-c(y.all,y[[m]])
    wt.all<-c(wt.all,wt)
  }

  
  browser()
  # Zou 2009 adaptive elastic net paper uses initial betas from enet, not unpenalized
  fit <- glm.fit(xmat.all,
                 y.all,
                 weights=wt.all,
                 family=binomial(link="logit"),
                 intercept=FALSE)

  coef.mags <- abs(coef(fit))
  nu<-log(p)/log(n)
  gam<-ceiling(2*nu/(1-nu))+1
  adWeight <- (coef.mags + 1/(n*M))^(-gam)

  #browser()
  
  afit <- cv.saenet(xmat, 
                    y, 
                    pf, 
                    adWeight, 
                    wt, 
                    family = "binomial",
                    alpha = alpha, 
                    nfolds = 5,
                    intercept=FALSE)
  
  acoef<-coef(afit)[-1]

  lpmat<-matrix(,length(cve),M)
  for(m in 1:M){
    xpred0 <-xpred[[m]]
    wpred0 <-xpred0[,1]
    for(j in 1:p){
      xpred0[,j]<-xpred0[,j]/wpred0
    }
    lpmat[,m]<- (xpred0 %*% acoef)
  }

  lp.mean<-apply(lpmat,1,mean)
  cate<-(exp(lp.mean/2)-1)/(exp(lp.mean/2)+1)

  phat<-exp(lp.mean)/(1+exp(lp.mean))
  
  #browser()
  list(cate, 
       acoef,
       coef(fit),
       phat)
  
}


ModCovLogitLASSO.CATE.xinter<-function(train, cvt, cve, M, y, w, x, pi, xinter=TRUE){
  ModCovLogitLASSO.CATE(train, cvt, cve, M, y, w, x, pi, xinter=TRUE)
}



logreg.risk<-function(train, cvt, cve, M, y, w, x, pi){
  
  iw<-which(colnames(train[[1]])==w)
  
  invlogit<-function(B){
    exp(B)/(1+exp(B))
  }
  
  pred.all<-c()
  for (m in 1:M){
    
    ds<-train[[m]]
    dst<-ds[cvt,]
    i0<-which(dst[,iw]==0)
    dst<-dst[i0,]
    dse<-ds[cve,]
    
    form<-as.formula(paste0(y,"~",paste(x,collapse="+")))
    
    fit0<-glm(form,
              data=dst,
              family="binomial")
    
    pred<-predict(fit0,newdata=dse,type="link")
    pred.all<-cbind(pred.all, pred)
    
  }
  
  mean.pred<-apply(pred.all,1,mean)
  
  list(invlogit(mean.pred),NA)
  
}



rf.risk<-function(train, cvt, cve, M, y, w, x, pi, ntrees=5000, hon=TRUE){
  
  iw<-which(colnames(train[[1]])==w)
  iy<-which(colnames(train[[1]])==y)
  pred.all<-imp.all<-c()
  
  for (m in 1:M) {
    ds<-train[[m]]
    dst<-ds[cvt,]
    dse<-ds[cve,]
    same<-identical(dst,dse)
    
    if (identical(ds[cvt,],ds[cve,])){
      perm<-order(dst[,w]) #to pop back, use order(perm)
    }
    
    dst.0<-dst[dst[,w]==0,]
    dst.1<-dst[dst[,w]==1,]
    
    # Step A: Fit rf models to each arm separately, and make predictions on the opposite arm
    
    #browser()
    
    rf.0<-regression_forest(X=as.matrix(dst.0[,x]),
                            Y=dst.0[,y],
                            num.trees=ntrees,
                            honesty=hon)
    
    if (!same){
      pred<-predict(rf.0,
                         newdata = as.matrix(dse[,x]))[,1] 
    } else {
      pred<-predict(rf.0)[,1]
    }
  
    pred.all<-cbind(pred.all, pred)
    
    imp.all<-cbind(imp.all,variable_importance(rf.0))
    
  }
  
  #browser()
  list(apply(pred.all,1,mean),
       imp.all)
  
}


rf.risk.10<-function(train, cvt, cve, M, y, w, x, pi){
  
  # Define the formula argument per Foster et al 2011
  covars<-x
  formstart<-paste(covars,collapse="+")
  form<-as.formula(paste0(y,"~",paste(formstart,sep="+")))
  
  iw<-which(colnames(train[[1]])==w)
  #browser()
  rf.pred.all<-rf.imp.all<-c()
  
  #browser()
  
  for (m in 1:M) {
    ds<-train[[m]]
    #i0<-which(ds[cvt,iw]==0)
    #ds[,y]<-as.factor(ds[,y])
    dst<-ds[cvt,]
    #i0<-which(dst[,iw]==0)
    #dst<-dst[i0,]
    dse<-ds[cve,]
    same<-identical(dst,dse)
    
    #browser()
    
    rf<-randomForest(form,
                     data=dst,
                     ntree=2000,
                     #mtry=floor(2*(length(x)/3)), 
                     importance=TRUE,
                     na.action="na.omit"
    )
    
    if (!same){
      #browser()
      pred<-predict(rf,
                    newdata = dse,
                    importance=FALSE,
                    na.action="na.omit",
                    var.used=FALSE,
                    #type="prob",
                    split.depth=FALSE)
      obs_pred<-pred 
    } else {
      pred<-predict(rf,
                    importance=F,
                    na.action="na.omit",
                    var.used=F,
                    #type="prob",
                    split.depth=F)
      obs_pred<-pred 
    }
    
    # RF predictions 
    
    rf.pred<-as.numeric(obs_pred)
    rf.pred.all<-rbind(rf.pred.all, rf.pred)
    
    rf.imp<-rf$importance[,1]
    rf.imp.all<-cbind(rf.imp.all, rf.imp)
    
  }
  
  #if(nrow(rf.imp.all)<length(x)) browser()
  list(apply(rf.pred.all,2,mean),
       apply(rf.imp.all,1,mean))
  
}




rf.CATE.slearn<-function(train, cvt, cve, M, y, w, x, pi=0.5,ntrees=1000){
  
  pi<-NA # Don't need this
  
  # Define the formula argument per Foster et al 2011
  covars<-x
  formstart<-paste(c(covars,w),collapse="+")
  interactions<-interactions2<-c()
  for(i in 1:totp){
    interactions[i]<-paste0("W*",covars[i])
    interactions2[i]<-paste0("W_rev*",covars[i])
  }
  formint<-paste(interactions,collapse="+")
  formint2<-paste(interactions2,collapse="+")
  form<-as.formula(paste0(y,"~",paste(formstart,formint,formint2,sep="+")))
  
  rf.pred.all<-c()
  
  iw<-which(colnames(train[[1]])==w)
  
  for (m in 1:M) {
    ds<-train[[m]]
    colnames(ds)[iw]<-"W"
    ds[,y]<-as.factor(ds[,y])
    dst<-ds[cvt,]
    dst$W_rev<-1-dst$W
    dse<-ds[cve,]
    dse$W_rev<-1-dse$W
    same<-identical(dst,dse)
    
    rf<-randomForest(form,
                     data=dst,
                     ntree=ntrees,
                     mtry=floor(sqrt(length(x)*3+1)), 
                     importance=TRUE,
                     na.action="na.omit"
    )
    
    if (!same){
      pred<-predict(rf,
                    newdata = dse,
                    importance=FALSE,
                    na.action="na.omit",
                    var.used=FALSE,
                    split.depth=FALSE)
      obs_pred<-pred 
    } else {
      pred<-predict(rf,
                    importance=F,
                    na.action="na.omit",
                    var.used=F,
                    split.depth=F)
      obs_pred<-pred 
    }
    dse_counter<-dse
    dse_counter$W<-1-dse$W
    dse_counter$W_rev<-1-dse_counter$W
    pred_counter<-predict(rf,
                          newdata = dse_counter,
                          importance=FALSE,
                          na.action="na.omit",
                          var.used=FALSE,
                          split.depth=FALSE)
    
    counter<-ifelse(dse$W==1,1,-1)
    
    # RF predictions 
    
    counter_pred<-pred_counter 
    rf.pred<-as.numeric((obs_pred-counter_pred)*counter)
    
    rf.pred.all<-rbind(rf.pred.all, rf.pred)
    browser()
    
  }
  
  list(apply(rf.pred.all,2,mean),
       rf$importance)
  
}



