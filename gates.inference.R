
library("ic.infer")
library("scales")
library("sandwich")

## General functions for use

cmat_maker<-function(num.q){
  mat<-matrix(0,num.q-1,num.q)
  for(r in 1:nrow(mat)){
    mat[r,r]<-(-1)
    mat[r,r+1]<-1
  }
  mat
}

# Should modify the test function to output (y,w,score,b,s), where b is rf prediction and s is cf prediction
pval_constrain<-function(gates,type,num.q,vcov.type){
  
  if (type==1){   # at least one strict inequality
    #browser()
    ui1<-cmat_maker(num.q)
    gatesest<-ic.est(coef(gates)[2:(num.q+1)],
                     ui=ui1,
                     Sigma=vcovHC(gates,vcov.type)[2:(num.q+1),2:(num.q+1)])
    ic.test(gatesest)$p.value
  } else if (type==2){   # inequality between ATEs in extreme bins
    ui2<-t(matrix(c(-1,1),2,1))
    
    gatesest2<-ic.est(coef(gates)[c(2,(num.q+1))],
                      ui=ui2,
                      Sigma=vcovHC(gates,vcov.type)[c(2,(num.q+1)),c(2,(num.q+1))])
    ic.test(gatesest2)$p.value
  } 
}



gates.fun<-function(master,pi,num.q,loo=FALSE,irev=NULL,lintest=TRUE,unknown.pi=FALSE,
                    alpha=0.05,q.breaks=NULL,vcov.type="HC3"){
  
  if(!is.null(q.breaks)){
    if(num.q != (length(q.breaks)-1)) stop("Make sure that length(q.breaks)-1 equals num.q.")
  }
  
  thetalist<-ulist<-llist<-matrix(,length(master),num.q)
  pval_monoton<-pval_extreme<-pval_lintest<-c()
  pval_lintest<-rep(0,length(master))
  
  bb_all<-matrix(,length(master),num.q+1)
  n0_all<-matrix(,length(master),num.q)
  N0_all<-matrix(,length(master),num.q)
  n1_all<-matrix(,length(master),num.q)
  N1_all<-matrix(,length(master),num.q)
  
  
  for (i in 1:length(master)){

    t<-master[[i]]
    
    if (i %in% irev) {
      t$score<- (t$score*-1)  #so that the lowest bins always have the people we think stand the most to gain
    }

    if(length(unique(quantile(t$score,probs=seq(0,1,1/num.q))))<(num.q+1)) {
      #browser()
      t$score<-rnorm(nrow(t),0,0.001) 
    }

    if(!is.null(q.breaks)){
      bb<-q.breaks
    } else {
      bb<-quantile(t$score,probs=seq(0,1,1/num.q))
      bb_all[i,]<-bb
    }

    t$e<-as.numeric(cut(t$score, 
                        breaks=bb,
                        names=FALSE,
                        include.lowest=TRUE))
    
    n0_all[i,]<-tapply(t$y[t$w==0],t$e[t$w==0],sum)
    N0_all[i,]<-tapply(t$y[t$w==0],t$e[t$w==0],length)
    n1_all[i,]<-tapply(t$y[t$w==1],t$e[t$w==1],sum)
    N1_all[i,]<-tapply(t$y[t$w==1],t$e[t$w==1],length)
    
    # If doing internal validation (loo=FALSE), then z quantile should be 2.28
    # Otherwise use the standard 1.96
    if(!loo){
      cov.coef<-qnorm(1-alpha/4)
    } else {
      cov.coef<-qnorm(1-alpha/2)
    }

    if(!unknown.pi){
      eis<-c()
      for(enum in 1:num.q){
        t$new.e<-ifelse(t$e==enum,1,0)*(t$w-pi)
        ei<-which(colnames(t)=="new.e")
        eis<-c(eis,ei)
        colnames(t)[ei]<-paste0("e",enum)
      }
      iy<-which(colnames(t)=="y")
      tg<-t[,c(iy,eis)]
      
      #Effects coding suggested in Chernozhukov 2018. Note slightly different betas 
      #   compared  with usual dummy coding. These differences completely go away if, 
      #   instead of subtracting 0.5 in the terms above, subtract the estimated propensity
      #   That said, the substantive results of the analysis do not change over the 100 splits.
      gates<-lm(y~.,data=tg)  
      #if(i %in% c(20,40,60,80)) browser()
      thetalist[i,]<-coef(gates)[2:(1+num.q)]
      ulist[i,]<-coef(gates)[2:(1+num.q)]+cov.coef*sqrt(diag(vcovHC(gates,vcov.type))[-1])
      llist[i,]<-coef(gates)[2:(1+num.q)]-cov.coef*sqrt(diag(vcovHC(gates,vcov.type))[-1])
      
      pval_monoton[i]<-pval_constrain(gates,1,num.q,vcov.type)
      pval_extreme[i]<-pval_constrain(gates,2,num.q,vcov.type)
      
    } else {
      #browser()
      ptlist<-lapply(1:num.q,function(x){
        dat<-t[t$e==x,]
        yy<-dat$y
        ww<-dat$w
        pt<-prop.test(c(sum(yy[ww==0]),sum(yy[ww==1])),
                  c(sum(ww==0),sum(ww==1)))
        list(pt$estimate[2]-pt$estimate[1],
             -pt$conf.int[2],
             -pt$conf.int[1])
      })
      
      mat<-matrix(unlist(ptlist),3,num.q)
      
      thetalist[i,]<-mat[1,]
      ulist[i,]<-mat[3,]
      llist[i,]<-mat[2,]
      
      # would need to modify the pval_constrain code to admit the above lm output
      pval_monoton[i]<-NA
      pval_extreme[i]<-NA
      
    }
    
    

    # Also do linear test of calibration?
    if(lintest){
      #browser()
      newt<-t #[!(t$e %in% c(2,3)),]
      score <- newt$score
      mean.score <- mean(score)
      # DF <- data.frame(y = newt$y-newt$b,  
      #                  main.score.term = (newt$w - pi) * mean.score, 
      #                  interact.score.term = (newt$w - pi) * (score - mean.score))
      
      # fit.blp <- lm(y ~ b + main.score.term + interact.score.term , 
      #               data = DF)
      DF <- data.frame(y = newt$y,
                       b = newt$b,
                       main.score.term = (newt$w - pi),
                       interact.score.term = (newt$w - pi) * (score-mean.score))
      
      fit.blp <- lm(y ~ b + main.score.term + interact.score.term , 
                    data = DF)
      
      #browser()
      #fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
      #fit.blp<-lm(y~w+score+w*score,data=t)
      blp.summary <- lmtest::coeftest(fit.blp, 
                                      vcov = sandwich::vcovCL, 
                                      type = "HC3")
      #browser()
      dimnames(blp.summary)[[2]][4] <- gsub("[|]", "", dimnames(blp.summary)[[2]][4])
      
      blp.summary[, 4] <- ifelse(blp.summary[, 3] < 0, 1 - blp.summary[,4]/2, blp.summary[,4]/2)
      #browser()
      pval_lintest[i]<-blp.summary[4,4]
      #coef_lintest[i]<-blp.summary[4,4]
      
    }
    
    #browser()
  }

  # The original vein result is predicated on having a set of p-values which are random conditional on the training set
  # if we generate a single p-value from the LOO estimates, we have a non-random p-value
  
  #browser()
  
  list("Theta Median" = apply(thetalist,2,median),
       "Theta Upper Median" = apply(ulist,2,median),
       "Theta Lower Median" =apply(llist,2,median),
       "Pval Monoton" = median(pval_monoton)*(2-loo*1),
       "Pval Extreme" = median(pval_extreme)*(2-loo*1),
       "Theta P25" = apply(thetalist,2,quantile,probs=0.25),
       "Theta P75" = apply(thetalist,2,quantile,probs=0.75),
       "Pval BLP test" = median(pval_lintest)*(2-loo*1),
       "Pval BLP coef intercept"= median(blp.summary[3,1]),
       "Pval BLP coef slope" = median(blp.summary[4,1]),
       "N quantile" = table(t$e),
       "QUantile Thresholds All" = bb_all,
       "Mean prediction in quantiles" = tapply(t$score,t$e,mean),
       "n0 quantile" = n0_all,
       "N0 quantile" = N0_all,
       "n1 quantile" = n1_all,
       "N1 quantile" = N1_all)
}


