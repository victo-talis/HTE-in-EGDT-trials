

# This function takes the cv_k datasets from oest.cv and post-processes for plotting

# It has 3 arguments
# 1. cv.obj is the output from oest.cv
# 2. stats: tells the function what summary statistics to use. Currently uses match.fun
# 3. scale: tells the function whether to summarize output by rank or original scale of the score; 
#           0 means rank; integer>0 means this many quantiles of the distribution of scores will be evaluated
# 4. cores: tells the function how many cores to pass to mclapply. If cores=1, uses lapply
# 5. direction: should we find the region where outcome is greater on treatment vs control? Less? Both?


# Define the function to obtain stats for binary outcomes

stats.bin<-function(data=NULL, 
                    ntest=nrow(data),
                    alternative="less",
                    names=FALSE){
  if (names==FALSE){
    #if(nrow(data)<20)browser()
    n.w1<-length(data[data$W==1,"Y"])
    n.w0<-length(data[data$W==0,"Y"])
    x.w1<-sum(data[data$W==1,"Y"])
    x.w0<-sum(data[data$W==0,"Y"])
    ss<-mean(data$score)
    
    if (n.w1>0 & n.w0>0){
      p1<-(x.w1/n.w1)
      p0<-(x.w0/n.w0)
      ad<-p1-p0
      pp<-(n.w1*p1 + n.w0*p0)/(n.w1+n.w0)
      denom<-sqrt(pp*(1-pp)*(n.w1^-1 + n.w0^-1))
      x2<-(ad/denom)^2
      pval<-pnorm(ad/denom,mean=0, sd=1,
                  lower.tail=ifelse(alternative=="less",TRUE,FALSE))
      #browser()
      c(ad, x2, pval,(pval<0.05)*1,ss)
        
    } else {
      rep(NA,5)
    }
  } else {
    c("ad","x2","pval","pval.lt0.05","mean.score")
  }
}


# Define the function to obtain stats for binary outcomes when pi unknown
# Estimates AIPW double robust

stats.bin.unkp<-function(data=NULL, 
                         ntest=nrow(data),
                         alternative="less",
                         names=FALSE){
  if (names==FALSE){
    ww<-data$W
    yy<-data$Y
    px<-data$px
    m1x<-data$m1x
    m0x<-data$m0x
    
    if (length(yy[ww==0])>0 & length(yy[ww==1])>0){
      
      p1i<-(ww*yy/px)-(m1x*(ww-px)/px)
      p0i<-((1-ww)*yy/(1-px))-m0x*((ww-px)/(1-px))
      ad<-mean(p1i-p0i)
      
      #browser()
      
      c(ad, NA, NA, NA)
      
    } else {
      rep(NA,4)
    }
  } else {
    c("ad","x2","pval","pval.lt0.05")
  }
}


# Define the function to obtain stats for normally distributed outcome

stats.norm<-function(data=NULL, 
                    ntest=nrow(data),
                    alternative="less",
                    names=FALSE){
  if (names==FALSE){
    y1<-data[data$W==1,"Y"]
    y0<-data[data$W==0,"Y"]
    n1<-length(y1)
    n0<-length(y0)
    x.w1<-sum(y1)
    x.w0<-sum(y0)

    if (n1>0 & n0>0){
      mu1<-(x.w1/n1)
      mu0<-(x.w0/n0)
      var1<-sum((y1-mu1)^2)/(n1-1)
      var0<-sum((y0-mu0)^2)/(n0-1)
      ad<-mu1-mu0
      var_satt<-var1/n1 + var0/n0
      df<-(var_satt^2)/((var1/n1)^2/(n1-1) + (var0/n0)^2/(n0-1))
      
      t<-ad/sqrt(var_satt)
      
      pval<-pt(t,df,lower.tail=ifelse(alternative=="less",TRUE,FALSE))
      
      c(ad, t, pval,(pval<0.05)*1)
      
    } else {
      rep(NA,4)
    }
  } else {
    c("ad","t","pval","pval.lt0.05")
  }
}



post.cv.fun<-function(cv.obj,
                      stats=c("stats.bin"),
                      scale=1000,
                      cores=1,
                      direction="both",
                      score="cate"){
  
  summary.list<-list()
  
  summarize.outer<-function(cv.ds,
                          drop.high.first=TRUE,
                          stats){
    
    dat<-cv.ds[[1]]
    #if(drop.high.first==FALSE) browser()
    dat.ord<-dat[order(dat$score, decreasing=drop.high.first),]   #least to most benefit with increasing row index
    
    alt<-ifelse((drop.high.first==TRUE & score=="cate") | 
                  (drop.high.first==FALSE & score=="risk"),"less","greater")
    
    new.list<-lapply(1:nrow(dat.ord), summarize.inner, dat.ord, alt, stats)

    new.df<-as.data.frame(t(matrix(unlist(new.list),,nrow(dat.ord))))

    names(new.df)[1]<-"q"
    names(new.df)[-1]<-match.fun(stats)(names=TRUE)
    
    cbind(dat.ord,new.df)
    #browser()
    
  }
  
  
  #takes an ordered dataset and pumps out statistics for it
  summarize.inner<-function(k, ordered.dat, alternative="less", stats){
    if (k>1){
      dat.sub<-ordered.dat[-c(1:(k-1)),]
    } else {
      dat.sub<-ordered.dat
    }
    #browser()
    stats<-match.fun(stats)(dat.sub,nrow(data),alternative)
    out<-c(1-((k-1)/nrow(ordered.dat)),stats)
  }
  

  
  if (scale>0){
    
    nte<-nrow(cv.obj[[1]][[1]])
    nte.trunc<-nte-30   # need at least ~15 people per treatment arm in the smallest group tested
    
    all.scores<-unlist(lapply(1:length(cv.obj),function(x){ #concatenate scores from all CV iterations
      cv.obj[[x]][[1]][31:nte.trunc,"score"]
    }))
    #scores.q<-quantile(all.scores,prob=seq(0.02,0.98,0.96/scale))[-(scale+1)]    #produce the score quantiles corresponding to (1/scale)*100 percentiles
    scores.q<-quantile(all.scores,prob=seq(0,1,1/scale)[-1])

    #browser()
    
    if (direction=="both"){
      # gives me a list, with augmented datasets, now ordered and with statistics
      summary.list.less<-lapply(cv.obj, summarize.outer, drop.high.first=TRUE,stats)
      summary.list.greater<-lapply(cv.obj, summarize.outer, drop.high.first=FALSE,stats)

      #browser()
      
      cols<-colnames(summary.list.less[[1]])
      remove<-which(cols %in% c("Y","W","score","q"))
      stats<-cols[-remove]

      #browser()
      
      outlist<-list()
      ratelist<-list()
      k=1
      #data.frame(scores.q)
      for (dir in c(1,-1)){
        if (dir==1){
          ldat<-summary.list.less
          sq<-rev(scores.q)
        } else {
          ldat<-summary.list.greater
          sq<-scores.q
        }
        outmat_<-data.frame(sq)
        for (st in stats){
          #browser()
          cv.mean<-unlist(lapply(sq,function(x){
            out<-unlist(lapply(ldat,function(y){
              if(sum(dir*y$score<=dir*x)>0){
                y[,st][dir*y$score==max(dir*y$score[dir*y$score<=dir*x])[1]]
              } else {
                NA
              }
            }))
            mean(out,na.rm=TRUE)
          }))
          #if(dir==1 & st=="ad") browser()
          outmat_<-cbind(outmat_,cv.mean)
        }
        outlist[[k]]<-outmat_

        #browser()
        # Insert here AUTOC and Adj Qini
        # Note here is that the RATE metrics are calculated
        # using unaggregated quantiles from each CV

        ratecalc<-lapply(ldat,function(x){
          # first 
          x$score_dir<-dir*x$score
          x$grp<-as.numeric(as.character(cut(x$score_dir,
                     quantile(x$score_dir,probs=seq(0,1,1/scale)),
                     include.lowest=TRUE,
                     labels=1:scale)))
          # several statistics per cumulative bin
          # several more for each mutually exclusive bin
          # need 
          sdat<-as.data.frame(matrix(numeric(),scale,9))
          colnames(sdat)<-c("sub","Y1_c","Y0_c","N1_c","N0_c",
                                  "Y1_m","Y0_m","N1_m","N0_m")
          sdat$sub<-1:scale
          for(ii in 1:scale){
            cumsub<-x[x$grp<=ii,]
            sdat$Y1_c[ii]<-sum(cumsub$Y[cumsub$W==1])
            sdat$Y0_c[ii]<-sum(cumsub$Y[cumsub$W==0])
            sdat$N1_c[ii]<-sum(cumsub$W==1)
            sdat$N0_c[ii]<-sum(cumsub$W==0)
            mutsub<-x[x$grp==ii,]
            sdat$Y1_m[ii]<-sum(mutsub$Y[mutsub$W==1])
            sdat$Y0_m[ii]<-sum(mutsub$Y[mutsub$W==0])
            sdat$N1_m[ii]<-sum(mutsub$W==1)
            sdat$N0_m[ii]<-sum(mutsub$W==0)
          }
          
          # Calculate AUQINI
          # per Melbahri arxiv paper
          nb<-scale
          phi_j<-sdat$N1_c/sdat$N1_c[nb]
          gphi_j<-(sdat$Y1_c - sdat$Y0_c*(sdat$N1_c/sdat$N0_c))/sdat$N1_c[nb]     # "incremental uplift"
          Qphi_j<-gphi_j-phi_j*gphi_j[nb]                          # qini curve
          qini_area<-(1/2)*sum((phi_j[-1]-phi_j[-nb])*(Qphi_j[-1]+Qphi_j[-nb]))
          
          # Calculate AUTOC
          # (per yadlowski arxiv paper)
          gate<-sdat$Y1_c/sdat$N1_c - sdat$Y0_c/sdat$N0_c
          u_j<-seq(1/scale,1,1/scale)
          autoc<-(1/2)*sum((u_j[-1]-u_j[-nb])*(gate[-1]+gate[-nb]))
          
          # Calculate Kendall tau
          tau<-cor(gate,1:scale,method="kendall")
          
          # Calculate adjusted AUQINI
          adjqini<-tau*qini_area

          c(qini_area,autoc,tau,adjqini)
          
          #browser()
          
        })
        
        #browser()
        ratelist[[k]]<-as.data.frame(t(matrix(unlist(ratecalc),4,length(ldat))))
        colnames(ratelist[[k]])<-c("AUQINI","AUTOC","Kendall Tau","Adj AUQINI")
        
        if(dir==1){
          names(ratelist)[k]<-"Drop High First"
        } else {
          names(ratelist)[k]<-"Drop Low First"
        }
        
        #browser()
        
        k=k+1
      }
      #browser()
      
      drop_high_names<-unlist(lapply(stats,function(x){
        paste0(x,"_DropHighFirst")
      }))
      drop_low_names<-unlist(lapply(stats,function(x){
        paste0(x,"_DropLowFirst")
      }))
      colnames(outlist[[1]])<-c("scores.q_DropHighFirst",drop_high_names)
      colnames(outlist[[2]])<-c("scores.q_DropLowFirst",drop_low_names)
      
      outmat<-cbind(outlist[[1]],outlist[[2]])

      
    } else if (direction=="greater"){
      summary.list.less<-NA
      summary.list.greater<-lapply(1:length(cv.obj), summarize.outer, drop.high.first=FALSE)
    } else {
      summary.list.less<-lapply(1:length(cv.obj), summarize.outer, drop.high.first=TRUE)
      summary.list.greater<-NA
    }

  } else if (scale==0){
    if (direction=="both"){
      summary.list.less<-lapply(cv.obj, summarize.outer, drop.high.first=TRUE,stats)
      summary.list.greater<-lapply(cv.obj, summarize.outer, drop.high.first=FALSE,stats)
      
      #browser()
      
      cols<-colnames(summary.list.less[[1]])
      remove<-which(cols %in% c("Y","W","score","q"))
      stats<-cols[-remove]
      
      outmat<-data.frame(q=summary.list.less[[1]][,"q"])
      for (dir in c(1,-1)){
        if (dir==1){
          ldat<-summary.list.less
          #sq<-rev(scores.q)
        } else {
          ldat<-summary.list.greater
          #sq<-scores.q
        }
        for (st in stats){
          alldat<-ldat[[1]]
          if (length(ldat)>1){
            for (i in 2:length(ldat)){
              #browser()
              alldat<-rbind(alldat,ldat[[i]])
            }
          }
          #browser()
          cv.mean<-tapply(alldat[,st],alldat[,"q"],mean,na.rm=TRUE)
          cv.mean<-cv.mean[order(names(cv.mean),decreasing=TRUE)]
          outmat<-cbind(outmat,cv.mean)
        }
      }
      
      benefit_names<-unlist(lapply(stats,function(x){
        paste0(x,"_DropHighFirst")
      }))
      harm_names<-unlist(lapply(stats,function(x){
        paste0(x,"_DropLowFirst")
      }))
      colnames(outmat)[-1]<-c(benefit_names,harm_names)
    }
    
    ratelist<-NA
    
  }
  
  list("outmat"=outmat,
       "ratelist"=ratelist)
  
}
