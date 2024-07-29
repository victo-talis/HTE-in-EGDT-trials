
# OEST test set confirmatory analysis for binary outcomes (via linear model analysis)
# Takes a list of training sets (e.g. from MICE), and a list of test sets (same)
# Also takes a pre-selected threshold value, cval, to apply to the test set
# If c_from_q=TRUE, then the threshold is determined from the supplied quantile q (optional)


confirm_analysis<-function(train_list,
                           test_list,
                           mod_ben="cf.CATE",
                           mod_harm="cf.CATE",
                           cval_ben,
                           cval_harm,
                           c_from_q=FALSE,
                           qval=NULL,
                           y="Y",
                           w="W",
                           xvars,
                           pi=0.5,
                           type="CATE",
                           out="preds",
                           ...){
  

  all_list<-list()
  for (i in 1:length(train_list)){
    all_list[[i]]<-rbind(train_list[[i]],test_list[[i]])
    rownames(all_list[[i]])<-1:(nrow(train_list[[1]])+nrow(test_list[[1]]))
  }
  
  ii<-1:nrow(train_list[[1]])
  jj<-nrow(train_list[[1]])+c(1:nrow(test_list[[1]]))
  
  #browser()
  ben_out<-match.fun(mod_ben)(all_list, ii, jj, 
                     M=length(train_list), y=y, w=w, x=xvars, pi=pi,
                     ...)
  if (mod_harm!=mod_ben){
    harm_out<-match.fun(mod_harm)(all_list, ii, jj, 
                                  M=length(train_list), y=y, w=w, x=xvars, pi=pi)
  } else {
    harm_out<-ben_out
  }
  
  
  if (out=="preds"){

    if(type %in% c("rf.risk.DropLowFirst")){
      benscore= -ben_out[[1]]
    } else {
      benscore= ben_out[[1]]
    }
    
    b<-match.fun("rf.risk.10")(all_list, ii, jj, 
                                M=length(train_list), y=y, w=w, x=xvars, pi=pi)
    stest<-match.fun("cf.CATE")(test_list, 1:nrow(test_list[[1]]), 1:nrow(test_list[[1]]), 
                            M=length(train_list), y=y, w=w, x=xvars, pi=pi)
    spitout<-data.frame(y=test_list[[1]][,y],
               w=test_list[[1]][,w],
               score=benscore,
               b=b[[1]],    # in case I want to control for baseline risk (and/or estimated CATEs) when I estimate GATEs
               s=stest[[1]])
    #browser()
    return(list(spitout,
                ben_out[[2]]))
  }
  
  if (!c_from_q){
    sub_ben<-ifelse(type=="CATE" & ben_out[[1]]<cval_ben,1,
                    ifelse(type=="risk" & ben_out[[1]]>cval_ben,1,0))
    sub_harm<-ifelse(type=="CATE" & harm_out[[1]]>cval_harm,1,
                    ifelse(type=="risk" & harm_out[[1]]<cval_harm,1,0))
    sub_indet<-ifelse( sub_ben==sub_harm,1,0)
    browser()
  }
  
  yw<-test_list[[1]][,c(y,w)]
  
  ben_pval<-prop.test(x=c(sum(sub_ben*yw[,2]*yw[,1]), sum(sub_ben*(1-yw[,2])*yw[,1])),
                      n=c(sum(sub_ben*yw[,2]), sum(sub_ben*(1-yw[,2]))),
                      alternative="less")$p.value
  
  
  indet_pval<-try(prop.test(x=c(sum(sub_indet*yw[,2]*yw[,1]), sum(sub_indet*(1-yw[,2])*yw[,1])),
                      n=c(sum(sub_indet*yw[,2]), sum(sub_indet*(1-yw[,2]))),
                      alternative="two.sided")$p.value,silent=TRUE)
  indet_pval<-ifelse(is.numeric(indet_pval),indet_pval,NA)
  harm_pval<-prop.test(x=c(sum(sub_harm*yw[,2]*yw[,1]), sum(sub_harm*(1-yw[,2])*yw[,1])),
                        n=c(sum(sub_harm*yw[,2]), sum(sub_harm*(1-yw[,2]))),
                        alternative="greater")$p.value
  
  #browser()
  
  resmat<-matrix(,3,7)
  resmat[1,]<-c(sum(sub_ben),sum(sub_ben*yw[,1]*yw[,2]),sum(sub_ben*yw[,1]*yw[,2])/sum(sub_ben*yw[,2]),
                sum(sub_ben*yw[,1]*(1-yw[,2])),sum(sub_ben*yw[,1]*(1-yw[,2]))/sum(sub_ben*(1-yw[,2])),
                sum(sub_ben*yw[,1]*yw[,2])/sum(sub_ben*yw[,2])-sum(sub_ben*yw[,1]*(1-yw[,2]))/sum(sub_ben*(1-yw[,2])),
                ben_pval)
                
  resmat[2,]<-c(sum(sub_indet),sum(sub_indet*yw[,1]*yw[,2]),sum(sub_indet*yw[,1]*yw[,2])/sum(sub_indet*yw[,2]),
                sum(sub_indet*yw[,1]*(1-yw[,2])),sum(sub_indet*yw[,1]*(1-yw[,2]))/sum(sub_indet*(1-yw[,2])),
                sum(sub_indet*yw[,1]*yw[,2])/sum(sub_indet*yw[,2])-sum(sub_indet*yw[,1]*(1-yw[,2]))/sum(sub_indet*(1-yw[,2])),
                indet_pval)
  
  resmat[3,]<-c(sum(sub_harm),sum(sub_harm*yw[,1]*yw[,2]),sum(sub_harm*yw[,1]*yw[,2])/sum(sub_harm*yw[,2]),
                sum(sub_harm*yw[,1]*(1-yw[,2])),sum(sub_harm*yw[,1]*(1-yw[,2]))/sum(sub_harm*(1-yw[,2])),
                sum(sub_harm*yw[,1]*yw[,2])/sum(sub_harm*yw[,2])-sum(sub_harm*yw[,1]*(1-yw[,2]))/sum(sub_harm*(1-yw[,2])),
                harm_pval)

  colnames(resmat)<-c("N","freq_death_W1","perc_death_W1",
                      "freq_death_W0","perc_death_W0",
                      "perc_diff","pval")
  rownames(resmat)<-c("benefit","indet","harm")
  
  resmat
}






