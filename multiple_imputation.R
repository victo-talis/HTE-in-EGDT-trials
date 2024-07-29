


devtools::load_all(paste0(machine,"_BoxMigration/HTE analysis programs/mice-master/mice-master"),recompile=TRUE)

mi_fun<-function(dat,             # data to be multiply imputed (ds000)
                 split_var_name,  # variable name to split (e.g. derivation/testing)
                 trset,           # value fo split_var_name to use as derivation
                 teset,           # value of split_var_name to use as validation
                 mi.use,          # number of MI datasets to generate
                 usey=FALSE,      #condition x models on outcome?
                 xvars){          # character vector of covariates to use in MICE
  
  xuse<-xvars
  
  ds<-dat[dat[,split_var_name] %in% c(trset,teset),c("y","w",xuse,split_var_name)]
  rownames(ds)<-1:nrow(ds)
  
  tr.ids<-which(ds[,split_var_name] %in% trset)
  #tr.ids<-1:nrow(ds)
  
  # Need to perform variable transformations
  
  ds$ln_bili<-log(ds$bili)
  ds$ln_sbp<-log(ds$sbp)
  ds$ln_map<-log(ds$map)
  
  ds$ln_cr<-log(ds$cr)
  
  ds$ln_lac<-log(ds$lac)
  ds$ln_plt<-log(ds$plt)
  ds$ln_sat<-log(ds$sat)
  ds$sqrt_wbc<-sqrt(ds$wbc)
  ds$ln_pao2<-log(ds$pao2)
  ds$sqrt_sofa<-sqrt(ds$sofa_total)
  
  
  lm<-c("age","rr","hr","hgb","ln_sbp","ln_bili","ln_map",
        "ln_bun","ln_cr","sqrt_wbc","apache",
        "ln_pao2","charlson",
        "gcs","sqrt_sofa",
        "albumin","ln_lac","ln_g","ln_plt","ln_sat",
        "temp")
  
  if(any(c("bun","g") %in% xuse)){
    ds$ln_g<-log(ds$g)
    ds$ln_bun<-log(ds$bun)
    
    lm<-c("age","rr","hr","hgb","ln_sbp","ln_bili","ln_map",
          "ln_bun","ln_cr","sqrt_wbc","apache",
          "ln_pao2","charlson",
          "gcs","sqrt_sofa",
          "albumin","ln_lac","ln_g","ln_plt","ln_sat",
          "temp")
  } else {
    lm<-c("age","rr","hr","ln_sbp","ln_bili","ln_map",
          "ln_cr","sqrt_wbc","apache",
          "ln_pao2",
          "gcs","sqrt_sofa",
          "ln_lac","ln_plt","ln_sat",
          "temp")
  }
  
  logreg<-c("y","w","male","site_lung","site_urine","site_abdom",
            "mvent","vaso")
  
  
  #rf<-c("lac","sat","ln_g","pao2","charlson","gcs","sofa_total")
  
  #xvars<-c(lm,logreg[-c(1:2)],rf)
  xvars<-c(lm,logreg[-c(1:2)])
  
  
  ds<-ds[,c("y","w",xvars)]
  n<-nrow(ds)
  
  
  methodds<-data.frame(var=names(ds),meth=NA)
  methodds$meth[which(methodds$var %in% lm)]<-"norm_mod"
  methodds$meth[which(methodds$var %in% logreg)]<-"logreg_mod"
  
  remove<-removen<-c()
  remov.mi.vars<-c()
  
  # training set T/F to apply to ds
  logic.tr<-rep(FALSE,n)
  logic.tr[tr.ids]<-TRUE
  # Training set imputations: incorporate Y and W
  logic.tr<-rep(FALSE,n)
  logic.tr[tr.ids]<-TRUE
  dsw<-ds[logic.tr,]
  pp<-order(dsw$w)
  dsw<-dsw[pp,]
  rownames(dsw)<-1:nrow(dsw)
  # logic.tr.w<-logic.tr[pp]
  dsw0<-dsw[dsw$w==0,]
  # logic.tr.w0<-logic.tr.w[dsw$w==0]
  dsw1<-dsw[dsw$w==1,]
  # logic.tr.w1<-logic.tr.w[dsw$w==1]
  rownames(dsw0)<-1:nrow(dsw0)
  rownames(dsw1)<-1:nrow(dsw1)
  mi0<-mice_mod(dsw0[,-c(2)],
                m=mi.use,
                method=methodds$meth[-c(2)],
                tr=1:nrow(dsw0))
  mi1<-mice_mod(dsw1[,-c(2)],
                m=mi.use,
                method=methodds$meth[-c(2)],
                tr=1:nrow(dsw1))
  
  if(teset==""){ #don't need to do fancy combining 
    midat_list0<-lapply(1:mi.use,function(x){
      xout<-rbind(complete(mi0,action=x),complete(mi1,action=x))[order(pp),-1]
      out<-cbind(ds[,c(1:2)],xout)
      for(i in logreg){
        if(length(removen)==0){
          out[,i]<-as.numeric(as.character(out[,i]))
        } else if (!(i %in% removen)){
          out[,i]<-as.numeric(as.character(out[,i]))
        }
      }
      rownames(out)<-1:nrow(out)
      out
    })
  } else { # need to make sure datasets combined in same way
    dsw<-ds
    if(usey){
      remy=2
    } else {
      remy=1:2
    }
    mib<-mice_mod(dsw[,-remy],
                  m=mi.use,
                  method=methodds$meth[-remy],
                  tr=logic.tr)
    
    midat_list0<-lapply(1:mi.use,function(x){
      xout_tr<-rbind(complete(mi0,action=x),complete(mi1,action=x))[order(pp),-1]
      if(usey){
        xout_te<-complete(mib,action=x)[-tr.ids,-1]
      } else {
        xout_te<-complete(mib,action=x)[-tr.ids,]
      }
      #browser()
      if(tr.ids[1]==1){ #standardize the order of patients in the output dataset
        xout<-rbind(xout_tr,xout_te)
      } else {
        xout<-rbind(xout_te,xout_tr)
      }
      out<-cbind(ds[,c(1:2)],xout)
      for(i in logreg){
        if(length(removen)==0){
          out[,i]<-as.numeric(as.character(out[,i]))
        } else if (!(i %in% removen)){
          out[,i]<-as.numeric(as.character(out[,i]))
        }
      }
      rownames(out)<-1:nrow(out)
      out
    })
  }
  
  midat_list0
}