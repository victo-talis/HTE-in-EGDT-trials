
# Table showing covariate values for each quintile from the model 
# fit to the combined sample (E + F)


# First: Construct ds000 dataset from:
#D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/Programs/ITE Score approach/EGDT trials/arise_process_v2_promise_randsplit_v4.R

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
rdat<-read.csv(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/arise_process_r6.csv"))
rdat<-rdat[!is.na(rdat$y),]
rdat<-rdat[order(rdat$trial),]
rdat$site2<-paste0(rdat$site_id,rdat$trial)

llrf.pred.e2f<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_eANDf_MIglue_min20_MI100.RDS"))
preds.eANDf_v1<-apply(llrf.pred.e2f,1,mean)

ds000<-ds000_cci

ds000$iarr<-(-preds.eANDf_v1)
ds000<-ds000[order(ds000$trial_id,ds000$site,ds000$map,ds000$age,ds000$g,ds000$temp,ds000$sbp),]
rdat<-rdat[order(rdat$trial,rdat$site_id,rdat$map,rdat$age,rdat$g,rdat$temp,rdat$sbp),]
colnames(rdat)[39]<-"trial_id"
rdat<-cbind(rdat,ds000[,c("iarr","trial","cirr","chf","dial","r6_dobut","r6_blood")])




ee<-as.numeric(cut(rdat$iarr,quantile(rdat$iarr,probs=seq(0,1,0.2)),include.lowest=TRUE,labels=1:5))
xvars<-c("age","male","temp","rr","hr","map","sbp",
         "apache","sofa_total","gcs","charlson",
         "albumin","hgb","pao2","bili","bun","cr",
         "g","lac","plt","sat","wbc",
         "site_lung","site_urine","site_abdom",
         "mvent","vaso","cirr","chf","dial",
         "r6_dobut","r6_blood")
bin<-c("male","site_lung","site_urine","site_abdom",
       "mvent","vaso","cirr","chf","dial",
       "r6_dobut","r6_blood")

resmat<-matrix(,length(xvars),5)
rownames(resmat)<-xvars
colnames(resmat)<-c("Q1 stats",
                    "Q2 stats",
                    "Q3 stats",
                    "Q4 stats",
                    "Q5 stats")
prnd<-1
mrnd1<-2
sdrnd1<-2

#ds<-rdat
ds<-rdat[rdat$w==0,]
#EEuse<-ee
EEuse<-ee[rdat$w==0]
k=1
for(x in xvars){
  qqq=1
  for(QQ in 1:5){
    aa<-ds[EEuse==QQ,x]
    if(x %in% bin){ # freq (%)
      stat<-paste0(sum(aa,na.rm=TRUE)," (",round(mean(aa,na.rm=TRUE)*100,prnd),")")
      mis<-sum(is.na(aa))
    } else {  # mean (SD); median (IQR)
      stat<-paste0(#signif(mean(aa,na.rm=TRUE),mrnd1)," (",signif(sd(aa,na.rm=TRUE),sdrnd1),"); ",
        signif(median(aa,na.rm=TRUE),mrnd1)," (",signif(quantile(aa,probs=0.25,na.rm=TRUE),mrnd1),",",
        signif(quantile(aa,probs=0.75,na.rm=TRUE),mrnd1),")")
      mis<-sum(is.na(aa))
    }
    #resmat[k,qqq:(qqq+1)]<-c(stat,mis)
    resmat[k,qqq]<-stat
    qqq=qqq+1
  }
  k=k+1
}

write.csv(resmat,paste0(machine,"CATES in Critical Care/EGDT manuscript 1/CCM review/quintiles_compare_cci_w0.csv"))


# to calculate deltas

delvars<-c("cirr","chf","dial","r6_dobut","r6_blood")
resmat<-matrix(,length(delvars),5)
rownames(resmat)<-delvars
colnames(resmat)<-c("Q1 stats",
                    "Q2 stats",
                    "Q3 stats",
                    "Q4 stats",
                    "Q5 stats")
prnd<-1
mrnd1<-2
sdrnd1<-2


ds<-rdat
#ds<-rdat[rdat$w==0,]
EEuse<-ee
#EEuse<-ee[rdat$w==0]
k=1
for(x in delvars){
  qqq=1
  for(QQ in 1:5){
    aa<-ds[EEuse==QQ & ds$w==1,x]
    bb<-ds[EEuse==QQ & ds$w==0,x]
    if(x %in% bin){ # freq (%)
      stat<-paste0(sum(aa,na.rm=TRUE)-sum(bb,na.rm=TRUE)," (",round(mean(aa,na.rm=TRUE)*100 - mean(bb,na.rm=TRUE)*100,prnd),")")
      mis<-sum(is.na(aa))
    } else {  # mean (SD); median (IQR)
      stat<-paste0(#signif(mean(aa,na.rm=TRUE),mrnd1)," (",signif(sd(aa,na.rm=TRUE),sdrnd1),"); ",
        signif(median(aa,na.rm=TRUE),mrnd1)," (",signif(quantile(aa,probs=0.25,na.rm=TRUE),mrnd1),",",
        signif(quantile(aa,probs=0.75,na.rm=TRUE),mrnd1),")")
      mis<-sum(is.na(aa))
    }
    #resmat[k,qqq:(qqq+1)]<-c(stat,mis)
    resmat[k,qqq]<-stat
    qqq=qqq+1
  }
  k=k+1
}

write.csv(resmat,paste0(machine,"CATES in Critical Care/EGDT manuscript 1/CCM review/quintiles_compare_cci_del.csv"))



# Did patients with lower predictions have higher odds of blood and/or dobutamine?

fit_cirr<-glm(cirr~iarr,
               family="binomial",
               data=ds)
summary(fit_cirr)
fit_chf<-glm(chf~iarr,
              family="binomial",
              data=ds)
summary(fit_chf)
fit_dial<-glm(dial~iarr,
              family="binomial",
              data=ds)
summary(fit_dial)

fit_dobut<-glm(r6_dobut~iarr,
               family="binomial",
               data=ds)
summary(fit_dobut)
fit_blood<-glm(r6_blood~iarr,
               family="binomial",
               data=ds)
summary(fit_blood)


# Doing the above except comparing Q1 vs rest

source("C:/Users/vbtal/OneDrive - University of Pittsburgh/_BoxMigration/HTE analysis programs/michael.q.R")

pfun<-function(var){
  vv<-ds[,var]
  ww<-ds[,"w"]
  ee0<-1*(ee==1)
  
  pmat<-matrix(c(mean(vv[ww==0 & ee0==0]),
                mean(vv[ww==1 & ee0==0]),
                     mean(vv[ww==0 & ee0==1]),
                          mean(vv[ww==1 & ee0==1])),2,2)
  
  nmat<-matrix(c(length(vv[ww==0 & ee0==0]),
                 length(vv[ww==1 & ee0==0]),
                 length(vv[ww==0 & ee0==1]),
                 length(vv[ww==1 & ee0==1])),2,2)
  
  list(pmat,nmat)

}

# cirr_<-pfun("cirr")
# chf_<-pfun("chf")
# dial_<-pfun("dial")
dobut_<-pfun("r6_dobut")
blood_<-pfun("r6_blood")

michaelq(dobut_[[1]],dobut_[[2]],0.05)
michaelq(blood_[[1]],blood_[[2]],0.05)

pfun2<-function(var){
  vv<-ds[,var]
  ww<-ds[,"w"]
  ee0<-1*(ee==1)
  
  pmat<-c(sum(vv[ee0==0]),sum(vv[ee0==1]))
  
  nmat<-c(length(vv[ee0==0]),length(vv[ee0==1]))
  
  list(pmat,nmat)
  
}

cirr_<-pfun2("cirr")
chf_<-pfun2("chf")
dial_<-pfun2("dial")

prop.test(cirr_[[1]],cirr_[[2]])
prop.test(chf_[[1]],chf_[[2]])
prop.test(dial_[[1]],dial_[[2]])

