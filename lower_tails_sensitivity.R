


llrf.pred.f2e<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_f2e_v2.RDS"))
iarre<- -apply(llrf.pred.f2e,1,mean)
llrf.pred.f2e<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_e2f_v2.RDS"))
iarrf<- -apply(llrf.pred.e2f,1,mean)

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
rdat<-read.csv(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/arise_process_r6.csv"))
rdat<-rdat[!is.na(rdat$y),]
rdat<-rdat[order(rdat$trial),]
rdat$site2<-paste0(rdat$site_id,rdat$trial)

ds000<-ds000_cci

ds000$iarre<-ds000$iarrf<-NA
ds000<-ds000[order(ds000$trial_id,ds000$site,ds000$map,ds000$age,ds000$g,ds000$temp,ds000$sbp),]
rdat<-rdat[order(rdat$trial,rdat$site_id,rdat$map,rdat$age,rdat$g,rdat$temp,rdat$sbp),]
colnames(rdat)[39]<-"trial_id"
rdat<-cbind(rdat,ds000[,c("iarre","iarrf","trial","cirr","chf","dial","r6_dobut","r6_blood")])
rdat$iarr[rdat$trial=="E"]<-iarre
rdat$iarr[rdat$trial=="F"]<-iarrf
rdat$low_grp<-NA
rdat$low_grp[rdat$trial=="E" & iarre<quantile(iarre,probs=0.1)]<-1
rdat$low_grp[rdat$trial=="F" & iarrf<quantile(iarrf,probs=0.1)]<-2


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

resmat<-matrix(,length(xvars),2)
rownames(resmat)<-xvars
colnames(resmat)<-c("Cohort A",
                    "Cohort B")
prnd<-1
mrnd1<-2
sdrnd1<-2

ds<-rdat
#EEuse<-ee
EEuse<-rdat$low_grp
k=1
for(x in xvars){
  qqq=1
  for(QQ in 1:2){
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

write.csv(resmat,paste0(machine,"CATES in Critical Care/EGDT manuscript 1/CCM review/lowtails_compare.csv"))
