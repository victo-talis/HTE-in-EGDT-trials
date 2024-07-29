

# Construct the main cohort for analysis


#machine<-"D:/Box Sync/Box Sync/Box Sync/"
#machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
prism.dir<-"_BoxMigration/Faraaz PROWESS/Programs/ITE Score approach/PRISM/Dataset/"
proc.dir<-"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/ProCESS/"

# Read in PRISM

dsp<-read.csv(paste0(machine,prism.dir,"prism_data_scvo2.csv"))
dsp<-dsp[dsp$trial_id %in% c("A","C"),]
dspc<-dsp[dsp$trial_id %in% c("C"),]
pris.xwalk<-read.csv(paste0(machine,prism.dir,"pris_xwalk.csv"))[,-1]




# Read in ProCESS

setwd(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/ProCESS/"))

vars.shu<-c("stnum","age","male","rr_pre","temp_pre","hr_pre","charlson",
            "albumin","chlor","hgb","pao2","gcs_pre","ln_bands","ln_bili",
            "ln_bun","ln_cr","ln_g","ln_lac","ln_plt","ln_sat","ln_sbp",
            "ln_wbc","na")
vars.pet1<-c("stnum","InHospCens60","BASELINE__APACHEII",
             "Cens90","Mort1Yr",
             "PREENROLLMENT__dbp","PREENROLLMENT__mechvent",
             "PREENROLLMENT__sbp",
             "PREENROLLMENT__urineout","BASELINE__SOFA_Total")

#             "BASELINE__SOFA_Cardiac","BASELINE__SOFA_CNS",
#             "BASELINE__SOFA_Coag","BASELINE__SOFA_Liver",
#             "BASELINE__SOFA_Renal",

vars.pet2<-c("stnum","arm","race","infectsite","Hosp_LOS",
             "rendis","achf","cirrhos")

ds.shu<-read.csv("OLD/ProCESS - Modeling Variables.csv")[,vars.shu]
ds.pet1<-read.csv("OLD/ProCESS-Petros_FromMasterDataFile2_August19.csv")[,vars.pet1]
ds.pet1$dbp<-as.numeric(ds.pet1$PREENROLLMENT__dbp)
ds.pet1$mvent<-as.numeric(ds.pet1$PREENROLLMENT__mechvent)
ds.pet1$mvent<-ifelse(ds.pet1$mvent==3,NA,ds.pet1$mvent)
ds.pet1$urine<-as.numeric(ds.pet1$PREENROLLMENT__urine)
ds.pet1$dial<-as.numeric(ds.pet2$rendis==1)
ds.pet2<-read.csv("OLD/ProCESS Data Request 10 20 2015.csv")[,vars.pet2]
ds.pet2$race<-ifelse(ds.pet2$race==10,NA,ds.pet2$race)
ds.pet2$race<-ifelse(ds.pet2$race==1,1,
                     ifelse(ds.pet2$race==2,2,
                            ifelse(is.na(ds.pet2$race),NA,3)))
ds.pet2$site_lung<-ifelse(ds.pet2$infectsite==1,1,0)
ds.pet2$site_abdom<-ifelse(ds.pet2$infectsite==2,1,0)
ds.pet2$site_urine<-ifelse(ds.pet2$infectsite==3,1,0)
ds.pet2$dial<-as.numeric(ds.pet2$rendis==1)
ds.pet2$chf<-as.numeric(ds.pet2$achf==1)
ds.pet2$cirr<-as.numeric(ds.pet2$cirrhos==1)

ds<-merge(ds.shu,ds.pet1,by="stnum")
ds<-merge(ds,ds.pet2,by="stnum")
ds$w<-ifelse(ds$arm==1,1,
             ifelse(ds$arm==2,-1,0))
ds$y<-ds$InHospCens60
ds<-ds[,-which(colnames(ds) %in% c("arm","InHospCens60"))]
ds<-ds[ds$w %in% c(0,1),]
rownames(ds)<-1:nrow(ds)



# read in Jason's .dta file

library(haven)
dta<-read_dta(paste0(machine,proc.dir,"Files from Jason June 2022/ProCESS - Full Clean - Final.dta"))


# Here starts vetted data

dsc<-read.csv(paste0(machine,proc.dir,"Files from BDMC/ProCESS_Master_Data_File_SAFE_HARBORED.csv"))
dsc$death_90<-ifelse(dsc$Mort90_NewNDI_09022015<90 & dsc$Cens90_NewNDI_09022015==1,1,
                     ifelse(dsc$Mort90_NewNDI_09022015==90 & dsc$Cens90_NewNDI_09022015==0,0,NA))

dsc2<-read.csv(paste0(machine,proc.dir,"Files from Jason June 2022/ProCESS NEJM TS2.csv"))
dsc2<-dsc2[dsc2$arm %in% c(1,3),c(1,3:19)]

dsc3<-merge(dsc,dsc2,by="stnum",all=FALSE)
dsc3$gcs<-apply(dsc3[,c(139:141)],1,function(x)sum(as.numeric(x)))

# use this for some labs,
dscl<-read.csv(paste0(machine,proc.dir,"Files from Jason June 2022/ProCESS_Cohort_All_Labs_Table.csv"))
dscl<-dscl[dscl$arm %in% c(1,3),]
dscl<-as.data.frame(lapply(dscl,as.numeric))
dscl<-dscl[order(dscl$stnum),]

# Notes: pao2 from shu's dataset contained values from non-baseline timepoints


# Need to get directly from dscl: albumin, bands, bili, lactate
# Also need to check the *.0 variables against my extraction method for dscl; these may exclude negative timestamps?

# which_lab_alb<-function(id,var,tol=0){
#   d0<-dta$albumin[dta$stnum==id]
#   dat<-na.omit(dscl[dscl$stnum==id,c("daidx",var)])
#   if(nrow(dat)==0 & is.na(d0)){
#     -77777
#   } else if (nrow(dat)==0){
#     -88888
#   } else if (nrow(dat)>0 & is.na(d0)){
#     browser()
#     -99999
#   } else {
#     wh<-which(dat[,2]>=d0-tol & dat[,2]<=d0+tol)
#     dat[wh,1]
#   }
# }
# 
# alb_hr<-unlist(lapply(unique(dscl$stnum),which_lab_alb,"albumin",0))
# dta$albumin[which(alb_hr==-99999)]

labs_bl<-function(id,var,min_incl=-24,max_incl=6){
  #browser()
  dat<-na.omit(dscl[dscl$stnum==id,c("daidx",var)])
  dat<-dat[dat[,1]<=max_incl & dat[,1]>=min_incl,]
  if(nrow(dat)==0){
    NA
  } else {
    best<-which.min(abs(dat[,1])) #should give the closest value to 0, perferring negatives
    as.numeric(dat[best,2])
  }
}

# best_lab_bl<-function(id,var){
#   dat<-na.omit(dscl[dscl$stnum==id,c("daidx",var)])
#   if(nrow(dat)==0){
#     NA
#   } else {
#     best<-which.min(abs(dat[,1])) #==min(abs(dat[,1])))
#     dat[best,1]
#   }
# }
# alb_best<-unlist(lapply(unique(dscl$stnum),best_lab_bl,"albumin"))
# table(cut(alb_best,breaks=c(-99999999,-12,-1,0,4,99999)),useNA="always")

albumin<-unlist(lapply(unique(dscl$stnum),labs_bl,"albumin",-12,6))
bands<-unlist(lapply(unique(dscl$stnum),labs_bl,"bands",-12,6))
bili<-unlist(lapply(unique(dscl$stnum),labs_bl,"bili",-12,6))
lactate<-unlist(lapply(unique(dscl$stnum),labs_bl,"lactate",-12,6))
creat<-unlist(lapply(unique(dscl$stnum),labs_bl,"creat",-12,6))

chlor<-unlist(lapply(unique(dscl$stnum),labs_bl,"chlor",-12,6))
hgb<-unlist(lapply(unique(dscl$stnum),labs_bl,"Hgb",-12,6))
pao2<-unlist(lapply(unique(dscl$stnum),labs_bl,"pao2",-12,6))
bun<-unlist(lapply(unique(dscl$stnum),labs_bl,"bun",-12,6))
gluc<-unlist(lapply(unique(dscl$stnum),labs_bl,"gluc",-12,6))
plt<-unlist(lapply(unique(dscl$stnum),labs_bl,"plate",-12,6))
wbc<-unlist(lapply(unique(dscl$stnum),labs_bl,"WBC",-12,6))

dscl0<-data.frame(stnum=unique(dscl$stnum),
                  albumin,bili,lac=lactate,
                  cr=creat,hgb,pao2,bun,g=gluc,
                  plt,wbc)

# 
# # usig a window of -12 to +6 hrs, we wind up with just a few fewer NAs for each variable
# summary(chlor); summary(ds$chlor)
# summary(creat); summary(exp(ds$ln_cr))
# summary(hgb); summary(exp(ds$hgb))
# summary(pao2); summary(exp(ds$pao2))
# summary(bun); summary(exp(ds$ln_bun))
# summary(gluc); summary(exp(ds$ln_g))
# summary(plt); summary(exp(ds$ln_plt))
# summary(wbc); summary(exp(ds$ln_wbc))


dsc4<-dsc3[dsc3$arm %in% c(1,3),c("stnum","death_90","arm","age_at_enrollment","SCREENING__sex",
                                  "PREENROLLMENT__resp","PREENROLLMENT__temp","PREENROLLMENT__hr",
                                  "BASELINE__Charleson","BASELINE__SOFA_Total",
                                  "BASELINE__APACHEII","gcs",
                                  "PREENROLLMENT__o2sat","PREENROLLMENT__sbp",
                                  "BASELINE__Intubated",
                                  "PREENROLLMENT__dbp","DCF__infectsite",
                                  "Hosp_LOS","DCF__achf","DCF__cirrhos","DCF__rendis",
                                  "Randomization_to_Hour6_Dobutamine",
                                  "Randomization_to_Hour6_Blood")]
dsc4$arm<-1*(dsc4$arm==3)
dsc4$SCREENING__sex<-1*(dsc4$SCREENING__sex==1)
dsc4$site_lung<-ifelse(dsc4$DCF__infectsite==1,1,
                       ifelse(dsc4$DCF__infectsite==9,NA,0))
dsc4$site_urine<-ifelse(dsc4$DCF__infectsite==3,1,
                        ifelse(dsc4$DCF__infectsite==9,NA,0))
dsc4$site_abdom<-ifelse(dsc4$DCF__infectsite==2,1,
                        ifelse(dsc4$DCF__infectsite==9,NA,0))


colnames(dsc4)[-c(24:26)]<-c("stnum","y","w","age","male","rr",
                             "temp","hr","charlson","sofa_total",
                             "apache","gcs","sat","sbp","mvent","dbp","infectsite","LOS",
                             "chf","cirr","dial","r6_dobut","r6_blood")
dsc4<-dsc4[,-17]
dsc4$cirr<-1*(dsc4$cirr=="1")


ds0v<-merge(dsc4,dscl0,by="stnum",all.x=TRUE,all.y=TRUE)
ds0<-cbind(ds0v[,1:17],ds0v[,23:35],ds0v[,18:22])

# Add prism variables
wv<-which(colnames(ds0) %in% c("age","male","y","Hosp_LOS","LOS","apache","sbp","dbp",
                               "lac","sofa_total","mvent","w"))

ds01<-ds0[,-wv]
colnames(ds01)[1]<-"pt_id_origtrial"

ds.procv<-merge(ds01,pris.xwalk[pris.xwalk$trial_id=="C",-2],by="pt_id_origtrial")
ds.proc<-cbind(ds.procv[,1:19],ds.procv[,25:ncol(ds.procv)],ds.procv[,20:24])
ds.proc2<-cbind(ds.proc[,c(20,23)],ds.proc[,-c(1,20,23,24,32)])
colnames(ds.proc2)[1:2]<-c("y","w")
colnames(ds.proc2)[22:29]<-c("male","apache","sbp","map","lac","sofa_total","mvent","vaso")

xvars_process<-c("age","rr","temp","hr","charlson","albumin",
                 "hgb","pao2","map","sbp","male",
                 "bili","bun","cr","g","lac","plt",
                 "sat","wbc","apache","sofa_total","gcs",
                 "site_lung","site_urine","site_abdom",
                 "mvent","vaso","cirr","chf","dial",
                 "r6_dobut","r6_blood")


ds_c<-ds.proc2[,c("y","w",xvars_process)]






# ARISE

#machine<-"D:/Box Sync/Box Sync/Box Sync/"
data.dir<-"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/ARISE/"
anal.dir<-"_BoxMigration/Faraaz PROWESS/Programs/ITE Score approach/ARISE/"


ds<-read.csv(paste0(machine,data.dir,"IPDMA_Arise_xtra.csv"))

ds$y0<-ds$D_DAY
ds$y0[is.na(ds$y0)]<-9999
ds$y<-1*(ds$y0<60)
ds$w<-1*(ds$GROUP=="E")
ds$age<-ds$AGE
ds$male<-1*(ds$SEX=="M")
ds$temp<-ds$core_temperature
ds$hr<-ds$heart_rate
ds$rr<-ds$respiratory_rate
ds$map<-ds$LAST_MAP           # is this pre-treatment? Don't have a DD entry
ds$sbp<-ds$systolic_blood_pressure
ds$dbp<-ds$diastolic_blood_pressure
ds$sat<-ds$spo2
ds$urine<-ds$urine_output
ds$apache<-ds$APIIPMH_CVS           # what to use for apache?
ds$sofa_total<-ds$SOFA_all5         # there is also LAST_SOFA...which is better?
ds$vaso<-ds$vasopressor_drugs_given # vasopressor drugs given 1 hr prior (assuming this is prior to enrollment)
#ds$ph<-ds$ph
#ds$pao2<-ds$pao2
ds$cr<-ds$creatinine
#ds$sodium
ds$bicarb<-ds$bicarbonate
ds$bun<-ds$urea
ds$k<-ds$potassium
#ds$albumin<-ds$albumin
ds$g<-ds$glucose
ds$lac<-ds$lactate
ds$bili<-ds$bilirubin
ds$plt<-ds$platelets
ds$hgb<-ds$haemoglobin
ds$wbc<-ds$wcc                    # is this white cell count
ds$pt<-ds$pt_ratio               # this is prothrombin; we should have for process too, but not yet
ds$cirr<-1*(ds$cirrhosis_with_portal_hyper==1 | ds$cirrhosis_without_portal_hyper==1)
ds$dial<-ds$dialysis
ds$chf<-ds$dyspnoea_with_therapy
ds$r6_dobut<-1*(ds$dobut==1)
ds$r6_dobut[is.na(ds$r6_dobut)]<-0
ds$r6_blood<-1*(ds$R6_BLOODVOL>0)


# Calculate Charlson

cage<-ifelse(ds$age<50,1,
             ifelse(ds$age<=60,2,
                    ifelse(ds$age<=70,3,
                           ifelse(ds$age<=80,4,5))))

cmi<-1*(ds$myocardial_infarct)
chf<-1*(ds$myocardial_infarct==1 | 
          ds$angina_dyspnoea_low_activity==1 | 
          ds$dyspnoea_with_therapy==1 |
          ds$untreated_aortic_aneurysm==1 |
          ds$pulmonary_hypertension==1)
pvd<-1*(ds$claudication_ischaemia_gangrene |
          ds$peripheral_bypass_surgery)
cereb_acc<-1*(ds$cerebral_vascular_accident)
dementia<-1*ds$cerebral_vascular_accident
copd<-1*(ds$dyspnoea_moderate_activity==1 |
           ds$dyspnoea_low_activity==1 |
           ds$chronic_hypoxia==1 |
           ds$chronic_co2_retention==1)
connect_tiss<-1*(ds$lupus_erythematosis==1 |
                   ds$polymyalgia_rheumatica==1 |
                   ds$rheumatoid_arthritis==1)
ulcer<-1*(ds$gastric_duodenal_ulcer)

liver_mild<-1*(ds$cirrhosis_without_portal_hyper==1)
liver_severe<-1*(ds$cirrhosis_with_portal_hyper==1 | ds$hepatic_failure==1)
liver<-ifelse(liver_severe==1,3,
              ifelse(liver_mild==1,1,0))
db_uncomp<-1*(ds$diabetes_oral==1 | ds$diabetes_insulin==1)
db_dmg<-1*(ds$diabetic_retino_neuro_nephro)
db<-ifelse(db_dmg==1,2,
           ifelse(db_uncomp==1,1,0))
hemiplegia<-1*(ds$cva_hemiplegia_paraplegia)
ckd_sev<-1*(ds$dialysis)
tumor_local<-1*(ds$tumour_without_metastases)
tumor_meta<-1*(ds$tumour_with_metastases)
tumor<-ifelse(tumor_meta==1,6,
              ifelse(tumor_local==1,2,0))
leukemia<-1*(ds$leukaemia)
lymphoma<-1*(ds$lymphoma)
aids<-1*(ds$aids)

ds$charlson<- cmi + chf + pvd + cereb_acc + dementia +
  copd + connect_tiss + ulcer + liver + db +
  2*hemiplegia + 2*ckd_sev + 2*leukemia + 2*lymphoma +
  tumor + 6*aids

ds$apache<-ds$LAST_APII 
ds$mvent0<-ds$ventilator_dependent

ds$mvent2a<-ds$invasive_mechanical_ventilation
ds$mvent2a[is.na(ds$mvent2a)]<-0
ds$mvent2b<-ds$noninvasive_mechanical_vent
ds$mvent2b[is.na(ds$mvent2b)]<-0
ds$mvent<-1*(ds$mvent2a==1 | ds$mvent2b==1)

ds$gcs<-ds$best_eye_response + ds$best_verbal_response + ds$best_motor_response 

ds$site_abdom<-1*(ds$E_INFECTION=="A")
ds$site_blood<-1*(ds$E_INFECTION=="B")
ds$site_cns<-1*(ds$E_INFECTION=="C")
ds$site_unk<-1*(ds$E_INFECTION=="U")
ds$site_lung<-1*(ds$E_INFECTION=="L")
ds$site_other<-1*(ds$E_INFECTION=="O")
ds$site_softis<-1*(ds$E_INFECTION=="S")
ds$site_urine<-1*(ds$E_INFECTION=="U")

ds$cr<-ds$cr*0.0113     #umol/L to mg/dL
ds$bun<-ds$bun*2.8     #mmol/L to mg/dL
ds$albumin<-ds$albumin*0.1     #g/L to g/dL
ds$g<-ds$g*18.0182         #mmol/L to mg/dL
ds$bili<-ds$bili*0.058     #umol/L to mg/dL
ds$hgb<-ds$hgb*0.1      #g/L to g/dL
#ds$wbc<-ds$wbc*1000     #10^9/L to count/mm3



xvars_arise<-c("age","rr","temp","hr","charlson","albumin",
               "hgb","pao2","dbp","sbp","male",
               "bili","bun","cr","g","lac","plt",
               "sat","wbc","apache","sofa_total","gcs",
               "site_lung","site_urine","site_abdom",
               "mvent","cirr","chf","dial",
               "r6_dobut","r6_blood")


ds_a<-ds[,c("patient_id","y","w",xvars_arise)]


# Add prism variables
wv<-which(colnames(ds_a) %in% c("age","male","y","apache","sbp","dbp",
                                "lac","sofa_total","mvent","w"))

ds_a<-ds_a[,-wv]
colnames(ds_a)[1]<-"pt_id_origtrial"

ds.arisev<-merge(ds_a,pris.xwalk[pris.xwalk$trial_id=="A",-2],by="pt_id_origtrial")
ds.arise<-cbind(ds.arisev[,1:19],ds.arisev[,25:ncol(ds.arisev)],ds.arisev[,20:24])
ds.arise2<-cbind(ds.arise[,c(20,23)],ds.arise[,-c(1,20,23,24,32)])
colnames(ds.arise2)[1:2]<-c("y","w")
colnames(ds.arise2)[22:29]<-c("male","apache","sbp","map","lac","sofa_total","mvent","vaso")

xvars_arise<-c("age","rr","temp","hr","charlson","albumin",
               "hgb","pao2","map","sbp","male",
               "bili","bun","cr","g","lac","plt",
               "sat","wbc","apache","sofa_total","gcs",
               "site_lung","site_urine","site_abdom",
               "mvent","vaso","cirr","chf","dial",
               "r6_dobut","r6_blood")


ds_a<-ds.arise2[,c("y","w",xvars_arise)]

ds0<-ds_c
ds0<-as.data.frame(lapply(ds0,as.numeric))
ds0$trial<-"C"
ds01<-ds_a
ds01<-as.data.frame(lapply(ds01,as.numeric))
ds01$trial<-"A"

ds00_cci<-rbind(ds0,ds01)


ds000<-ds00_cci
ds000<-ds000[!is.na(ds000$y),]
rownames(ds000)<-1:nrow(ds000)

ds000<-ds000[ds000$trial %in% c("C","A"),]


# Now merge with rdat to get site_id

rdat<-read.csv(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/arise_process_r6.csv"))
rdat<-rdat[!is.na(rdat$y),]
ds000$site<-rdat$site_id


# Now randomly allocate sites into two "new" trials, "E" and "F"

sites<-unique(ds000$site)
set.seed(2152023)
s1<-sample(sites,30,replace=FALSE)
ds000$trial2<-ifelse(ds000$site %in% s1,"E","F")
rdat$trial2<-ifelse(rdat$site_id %in% s1,"E","F")
ds000$trial_id<-ds000$trial

table(ds000$trial2,ds000$trial)

tapply(rdat$r6_fluidvol[rdat$w==0],rdat$trial2[rdat$w==0],mean,na.rm=TRUE)
tapply(rdat$r6_fluidvol[rdat$w==1],rdat$trial2[rdat$w==1],mean,na.rm=TRUE)

tapply(rdat$r6_vaso[rdat$w==0],rdat$trial2[rdat$w==0],mean,na.rm=TRUE)
tapply(rdat$r6_vaso[rdat$w==1],rdat$trial2[rdat$w==1],mean,na.rm=TRUE)

ds000$trial<-ds000$trial2
ds000<-ds000[,-which(names(ds000) %in% c("trial2"))]

trialord<-order(ds000$trial)
ds000<-ds000[trialord,]

ds000_cci<-ds000


