

machine<-"D:/Box Sync/Box Sync/Box Sync/"
ds000<-readRDS(paste0(machine,"Faraaz PROWESS/datasets/EGDT Trials/arise_process_ds000.RDS"))


# Explore proclivity to treat

prism.dir<-"Faraaz PROWESS/Programs/ITE Score approach/PRISM/Dataset/"
dsp<-read.csv(paste0(machine,prism.dir,"prism_data_scvo2.csv"))
dsp$site2<-paste0(dsp$site_id,dsp$trial_id)
dsp$w<-1
dsp$w[dsp$treatment=="Usual Care"]<-0

# Set up variables in expected model(s)
# mean impute 
dsp$a<-dsp$age
dsp$a[is.na(dsp$a)]<-mean(dsp$a,na.rm=TRUE)
dsp$n<-dsp$nhome
dsp$n[is.na(dsp$n)]<-0
dsp$cc<-1*(dsp$apiipmh_cvs+
             dsp$apiipmh_immuno+
             dsp$apiipmh_liver+
             dsp$apiipmh_renal+
             dsp$apiipmh_resp>0)
dsp$cc[is.na(dsp$cc)]<-0
dsp$s<-1*(dsp$sex=="M")

dsp$site_soft<-dsp$e_infection_s
dsp$site_lung<-dsp$e_infection_l
dsp$site_abdom<-dsp$e_infection_a
dsp$site_blood<-dsp$e_infection_b
dsp$site_cns<-dsp$e_infection_c
dsp$site_urin<-dsp$e_infection_u
dsp$site_oth<-dsp$e_infection_o
dsp$site_unkn<-dsp$e_infection_k
dsp$site_none<-dsp$e_infection_n

dsp$crit_rh<-dsp$e_rh_only
dsp$crit_hlact<-dsp$e_hl_only
dsp$crit_both<-dsp$e_both

dsp$ap<-dsp$last_apsapii

dsp$rmv<-dsp$r_imv
dsp$rvaso<-dsp$r_vaso
dsp$day1<-dsp$r_dow1
dsp$day2<-dsp$r_dow2
dsp$day3<-dsp$r_dow3
dsp$day4<-dsp$r_dow4
dsp$day5<-dsp$r_dow5
dsp$day6<-dsp$r_dow6
dsp$day7<-dsp$r_dow7

dsp$time_of_day<-dsp$r_tod
dsp$pre_fluidvol_liters<-dsp$prer_fluidvol/1000
dsp$r6_fluidvol_liters<-dsp$r6_fluidvol/1000
dsp$r6_bloodvol_liters<-dsp$r6_bloodvol/1000
dsp$prer_r6_fluidvol<-(dsp$prer_fluidvol + dsp$r6_fluidvol)/1000

dsp$wt<-dsp$bodywt
dsp$time_ed2rand<-dsp$time_prer

vaso_dat<-dsp[dsp$w==0,c("r6_vaso",
                         "a","s","n","cc","site_soft",
                         "site_lung","site_abdom",
                         "site_blood","site_cns","site_urin",
                         "site_oth","site_unkn","site_none",
                         "crit_hlact","crit_both",
                         "ap","rmv","rvaso",
                         "day1","day2","day3","day4",
                         "day5","day6","day7",
                         "time_of_day",
                         "pre_fluidvol_liters",
                         "r6_fluidvol_liters",
                         "r6_bloodvol_liters")]

vaso_dat<-na.omit(vaso_dat)

vaso_expect_fit<-glm(r6_vaso~a+s+n+cc+
                       #site_soft+
                       site_lung+site_abdom+
                       site_blood+site_cns+site_urin+
                       site_oth+site_unkn+site_none+
                       #crit_rh+
                       crit_hlact+crit_both+
                       ap+rmv+rvaso+
                       day1+
                       day2+day3+day4+
                       day5+day6+
                       #day7+
                       time_of_day+
                       pre_fluidvol_liters+r6_fluidvol_liters+
                       r6_bloodvol_liters,
                     data=vaso_dat,
                     family=binomial(link="logit"))

fit_sum<-as.data.frame(summary(vaso_expect_fit)$coefficients)
fit_sum$OR<-round(exp(fit_sum[,1]),2)
fit_sum$OR_nejm<-c(NA,1.01,1.03,0.95,1.00,
                   1.16,1.35,1.5,2.59,1.13,1.10,
                   0.82,0.61,0.65,1.23,1.07,1.37,
                   1.24,1.09,1.22,1.08,1.23,1.14,
                   0.7,1.39,1.3,1.14,1.22)
fit_sum$OR_diff<-fit_sum$OR-fit_sum$OR_nejm

#write.csv(fit_sum,paste0(machine,prism.dir,"vaso_fit_sum.csv"))


fluid_dat<-dsp[dsp$w==0,c("r6_fluidvol_liters",
                          "a","s","wt","n","cc",
                          #"site_soft",
                          "site_lung","site_abdom",
                          "site_blood","site_cns","site_urin",
                          "site_oth","site_unkn","site_none",
                          #"crit_rh",
                          "crit_hlact","crit_both",
                          "ap","rmv","rvaso",
                          "time_ed2rand",
                          "day1",
                          "day2","day3","day4",
                          "day5","day6",
                          #"day7",
                          "time_of_day",
                          "pre_fluidvol_liters",
                          "r6_bloodvol_liters",
                          "r6_vaso")]
fluid_dat<-na.omit(fluid_dat)

fluidvol_expect_fit<-lm(r6_fluidvol_liters~.,
                        data=fluid_dat)

fit_sum<-as.data.frame(summary(fluidvol_expect_fit)$coefficients)
fit_sum$BETA_nejm<-c(NA,-0.008,0.053,0.001,0.121,-0.007,-0.173,0.161,
                     -0.108,-0.160,0.042,0.253,-0.061,0.353,
                     -0.161,0.322,0.023,-0.120,-0.207,-0.125,
                     0.188,0.296,0.259,0.244,0.275,0.192,0.191,
                     -0.034,0.374,0.250)
fit_sum[,1]<-round(fit_sum[,1],3)
fit_sum$BETA_diff<-fit_sum$Estimate-fit_sum$BETA_nejm
fit_sum
# 
# site1A<-dsp[dsp$w==0 & dsp$site2=="1A",]
# #nrow(site1A) # 71 usual care patients
# summary(site1A[,c("fluid_actual","vaso_actual")])
# summary(site1A[,c("r6_fluidvol","r6_vaso")])


#write.csv(fit_sum,paste0(machine,prism.dir,"fluids_fit_sum.csv"))

dsp$pred_r6_vaso<-predict(vaso_expect_fit,
                          newdata=dsp,
                          type="response")
vaso_pred<-tapply(dsp[dsp$w==0,"pred_r6_vaso"],dsp[dsp$w==0,"site2"],mean,na.rm=TRUE)

dsp$pred_r6_fluidvol<-predict(fluidvol_expect_fit,
                              newdata=dsp)
fluidvol_pred<-tapply(dsp[dsp$w==0,"pred_r6_fluidvol"],dsp[dsp$w==0,"site2"],mean,na.rm=TRUE)



# Select sites with 3 or more UC patients

uc_sites<-tapply(dsp[dsp$w==0,"vaso_actual"],dsp[dsp$w==0,"site2"],length)
remove<-names(which(uc_sites<4))

dsp0<-dsp[! dsp$site2 %in% remove,]
dsp0<-dsp0[dsp0$w==0,]

vaso_obssite<-tapply(dsp0$r6_vaso,dsp0$site2,mean,na.rm=TRUE)
vaso_predsite<-tapply(dsp0$pred_r6_vaso,dsp0$site2,mean,na.rm=TRUE)
vaso_myratio<-vaso_obssite/vaso_predsite
vaso_actual<-tapply(dsp0$vaso_actual,dsp0$site2,mean,na.rm=TRUE)
vaso_expected<-tapply(dsp0$vaso_expected,dsp0$site2,mean,na.rm=TRUE)
vaso_davidratio<-tapply(dsp0$vaso_ratio,dsp0$site2,mean,na.rm=TRUE)
# vaso_myratio is nearly perfectly correlated with vaso_davidratio
cor(vaso_myratio,vaso_davidratio,use="complete.obs")


fluidvol_obssite<-tapply(dsp0$r6_fluidvol/1000,dsp0$site2,mean,na.rm=TRUE)
fluidvol_obssite2<-tapply((dsp0$r6_fluidvol+dsp0$prer_fluidvol)/1000,dsp0$site2,mean,na.rm=TRUE)
fluidvol_predsite<-tapply(dsp0$pred_r6_fluidvol,dsp0$site2,mean,na.rm=TRUE)
fluidvol_myratio<-fluidvol_obssite/fluidvol_predsite
fluidvol_myratio2<-fluidvol_obssite2/fluidvol_predsite
fluidvol_actual<-tapply(dsp0$fluid_actual,dsp0$site2,mean,na.rm=TRUE)
fluidvol_expected<-tapply(dsp0$fluid_expected,dsp0$site2,mean,na.rm=TRUE)
fluidvol_davidratio<-tapply(dsp0$fluid_ratio,dsp0$site2,mean,na.rm=TRUE)
cor(fluidvol_myratio,fluidvol_davidratio,use="complete.obs")


dev.new(height=6,width=12,noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(vaso_myratio,vaso_davidratio,main="Vasopressor")
plot(fluidvol_myratio,fluidvol_davidratio,main="Fluid Vol")


sites_mat<-data.frame(site2=names(fluidvol_obssite),
                      num_pat_uc=uc_sites[-which(names(uc_sites) %in% remove)],
                      fluidvol_myratio,
                      fluidvol_davidratio,
                      vaso_myratio,
                      vaso_davidratio,
                      fluidvol_obssite,
                      vaso_obssite)
sites_mat$fluid_diff<-1*((sites_mat$fluidvol_myratio<1) & (sites_mat$fluidvol_davidratio>1)) | 
  ((sites_mat$fluidvol_myratio>1) & (sites_mat$fluidvol_davidratio<1))
sites_mat$vaso_diff<-1*((sites_mat$vaso_myratio<1) & (sites_mat$vaso_davidratio>1)) | 
  ((sites_mat$vaso_myratio>1) & (sites_mat$vaso_davidratio<1))
sites_mat$trial<-gsub('[[:digit:]]+', '', sites_mat$site2)



# Connect model data with my O/E ratios and other useful info in sites_mat

rdat<-read.csv(paste0(machine,"Faraaz PROWESS/datasets/EGDT Trials/arise_process_r6.csv"))
rdat<-rdat[!is.na(rdat$y),]
rdat<-rdat[order(rdat$trial),]
rdat$site2<-paste0(rdat$site_id,rdat$trial)

smac<-sites_mat[sites_mat$trial %in% c("A","C"),]

rdat<-merge(rdat,smac[,-which(names(smac)=="trial")],by="site2",all.x=TRUE,all.y=FALSE)

num_pat_tot<-data.frame(num_pat_tot=tapply(rdat$X,rdat$site2,length),
                        site2=names(tapply(rdat$X,rdat$site2,length)),
                        num_pat_uc=tapply(rdat$num_pat_uc,rdat$site2,mean))
rdat<-merge(rdat,num_pat_tot[,-3],by="site2",all=TRUE)


# Plot quantities for 

sites<-data.frame(site2=names(tapply(rdat$X,rdat$site2,length)),
                  num_pat_tot=tapply(rdat$X,rdat$site2,length),
                  num_pat_uc=tapply(rdat$num_pat_uc,rdat$site2,mean),
                  fluidvol_myratio=tapply(rdat$fluidvol_myratio,rdat$site2,mean),
                  vaso_myratio=tapply(rdat$vaso_myratio,rdat$site2,mean),
                  fluidvol_obssite=tapply(rdat$fluidvol_obssite,rdat$site2,mean),
                  vaso_obssite=tapply(rdat$vaso_obssite,rdat$site2,mean))

sites$site_size_cat<-ifelse(sites$num_pat_uc<11,1,
                            ifelse(sites$num_pat_uc<21,2,3))


# Generate O/E groups

sites0<-sites[!is.na(sites$num_pat_uc),]

hi<-1.1
lo<-0.9

sites0$less<-1*(sites0$fluidvol_myratio<lo & sites0$vaso_myratio<lo)
sites0$dry<-1*(sites0$fluidvol_myratio<lo & sites0$vaso_myratio>hi)
sites0$wet<-1*(sites0$fluidvol_myratio>hi & sites0$vaso_myratio<lo)
sites0$max<-1*(sites0$fluidvol_myratio>hi & sites0$vaso_myratio>hi)
sites0$avg<-1*(sites0$fluidvol_myratio<hi & sites0$fluidvol_myratio>lo & 
                 sites0$vaso_myratio<hi & sites0$vaso_myratio>lo)

sites0$none<-apply(sites0[,c("less","dry","max","wet","avg")],1,function(x) sum(x)==0)

sum(sites0$num_pat_tot[sites0$less==1])
sum(sites0$num_pat_tot[sites0$dry==1])
sum(sites0$num_pat_tot[sites0$wet==1])
sum(sites0$num_pat_tot[sites0$max==1])
sum(sites0$num_pat_tot[sites0$avg==1])
sum(sites0$num_pat_tot[sites0$none==1])



dev.new(height=6,width=10,noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(sites$num_pat_uc,
     sites$fluidvol_myratio,
     xlab="Number patients UC arm",
     ylab="Fluid volume O/E ratio")
plot(sites$num_pat_uc,
     sites$vaso_myratio,
     xlab="Number patients UC arm",
     ylab="Vasopressor O/E ratio")

dev.new()
plot(sites$fluidvol_myratio,
     sites$vaso_myratio,
     ylab="Vasopressor O/E ratio",
     xlab="Fluid volume O/E ratio")
abline(h=1,lty=2)
abline(v=1,lty=2)
symbols(x=sites$fluidvol_myratio,
        y=sites$vaso_myratio, 
        circles=sites$site_size_cat, 
        inches=1/10,
        add=TRUE, 
        bg="steelblue2", 
        fg=NULL)


dev.new(height=6,width=10,noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(sites$num_pat_uc,
     sites$fluidvol_myratio,
     xlab="Number patients UC arm",
     ylab="Fluid volume O/E ratio")
plot(sites$num_pat_uc,
     sites$vaso_myratio,
     xlab="Number patients UC arm",
     ylab="Vasopressor O/E ratio")

dev.new()
plot(sites$fluidvol_myratio,
     sites$vaso_myratio,
     ylab="Vasopressor O/E ratio",
     xlab="Fluid volume O/E ratio")
abline(h=1.1,lty=2,col=2)
abline(v=1.1,lty=2,col=2)
abline(h=0.9,lty=2,col=2)
abline(v=0.9,lty=2,col=2)
symbols(x=sites$fluidvol_myratio,
        y=sites$vaso_myratio, 
        circles=sites$site_size_cat, 
        inches=1/10,
        add=TRUE, 
        bg="steelblue2", 
        fg=NULL)


dev.new()
plot(sites$fluidvol_myratio,
     sites$vaso_myratio,
     ylab="Vasopressor O/E ratio",
     xlab="Fluid volume O/E ratio")
segments(x0=-99,x1=0.9,y0=1.1,y1=1.1,col=2,lty=2,lwd=3)
segments(x1=99,x0=1.1,y0=0.9,y1=0.9,col=2,lty=2,lwd=3)
segments(x0=0.9,x1=0.9,y0=99,y1=1.1,col=2,lty=2,lwd=3)
segments(x0=1.1,x1=1.1,y1=-99,y0=0.9,col=2,lty=2,lwd=3)
# 
# abline(h=1.1,lty=2,col=2)
# abline(v=1.1,lty=2,col=2)
# abline(h=0.9,lty=2,col=2)
# abline(v=0.9,lty=2,col=2)
symbols(x=sites$fluidvol_myratio,
        y=sites$vaso_myratio, 
        circles=sites$site_size_cat, 
        inches=1/10,
        add=TRUE, 
        bg="steelblue2", 
        fg=NULL)







# weighted correlation

library(wCorr)
weightedCorr(sites$fluidvol_myratio[!is.na(sites$fluidvol_myratio) & !is.na(sites$vaso_myratio)],
             sites$vaso_myratio[!is.na(sites$fluidvol_myratio) & !is.na(sites$vaso_myratio)],
             weights=sites$num_pat_uc[!is.na(sites$fluidvol_myratio) & !is.na(sites$vaso_myratio)],
             method="Spearman")


# Relationship between Observed and O/E

dev.new(height=6,width=10,noRStudioGD = TRUE)
par(mfrow=c(1,2))
plot(sites$fluidvol_obssite,
     sites$fluidvol_myratio,
     xlab="Observed fluid volume",
     ylab="Fluid volume O/E ratio")
plot(sites$vaso_obssite,
     sites$vaso_myratio,
     xlab="Observed % vasopressors",
     ylab="Vasopressor O/E ratio")




# Histograms of patients at each site
dev.new()
hist(sites$num_pat_uc,
     xlab="Number of UC patients",
     main="",
     breaks=25)


# (May need to re-do all below depending on feedback)


# Now look at CATE predictions for the patients in the red and green boxes

# Need to read in files with CATES on them first

# These files are in teh order of the patietns in the ORIGINAL ds000 
# llrf.pred.e2f<-readRDS(paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_e2f_v2.RDS"))
# llrf.pred.f2e<-readRDS(paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_f2e_v2.RDS"))
# 
# pred.e2f<-apply(llrf.pred.e2f,1,mean)
# pred.f2e<-apply(llrf.pred.f2e,1,mean)
# 
# preds.EthenF<-c(pred.f2e,pred.e2f)

llrf.pred.ef<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/llrf_pred_eANDf_MIglue_min20_MI100.RDS"))
pred.ef<-apply(llrf.pred.ef,1,mean)
ds000$iarr<-(-pred.ef)

ord000<-order(ds000$trial_id,ds000$site,ds000$map,ds000$age,ds000$g,ds000$temp,ds000$sbp)

ds000<-ds000[order(ds000$trial_id,ds000$site,ds000$map,ds000$age,ds000$g,ds000$temp,ds000$sbp),]
rdat<-rdat[order(rdat$trial,rdat$site_id,rdat$map,rdat$age,rdat$g,rdat$temp,rdat$sbp),]
rdat<-cbind(rdat,ds000[,c("iarr","trial")])
colnames(rdat)[ncol(rdat)]<-"trial"

#saveRDS(ord000,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/ord000.RDS")
#saveRDS(ds000,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/ds000_reord.RDS")

lo<-0.9
hi<-1.1

rdat$ratio_grp<-ifelse(is.na(rdat$fluidvol_myratio) | is.na(rdat$vaso_myratio),NA,
                       ifelse(rdat$fluidvol_myratio<lo & rdat$vaso_myratio<lo,"less",
                              ifelse(rdat$fluidvol_myratio<lo & rdat$vaso_myratio>hi,"dry",
                                     ifelse(rdat$fluidvol_myratio>hi & rdat$vaso_myratio<lo,"wet",
                                            ifelse(rdat$fluidvol_myratio>hi & rdat$vaso_myratio>hi,"max",
                                                   ifelse(rdat$fluidvol_myratio<hi & rdat$fluidvol_myratio>lo & 
                                                            rdat$vaso_myratio<hi & rdat$vaso_myratio>lo,"avg","none"))))))


# Patient-level fluids and pressors by ratio group

dev.new(height=6,width=14,noRStudioGD = TRUE)
par(mfrow=c(1,3))
boxplot(r6_fluidvol/1000~ratio_grp,
        data=rdat[rdat$w==0,],
        ylab="Fluid volume (L)",
        xlab="Ratio group",
        main="Usual Care",
        ylim=c(0,10))
boxplot(r6_fluidvol/1000~ratio_grp,
        data=rdat[rdat$w==1,],
        ylab="Fluid volume (L)",
        xlab="Ratio group",
        main="EGDT",
        ylim=c(0,10))
del<-sapply(unique(rdat$ratio_grp),function(x){
  dat<-rdat[rdat$ratio_grp==x,]
  yy<-dat$r6_fluidvol/1000
  ww<-dat$w
  mean(yy[ww==1],na.rm=TRUE)-mean(yy[ww==0],na.rm=TRUE)
})
bardat<-data.frame(delta=na.omit(del),
                   grp=na.omit(unique(rdat$ratio_grp)))
barplot(delta~grp,
        bardat,
        ylab="Fluid vol difference (L)",
        xlab="Ratio group",
        main="Difference (EGDT-UC)",
        ylim=c(0,1))


# Patient-level fluids and pressors by ratio group (3 groups)


rdat$ratio_grp3<-ifelse(is.na(rdat$fluidvol_myratio) | is.na(rdat$vaso_myratio),NA,
                       ifelse(rdat$fluidvol_myratio<lo & rdat$vaso_myratio>hi,"Restrictive",
                              ifelse(rdat$fluidvol_myratio>hi & rdat$vaso_myratio<lo,"Liberal","Other")))
rdat$ratio_grp3<-factor(rdat$ratio_grp3,levels=c("Restrictive","Liberal","Other"))
#levels(rdat$ratio_grp3)<-c("Restrictive","Liberal","Other")

#saveRDS(rdat,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/rdat_with_grp3.RDS")


dev.new(height=6,width=14,noRStudioGD = TRUE)
par(mfrow=c(1,3))
boxplot(r6_fluidvol/1000~ratio_grp3,
        data=rdat[rdat$w==0,],
        ylab="Fluid volume (L)",
        xlab="Site group",
        main="Usual Care",
        ylim=c(0,10))
boxplot(r6_fluidvol/1000~ratio_grp3,
        data=rdat[rdat$w==1,],
        ylab="Fluid volume (L)",
        xlab="Site group",
        main="EGDT",
        ylim=c(0,10))
del<-sapply(unique(rdat$ratio_grp3),function(x){
  dat<-rdat[rdat$ratio_grp3==x,]
  yy<-dat$r6_fluidvol/1000
  ww<-dat$w
  mean(yy[ww==1],na.rm=TRUE)-mean(yy[ww==0],na.rm=TRUE)
})
bardat<-data.frame(delta=na.omit(del),
                   grp=na.omit(unique(rdat$ratio_grp3)))
barplot(delta~grp,
        bardat,
        ylab="Fluid vol difference (L)",
        xlab="Site group",
        main="Difference (EGDT-UC)",
        ylim=c(0,1))






dev.new(height=6,width=14,noRStudioGD = TRUE)
par(mfrow=c(1,3))
uc<-sapply(unique(rdat$ratio_grp3),function(x){
  dat<-rdat[rdat$ratio_grp3==x,]
  yy<-dat$r6_vaso
  ww<-dat$w
  mean(yy[ww==0],na.rm=TRUE)
})
bardat<-data.frame(delta=na.omit(uc),
                   grp=na.omit(unique(rdat$ratio_grp3)))
barplot(delta~grp,
        bardat,
        ylab="Vasopressor use (%)",
        xlab="Ratio group",
        main="Usual Care",
        ylim=c(0,0.6))
egdt<-sapply(unique(rdat$ratio_grp3),function(x){
  dat<-rdat[rdat$ratio_grp3==x,]
  yy<-dat$r6_vaso
  ww<-dat$w
  mean(yy[ww==1],na.rm=TRUE)
})
bardat<-data.frame(delta=na.omit(egdt),
                   grp=na.omit(unique(rdat$ratio_grp3)))
barplot(delta~grp,
        bardat,
        ylab="Vasopressor use (%)",
        xlab="Ratio group",
        main="EGDT",
        ylim=c(0,0.6))
del<-sapply(unique(rdat$ratio_grp3),function(x){
  dat<-rdat[rdat$ratio_grp3==x,]
  yy<-dat$r6_vaso
  ww<-dat$w
  mean(yy[ww==1],na.rm=TRUE)-mean(yy[ww==0],na.rm=TRUE)
})
bardat<-data.frame(delta=na.omit(del),
                   grp=na.omit(unique(rdat$ratio_grp3)))
barplot(delta~grp,
        bardat,
        ylab="Vasopressor use (%)",
        xlab="Ratio group",
        main="% Difference (EGDT-UC)",
        ylim=c(0,0.3))




# Use gates function to get tertiles by group

gates.dat<-rdat[,c("y","w","iarr","ratio_grp","ratio_grp3")]
names(gates.dat)[3]<-"score"
gates.dat$score<-(-gates.dat$score)

b<-regression_forest(X=as.matrix(rdat[,-c(1:4,32:54)]),
                     Y=rdat[,"y"])$predictions
gates.dat$b<-as.numeric(b)

# Check overall
gates<-gates.fun(list(gates.dat),
                 pi=0.5,
                 num.q=5,
                 lintest=TRUE,
                 loo=TRUE)

# Make function to tabulate for me

gates.dat$ratio_grp[is.na(gates.dat$ratio_grp)]<-"NA"

gates.dat$ratio_grp3<-as.character(gates.dat$ratio_grp3)
gates.dat$ratio_grp3[is.na(gates.dat$ratio_grp3)]<-"NA"

gates.tab<-function(sub,num.q=3){
  gates<-gates.fun(list(gates.dat[sub,]),
                   pi=0.5,
                   num.q=num.q,
                   lintest=TRUE,
                   loo=TRUE)
  
  nn<-rev(gates$`N quantile`)
  pred<-rev(round(-gates$`Mean prediction in quantiles`*100,1))
  obs<-paste0(-rev(round(gates$`Theta Median`*100,1))," (",
              -rev(round(gates$`Theta Upper Median`*100,1)),",",
              -rev(round(gates$`Theta Lower Median`*100,1)),")")
  matrix(c(nn,pred,obs),num.q,3)
}

all<-gates.tab(sub=rep(TRUE,nrow(gates.dat)),3)
dry<-gates.tab(sub=(gates.dat$ratio_grp=="dry"),3)
wet<-gates.tab(sub=(gates.dat$ratio_grp=="wet"),3)
max<-gates.tab(sub=(gates.dat$ratio_grp=="max"),3)
less<-gates.tab(sub=(gates.dat$ratio_grp=="less"),3)
avg<-gates.tab(sub=(gates.dat$ratio_grp=="avg"),3)
write.csv(all,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/all_oe.csv")
write.csv(dry,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/dry_oe.csv")
write.csv(wet,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/wet_oe.csv")
write.csv(max,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/max_oe.csv")
write.csv(less,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/less_oe.csv")
write.csv(avg,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/avg_oe.csv")

lib<-gates.tab(sub=(gates.dat$ratio_grp3=="Liberal"),3)
rest<-gates.tab(sub=(gates.dat$ratio_grp3=="Restrictive"),3)
oth<-gates.tab(sub=(gates.dat$ratio_grp3=="Other"),3)
write.csv(lib,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/liberal_oe.csv")
write.csv(rest,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/restrictive_oe.csv")
write.csv(oth,"D:/Box Sync/Box Sync/Box Sync/Faraaz PROWESS/HTE Manuscript/Results 2023/other_oe.csv")


# Get p-value comparing HTE between two groups

newt<-gates.dat[gates.dat$ratio_grp3 %in% c("Restrictive","Liberal"),]
score <- newt$score
compare<- 1*(newt$ratio_grp3=="Liberal")
mean.score <- mean(score)
mean.compare <- mean(compare)
DF <- data.frame(y = newt$y-newt$b,  
                 main.score.term = (newt$w - 0.5) * mean.score, 
                 interact.score.term = (newt$w - 0.5) * (score - mean.score),
                 compare.score.term = (newt$w - 0.5) * (score - mean.score) * compare)

fit.blp <- lm(y ~ main.score.term + interact.score.term + compare.score.term + 0, 
              data = DF)
#browser()
#fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
#fit.blp<-lm(y~w+score+w*score,data=t)
blp.summary <- lmtest::coeftest(fit.blp, 
                                vcov = sandwich::vcovCL, 
                                type = "HC3")

# Get p-value comparing HTE between two groups

newt<-gates.dat[gates.dat$ratio_grp3 %in% c("Restrictive"),]
score <- newt$score
compare<- 1*(newt$ratio_grp3=="Liberal")
mean.score <- mean(score)
mean.compare <- mean(compare)
DF <- data.frame(y = newt$y-newt$b,  
                 main.score.term = (newt$w - 0.5) * mean.score, 
                 interact.score.term = (newt$w - 0.5) * (score - mean.score))

fit.blp <- lm(y ~ main.score.term + interact.score.term + 0, 
              data = DF)
#browser()
#fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
#fit.blp<-lm(y~w+score+w*score,data=t)
blp.summary <- lmtest::coeftest(fit.blp, 
                                vcov = sandwich::vcovCL, 
                                type = "HC3")

newt<-gates.dat[gates.dat$ratio_grp3 %in% c("Liberal"),]
score <- newt$score
compare<- 1*(newt$ratio_grp3=="Liberal")
mean.score <- mean(score)
mean.compare <- mean(compare)
DF <- data.frame(y = newt$y-newt$b,  
                 main.score.term = (newt$w - 0.5) * mean.score, 
                 interact.score.term = (newt$w - 0.5) * (score - mean.score))

fit.blp <- lm(y ~ main.score.term + interact.score.term + 0, 
              data = DF)
#browser()
#fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
#fit.blp<-lm(y~w+score+w*score,data=t)
blp.summary <- lmtest::coeftest(fit.blp, 
                                vcov = sandwich::vcovCL, 
                                type = "HC3")

newt<-gates.dat[gates.dat$ratio_grp3 %in% c("Other"),]
score <- newt$score
compare<- 1*(newt$ratio_grp3=="Liberal")
mean.score <- mean(score)
mean.compare <- mean(compare)
DF <- data.frame(y = newt$y-newt$b,  
                 main.score.term = (newt$w - 0.5) * mean.score, 
                 interact.score.term = (newt$w - 0.5) * (score - mean.score))

fit.blp <- lm(y ~ main.score.term + interact.score.term + 0, 
              data = DF)
#browser()
#fit.blp<-glm(y~w+score+w*score,data=t,family="binomial")
#fit.blp<-lm(y~w+score+w*score,data=t)
blp.summary <- lmtest::coeftest(fit.blp, 
                                vcov = sandwich::vcovCL, 
                                type = "HC3")





# Get BLP p-values
gates.fun(list(gates.dat[(gates.dat$ratio_grp3=="Liberal"),]),
          pi=0.5,
          num.q=3,
          lintest=TRUE,
          loo=TRUE)$`Pval BLP test`