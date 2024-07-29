
##############################################
## Sensitivity analysis: Splitting by trial ##
##############################################
# Step 0: load dataset (created via create_ds000.R)

#source("C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")
#source("C:/Users/VIT13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")
source("C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")

# Step 1: Multiple imputation

#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
source(paste0(machine,src.dir,"multiple_imputation.R"))

xvars<-c("age","temp","rr","hr","charlson", 
         "albumin","g","lac","plt",
         "hgb","pao2","map","sbp","male",
         "bili","bun","cr",
         "sat","wbc","apache","sofa_total","gcs",
         "site_lung","site_urine","site_abdom",
         "mvent","vaso")

midat_list_a2c<-mi_fun(dat=ds000,
                         split_var_name="trial_id",
                         trset="A",
                         teset="C",
                         mi.use=100,
                         xvars=xvars)

midat_list_c2a<-mi_fun(dat=ds000,
                       split_var_name="trial_id",
                       trset="C",
                       teset="A",
                       mi.use=100,
                       xvars=xvars)

saveRDS(midat_list_c2a,paste0(machine,src.dir,"midat_list_c2a.RDS"))
saveRDS(midat_list_a2c,paste0(machine,src.dir,"midat_list_a2c.RDS"))

# midat_list_a2c<-readRDS(paste0(machine,src.dir,"midat_list_a2c.RDS"))

# Step 2: model training comparison

# Compare model performance in the training set(s)

# K-fold CV for R-loss, but what about RATE?

# machine<-"D:/Box Sync/Box Sync/Box Sync/"
#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
machine<-"C:/Users/vit13/OneDrive - University of Pittsburgh/"
# machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/_BoxMigration/"
sav.dir<-"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/"
hte.dir<-"_BoxMigration/HTE analysis programs/"
source(paste0(machine,hte.dir,"scoring.methods.R"))
source(paste0(machine,hte.dir,"cross.validation.R"))
source(paste0(machine,hte.dir,"post.process_v2.R"))


###########
# ProCESS #

# For RATE, use 20 random train/test splits with 5 imputations each time

midat_list0<-readRDS(paste0(machine,src.dir,"midat_list_c2a.RDS"))
ds<-ds000[ds000$trial_id %in% c("A","C"),c("y","w","trial_id")]
rownames(ds)<-1:nrow(ds)

tr.ids<-which(ds$trial_id=="C")
setA_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cv<-oest.cv.fast(setA_list[1:100],
                 cv.tr.pr=0.5,
                 cv_k=100,
                 y="y",
                 w="w",
                 x=names(setA_list[[1]])[-c(1:2)],
                 pi=0.5,
                 mods=c("rf.risk","cf.CATE","llrf.CATE.rlearn"),
                 cores=1,
                 seed=15223)

#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
#saveRDS(cv,paste0(machine,src.dir,"cv_c2a.RDS"))
cv<-readRDS(paste0(machine,src.dir,"cv_c2a.RDS"))

rf.risk.TOC<-post.cv.fun(cv$rf.risk,
                         stats=c("stats.bin"),
                         scale=20,
                         cores=1,
                         direction="both",
                         score="risk")
cf.risk.TOC<-post.cv.fun(cv$cf.CATE,
                         stats=c("stats.bin"),
                         scale=20,
                         cores=1,
                         direction="both",
                         score="cate")
llrf.risk.TOC<-post.cv.fun(cv$llrf.CATE.rlearn,
                           stats=c("stats.bin"),
                           scale=20,
                           cores=1,
                           direction="both",
                           score="cate")

rfA_rates<-colMeans(rf.risk.TOC$ratelist[[2]])
rfB_rates<-colMeans(rf.risk.TOC$ratelist[[1]])
cf_rates<-colMeans(cf.risk.TOC$ratelist[[1]])
llrf_rates<-colMeans(llrf.risk.TOC$ratelist[[1]])
resmat<-matrix(c(rfA_rates,rfB_rates,cf_rates,llrf_rates),4,4)
resmat[c(1,2,4),]<- round(-resmat[c(1,2,4),]*100,3)
rownames(resmat)<-names(rfA_rates)
colnames(resmat)<-c("rfa","rfb","cf","llrf")


# Draw TOC curves

mean.ate<-mean(sapply(cv[[1]],function(x){
  dat<-x[[1]]
  yy<-dat$Y
  ww<-dat$W
  mean(yy[ww==1])-mean(yy[ww==0])
}))


dev.new()
qq<-seq(0,1,0.05)[-1]
qq<-qq[-length(qq)]
qq<-c(qq,0.995,1)
plot(x=qq,
     y=-c(rev(rf.risk.TOC$outmat$ad_DropLowFirst),mean.ate)*100,
     type="l",
     ylim=c(-0.05,0.15)*100,
     lwd=2,
     bty="l",
     ylab="Absolute Risk Difference (%)",
     xlab="Proportion treated",
     cex.axis=1.3,
     cex.lab=1.3)
points(x=qq,
       y=-c(rev(rf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=1,
       lty=2,
       lwd=2)
points(x=qq,
       y=-c(rev(cf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=2,
       lwd=2)
points(x=qq,
       y=-c(rev(llrf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=3,
       lwd=2)
abline(h=-mean.ate,lty=2)

legend("topright",col=c(1,1,2,3),lty=c(1,2,1,1),lwd=3,legend=c("RF Baseline Risk A","RF Baseline Risk B",
                                         "Causal Forest",
                                         "Local Linear RF R-learner"),
       bty="n")


# Generate LOO R-loss

midat_list0<-readRDS(paste0(machine,src.dir,"midat_list_c2a.RDS"))
ds<-ds000[ds000$trial_id %in% c("A","C"),c("y","w","trial_id")]
rownames(ds)<-1:nrow(ds)

# Pull out the patients in ProCESS
tr.ids<-which(ds$trial=="C")
setA_list<-lapply(1:20,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

# Use all of Trial C to generate LOO predictions from m(x) function
mx.dst<-matrix(NA,length(tr.ids),20)
for (mm in 1:20){
  dst<-setA_list[[mm]]
  mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
                          Y=dst[,"y"])
  mx.dst[,mm]<-predict(mfit)$predictions
}


# Do 5-fold cross-validated R-loss

nn<-nrow(dst)
ii<-1:nn
kmat<-matrix(,2,5)
kmat[,1]<-c(1,floor(nn/5))
kmat[,2]<-c(floor(nn/5)+1,floor(2*nn/5))
kmat[,3]<-c(floor(2*nn/5)+1,floor(3*nn/5))
kmat[,4]<-c(floor(3*nn/5)+1,floor(4*nn/5))
kmat[,5]<-c(floor(4*nn/5)+1,nn)

cf.cvpreds<-llrf.cvpreds<-matrix(,nn,5)

for(kk in 1:5){
  
  kv<-kmat[,kk]
  cvte<-kv[1]:kv[2]
  cvtr<-c(1:nn)[-cvte]
  
  for(mm in 1:5){
    
    dst<-setA_list[[mm]]
    
    cftr<-causal_forest(X=as.matrix(dst[cvtr,-c(1:2)]),
                        Y=dst[cvtr,"y"],
                        W=dst[cvtr,"w"],
                        W.hat=0.5,
                        tune.parameters="none",
                        num.trees=2500,
                        honesty=TRUE,   
                        honesty.fraction=0.5)
    
    mxtr<-mx.dst[cvtr,mm]
    yst<-(dst[cvtr,"y"]-mxtr)/(dst[cvtr,"w"]-0.5)
    llrftr<-ll_regression_forest(X=as.matrix(dst[cvtr,-c(1:2)]),
                                 Y=yst,
                                 enable.ll.split=TRUE,
                                 ll.split.weight.penalty=TRUE,
                                 ll.split.lambda=0.1,
                                 num.trees=2500)
    
    cf.cvpreds[cvte,mm]<-predict(cftr, 
                                 newdata = dst[cvte,-c(1:2)])$predictions
    llrf.cvpreds[cvte,mm]<-predict(llrftr, 
                                   newdata = dst[cvte,-c(1:2)],
                                   linear.correction.variables = c(1:21))$predictions
    
  }
  
}


# Calculate R-loss
mx<-apply(mx.dst,1,mean)
dst1<-setA_list[[1]]
yst<-(dst1[,"y"]-mx)/(dst1[,"w"]-0.5)
cf.final<-apply(cf.cvpreds,1,mean)
llrf.final<-apply(llrf.cvpreds,1,mean)

rloss.cfe<-mean((yst-cf.final)^2)
rloss.llrfe<-mean((yst-llrf.final)^2)





###########
## ARISE ##

midat_list0<-readRDS(paste0(machine,src.dir,"midat_list_a2c.RDS"))
ds<-ds000[ds000$trial_id %in% c("A","C"),c("y","w","trial_id")]
rownames(ds)<-1:nrow(ds)

tr.ids<-which(ds$trial_id=="A")
setA_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cv<-oest.cv.fast(setA_list[1:100],
                 cv.tr.pr=0.5,
                 cv_k=100,
                 y="y",
                 w="w",
                 x=names(setA_list[[1]])[-c(1:2)],
                 pi=0.5,
                 mods=c("rf.risk","cf.CATE","llrf.CATE.rlearn"),
                 cores=1,
                 seed=15223)

#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
#machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
machine<-"C:/Users/vit13/OneDrive - University of Pittsburgh/"
#saveRDS(cv,paste0(machine,src.dir,"cv_a2c.RDS"))
cv<-readRDS(paste0(machine,src.dir,"cv_a2c.RDS"))

rf.risk.TOC<-post.cv.fun(cv$rf.risk,
                         stats=c("stats.bin"),
                         scale=20,
                         cores=1,
                         direction="both",
                         score="risk")
cf.risk.TOC<-post.cv.fun(cv$cf.CATE,
                         stats=c("stats.bin"),
                         scale=20,
                         cores=1,
                         direction="both",
                         score="cate")
llrf.risk.TOC<-post.cv.fun(cv$llrf.CATE.rlearn,
                           stats=c("stats.bin"),
                           scale=20,
                           cores=1,
                           direction="both",
                           score="cate")

rfA_rates<-colMeans(rf.risk.TOC$ratelist[[2]])
rfB_rates<-colMeans(rf.risk.TOC$ratelist[[1]])
cf_rates<-colMeans(cf.risk.TOC$ratelist[[1]])
llrf_rates<-colMeans(llrf.risk.TOC$ratelist[[1]])
resmat<-matrix(c(rfA_rates,rfB_rates,cf_rates,llrf_rates),4,4)
resmat[c(1,2,4),]<- round(-resmat[c(1,2,4),]*100,3)


# Draw TOC curves

mean.ate<-mean(sapply(cv[[1]],function(x){
  dat<-x[[1]]
  yy<-dat$Y
  ww<-dat$W
  mean(yy[ww==1])-mean(yy[ww==0])
}))


dev.new()
qq<-seq(0,1,0.05)[-1]
qq<-qq[-length(qq)]
qq<-c(qq,0.995,1)
plot(x=qq,
     y=-c(rev(rf.risk.TOC$outmat$ad_DropLowFirst),mean.ate)*100,
     type="l",
     ylim=c(-0.05,0.15)*100,
     lwd=2,
     bty="l",
     ylab="Absolute Risk Difference (%)",
     xlab="Proportion treated",
     cex.axis=1.3,
     cex.lab=1.3)
points(x=qq,
       y=-c(rev(rf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=1,
       lty=2,
       lwd=2)
points(x=qq,
       y=-c(rev(cf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=2,
       lwd=2)
points(x=qq,
       y=-c(rev(llrf.risk.TOC$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=3,
       lwd=2)
abline(h=-mean.ate,lty=2)

legend("topright",col=c(1,1,2,3),lty=c(1,2,1,1),lwd=3,legend=c("RF Baseline Risk A","RF Baseline Risk B",
                                                               "Causal Forest",
                                                               "Local Linear RF R-learner"),
       bty="n")



# Use all of Trial C to generate LOO predictions from m(x) function
mx.dst<-matrix(NA,length(tr.ids),20)
for (mm in 1:20){
  dst<-setA_list[[mm]]
  mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
                          Y=dst[,"y"])
  mx.dst[,mm]<-predict(mfit)$predictions
}

# Do 5-fold cross-validated R-loss

nn<-nrow(dst)
ii<-1:nn
kmat<-matrix(,2,5)
kmat[,1]<-c(1,floor(nn/5))
kmat[,2]<-c(floor(nn/5)+1,floor(2*nn/5))
kmat[,3]<-c(floor(2*nn/5)+1,floor(3*nn/5))
kmat[,4]<-c(floor(3*nn/5)+1,floor(4*nn/5))
kmat[,5]<-c(floor(4*nn/5)+1,nn)

cf.cvpreds<-llrf.cvpreds<-matrix(,nn,5)

for(kk in 1:5){
  
  kv<-kmat[,kk]
  cvte<-kv[1]:kv[2]
  cvtr<-c(1:nn)[-cvte]
  
  for(mm in 1:5){
    
    dst<-setA_list[[mm]]
    
    cftr<-causal_forest(X=as.matrix(dst[cvtr,-c(1:2)]),
                        Y=dst[cvtr,"y"],
                        W=dst[cvtr,"w"],
                        W.hat=0.5,
                        tune.parameters="none",
                        num.trees=2500,
                        honesty=TRUE,   
                        honesty.fraction=0.5)
    
    mxtr<-mx.dst[cvtr,mm]
    yst<-(dst[cvtr,"y"]-mxtr)/(dst[cvtr,"w"]-0.5)
    llrftr<-ll_regression_forest(X=as.matrix(dst[cvtr,-c(1:2)]),
                                 Y=yst,
                                 enable.ll.split=TRUE,
                                 ll.split.weight.penalty=TRUE,
                                 ll.split.lambda=0.1,
                                 num.trees=2500)
    
    cf.cvpreds[cvte,mm]<-predict(cftr, 
                                 newdata = dst[cvte,-c(1:2)])$predictions
    llrf.cvpreds[cvte,mm]<-predict(llrftr, 
                                   newdata = dst[cvte,-c(1:2)],
                                   linear.correction.variables = c(1:21))$predictions
    
  }
  
}


# Calculate R-loss
mx<-apply(mx.dst,1,mean)
dst1<-setA_list[[1]]
yst<-(dst1[,"y"]-mx)/(dst1[,"w"]-0.5)
cf.final<-apply(cf.cvpreds,1,mean)
llrf.final<-apply(llrf.cvpreds,1,mean)

rloss.cfe<-mean((yst-cf.final)^2)
rloss.llrfe<-mean((yst-llrf.final)^2)



