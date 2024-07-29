

# Model validation analyses


machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
src.dir<-"_BoxMigration/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
oest.dir<-"_BoxMigration/HTE analysis programs/"
#machine2<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/_BoxMigration/"
source(paste0(machine,oest.dir,"scoring.methods.R"))
source(paste0(machine,oest.dir,"confirm.analysis.R"))
source(paste0(machine,oest.dir,"gates.inference.R"))
source(paste0(machine,oest.dir,"cross.validation.R"))
source(paste0(machine,oest.dir,"post.process_v2.R"))



# Load ds000
#source("C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")
source("C:/Users/VIT13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")



########################################
# Derive in cohort A, validate in cohort B #
########################################

src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
midat_list0<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_e2f_MI100.RDS"))
ds<-ds000[ds000$trial %in% c("E","F"),c("y","w","trial")]
rownames(ds)<-1:nrow(ds)

MI<-10

tr.ids<-which(ds$trial=="E")
setA_list<-lapply(1:MI,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})
setB_list<-lapply(1:MI,function(x){
  out<-midat_list0[[x]][-tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cf.pred.dse<-llrf.pred.dse<-matrix(NA,sum(ds$trial_id=="F"),MI)
cf.pred.dst<-llrf.pred.dst<-matrix(NA,sum(ds$trial_id=="E"),MI)

for (mm in 1:MI){
  
  dst<-setA_list[[mm]]
  dse<-setB_list[[mm]]
  
  cf<-causal_forest(X=as.matrix(dst[,-c(1:2)]),
                    Y=dst[,"y"],
                    W=dst[,"w"],
                    W.hat=0.5,
                    tune.parameters="none",
                    num.trees=5000,
                    honesty=TRUE,   
                    honesty.fraction=0.5)
  
  mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
                          Y=dst[,"y"])
  mx<-mfit$predictions[,1]
  yst<-(dst$y-mx)/(dst$w-0.5)
  rtest_ll<-ll_regression_forest(X=as.matrix(dst[,-c(1:2)]),
                                 Y=yst,
                                 enable.ll.split=TRUE,
                                 ll.split.weight.penalty=TRUE,
                                 ll.split.lambda=0.1,
                                 num.trees=5000)
  
  cf.pred.dse[,mm]<-predict(cf, 
                            newdata = dse[,-c(1:2)])$predictions
  llrf.pred.dse[,mm]<-predict(rtest_ll, 
                              newdata = dse[,-c(1:2)],
                              linear.correction.variables = c(1:21))$predictions
  cf.pred.dst[,mm]<-predict(cf)$predictions
  llrf.pred.dst[,mm]<-predict(rtest_ll,
                              linear.correction.variables = c(1:21))$predictions
  
}

val_e2f<-list(cf.pred.dst,cf.pred.dse,llrf.pred.dst,llrf.pred.dse)

#saveRDS(val_e2f,paste0(machine,src.dir,"val_e2f.RDS"))
val_e2f<-readRDS(paste0(machine,src.dir,"val_e2f.RDS"))



# need to insert calculation of AUTOC, Adj Qini, and R-loss here

up<-oest.cv.fun(midat_list0[1:10],
                 cv.tr.pr=0.5,
                 cv_k=1,
                 trids0=tr.ids,
                 y="y",
                 w="w",
                 x=names(setA_list[[1]])[-c(1:2)],
                 pi=0.5,
                 mods=c("llrf.CATE.rlearn"),
                 cores=1,
                 seed=15223)

llrf.risk.TOC<-post.cv.fun(up$llrf.CATE.rlearn,
                           stats=c("stats.bin"),
                           scale=20,
                           cores=1,
                           direction="both",
                           score="cate")
llrf_rates<-colMeans(llrf.risk.TOC$ratelist[[1]])



#


########################################
# Derive in cohort B, validate in cohort A #
########################################

src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
midat_list0<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_f2e_MI100.RDS"))
ds<-ds000[ds000$trial %in% c("E","F"),c("y","w","trial")]
rownames(ds)<-1:nrow(ds)

MI<-10

tr.ids<-which(ds$trial=="F")
setA_list<-lapply(1:MI,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})
setB_list<-lapply(1:MI,function(x){
  out<-midat_list0[[x]][-tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cf.pred.dse<-llrf.pred.dse<-matrix(NA,sum(ds$trial_id=="F"),MI)
cf.pred.dst<-llrf.pred.dst<-matrix(NA,sum(ds$trial_id=="E"),MI)

for (mm in 1:MI){
  
  dst<-setA_list[[mm]]
  dse<-setB_list[[mm]]
  
  cf<-causal_forest(X=as.matrix(dst[,-c(1:2)]),
                    Y=dst[,"y"],
                    W=dst[,"w"],
                    W.hat=0.5,
                    tune.parameters="none",
                    num.trees=5000,
                    honesty=TRUE,   
                    honesty.fraction=0.5)
  
  mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
                          Y=dst[,"y"])
  mx<-mfit$predictions[,1]
  yst<-(dst$y-mx)/(dst$w-0.5)
  rtest_ll<-ll_regression_forest(X=as.matrix(dst[,-c(1:2)]),
                                 Y=yst,
                                 enable.ll.split=TRUE,
                                 ll.split.weight.penalty=TRUE,
                                 ll.split.lambda=0.1,
                                 num.trees=5000)
  
  cf.pred.dse[,mm]<-predict(cf, 
                            newdata = dse[,-c(1:2)])$predictions
  llrf.pred.dse[,mm]<-predict(rtest_ll, 
                              newdata = dse[,-c(1:2)],
                              linear.correction.variables = c(1:21))$predictions
  cf.pred.dst[,mm]<-predict(cf)$predictions
  llrf.pred.dst[,mm]<-predict(rtest_ll,
                              linear.correction.variables = c(1:21))$predictions
  
}

val_e2f<-list(cf.pred.dst,cf.pred.dse,llrf.pred.dst,llrf.pred.dse)

#saveRDS(val_f2e,paste0(machine,src.dir,"val_f2e.RDS"))
val_f2e<-readRDS(paste0(machine,src.dir,"val_f2e.RDS"))



# need to insert calculation of AUTOC, Adj Qini, and R-loss here

up<-oest.cv.fun(midat_list0[1:10],
                cv.tr.pr=0.5,
                cv_k=1,
                trids0=tr.ids,
                y="y",
                w="w",
                x=names(setA_list[[1]])[-c(1:2)],
                pi=0.5,
                mods=c("llrf.CATE.rlearn"),
                cores=1,
                seed=15223)

llrf.risk.TOC<-post.cv.fun(up$llrf.CATE.rlearn,
                           stats=c("stats.bin"),
                           scale=20,
                           cores=1,
                           direction="both",
                           score="cate")
llrf_rates<-colMeans(llrf.risk.TOC$ratelist[[1]])