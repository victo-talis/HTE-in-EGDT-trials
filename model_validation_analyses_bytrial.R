

# Model validation analyses


#machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
oest.dir<-"HTE analysis programs/"
machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/_BoxMigration/"
source(paste0(machine2,oest.dir,"scoring.methods.R"))
source(paste0(machine2,oest.dir,"confirm.analysis.R"))
source(paste0(machine2,oest.dir,"gates.inference.R"))

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"


# Load ds000
source("C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/create_ds000.R")



########################################
# Derive in process, validate in arise #
########################################

machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
midat_list0<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_c2a_MI100.RDS"))
ds<-ds000[ds000$trial_id %in% c("C","A"),c("y","w","trial_id")]
rownames(ds)<-1:nrow(ds)

tr.ids<-which(ds$trial_id=="C")
setA_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})
setB_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][-tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cf.pred.dse<-matrix(NA,sum(ds$trial_id=="A"),100)
cf.pred.dst<-matrix(NA,sum(ds$trial_id=="C"),100)

for (mm in 1:100){
  
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
  
  # mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
  #                         Y=dst[,"y"])
  # mx<-mfit$predictions[,1]
  # yst<-(dst$y-mx)/(dst$w-0.5)
  # rtest_ll<-ll_regression_forest(X=as.matrix(dst[,-c(1:2)]),
  #                                Y=yst,
  #                                enable.ll.split=TRUE,
  #                                ll.split.weight.penalty=TRUE,
  #                                ll.split.lambda=0.1,
  #                                num.trees=5000)
  
  cf.pred.dse[,mm]<-predict(cf, 
                            newdata = dse[,-c(1:2)])$predictions
  # llrf.pred.dse[,mm]<-predict(rtest_ll, 
  #                             newdata = dse[,-c(1:2)],
  #                             linear.correction.variables = c(1:21))$predictions
  cf.pred.dst[,mm]<-predict(cf)$predictions
  # llrf.pred.dst[,mm]<-predict(rtest_ll,
  #                             linear.correction.variables = c(1:21))$predictions
  
}

val_c2a<-list(cf.pred.dst,cf.pred.dse)
              #llrf.pred.dst,llrf.pred.dse)

#saveRDS(val_c2a,paste0(machine,src.dir,"val_c2a.RDS"))
val_c2a<-readRDS(paste0(machine,src.dir,"val_c2a.RDS"))


# Evaluation analyses in trial F (model developed trial E)
# Best model was causal forest

cf.pred.dse<-readRDS(paste0(machine,src.dir,"val_c2a.RDS"))[[2]]

# Construct 1 test set
ydat<-setB_list[[1]][,c("y","w")]
score<-apply(cf.pred.dse,1,mean)
bfitdat<-lapply(1:100,function(x){
  regression_forest(X=as.matrix(setB_list[[x]][,-c(1:2)]),
                    Y=setB_list[[x]][,1])$predictions
})
bfit<-apply(matrix(unlist(bfitdat),length(score),100),1,mean)
mast_cftest<-as.data.frame(cbind(y=ydat[,1],
                                  w=ydat[,2],
                                  score=score,
                                  b=bfit))
#saveRDS(mast_llrtest,paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/mast_llrtest_e2f.RDS"))
#mast_llrtest<-readRDS(paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/mast_llrtest_e2f.RDS"))


# Calculate validation statistics
gates_cf_c2a<-gates.fun(list(mast_cftest),0.5,5,lintest=TRUE,loo=TRUE,unknown.pi=TRUE)


# Calculate relative risk for each quintile
bb<-quantile(score,probs=seq(0,1,1/5))
qq<-as.numeric(cut(score, 
                   breaks=bb,
                   names=FALSE,
                   include.lowest=TRUE))
rrdat<-mast_cftest
rrdat$qq<-qq
rrfun<-function(yy,ww){
  #browser()
  x1<-sum(yy[ww==1])
  n1<-sum(ww==1)
  x0<-sum(yy[ww==0])
  n0<-sum(ww==0)
  rr<-mean(yy[ww==1])/mean(yy[ww==0])
  var<- (n1-x1)/(x1*n1) + (n0-x0)/(x0*n0)
  low<-exp(log(rr)-1.96*sqrt(var))
  hi<-exp(log(rr)+1.96*sqrt(var))
  c(x1,n1,x0,n0,rr,low,hi)
}
rrv<-c(rrfun(rrdat[rrdat$qq==5,"y"],rrdat[rrdat$qq==5,"w"]),
       rrfun(rrdat[rrdat$qq==4,"y"],rrdat[rrdat$qq==4,"w"]),
       rrfun(rrdat[rrdat$qq==3,"y"],rrdat[rrdat$qq==3,"w"]),
       rrfun(rrdat[rrdat$qq==2,"y"],rrdat[rrdat$qq==2,"w"]),
       rrfun(rrdat[rrdat$qq==1,"y"],rrdat[rrdat$qq==1,"w"]))
rrmat<-data.frame(t(matrix(rrv,7,5)))
colnames(rrmat)<-c("x1","n1","x0","n0","rr","rr.lo","rr.hi")

quantile_mat<-data.frame(qq=1:5,
                         mean_iarr=-rev(gates_cf_c2a$`Mean prediction in quantiles`),
                         arr=-rev(gates_cf_c2a$`Theta Median`),
                         arr.lo=-rev(gates_cf_c2a$`Theta Lower Median`),
                         arr.hi=-rev(gates_cf_c2a$`Theta Upper Median`))
quantile_mat<-cbind(quantile_mat,rrmat)

quantile_mat$blp_coef<-c(gates_cf_c2a$`Pval BLP coef slope`,NA,NA,NA,NA)
quantile_mat$blp_pval<-c(gates_cf_c2a$`Pval BLP test`,NA,NA,NA,NA)

#saveRDS(quantile_mat,paste0(machine,src.dir,"quantile_mat_c2a.RDS"))
















########################################
# Derive in arise, validate in process #
########################################

machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
midat_list0<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_a2c_MI100.RDS"))
ds<-ds000[ds000$trial_id %in% c("C","A"),c("y","w","trial_id")]
rownames(ds)<-1:nrow(ds)

tr.ids<-which(ds$trial_id=="A")
setA_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})
setB_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][-tr.ids,]
  out<-na.omit(out)
  rownames(out)<-1:nrow(out)
  out
})

cf.pred.dse<-matrix(NA,sum(ds$trial_id=="C"),100)
cf.pred.dst<-matrix(NA,sum(ds$trial_id=="A"),100)

for (mm in 1:100){
  
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
  
  # mfit<-regression_forest(X=as.matrix(dst[,-c(1:2)]),
  #                         Y=dst[,"y"])
  # mx<-mfit$predictions[,1]
  # yst<-(dst$y-mx)/(dst$w-0.5)
  # rtest_ll<-ll_regression_forest(X=as.matrix(dst[,-c(1:2)]),
  #                                Y=yst,
  #                                enable.ll.split=TRUE,
  #                                ll.split.weight.penalty=TRUE,
  #                                ll.split.lambda=0.1,
  #                                num.trees=5000)
  
  cf.pred.dse[,mm]<-predict(cf, 
                            newdata = dse[,-c(1:2)])$predictions
  # llrf.pred.dse[,mm]<-predict(rtest_ll, 
  #                             newdata = dse[,-c(1:2)],
  #                             linear.correction.variables = c(1:21))$predictions
  cf.pred.dst[,mm]<-predict(cf)$predictions
  # llrf.pred.dst[,mm]<-predict(rtest_ll,
  #                             linear.correction.variables = c(1:21))$predictions
  
}

val_a2c<-list(cf.pred.dst,cf.pred.dse)
              # ,llrf.pred.dst,llrf.pred.dse)

#saveRDS(val_a2c,paste0(machine,src.dir,"val_a2c.RDS"))



# Build table for plotting
# Best model was causal forest

cf.pred.dse<-readRDS(paste0(machine,src.dir,"val_a2c.RDS"))[[2]]

# Construct 1 test set
ydat<-setB_list[[1]][,c("y","w")]
score<-apply(cf.pred.dse,1,mean)
bfitdat<-lapply(1:100,function(x){
  regression_forest(X=as.matrix(setB_list[[x]][,-c(1:2)]),
                    Y=setB_list[[x]][,1])$predictions
})
bfit<-apply(matrix(unlist(bfitdat),length(score),100),1,mean)
mast_cftest<-as.data.frame(cbind(y=ydat[,1],
                                 w=ydat[,2],
                                 score=score,
                                 b=bfit))
#saveRDS(mast_llrtest,paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/mast_llrtest_e2f.RDS"))
#mast_llrtest<-readRDS(paste0(machine,"Faraaz PROWESS/HTE Manuscript/Results 2023/Calibration/mast_llrtest_e2f.RDS"))


# Calculate validation statistics
gates_cf_a2c<-gates.fun(list(mast_cftest),0.5,5,lintest=TRUE,loo=TRUE,unknown.pi=TRUE)


# Calculate relative risk for each quintile
bb<-quantile(score,probs=seq(0,1,1/5))
qq<-as.numeric(cut(score, 
                   breaks=bb,
                   names=FALSE,
                   include.lowest=TRUE))
rrdat<-mast_cftest
rrdat$qq<-qq
rrfun<-function(yy,ww){
  #browser()
  x1<-sum(yy[ww==1])
  n1<-sum(ww==1)
  x0<-sum(yy[ww==0])
  n0<-sum(ww==0)
  rr<-mean(yy[ww==1])/mean(yy[ww==0])
  var<- (n1-x1)/(x1*n1) + (n0-x0)/(x0*n0)
  low<-exp(log(rr)-1.96*sqrt(var))
  hi<-exp(log(rr)+1.96*sqrt(var))
  c(x1,n1,x0,n0,rr,low,hi)
}
rrv<-c(rrfun(rrdat[rrdat$qq==5,"y"],rrdat[rrdat$qq==5,"w"]),
       rrfun(rrdat[rrdat$qq==4,"y"],rrdat[rrdat$qq==4,"w"]),
       rrfun(rrdat[rrdat$qq==3,"y"],rrdat[rrdat$qq==3,"w"]),
       rrfun(rrdat[rrdat$qq==2,"y"],rrdat[rrdat$qq==2,"w"]),
       rrfun(rrdat[rrdat$qq==1,"y"],rrdat[rrdat$qq==1,"w"]))
rrmat<-data.frame(t(matrix(rrv,7,5)))
colnames(rrmat)<-c("x1","n1","x0","n0","rr","rr.lo","rr.hi")

quantile_mat<-data.frame(qq=1:5,
                         mean_iarr=-rev(gates_cf_a2c$`Mean prediction in quantiles`),
                         arr=-rev(gates_cf_a2c$`Theta Median`),
                         arr.lo=-rev(gates_cf_a2c$`Theta Lower Median`),
                         arr.hi=-rev(gates_cf_a2c$`Theta Upper Median`))
quantile_mat<-cbind(quantile_mat,rrmat)

quantile_mat$blp_coef<-c(gates_cf_a2c$`Pval BLP coef slope`,NA,NA,NA,NA)
quantile_mat$blp_pval<-c(gates_cf_a2c$`Pval BLP test`,NA,NA,NA,NA)

#saveRDS(quantile_mat,paste0(machine,src.dir,"quantile_mat_a2c.RDS"))

