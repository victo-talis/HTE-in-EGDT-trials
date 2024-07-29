

# Setup for SHAP calculations

library(grf)

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
#machine<-"C:/Users/vit13/OneDrive - University of Pittsburgh/"


midat_eANDf<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_eANDf.RDS"))
midat_e2f<-readRDS(paste0(machine,"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/midat_list_e2f.RDS"))
midat_f2e<-readRDS(paste0(machine,"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/midat_list_f2e.RDS"))
midat_a2c<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_a2c_MI100.RDS"))
midat_c2a<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_c2a_MI100.RDS"))

ds<-ds000[,c("y","w","trial","trial_id")]
rownames(ds)<-1:nrow(ds)


# All patients: best model is LLRF
dsmALL<-midat_eANDf[[1]]
mfit<-regression_forest(X=as.matrix(dsmALL[,-c(1:2)]),
                        Y=dsmALL[,"y"])
mx<-mfit$predictions[,1]
yst<-(dsmALL$y-mx)/(dsmALL$w-0.5)
dsmALL<-cbind(yst,dsmALL[,-c(1:2)])
saveRDS(dsmALL,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/dsmALL_eANDf.RDS"))

modEFtrain<-grf::ll_regression_forest(X=as.matrix(dsmALL[,-1]),
                                      Y=dsmALL$yst,
                                      enable.ll.split=TRUE,
                                      ll.split.weight.penalty=TRUE,
                                      ll.split.lambda=0.1,
                                      num.trees=1000)
saveRDS(modEFtrain,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/modEFtrain.RDS"))


# Cohort A: best model is LLRF
dsmALL<-midat_e2f[[1]]
mfit<-regression_forest(X=as.matrix(dsmALL[,-c(1:2)]),
                        Y=dsmALL[,"y"])
mx<-mfit$predictions[,1]
yst<-(dsmALL$y-mx)/(dsmALL$w-0.5)
dsmALL<-cbind(yst,dsmALL[,-c(1:2)])
saveRDS(dsmALL,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/dsmALL_e2f.RDS"))

modE2Ftrain<-grf::ll_regression_forest(X=as.matrix(dsmALL[ds$trial=="E",-1]),
                                      Y=dsmALL$yst[ds$trial=="E"],
                                      enable.ll.split=TRUE,
                                      ll.split.weight.penalty=TRUE,
                                      ll.split.lambda=0.1,
                                      num.trees=1000)
saveRDS(modE2Ftrain,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/modE2Ftrain.RDS"))


# Cohort B: best model is LLRF
dsmALL<-midat_f2e[[1]]
mfit<-regression_forest(X=as.matrix(dsmALL[,-c(1:2)]),
                        Y=dsmALL[,"y"])
mx<-mfit$predictions[,1]
yst<-(dsmALL$y-mx)/(dsmALL$w-0.5)
dsmALL<-cbind(yst,dsmALL[,-c(1:2)])
saveRDS(dsmALL,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/dsmALL_f2e.RDS"))

modF2Etrain<-grf::ll_regression_forest(X=as.matrix(dsmALL[ds$trial=="F",-1]),
                                      Y=dsmALL$yst[ds$trial=="F"],
                                      enable.ll.split=TRUE,
                                      ll.split.weight.penalty=TRUE,
                                      ll.split.lambda=0.1,
                                      num.trees=1000)
saveRDS(modF2Etrain,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/modF2Etrain.RDS"))


# ARISE: best model is LLRF
dsmALL<-midat_a2c[[1]]
mfit<-regression_forest(X=as.matrix(dsmALL[,-c(1:2)]),
                        Y=dsmALL[,"y"])
mx<-mfit$predictions[,1]
yst<-(dsmALL$y-mx)/(dsmALL$w-0.5)
dsmALL<-cbind(yst,dsmALL[,-c(1:2)])
saveRDS(dsmALL,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/dsmALL_a2c.RDS"))

modA2Ctrain<-grf::ll_regression_forest(X=as.matrix(dsmALL[ds$trial_id=="A",-1]),
                                       Y=dsmALL$yst[ds$trial_id=="A"],
                                       enable.ll.split=TRUE,
                                       ll.split.weight.penalty=TRUE,
                                       ll.split.lambda=0.1,
                                       num.trees=1000)
saveRDS(modA2Ctrain,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/modA2Ctrain.RDS"))


# ProCESS: best model is Causal Forest
dsmALL<-midat_c2a[[1]]
saveRDS(dsmALL,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/dsmALL_c2a.RDS"))

modC2Atrain<-grf::causal_forest(X=as.matrix(dsmALL[ds$trial_id=="C",-c(1:2)]),
                                Y=dsmALL$y[ds$trial_id=="C"],
                                W=dsmALL$w[ds$trial_id=="C"],
                                W.hat=0.5,
                                num.trees=1000)
saveRDS(modC2Atrain,paste0(machine,"_BoxMigration/Faraaz PROWESS/HTE Manuscript/Results 2023/modC2Atrain.RDS"))





