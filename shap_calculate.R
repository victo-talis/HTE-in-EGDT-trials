# args = c(2,1,200,100,100,
#          1,2,1)
args =  commandArgs(trailingOnly=TRUE)

#options(error=browser)

MC<-as.numeric(args[1])
start<-as.numeric(args[2])
nb<-as.numeric(args[3])            # n_batch
ns<-as.numeric(args[4])            # n_samples
nc<-as.numeric(args[5])            # n_combinations
tr<-as.numeric(args[6])            # 0=all patients; 1=cohort a; 2=cohort b; 3=arise; 4=process

MC
start
nb
ns
nc
tr

return(args)

#proj.dir<-"C:/Users/vit13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/shap_cluster_test_2024/"
proj.dir<-"/ihome/vtalisa/vit13/CCM/Faraaz/egdt_shapr_2024/"

library(grf)
library(shapr)

if(tr==0){
  dsm<-readRDS(paste0(proj.dir,"dsmALL_eANDf.RDS"))
  modtrain<-readRDS(paste0(proj.dir,"modEFtrain.RDS"))
} else if (tr==1){
  dsm<-readRDS(paste0(proj.dir,"dsmALL_e2f.RDS"))
  modtrain<-readRDS(paste0(proj.dir,"modE2Ftrain.RDS"))
} else if (tr==2){
  dsm<-readRDS(paste0(proj.dir,"dsmALL_f2e.RDS"))
  modtrain<-readRDS(paste0(proj.dir,"modF2Etrain.RDS"))
}else if (tr==3){
  dsm<-readRDS(paste0(proj.dir,"dsmALL_a2c.RDS"))
  modtrain<-readRDS(paste0(proj.dir,"modA2Ctrain.RDS"))
} else if (tr==4){
  dsm<-readRDS(paste0(proj.dir,"dsmALL_c2a.RDS"))
  modtrain<-readRDS(paste0(proj.dir,"modC2Atrain.RDS"))
}


# Utility functions
if(tr %in% c(0,1,2,3)){
  predict_model.ll_regression_forest<-function(x, newdata){
    #browser()
    predict(x, 
            newdata = newdata,
            linear.correction.variables = c(1:27))$predictions
  }
  get_model_specs.ll_regression_forest <- function(x){
    feature_list = list()
    feature_list$labels <- colnames(x$X.orig)
    m <- length(feature_list$labels)
    feature_list$classes <- class(x$X.orig[,1])
    feature_list$factor_levels <- setNames(vector("list", m), feature_list$labels)
    feature_list$factor_levels[feature_list$classes=="factor"] <- NA # the model object doesn't contain factor levels info
    
    return(feature_list)
  }
  p0<-mean(dsm$yst)
  
} else {
  
  predict_model.causal_forest<-function(x, newdata){
    predict(x, 
            newdata = newdata)$predictions
  }
  get_model_specs.causal_forest <- function(x){
    feature_list = list()
    feature_list$labels <- colnames(x$X.orig)
    m <- length(feature_list$labels)
    feature_list$classes <- class(x$X.orig[,1])
    feature_list$factor_levels <- setNames(vector("list", m), feature_list$labels)
    feature_list$factor_levels[feature_list$classes=="factor"] <- NA # the model object doesn't contain factor levels info
    
    return(feature_list)
  }
  
  yy<-dsm$y
  ww<-dsm$w
  p0<-mean(yy[ww==1])-mean(yy[ww==0])
}




if ((start+MC-1)<=2480){
  
  explainer<-shapr(dsm, 
                   modtrain,
                   n_combinations=nc)  # approx exponential in n_combinations; shoot for 10000
  explanation<-explain(dsm[start:(start+MC-1),],
                     explainer,
                     approach="ctree",
                     n_batches=nb,  # sampling and prediction can be done in batches
                     prediction_zero=p0,
                     mincriterion = 0.95,
                     minsplit = 20,
                     minbucket = 7,
                     sample = TRUE,
                     n_samples=ns)  # approx linear in n_sample Monte Carlo; default=1000

  saveRDS(explanation,paste0(proj.dir,"output/output_tr",tr,"_START",start,"_END",start+MC-1))

} else if (start<2481 & (start+MC-1)>2480) {

  explainer<-shapr(dsm, 
                   modtrain,
                   n_combinations=nc)  # approx exponential in n_combinations; shoot for 10000
  explanation<-explain(dsm[start:2480,],
                       explainer,
                       approach="ctree",
                       n_batches=nb,  # sampling and prediction can be done in batches
                       prediction_zero=p0,
                       mincriterion = 0.95,
                       minsplit = 20,
                       minbucket = 7,
                       sample = TRUE,
                       n_samples=ns)  # approx linear in n_sample Monte Carlo; default=1000
  
  saveRDS(explanation,paste0(proj.dir,"output/output_tr",tr,"_START",start,"_END",2480))
}
