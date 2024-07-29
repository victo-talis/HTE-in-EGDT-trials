

machine<-"C:/Users/vit13/OneDrive - University of Pittsburgh/"
#machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
output.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/shap_raw/"
dat.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/shap_cluster_test_2024/"
code.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"

# Figure for tr0 (model on 100% sample)

tr<-"tr0/"
dsm_name<-"dsmALL_eANDf.RDS"
#dsm_name<-"dsmALL_e2f.RDS"
#dsm_name<-"dsmALL_f2e.RDS"
#dsm_name<-"dsmALL_a2c.RDS"
#dsm_name<-"dsmALL_c2a.RDS"

# Consolidate the SHAPS output
files<-list.files(paste0(machine,output.dir,tr))
startnum<-sapply(files,function(x){
  as.numeric(gsub("START", "", unlist(strsplit(x,"_"))[3]))
})
ord<-order(startnum)
files<-files[ord]

rds<-readRDS(paste0(machine,output.dir,tr,files[1]))
shaps<-do.call("rbind",rds)
len<-c(length(rds))
for(ff in 2:length(files)){
  rds<-readRDS(paste0(machine,output.dir,tr,files[ff]))
  len[ff]<-length(rds)
  shaps<-rbind(shaps,do.call("rbind",rds))
}

# Get the patient data
x_test<-readRDS(paste0(machine,dat.dir,dsm_name))[,-1]

shaps<-(-shaps*100)

library(grf)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(gridExtra)
source(paste0(machine,code.dir,"shap_plot_fun.R"))

windowsFonts("Arial" = windowsFont("Arial"))

plot_var_names<-c("Age","Respiratory Rate","Heart rate","Hemoglobin","Systolic BP","Bilirubin",
             "Mean Arterial Pressure","Blood Urea Nitrogen","Creatinine","White blood cell count",
             "APACHE II","PaO2","Charlson","Glasgow Coma Scale","SOFA","Albumin","Lactate",
             "Glucose","Platelets","O2 Saturation","Temperature","Male","Lung infection",
            "Urinary tract infection","Abdominal infection","Mechanical Vent","Vasopressor use")


# Importance + Beeswarm plot

mean_shap<-plot_feature_importance(-shaps[,-1],
                                  desc_sorting = TRUE,
                                  max_vars = ncol(shaps[,-1]),
                                  title = "",
                                  subtitle = NULL,
                                  custom_var_names = plot_var_names,
                                  tsize=12,
                                  psize=2,
                                  ljust=c(50,0),
                                  BREAKS=c(0,0.5,1.0,1.5,2.0,2.5))

ind_shap<-plot_feature_beeswarm2(shaps[,-1], 
                       x_test, 
                       title = "", 
                       subtitle = NULL,
                       custom_var_names = plot_var_names,
                       width=1,
                       wh_bin=22:27,
                       xlim=c(-0.25,0.15)*100,
                       tsize=12,psize=2,ljust=c(0.5,0),
                       max_vars=NULL) 

dev.new(width=12,height=10,noRStudioGD = TRUE)
grid.arrange(grobs=list(mean_shap,ind_shap),
             width=c(1,1,1,1,1),
             layout_matrix=rbind(c(1,1,2,2,2)))




dev.new(height=12,width=15,noRStudioGD = TRUE)
par(mfrow=c(2,3))
scatter.smooth(x=x_test$albumin,
               y=shaps$albumin,
               main="Albumin",
               ylab="Albumin SHAP",
               xlab="Albumin value",
               xlim=c(0,6.2),
               lpars = list(col = "red", lwd = 3, lty = 3))
scatter.smooth(x=x_test$temp,
               y=shaps$temp,
               main="Body Temperature",
               ylab="Temperature SHAP",
               xlab="Temperature value",
               xlim=c(28,43),
               lpars = list(col = "red", lwd = 3, lty = 3))
scatter.smooth(x=x_test$hr,
               y=shaps$hr,
               main="Heart Rate",
               ylab="Heart Rate SHAP",
               xlab="Heart Rate value",
               xlim=c(30,200),
               lpars = list(col = "red", lwd = 3, lty = 3))
scatter.smooth(x=x_test$age,
     y=shaps$age,
     main="Age",
     ylab="Age SHAP",
     xlab="Age value",
     lpars = list(col = "red", lwd = 3, lty = 3))
scatter.smooth(x=x_test$ln_g,
     y=shaps$ln_g,
     main="Log Glucose",
     ylab="Log Glucose SHAP",
     xlab="Log Glucose value",
     xlim=c(2.5,7.2),
     lpars = list(col = "red", lwd = 3, lty = 3))
scatter.smooth(x=x_test$ln_bun,
     y=shaps$temp,
     main="Log BUN",
     ylab="Log BUN SHAP",
     xlab="Log BUN value",
     lpars = list(col = "red", lwd = 3, lty = 3))



# Individual SHAP plot

# Need the predictions
proj.dir<-"C:/Users/vit13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/shap_cluster_test_2024/"
# dsm<-readRDS(paste0(proj.dir,"dsmALL_eANDf.RDS"))
# p0<- -mean(dsm$yst)*100
modtrain<-readRDS(paste0(proj.dir,"modEFtrain.RDS"))
preds<-rowSums(shaps)


# Dominated by albumin
wh1<-which(shaps$albumin< -13 & shaps$albumin> -17 &
             preds< -14 & preds>-17) #3
harm1<-plot_contribution(shaps,
                         x_test,
                         obs=wh1[3],
                         title="Patient 1" ,
                         axtxtsz=12,
                         axttlsz=12,
                         pltttlsz=14,
                         max_vars=5,
                         manual_yaxis_vals=c("Intercept","Albumin = 1.5",
                                            "BUN = 87.4","Temperature = 36.1","Heart rate = 96","SOFA = 3.0",
                                            "+ All others","Final prediction"),
                         min_max=c(-16,3))


# Albumin and age

wh5<-which(shaps$albumin< 0 & 
             shaps$age< -5 & shaps$age> -7 &
           preds< -10)
harm3<-plot_contribution(shaps,
                         x_test,
                         obs=wh5[2],
                         title="Patient 1",
                         axtxtsz=12,
                         axttlsz=12,
                         pltttlsz=14,
                         max_vars=5,
                         manual_yaxis_vals=c("Intercept","Albumin = 1.4",
                                             "Age = 25","Glucose = 88","Temperature = 37.1","Platelets = 119",
                                             "+ All others","Final prediction"),
                         min_max=c(-32,3))


# Benefit: Albumin and age

wh5<-which(shaps$albumin>0 & 
             shaps$age> 5 &
             x_test$age>60 &
             #shaps$hr>5 &
             preds> 5)
ben1<-plot_contribution(shaps,
                         x_test,
                         obs=wh5[1],
                         title="Patient 1",
                         axtxtsz=12,
                         axttlsz=12,
                         pltttlsz=14,
                         max_vars=5,
                         manual_yaxis_vals=c("Intercept","Age = 88",
                                             "SBP = 166","Respiratory rate = 45","APACHE = 27","Heart rate = 134",
                                             "+ All others","Final prediction"),
                         min_max=c(-1,27))

#wh5[7]

ben3<-plot_contribution(shaps,
                        x_test,
                        obs=wh5[2],
                        title="Patient 1",
                        axtxtsz=12,
                        axttlsz=12,
                        pltttlsz=14,
                        max_vars=5,
                        manual_yaxis_vals=c("Intercept","Age = 90",
                                            "Albumin = 3.5","BUN = 79","SOFA = 11","Temperature = 37.7",
                                            "+ All others","Final prediction"),
                        min_max=c(-3,25))


dev.new(width=12,height=6.6,noRStudioGD = TRUE)
grid.arrange(grobs=list(harm1,harm3,
                        ben1,ben3),
             width=c(1,1,1,1),
             layout_matrix=matrix(1:4,2,2))


