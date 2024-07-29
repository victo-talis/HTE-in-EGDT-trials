
# Sensitivity analysis: Comparison of methods with and without 

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
hte.dir<-"_BoxMigration/HTE analysis programs/"
source(paste0(machine,hte.dir,"cross.validation.R"))
source(paste0(machine,hte.dir,"scoring.methods.R"))
source(paste0(machine,hte.dir,"post.process_v2.R"))

# CV in the 100% sample, comparing LLRF with and without clustering

site<-ds000$site
trial_id<-ds000$trial_id
trial_id<-ifelse(trial_id=="A","0001","0002")

site0<-as.numeric(paste0(site,trial_id))

trids<-sample(1:length(site0),length(site0)/2)
teids<-c(1:length(site0))[-trids]


# test<-regression_forest(X=as.matrix(setA_list[[1]][,3:29]),
#                         Y=setA_list[[1]][,"y"],
#                         cluster=setA_list[[1]][,"site"])
# 
# 
# test<-llrf.CATE.rlearn(train=setA_list[1], 
#                                  cvt=trids, 
#                                  cve=teids, 
#                                  M=1, 
#                                  y="y", 
#                                  w="w", 
#                                  x=names(setA_list[[1]])[-c(1:2,30)], 
#                                  pi=0.5, 
#                                  ntrees=5000,
#                                 cluster="site")



midat_list0<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_eANDf.RDS"))
ds<-ds000[ds000$trial %in% c("E","F"),c("y","w","trial")]
rownames(ds)<-1:nrow(ds)

tr.ids<-1:nrow(midat_list0[[1]])
setA_list<-lapply(1:100,function(x){
  out<-midat_list0[[x]][tr.ids,]
  out<-na.omit(out)
  out$site<-site0
  rownames(out)<-1:nrow(out)
  out
})

cv.noclust<-oest.cv.fast(setA_list[1:20],
                         cv.tr.pr=0.5,
                         cv_k=20,
                         y="y",
                         w="w",
                         x=names(setA_list[[1]])[-c(1:2)],
                         pi=0.5,
                         mods=c("llrf.CATE.rlearn"),
                         cores=1,
                         seed=15223)

cv.clust<-oest.cv.fast(setA_list[1:20],
                         cv.tr.pr=0.5,
                         cv_k=20,
                         y="y",
                         w="w",
                         x=names(setA_list[[1]])[-c(1:2)],
                         pi=0.5,
                         mods=c("llrf.CATE.rlearn"),
                         cores=1,
                         seed=15223,
                         cluster="site")



llrf.noclust<-post.cv.fun(cv.noclust$llrf.CATE.rlearn,
                         stats=c("stats.bin"),
                         scale=10,
                         cores=1,
                         direction="both",
                         score="cate")
llrf.clust<-post.cv.fun(cv.clust$llrf.CATE.rlearn,
                          stats=c("stats.bin"),
                          scale=10,
                          cores=1,
                          direction="both",
                          score="cate")

llrf_noclust_rates<-colMeans(llrf.noclust$ratelist[[1]])
llrf_clust_rates<-colMeans(llrf.clust$ratelist[[1]])
resmat<-matrix(c(llrf_noclust_rates,llrf_clust_rates,cf_rates,llrf_rates),4,4)
resmat[c(1,2,4),]<- round(-resmat[c(1,2,4),]*100,3)





mean.ate<-mean(sapply(cv.clust[[1]],function(x){
  dat<-x[[1]]
  yy<-dat$Y
  ww<-dat$W
  mean(yy[ww==1])-mean(yy[ww==0])
}))

dev.new()
qq<-seq(0,1,0.1)[-1]
qq<-qq[-length(qq)]
qq<-c(qq,0.995,1)
plot(x=qq,
     y=-c(rev(llrf.noclust$outmat$ad_DropHighFirst),mean.ate)*100,
     type="l",
     ylim=c(-0.05,0.15)*100,
     lwd=2,
     bty="l",
     ylab="Absolute Risk Difference (%)",
     xlab="Proportion treated",
     cex.axis=1.3,
     cex.lab=1.3)
points(x=qq,
       y=-c(rev(llrf.clust$outmat$ad_DropHighFirst),mean.ate)*100,
       type="l",
       col=1,
       lty=2,
       lwd=2)
abline(h=-mean.ate,lty=2)

legend("topright",col=c(1,1,2,3),lty=c(1,2,1,1),lwd=3,legend=c("Vanilla LLRF R-learner w/o site",
                                                               "Cluster-robust LLRF R-learner"),
       bty="n")


