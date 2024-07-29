

# Traditional analyses

machine<-"C:/Users/vbtal/OneDrive - University of Pittsburgh/"
midat_eANDf<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/midat_list_eANDf.RDS"))

alb<-unlist(lapply(midat_eANDf,function(x){
  x$albumin
}))
age<-unlist(lapply(midat_eANDf,function(x){
  x$age
}))

qsalb<-quantile(alb,probs=c(0.33,0.66))
qsage<-quantile(age,probs=c(0.33,0.66))

mimod<-lapply(midat_eANDf[1:20],function(x){
  #browser()
  x$ee_alb<-as.numeric(as.character(cut(x$albumin,c(0,qsalb,100),include.lowest=TRUE,labels=1:3)))
  x$ee_age<-as.numeric(as.character(cut(x$age,c(0,qsage,100),include.lowest=TRUE,labels=1:3)))
  x
})

fit_mialb<-lapply(mimod,function(x){
  glm(y~w*albumin,
      family="binomial",
      data=x)
})
summary(pool(fit_mialb))

fit_miage<-lapply(mimod,function(x){
  glm(y~w*age,
      family="binomial",
      data=x)
})
summary(pool(fit_miage))

cicalc<-function(s){
  p<-round(exp(s[2,2]),2)
  hi<-round(exp(s[2,2]+1.96*s[2,3]),2)
  lo<-round(exp(s[2,2]-1.96*s[2,3]),2)
  paste0(p," (",lo,",",hi,")")
}

fit_mialbee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_alb==1,])
})
cicalc(summary(pool(fit_mialbee)))

fit_mialbee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_alb==2,])
})
cicalc(summary(pool(fit_mialbee)))
fit_mialbee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_alb==3,])
})
cicalc(summary(pool(fit_mialbee)))



fit_miageee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_age==1,])
})
cicalc(summary(pool(fit_miageee)))

fit_miageee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_age==2,])
})
cicalc(summary(pool(fit_miageee)))
fit_miageee<-lapply(mimod,function(x){
  glm(y~w,
      family="binomial",
      data=x[x$ee_age==3,])
})
cicalc(summary(pool(fit_miageee)))
