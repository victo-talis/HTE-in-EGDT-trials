
# Make a fast uplift curve generator for ARR

uplift.fast<-function(preds,y,w,wABUC=FALSE){
  
  ord<-order(preds,decreasing=TRUE)
  y<-y[ord]
  w<-w[ord]
  
  x1<-sum(y[w==1])
  x0<-sum(y[w==0])
  n1<-sum(w==1)
  n0<-sum(w==0)
  
  up<-(x1/n1)-(x0/n0)
  
  i<-1
  while(n1>0 & n0>0){
    yi<-y[i]
    wi<-w[i]
    
    x1<-x1-wi*yi
    n1<-n1-wi
    x0<-x0-(1-wi)*yi
    n0<-n0-(1-wi)
    
    up<-c(up,(x1/n1)-(x0/n0))
    
    i<-i+1
  }
  
  if(i<length(preds)){
    up<-c(up,rep(NA,length(preds)-i))
  }
  
  if(wABUC){
    d<-up[1:(length(up)-sum(is.na(up)))]
    q<-rev(c(1:length(d))/length(up))
    mean((1-q)*(d-d[1]))
  } else {
    up
  }

}