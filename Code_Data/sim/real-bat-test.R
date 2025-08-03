

lge13<-function(theta){ 
  
  beta0<-theta[1]
  beta1<-theta[2]
  beta2<-theta[3]
  r1=(y==0)
  d1=(y>=cut)
  y=dat$y;
  x<-as.matrix(cbind(data[,11:12]))
  cov=cbind(1,x);
  betatr=c(beta0,beta1,beta2);
  pro=1/(1+exp(cov%*%betatr))
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  l1=(1-d1)*log(dgeom(y,prob=pro))+d1*log(leiji);
  ll=-sum(l1);
  ll;
}



lzige13<-function(p0,beta0,beta1,beta2){
  y=dat$y;
  x<-as.matrix(cbind(data[,11:12]))
  cov=cbind(1,x);
  
  r1=(y==0)
  d1=(y>=cut)
  betatr=c(beta0,beta1,beta2);
  pro=1/(1+exp(cov%*%betatr))
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  l1<-((1-p0)*dgeom(y,prob=pro));
  l2<-(1-p0)*leiji
  l3=p0+(1-p0)*pro;
  
  ll=-sum((1-r1)*((1-d1)*log(l1)+d1*log(l2))+r1*(1-d1)*log(l3));
  ll;
  
}


