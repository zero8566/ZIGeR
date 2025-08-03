
library(stats4)
library(pscl)
library(MASS)
library(gee)
library(Matrix)
library(xtable)
library(stats4)
library(pscl)
library(MASS)
library(gee)
library(Matrix)
library(xtable)
library("stats")
library("MASS")
library(splines)
library(VGAM)


lzige0<-function(p0,beta0){
  cov=matrix(rep(1,n),ncol=1)
  y=dat$y;

  r1=(y==0)
  d1=(y>=cut)
  pro=1/(1+exp(cov%*%beta0))
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  l1<-((1-p0)*dgeom(y,prob=pro));
  l2<-(1-p0)*leiji
  l3=p0+(1-p0)*pro;

  ll=-sum((1-r1)*((1-d1)*log(l1)+d1*log(l2))+r1*(1-d1)*log(l3));
  ll;
}
lge0<-function(beta0){
  
  cov=matrix(rep(1,n),ncol=1)
  y=dat$y;
  
  pro=1/(1+exp(cov%*%beta0))
  r1=(y==0)
  d1=(y>=cut)
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  # pro=plogis(beta);
  l1=(1-d1)*log(dgeom(y,prob=pro))+d1*log(leiji);
  ll=-sum(l1);
}



lzige1<-function(p0,beta0,beta1){
  y=dat$y;
  x=dat$x;
  cov=cbind(1,x);
  
  r1=(y==0)
  d1=(y>=cut)
  betatr=c(beta0,beta1);
  pro=1/(1+exp(cov%*%betatr))
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  l1<-((1-p0)*dgeom(y,prob=pro));
  l2<-(1-p0)*leiji
  l3=p0+(1-p0)*pro;
  
  ll=-sum((1-r1)*((1-d1)*log(l1)+d1*log(l2))+r1*(1-d1)*log(l3));
  ll;
 
}


lge1<-function(theta){ 
    beta0<-theta[1]
    beta1<-theta[2]
    r1=(y==0)
    d1=(y>=cut)
  y=dat$y;
  x=dat$x;
  cov=cbind(1,x);
  betatr=c(beta0,beta1);
  pro=1/(1+exp(cov%*%betatr))
  leiji<-pgeom(cut-1,pro,lower.tail = F)
  l1=(1-d1)*log(dgeom(y,prob=pro))+d1*log(leiji);
  ll=-sum(l1);
  ll;
}

  
  test.score<- function( beta,outcome,cov,cut){
    y=outcome
    d1=(y>=cut)
    n=length(outcome);
    ri=(outcome==0);
    
    mu=exp(cov %*% beta);
    dim(cov %*% beta)
    c1=1/(1+mu)
    
    p=length(beta);
    
    dU=matrix(rep(0, n^2), nrow=n);##n*n矩阵=50*50(无协变量时)
    dim(dU)##50*50
    diag(dU)=(1-d1)*(c1+y*c1-1)+d1*((1-c1)^cut)*c1*(1+cut-(1/c1))
    
    uWT1=t(cov)%*%dU;
    uWT2=(ri-c1)/c1;
    dim(uWT2)
    
    L11=uWT1%*%t(uWT1);
    
    L12=uWT1%*%(uWT2)
    L21=t(L12)
    L22=t(uWT2)%*%(uWT2)
    LL=rbind(cbind(L11,L12),cbind(L21,L22));
    dim(LL)
    
    score<-sum((ri-c1)/c1)*sqrt(solve(LL)[p+1,p+1])
    
    
    if(det(LL)!=0 && solve(LL)[p+1,p+1]>=0){
      try(score<-sum((ri-c1)/c1)*sqrt(solve(LL)[p+1,p+1]));
    }
    if(det(LL)!=0 && solve(LL)[p+1,p+1]<0){
      score<-88888;
    }
    if(det(LL)==0 ){
      score<-99999
    }
    
    return(score);
    
  }
  

  

test.new<- function( beta,outcome,cov,cut) {
   

    y=outcome
    d1=(y>=cut)
    n=length(outcome);
    ri=(outcome==0);
   
    mu=exp(cov %*% beta);
    c1=1/(1+mu)
    
    d1=(outcome>=cut)
    p=length(beta);
    
    ri= (y==0);
    pi=exp(-mu);
    
   
    dU=matrix(rep(0, n^2), nrow=n);##n*n
    dim(dU)
    
    diag(dU)=(1-d1)*(c1+y*c1-1)+d1*((1-c1)^cut)*c1*(1+cut-(1/c1))##SCORE一样
   
    dim(dU)
    uWT1=ri-c1;
    uWT2=t(cov)%*%dU;##1*n%*%n*n=1*n,,COV是协变量x加一个截距项1
    
    dim(cov)
    dim(uWT2)
    
    A11=-n;
    A12=t(uWT1)%*%t(uWT2);
    A21=matrix(rep(0,p),ncol=1);
    # p=1
    A22=uWT2%*%t(uWT2);
    dim(A21)
    AA=rbind(cbind(A11,-A12),cbind(A21,-A22));
    A=AA/n;
    
    
    B11=t(uWT1)%*%uWT1;
    B12=t(uWT1)%*%t(uWT2);
    B22=uWT2%*%t(uWT2);###t(matrix(rep(1,n),ncol=1))%*%dU%*%t(matrix(rep(1,n),ncol=1))%*%dU===
    # dim(B22)
    B21=t(B12);
    BB=rbind(cbind(B11,B12),cbind(B21,B22));
    B=BB/n;
    
    
    TT=solve(A)%*%B%*%t(solve(A));
    if(det(A)!=0 && (TT)[1,1]>=0){
      try(snew<-(sqrt(n)*mean(ri-c1))/sqrt((TT)[1,1]));
      # di<-(y>=cut)
      # ri*di
    }
    if(det(A)!=0 && (TT)[1,1]<0){
      snew<-88888;
    }
    if(det(A)==0){
      snew<-99999;
    }
    return(snew);
  }
  
 

