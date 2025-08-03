rm(list=ls())
setwd("D:/0Adas/Code_and_Data/")

source("./sim/real-likehood.R")
# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/")

setwd("D:/0Adas/results/30%删失/")
#########################################################################
#########################################################################
#   Section 3.1.1 no covariate 
# compare the four tests on controling  type I errors
###########################################################################


nsim=1000; 
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000
tp=seq(0.5,2,0.5)
# tp=seq(1.5,3,0.5)
set.seed(2024)

ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(tp)));
beta_mle=ymin=ymax=w_hat=array(dim=c(nsim, length(ns),length(tp)));

rcpro=0.3


for(jj in 1:length(tp)){k=0;
for( n in ns){k=k+1;
for (i in 1:nsim){
  cov<-matrix(rep(1,n),ncol=1)
  pro=1/(1+exp(tp[jj]));  #the success probability rho  #pro=1-plogis(tp[jj])
  yy=rgeom(n,prob=pro);
  rightcensor<-n-(length(yy)*rcpro)
  yy_sorted <- sort(yy)
  cut<- yy_sorted[rightcensor]
  y<- ifelse(yy> cut, cut, yy)
  dat=data.frame(y=y);
  a<-glm(y ~ 1,family=negative.binomial(1))
  beta=as.vector(coef(a));
  
  res<-optim(beta,lge0)##L1
  zz<-res$par
  beta=as.vector(zz);
  beta_mle[i,k,jj]= beta;
  
  try(loglike1<--res$value);
  
 
  try(z<-test.new(beta,y,cov,cut));
  try(ps[i,k,jj]<-pnorm(z)); #### the  p value of new test based on N(0,1)
 
  try(c<-test.score(beta,y,cov,cut));
  try(pc[i,k,jj]<-pnorm(c)); 
  
  

  res2<-try(fit2<-mle(lzige0,start =c(0,beta),nobs =NROW(dat),control=list(maxit=1000)));

  try(conv[i,k,jj]<-fit2@details$convergence);  ##find non-convergence samples

  
  if(class(res2)!= "try-error" && (conv[i,k,jj]==0))
  { 
    
    w_hat[i,k,jj]=coef(fit2)[1];  ## estimate the zero-inflated parameter omega
    w=coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1];  ##the  Wald  statistic
    try(pw[i,k,jj]<-pnorm(w));  ##  the p value of Wald test based on N(0,1)
  
    try(loglike2<-logLik(fit2));
   
    try(pl[i,k,jj]<-pchisq(2*(loglike2-loglike1),df =1));
   
    
  }
 
  
  remove(fit2);
  
}
}
}

save.image("./results/TypeI_1.RData");


#########################################################################
#   Section 3.1.3 covariate x ~U(0,2) 
# compare the four tests on controling  type I errors
###########################################################################

#########################################################################
#########################################################################
##runif(n,0,2);
nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000


tp=seq(0.5,2,0.5)
# tp=seq(1.5,3,0.5)
ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(tp)));
beta_mle=w_hat=array(dim=c(nsim, length(ns),length(tp)));
set.seed(2024)
for(jj in 1:length(tp)){k=0;
for( n in ns){k=k+1;
for (i in 1:nsim){
  
  x=runif(n,0,2);
  cov=cbind(1,x);
  betatr=c(tp[jj],1);
  
  
  pro=1/(1+exp(cov %*% betatr))
  
  yy=rgeom(n,prob=pro);                       
  
  rightcensor<-n-(length(yy)*rcpro)
  yy_sorted <- sort(yy)
  cut<- yy_sorted[rightcensor]
  
  
  y<- ifelse(yy> cut, cut, yy)
  dat=data.frame(x=x,y=y);
  
  
  a<-glm(y ~ x,family=negative.binomial(1))
  beta=as.vector(coef(a));
  res<-optim(beta,lge1)##L1似然估计值--
  zz<-res$par
  beta=as.vector(zz);
  beta_mle[i,k,jj]= beta[2];
  
  try(loglike1<--res$value);
  
  
  try(z<-test.new(beta,y,cov,cut));
  try(ps[i,k,jj]<-pnorm(z)); #### the  p value of new test based on N(0,1)
  
  try(c<-test.score(beta,y,cov,cut));
  try(pc[i,k,jj]<-pnorm(c)); 
  
  
  res2<-try(fit2<-mle(lzige1,start=c(0,beta),nobs=NROW(dat),control=list(maxit=1000)));
  
  try(conv[i,k,jj]<-fit2@details$convergence);  ##find non-convergence samples
  
  
  if(class(res2) != "try-error" && (conv[i,k,jj]==0))
  { 
    w_hat[i,k,jj]=coef(fit2)[1];  ## estimate the zero-inflated parameter omega
    w=coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]; ##the  Wald  statistic
    try(pw[i,k,jj]<-pnorm(w));  ## the p value of Wald test based on N(0,1)
    
    try(loglike2<-logLik(fit2));
    ##LR检验
    try(pl[i,k,jj]<-pchisq(2*(loglike2-loglike1),df =1));
  }
  
  remove(fit2);
  
}
}
}

# 
# 
# for(jj in 1:length(tp)){k=0;
# for( n in ns){k=k+1;
# for (i in 1:nsim){
#   
#   x=runif(n,0,2);
#   cov=cbind(1,x);
#   betatr=c(tp[jj],1);
#   
# 
#   pro=1/(1+exp(cov %*% betatr))
#  
#   yy=rgeom(n,prob=pro);                       
# 
#   rightcensor<-n-(length(yy)*rcpro)
#   yy_sorted <- sort(yy)
#   cut<- yy_sorted[rightcensor]
# 
#  
#   y<- ifelse(yy> cut, cut, yy)
#   dat=data.frame(x=x,y=y);
#   
#   
#   a<-glm(y ~ x,family=negative.binomial(1))
#   beta=as.vector(coef(a));
#   res<-optim(beta,lge1)##L1似然估计值--
#   zz<-res$par
#   beta=as.vector(zz);
#   beta_mle[i,k,jj]= beta[2];
#   
#   try(loglike1<--res$value);
#   
#  
#   try(z<-test.new(beta,y,cov,cut));
#   try(ps[i,k,jj]<-pnorm(z)); #### the  p value of new test based on N(0,1)
# 
#   try(c<-test.score(beta,y,cov,cut));
#   try(pc[i,k,jj]<-pnorm(c)); 
#   
#   theta<-c(0.1,beta)
#   res2<-optim(theta,lzig1)
#   res2 <- optim(theta, lzig1, hessian = TRUE)
#   hessian_matrix <- res2$hessian[2:3, 2:3]
# 
#   cov_matrix <- solve(-hessian_matrix)
#   std_errors <- sqrt(diag(cov_matrix))
#   
#   
#   
#   eigen_values <- eigen(hessian_matrix)$values
#   cat("Hessian矩阵的特征值：\n")
#   print(eigen_values)
#   
#   # 2. 验证Hessian矩阵是否负定（特征值全为负）
#   if(all(eigen_values < 0)) {
#     cat("Hessian矩阵是负定的，可以计算协方差矩阵\n")
#   } else {
#     cat("警告：Hessian矩阵非负定，需要修正\n")
#     
#     # 3. 使用特征值分解修正Hessian矩阵（Tikhonov正则化）
#     min_eigen <- min(eigen_values)
#     if(min_eigen >= 0) {
#       # 完全非负定矩阵，添加足够大的正则化项
#       hessian_fixed <- hessian_matrix - diag(rep(min_eigen + 1e-6, length(eigen_values)))
#     } else {
#       # 部分特征值为负，调整到安全范围
#       hessian_fixed <- hessian_matrix - diag(rep(min_eigen - 1e-6, length(eigen_values)))
#     }
#     
#     # 4. 重新计算协方差矩阵
#     cov_matrix <- solve(-hessian_fixed)
#     cat("修正后的协方差矩阵：\n")
#     print(cov_matrix)
#     cov_eigen <- eigen(cov_matrix)$values
#     if(all(cov_eigen > 0))cat("修正后的协方差矩阵是正定的\n")
#     
#     
#     
#   res2 <- optim(theta, lzig1_safe, hessian = TRUE) 
#   
#   res2<-try(fit2<-mle(lzige1_fixed,start=c(0.14,beta),nobs=NROW(dat),control=list(maxit=1000)));
#   
#   res2
#   res2<-try(fit2<-mle(lzige1,start=c(0.1,beta),nobs=NROW(dat),control=list(maxit=1000)));
# 
#   try(conv[i,k,jj]<-res2$convergence);  ##find non-convergence samples
#   
#   
#   if(class(res2) != "try-error" && (conv[i,k,jj]==0))
#   { 
#     w_hat[i,k,jj]=res2$par[1]  coef(fit2)[1];  ## estimate the zero-inflated parameter omega
#     w=coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]; ##the  Wald  statistic
#     try(pw[i,k,jj]<-pnorm(w));  ## the p value of Wald test based on N(0,1)
#     
#     try(loglike2<-logLik(fit2));
#     ##LR检验
#     try(pl[i,k,jj]<-pchisq(2*(loglike2-loglike1),df =1));
#   }
#   
#   remove(fit2);
#   
# }
# }
# }

save.image("./results/TypeI_2.RData");



#########################################################################
#########################################################################
#   Section 3.1.2  covariate x ~N(0,1) 
# compare the four tests on controling  type I errors
###########################################################################
nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000
#tp=seq(-1.5,1.5,1);

tp=seq(0.5,2,0.5)
#tp=seq(1.5,3,0.5)

#tp=seq(1,2.5,0.5)
#tp=c(-2.5,-1.5,-0.5,0.5)
#tp=c(2.5 ,3.0, 3.5, 4.0);
#tp=c(2.5,1.5,-1.5,-2.5);
ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(tp)));
beta_mle=ymin=ymax=w_hat=array(dim=c(nsim, length(ns),length(tp)));
# cut=5
set.seed(2024)
for(jj in 1:length(tp)){k=0;
for( n in ns){k=k+1;
for (i in 1:nsim){
  
  x=rnorm(n,0,1);
  betatr=c(tp[jj],1);
  cov=cbind(1,x);
  pro=1/(1+exp(cov %*% betatr))

  yy=rgeom(n,prob=pro); 
  rightcensor<-n-(length(yy)*rcpro)
  yy_sorted <- sort(yy)
  cut<- yy_sorted[rightcensor]
  y<- ifelse(yy> cut, cut, yy)
  dat=data.frame(x=x,y=y);
  
  a<-glm(y ~ x,family=negative.binomial(1))
  beta=as.vector(coef(a));
  res<-optim(beta,lge1)##L1似然估计值--
  zz<-res$par
  beta=as.vector(zz);
  beta_mle[i,k,jj]= beta[2];
  
  try(loglike1<--res$value);
  
  try(z<-test.new(beta,y,cov,cut));
  try(ps[i,k,jj]<-pnorm(z)); #### the  p value of new test based on N(0,1)
  # 2score统计量
  try(c<-test.score(beta,y,cov,cut));
  try(pc[i,k,jj]<-pnorm(c)); 
  
  res2<-try(fit2<-mle(lzige1,start=c(0,beta),nobs=NROW(dat),control=list(maxit=1000)));
  res2
  
  try(conv[i,k,jj]<-fit2@details$convergence);
  
  if(class(res2) != "try-error" && (conv[i,k,jj]==0))
  { 
    w_hat[i,k,jj]=coef(fit2)[1];  ## estimate the zero-inflated parameter omega
    w=coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]; ##the  Wald  statistic
    try(pw[i,k,jj]<-pnorm(w));  ## the p value of  Wald test
    
    try(loglike2<-logLik(fit2));
    #print(loglike2);
    
    try(pl[i,k,jj]<-pchisq(2*(loglike2-loglike1),df =1));
    
  }
  remove(fit2);
  
  
}
}
}


save.image("./results/TypeI_3.RData")
########################################################



##########################################################
###   Section 3.2  examine the powers based on ZICGE Responses 
##########################################################

#########################################################################
##  Section 3.2.1 No covariate for the ZICG  data
# compare the four tests on powers

nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000

tp=seq(0.5,2,0.5)
#tp=seq(1.5,3,0.5)

lp=0.01*(seq(5, 30, 5))  ### the probability of structural zeros

ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(lp),length(tp)));
beta_mle=omega_hat=omega_se=array(dim=c(nsim, length(ns),length(lp),length(tp)));

set.seed(2024)
for(tt in 1:length(tp)){
  for(jj in 1:length(lp)){k=0;
  for( n in ns){k=k+1;
  for (i in 1:nsim){
  
    cov<-matrix(rep(1,n),ncol=1)
    pro=1/(1+exp(tp[tt]));
    yy<-rzigeom(n,pro=pro,pstr0=lp[jj])
    rightcensor<-n-(length(yy)*rcpro)
    yy_sorted <- sort(yy)
    cut<- yy_sorted[rightcensor]
    y<- ifelse(yy> cut, cut, yy)
    dat=data.frame(y=y);
    
    
    a<-glm(y ~ 1,family=negative.binomial(1))
    beta=as.vector(coef(a));
    res<-optim(beta,lge0)##L1似然估计值--lge0删失几何回归模型
    zz<-res$par
    beta=as.vector(zz);
    beta_mle[i,k,jj,tt]=beta;
    try(loglike1<--res$value);
    
    try(ps[i,k,jj,tt]<-test.new(beta,y,cov,cut));  # the new statistic; 
    try(pc[i,k,jj,tt]<-test.score(beta,y,cov,cut)); #the score statistic; 
    
    ####fitting the data by zib model
    res2<-try(fit2<-mle(lzige0,start =c(0.1,beta),nobs=NROW(dat),control=list(maxit=1000)));
    
    try(conv[i,k,jj,tt]<-fit2@details$convergence); 
    if(class(res2) != "try-error"&& (conv[i,k,jj,tt]==0))
    { 
      omega_hat[i,k,jj,tt]=coef(fit2)[1];  ## estimate the zero-inflated parameter omega
      omega_se[i,k,jj,tt]=sqrt(diag(vcov(fit2)))[1];
      try(pw[i,k,jj,tt]<-coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]);  ## Wald statistic
      
      #try(loglike2<-logLik(fit2));
      try(pl[i,k,jj,tt]<-2*(logLik(fit2)-loglike1));### LR statistic
      
    }
    remove(fit2);
    
  }
  }
  } 
}
save.image("./results/power_1.RData"); 

#########################################################################
##  Section 3.2.3 (1)covariate x~U(0,2)  for the ZIB  data
# compare the four tests on powers

nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000
#tp=seq(1.5,3,0.5)
tp=seq(0.5,2,0.5)
lp=0.01*(seq(5, 30, 5))  ### the probability of structural zeros

ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(lp),length(tp)));
beta_mle=omega_hat=omega_se=array(dim=c(nsim, length(ns),length(lp),length(tp)));

set.seed(2024)
for(tt in 1:length(tp)){
  for(jj in 1:length(lp)){k=0;
  for( n in ns){k=k+1;
  for (i in 1:nsim){
    
    x=runif(n,0,2);
    cov=cbind(1,x);
    betatr=c(tp[tt],1);
    pro=1/(1+exp(cov %*% betatr))
    yy<-rzigeom(n,pro=pro,pstr0=lp[jj])
    ########################################################
    rightcensor<-n-(length(yy)*rcpro)
    yy_sorted <- sort(yy)
    cut<- yy_sorted[rightcensor]
    y<- ifelse(yy> cut, cut, yy)
    dat=data.frame(x=x,y=y);
    a<-glm(y ~ x,family=negative.binomial(1))
    beta=as.vector(coef(a));
    res<-optim(beta,lge1)
    
    zz<-res$par
    beta=as.vector(zz);
    
    beta_mle[i,k,jj,tt]=beta[2];
    
    try(loglike1<--res$value);
    try(loglike1<--res$value);
    try(ps[i,k,jj,tt]<-test.new(beta,y,cov,cut));  # the new statistic; 
    
    try(pc[i,k,jj,tt]<-test.score(beta,y,cov,cut)); # the score statistic; 
    res2<-try(fit2<-mle(lzige1,start=c(0,beta),nobs=NROW(dat),control=list(maxit=1000)));
  
    try(conv[i,k,jj,tt]<-fit2@details$convergence);
    if(class(res2) != "try-error" && (conv[i,k,jj,tt]==0))
    { 
      omega_hat[i,k,jj,tt]=coef(fit2)[1]; ## estimate the zero-inflated parameter omega
      omega_se[i,k,jj,tt]=sqrt(diag(vcov(fit2)))[1];
      try(pw[i,k,jj,tt]<-coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]);  ## Wald statistic
      
      #try(loglike2<-logLik(fit2)); 
      try(pl[i,k,jj,tt]<-2*(logLik(fit2)-loglike1));### LR statistic
      
    }
    remove(fit2);
    
  }
  }
  }
}

save.image("./results/power_2.RData");



#########################################################################
##  Section 3.2.3 (2)covariate x~U(0,2)  for the ZIB  data
# compare the four tests on powers
nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000

tp=seq(0.5,2,0.5)
#tp=seq(1.5,3,0.5)
lp=seq(-3,-0.5, 0.5)
ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(lp),length(tp)));
beta_mle=omega_hat=omega_se=array(dim=c(nsim, length(ns),length(lp),length(tp)));

set.seed(2024)
for(tt in 1:length(tp)){
  for(jj in 1:length(lp)){k=0;
  for( n in ns){k=k+1;
  for (i in 1:nsim){
    
    x=runif(n,0,2);
    cov=cbind(1,x);
    betatr=c(tp[tt],1);
    
    
    pro=1/(1+exp(cov %*% betatr))
    omega=plogis(lp[jj]+x);  ####the probabilty of structutal zeros
    yy<-rzigeom(n,pro=pro,pstr0=omega)
    rightcensor<-n-(length(yy)*rcpro)
    yy_sorted <- sort(yy)
    cut<- yy_sorted[rightcensor]
    y<- ifelse(yy> cut, cut, yy)
    dat=data.frame(x=x,y=y);
    
    a<-glm(y ~ x,family=negative.binomial(1))
    beta=as.vector(coef(a));
    res<-optim(beta,lge1)
    
    zz<-res$par
    beta=as.vector(zz);
    
    beta_mle[i,k,jj,tt]=beta[2];
    
    
    try(loglike1<--res$value);
    
    try(ps[i,k,jj,tt]<-test.new(beta,y,cov,cut)); #the new test; 
    
    try(pc[i,k,jj,tt]<-test.score(beta,y,cov,cut));# the  score test; 
    
    
    res2<-try(fit2<-mle(lzige1,start=c(0.1,beta),nobs=NROW(dat),control=list(maxit=1000)));

    try(conv[i,k,jj,tt]<-fit2@details$convergence);
    if(class(res2) != "try-error" && (conv[i,k,jj,tt]==0))
    { 
      omega_hat[i,k,jj,tt]=coef(fit2)[1]; ## estimate the zero-inflated parameter omega
      omega_se[i,k,jj,tt]=sqrt(diag(vcov(fit2)))[1];
      
      try(pw[i,k,jj,tt]<-coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]);  ## Wald statistic
      
      try(pl[i,k,jj,tt]<-2*(logLik(fit2)-loglike1));### LR statistic
      
    }
    remove(fit2);
    
  }
  }
  }
}

save.image("./results/power_3.RData");




########################################################################
#  Section 3.2.2  (1)covariate x~N(0,1)  
#compare the four tests on powers
nsim=1000;  # 1000 Monte Carlo replications
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000
#tp=seq(1.5,3,0.5)

tp=seq(0.5,2,0.5)
lp=0.01*(seq(5, 30, 5))  ### the probability of structural zeros


ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(lp),length(tp)));
beta_mle=omega_hat=omega_se=array(dim=c(nsim, length(ns),length(lp),length(tp)));

set.seed(2024)
for(tt in 1:length(tp)){
  for(jj in 1:length(lp)){k=0;
  for( n in ns){k=k+1;
  for (i in 1:nsim){

    x=rnorm(n,0,1);
    cov=cbind(1,x);
   
    betatr=c(tp[tt],1);
   
    pro=1/(1+exp(cov %*% betatr))
    
    
    yy<-rzigeom(n,pro=pro,pstr0=lp[jj])
    rightcensor<-n-(length(yy)*rcpro)
    yy_sorted <- sort(yy)
    cut<- yy_sorted[rightcensor]
    y<- ifelse(yy> cut, cut, yy)
    dat=data.frame(x=x,y=y);


    a<-glm(y ~ x,family=negative.binomial(1))
    beta=as.vector(coef(a));
    res<-optim(beta,lge1)##L1似然估计值--
    zz<-res$par
    beta=as.vector(zz);

    beta_mle[i,k,jj,tt]=beta[2];
    try(loglike1<--res$value);

    try(ps[i,k,jj,tt]<-test.new(beta,y,cov,cut));  # the new statistic;

    try(pc[i,k,jj,tt]<-test.score(beta,y,cov,cut)); # the score statistic;

    res2<-try(fit2<-mle(lzige1,start=c(0.1,beta),nobs=NROW(dat),control=list(maxit=1000)));
    
    
    try(conv[i,k,jj,tt]<-fit2@details$convergence);
    if(class(res2) != "try-error" && (conv[i,k,jj,tt]==0))
    {
      omega_hat[i,k,jj,tt]=coef(fit2)[1]; ## estimate the zero-inflated parameter omega
      omega_se[i,k,jj,tt]=sqrt(diag(vcov(fit2)))[1];
      try(pw[i,k,jj,tt]<-coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]);  ## Wald statistic

     

      try(pl[i,k,jj,tt]<-2*(logLik(fit2)-loglike1));### LR statistic

    }
    remove(fit2);

  }
  }
  }
}

save.image("./results/power_4.RData");



##################################################################################

########################################################################
#  Section 3.2.2 (2)covariate x~N(0,1)  
nsim=1000;  # 1000 Monte Carlo replications   
ns=c(50,100,200,500,1000); #sample size 50,100,200,500,1000

tp=seq(0.5,2,0.5)
#tp=seq(1.5,3,0.5)
lp=seq(-3,-0.5, 0.5)

ps=pc=pl=pw=conv=array(dim=c(nsim, length(ns),length(lp),length(tp)));
beta_mle=omega_hat=omega_se=array(dim=c(nsim, length(ns),length(lp),length(tp)));

set.seed(2024)
for(tt in 1:length(tp)){
  for(jj in 1:length(lp)){k=0;
  for( n in ns){k=k+1;
  for (i in 1:nsim){
    
    x=rnorm(n,0,1);
    cov=cbind(1,x);
    betatr=c(tp[tt],1);
    pro=1/(1+exp(cov %*% betatr))
    omega=plogis(lp[jj]+x); 
    yy<-rzigeom(n,pro=pro,pstr0=omega)
    rightcensor<-n-(length(yy)*rcpro)
    yy_sorted <- sort(yy)
    cut<- yy_sorted[rightcensor]
    y<- ifelse(yy> cut, cut, yy)
    dat=data.frame(x=x,y=y);
    
  
    a<-glm(y ~ x,family=negative.binomial(1))
    beta=as.vector(coef(a));
    res<-optim(beta,lge1)##L1
    zz<-res$par
    beta=as.vector(zz);
 
    beta_mle[i,k,jj,tt]= beta[2];

    try(loglike1<--res$value);
    try(ps[i,k,jj,tt]<-test.new(beta,y,cov,cut)); #the new test; 
    
    try(pc[i,k,jj,tt]<-test.score(beta,y,cov,cut));# the  score test; 
    
    res2<-try(fit2<-mle(lzige1,start=c(0,beta),nobs=NROW(dat),control=list(maxit=1000)));
   
    try(loglike1<--res$value);
    
    try(conv[i,k,jj,tt]<-fit2@details$convergence);
    if(class(res2) != "try-error" && (conv[i,k,jj,tt]==0)){ 
      omega_hat[i,k,jj,tt]=coef(fit2)[1]; ## estimate the zero-inflated parameter omega
      omega_se[i,k,jj,tt]=sqrt(diag(vcov(fit2)))[1];
      
      try(pw[i,k,jj,tt]<-coef(fit2)[1]/sqrt(diag(vcov(fit2)))[1]);  ## Wald statistic
      
      try(pl[i,k,jj,tt]<-2*(logLik(fit2)-loglike1));### LR statistic
      
    }
    remove(fit2);
    
  }
  }
  }
}

save.image("./results/power_5.RData")








