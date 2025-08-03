##########################################################
###   Section 3.2  examine the powers based on ZIB Responses 
##########################################################

# Plot of power for Section 3.2.1 No covariate for the ZIB  data
### constant trial size m0 
# sample size 50, 100, 200,500, 1000; 
# tp=seq(-2.5,-1,0.5);  
# lp= 0.01*(seq(5, 30, 5));
##[i,k,jj,tt]=[nsim, ns, lp, tp]??
setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/results/")



# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/results")

#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")
load("../results/power_1.RData");
# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/")


# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第2次修改/result2")
## check convergency;
apply(conv,2:4, sum); 
###########the zero-inflated parameter
omegamean=apply(omega_hat, 2:4,mean,na.rm=T);## 
omegasd=apply(omega_hat,2:4,sd,na.rm=T); ###
omegase=apply(omega_se, 2:4,mean,na.rm=T);  ###


##########calculate non-convergence samples
b1=apply(is.na(pl),2:4, sum);###the LR test
b2=apply(is.na(pw),2:4, sum);  ### the Wald test
b3=apply(is.na(pc),2:4, sum);  ### the score test
b4=apply(is.na(ps),2:4, sum);  ### the new test

######## empirical power of all test for two-sided test
## WebFigure S4
pl2=apply(pchisq(pl,1)>.95, 2:4,mean,na.rm=T);  ###the LR test
pw2=apply(pnorm(abs(pw))>.975, 2:4,mean,na.rm=T);  ### the Wald test
pc2=apply(pnorm(abs(pc))>.975, 2:4,mean,na.rm=T);  ### the score test
ps2=apply(pnorm(abs(ps))>.975, 2:4,mean,na.rm=T);  ### the New test

######## empirical power of all test for one-sided test
## Figure 4
pl1<- apply((pchisq(pl,1)>.90)&(omega_hat>0), 2:4, mean,na.rm=T);  #the LR test
pw1<- apply(pnorm(pw)>.95, 2:4, mean,na.rm=T); ### the Wald test
pc1<- apply(pnorm(pc)>.95, 2:4, mean,na.rm=T);###the score test
ps1<- apply(pnorm(ps)>.95, 2:4, mean,na.rm=T);### the new test


##############
## Figure 4
pdf('power_const_m0.pdf')

op<-par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200,500, 1000);
#gamma=seq(-2.5,-1,0.5);  ###4 values
#gamma=c(-1, 0.5, 1.5,2.5)##均值为正的结果power好一些
#gamma=c(-0.5, 0.5, 1.5,2.5);
#gamma=seq(-1.5,1.5,1)
#gamma=seq(-0.5,2.5,1)

gamma=c(0.5,1,1.5,2)
lp=0.01*(seq(5, 30, 5));  ### 6 values
###########################################################


kk=1
for (k in 1:4){ 
  txt= bquote(gamma == .(gamma[k]))
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black", xaxt='n',xlab="", ylab="", ylim=c(0,1),
                main = txt)}else{plot(lp,ps1[kk, ,k] ,type="l",xaxt='n',yaxt='n',xlab="", ylab="", ylim=c(0,1),
                                      main = txt
                )}
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green") 
  text( 0.25, .03, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:3){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.25, .03, paste("n = ", ns[kk]));
  }
}



for( kk in 4){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0.3,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0.3,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.25, .35, paste("n = ", ns[kk]));
  }
}

########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0.35,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0.35,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .4,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-0.6,-0.05,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T) 

par(op);##fig4

dev.off()


######## Plot of power for Section 3.2.2  covariate x~U(0,1)  for the ZIB  data

###No covariate for probability of structural zeros
### constant trial size m0=10 
# sample size 50, 100, 200,500, 1000; 
# tp=c(-3.5,-2.5, -2, -1.5)###4 values
# lp= 0.01*(seq(5, 30, 5))### 6 values
##[i,k,jj,tt]=[nsim, ns, lp, tp]??

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/results")

#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")
load("../results/power_2.RData");

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/")

# 
# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第2次修改/result2")


## check convergency;
apply(conv,2:4, sum); 
###########the zero-inflated parameter
omegamean=apply(omega_hat, 2:4,mean,na.rm=T);## 
omegasd=apply(omega_hat,2:4,sd,na.rm=T); ###
omegase=apply(omega_se, 2:4,mean,na.rm=T);  ###

##########calculate non-convergence samples
b1=apply(is.na(pl),2:4, sum);###the LR test
b2=apply(is.na(pw),2:4, sum);  ### the Wald test
b3=apply(is.na(pc),2:4, sum);  ### the score test
b4=apply(is.na(ps),2:4, sum);  ### the new test

######## empirical power of all test for two-sided test with gamma0=0.05
## WebFigure S7
pl2=apply(pchisq(pl,1)>.95, 2:4,mean,na.rm=T);  ###the LR test
pw2=apply(pnorm(abs(pw))>.975, 2:4,mean,na.rm=T);  ### the Wald test
pc2=apply(pnorm(abs(pc))>.975, 2:4,mean,na.rm=T);  ### the score test
ps2=apply(pnorm(abs(ps))>.975, 2:4,mean,na.rm=T);  ### the New test

######## empirical power of all test for one-sided test with gamma0=0.05
## Figure 5
pl1<- apply((pchisq(pl,1)>.90)&(omega_hat>0), 2:4, mean,na.rm=T);  #the LR test
pw1<- apply(pnorm(pw)>.95, 2:4, mean,na.rm=T); ### the Wald test
pc1<- apply(pnorm(pc)>.95, 2:4, mean,na.rm=T);###the score test
ps1<- apply(pnorm(ps)>.95, 2:4, mean,na.rm=T);### the new test

##############
## Figure 5
pdf('power_unif_m02.pdf')

op <- par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200, 500,1000);
#gamma=seq(-3,-1.5,0.5);  ###4 values
gamma=tp
#gamma=seq(-0.5,2.5,1)
#gamma=c(-0.5,0.5,1.5,2.5); 
lp=lp;  ### 6 values

kk=1
for (k in 1:4){ 
  txt= bquote(gamma == .(gamma[k]))
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black", xaxt='n',xlab="", ylab="", ylim=c(0,1),
                main = txt)}else{plot(lp,ps1[kk, ,k] ,type="l",xaxt='n',yaxt='n',xlab="", ylab="", ylim=c(0,1),
                                      main = txt
                )}
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green") 
  text( 0.25, .03, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:3){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.25, .03, paste("n = ", ns[kk]));
  }
}

kk=4
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,xaxt='n',ylab="", ylim=c(0.5,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,xaxt='n',ylab="", ylim=c(0.5,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .55,  paste("n = ", ns[kk]));
}



########################
## 5th row
#######################

# 
kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0.9,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0.9,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25,.91,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-0.63,.835,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T); 

par(op);
dev.off()####fig5


###########################################
###Covariate for probability of structural zeros
### constant trial size m0=10 
# sample size 50, 100, 200,500, 1000
# tp=c(-3.5,-2.5, -2, -1.5)
# the probability of structural zero: logit(w)=lp+x
# lp=seq(-3,-0.5, 0.5);  
# the probabilities of structural zero are 0.08, 0.12, 0.19, 0.27, 0.38 and 0.50
###[i,k,jj,tt]=[nsim, ns, lp, tp]??
########################################################### 


# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/results")

#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")

load("../results/power_3.RData");



## check convergency;
apply(conv,2:4, sum); 
###########the zero-inflated parameter
omegamean=apply(omega_hat, 2:4,mean,na.rm=T);## 
omegasd=apply(omega_hat,2:4,sd,na.rm=T); ###
omegase=apply(omega_se, 2:4,mean,na.rm=T);  ###

##########calculate non-convergence samples
b1=apply(is.na(pl),2:4, sum);###the LR test
b2=apply(is.na(pw),2:4, sum);  ### the Wald test
b3=apply(is.na(pc),2:4, sum);  ### the score test
b4=apply(is.na(ps),2:4, sum);  ### the new test

######## empirical power of all test for two-sided test with gamma0=0.05
## WebFigure S10
pl2=apply(pchisq(pl,1)>.95, 2:4,mean,na.rm=T);  ###the LR test
pw2=apply(pnorm(abs(pw))>.975, 2:4,mean,na.rm=T);  ### the Wald test
pc2=apply(pnorm(abs(pc))>.975, 2:4,mean,na.rm=T);  ### the score test
ps2=apply(pnorm(abs(ps))>.975, 2:4,mean,na.rm=T);  ### the New test

######## empirical power of all test for one-sided test with gamma0=0.05
## Figure 6
pl1<- apply((pchisq(pl,1)>.90)&(omega_hat>0), 2:4, mean,na.rm=T);  #the LR test
pw1<- apply(pnorm(pw)>.95, 2:4, mean,na.rm=T); ### the Wald test
pc1<- apply(pnorm(pc)>.95, 2:4, mean,na.rm=T);###the score test
ps1<- apply(pnorm(ps)>.95, 2:4, mean,na.rm=T);### the new test


##############
## Figure 6
pdf('power_unifboth_m01(lp=seq(-3,-0.5, 0.5)).pdf')

op <- par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200,500, 1000);
#gamma=seq(-3,-1.5,0.5);  ###4 values
#gamma=c(2.5 ,3.0, 3.5, 4.0)
gamma=c(0.5,1,1.5,2) ##
lp=lp;  ### 6 values
#lp=log(1+exp(1+lp))-log(1+exp(lp)); ### the proportion of structural zeros
lp=plogis(lp+0.5)### the proportion of structural zeros

kk=1
for (k in 1:4){ 
  txt= bquote(gamma == .(gamma[k]))
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black", xaxt='n',xlab="", ylab="", ylim=c(0.5,1),
                main = txt)}else{plot(lp,ps1[kk, ,k] ,type="l",xaxt='n',yaxt='n',xlab="", ylab="", ylim=c(0.5,1),
                                      main = txt
                )}
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green") 
  text( 0.3, .53, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:3){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0.8,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0.8,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.3, .83, paste("n = ", ns[kk]));
  }
}

kk=4
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,xaxt='n',ylab="", ylim=c(0.8,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,xaxt='n',ylab="", ylim=c(0.8,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.3, .82,  paste("n = ", ns[kk]));
}


kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0.8,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0.8,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.31, .82,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-1.08,.67,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T);

par(op);

dev.off()##fig6






################################################
##################################################

# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/results")

#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")
load("../results/power_4.RData");
# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/")

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第2次修改/result2")

# 
# 
# 
# load("../results/power_4.RData")
## check convergency;
apply(conv,2:4, sum); 
###########the zero-inflated parameter
omegamean=apply(omega_hat, 2:4,mean,na.rm=T);## 
omegasd=apply(omega_hat,2:4,sd,na.rm=T); ###
omegase=apply(omega_se, 2:4,mean,na.rm=T);  ###

##########calculate non-convergence samples
b1=apply(is.na(pl),2:4, sum);###the LR test
b2=apply(is.na(pw),2:4, sum);  ### the Wald test
b3=apply(is.na(pc),2:4, sum);  ### the score test
b4=apply(is.na(ps),2:4, sum);  ### the new test

######## empirical power of all test for two-sided test with gamma0=0.05
## WebFigure S13
pl2=apply(pchisq(pl,1)>.95, 2:4,mean,na.rm=T);  ###the LR test
pw2=apply(pnorm(abs(pw))>.975, 2:4,mean,na.rm=T);  ### the Wald test
pc2=apply(pnorm(abs(pc))>.975, 2:4,mean,na.rm=T);  ### the score test
ps2=apply(pnorm(abs(ps))>.975, 2:4,mean,na.rm=T);  ### the new test

######## empirical power of all test for one-sided test with gamma0=0.05
## Figure 7
pl1<- apply((pchisq(pl,1)>.90)&(omega_hat>0), 2:4, mean,na.rm=T);  #the LR test
pw1<- apply(pnorm(pw)>.95, 2:4, mean,na.rm=T); ### the Wald test
pc1<- apply(pnorm(pc)>.95, 2:4, mean,na.rm=T);###the score test
ps1<- apply(pnorm(ps)>.95, 2:4, mean,na.rm=T);### the new test



##############
## Figure 7
pdf('power_normal_m04.pdf')

op <- par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200 ,500,1000);

#gamma=seq(-0.5,2.5,1)
gamma=c(0.5,1,1.5,2)
#gamma=seq(-2.5,-1,0.5);  ###4 values
lp=0.01*(seq(5, 30, 5));  ### 6 values

kk=1
for (k in 1:4){ 
  txt= bquote(gamma == .(gamma[k]))
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black", xaxt='n',xlab="", ylab="", ylim=c(0,1),
                main = txt)}else{plot(lp,ps1[kk, ,k] ,type="l",xaxt='n',yaxt='n',xlab="", ylab="", ylim=c(0,1),
                                      main = txt
                )}
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green") 
  text( 0.25, .03, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:3){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.25, .03, paste("n = ", ns[kk]));
  }
}


kk=4
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,xaxt='n',ylab="", ylim=c(0.4,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,xaxt='n',ylab="", ylim=c(0.4,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .45,  paste("n = ", ns[kk]));
}


########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0.5,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0.5,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .53,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line =3)
par(xpd=NA)
legend(-0.595,.2,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T) 

par(op);##fig7

dev.off()

##################################################
###No covariate for probability of structural zeros
###varying trial size mi
# sample size 50, 100, 200,500, 1000; 
# tp=seq(-2.5,-1,0.5) ### 4 values
# lp= 0.01*(seq(5, 30, 5)) ### 6 values
##[i,k,jj,tt]=[nsim, ns, lp, tp]??


###########################################
###Covariate for probability of structural zeros
### constant trial size m0
# sample size 50, 100, 200,500, 1000
# tp=c(-2.5, -2, -1.5,-1)### 4 values
# the probability of structural zero: logit(w)=lp+x
# lp=seq(-3,-0.5, 0.5);  
# the probabilities of structural zeros  are:
# 0.07, 0.11, 0.16, 0.22, 0.30 and 0.40 
###[i,k,jj,tt]=[nsim, ns, lp, tp]??

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/results")

#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")
load("../results/power_5.RData");

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/")

#setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第2次修改/result2")

# 
# 
# 
# load("../results/power_5(lp=seq(-3,-0.5, 0.5)).RData")


## check convergency;
apply(conv,2:4, sum); 
###########the zero-inflated parameter
omegamean=apply(omega_hat, 2:4,mean,na.rm=T);## 
omegasd=apply(omega_hat,2:4,sd,na.rm=T); ###
omegase=apply(omega_se, 2:4,mean,na.rm=T);  ###

##########calculate non-convergence samples
b1=apply(is.na(pl),2:4, sum);###the LR test
b2=apply(is.na(pw),2:4, sum);  ### the Wald test
b3=apply(is.na(pc),2:4, sum);  ### the score test
b4=apply(is.na(ps),2:4, sum);  ### the new test

######## empirical power of all test for two-sided test with gamma0=0.05
## WebFigure S16
pl2=apply(pchisq(pl,1)>.95, 2:4,mean,na.rm=T);  ###the LR test
pw2=apply(pnorm(abs(pw))>.975, 2:4,mean,na.rm=T);  ### the Wald test
pc2=apply(pnorm(abs(pc))>.975, 2:4,mean,na.rm=T);  ### the score test
ps2=apply(pnorm(abs(ps))>.975, 2:4,mean,na.rm=T);  ### the New test

######## empirical power of all test for one-sided test with gamma0=0.05
## Figure 8
pl1<- apply((pchisq(pl,1)>.90)&(omega_hat>0), 2:4, mean,na.rm=T);  #the LR test
pw1<- apply(pnorm(pw)>.95, 2:4, mean,na.rm=T); ### the Wald test
pc1<- apply(pnorm(pc)>.95, 2:4, mean,na.rm=T);###the score test
ps1<- apply(pnorm(ps)>.95, 2:4, mean,na.rm=T);### the new test

##############
## Figure 8
pdf('power_normalboth_m0lp=seq(-3,-0.5, 0.5).pdf')

op <- par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200,500, 1000);

#tp=seq(-2,-0.5, 0.5)lp=seq(-4.5,0.5, 1)


gamma=tp
#gamma=seq(-0.5,2.5,1)
#gamma=c(0.5,1.5,2.5,3.5); ###4 values
lp=lp; ### 6 values
lp=plogis(lp+0)


kk=1
for (k in 1:4){ 
  txt= bquote(gamma == .(gamma[k]))
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black", xaxt='n',xlab="", ylab="", ylim=c(0,1),
                main = txt)}else{plot(lp,ps1[kk, ,k] ,type="l",xaxt='n',yaxt='n',xlab="", ylab="", ylim=c(0,1),
                                      main = txt
                )}
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green") 
  text( 0.25, .03, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:3){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.25, .03, paste("n = ", ns[kk]));
  }
}


kk=4
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,xaxt='n',ylab="", ylim=c(0.8,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,xaxt='n',ylab="", ylim=c(0.8,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red")
  lines(lp,pw1[kk, ,k],lty=3, col="blue")
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.26, .82,  paste("n = ", ns[kk]));
}


########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0.8,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0.8,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.26, .82,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
# legend(-0.720,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
#        col=c("black","red","blue","green"),  horiz = T);
legend(-0.8,.67,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T);


par(op);##fig8

dev.off()



#




