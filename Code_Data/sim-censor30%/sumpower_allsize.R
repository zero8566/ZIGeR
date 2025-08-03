# setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/提交code/Code_and_Data/sim_xiu/0.3censored/")
setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/results/30%删失")

# setwd("D:/0Adas/codes/results")
#setwd("D:/0Adas/01零膨胀/zero-test-weiyahui/bio-/wodeCode/test_result/results/")
load("./results/power_1.RData");

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

pdf('power_const_m0.pdf')

op<-par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200,500, 1000);

gamma=tp
#c(0.5,1,1.5,2)
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
for( kk in 2:4){
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
########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .03,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-0.565,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T) 

par(op);##fig4

dev.off()













load("./results/power_2.RData")

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
for( kk in 2:4){
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
########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .03,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-0.60,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T); 

par(op);
dev.off()####fig5



















########################################################### 

#load("./results/power_3.RData");

load("./results/power_3(lp=seq(-3,-0.5, 0.5)).RData");
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
pdf('power_unifboth_m0.pdf')

op <- par(mfrow = c(5,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1) 

###############################################3
# first row
# k for different gamma
ns=c(50, 100, 200,500, 1000);
#gamma=seq(-3,-1.5,0.5);  ###4 values
#gamma=c(2.5 ,3.0, 3.5, 4.0)
gamma=tp
  #c(0.5,1,1.5,2) ##
lp=lp;  ### 6 values
#lp=log(1+exp(1+lp))-log(1+exp(lp)); ### the proportion of structural zeros
lp=plogis(lp+0.5)### the proportion of structural zeros

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
  text( 0.3, .03, paste("n = ", ns[kk]));
}

########################
## 2-4 row
#######################
for( kk in 2:4){
  for (k in 1:4)  #k for different gamma
  {if(k==1){
    plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="", xaxt='n',ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",xlab="", xaxt='n',yaxt='n',ylab="", ylim=c(0,1))}
    
    lines(lp,pc1[kk, ,k],lty=2,col="red") 
    lines(lp,pw1[kk, ,k],lty=3, col="blue") 
    lines(lp,pl1[kk, ,k],lty=4, col="green") 
    text( 0.3, .03, paste("n = ", ns[kk]));
  }
}
########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.3, .03,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-1,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T);

par(op);

dev.off()##fig6


















################################################
##################################################


load("./results/power_4.RData")
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
# gamma=c(0.5,1,1.5,2)
gamma=tp
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
for( kk in 2:4){
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
########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .03,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
legend(-0.595,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T) 

par(op);##fig7

dev.off()


















load("./results/power_5.RData")


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
pdf('power_normalboth_m0.pdf')

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
for( kk in 2:4){
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
########################
## 5th row
#######################

kk=5
for (k in 1:4)  #k for different gamma
{
  if(k==1){plot(lp,ps1[kk, ,k] ,type="l", col="black",xlab="" ,ylab="", ylim=c(0,1))
  }else{plot(lp,ps1[kk, ,k] ,type="l",yaxt='n',xlab="" ,ylab="", ylim=c(0,1))
  }
  lines(lp,pc1[kk, ,k],lty=2,col="red") 
  lines(lp,pw1[kk, ,k],lty=3, col="blue") 
  lines(lp,pl1[kk, ,k],lty=4, col="green")
  text( 0.25, .03,  paste("n = ", ns[kk]));
}
title(xlab = "Probabilities of Structural Zeros", ylab = "Power", outer = TRUE, line = 3)
par(xpd=NA)
# legend(-0.720,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
#        col=c("black","red","blue","green"),  horiz = T);
legend(-0.8,-.60,legend=c("New Test","Score Test","Wald Test","LR Test"), lty = c(1,2,3,4), 
       col=c("black","red","blue","green"),  horiz = T);


par(op);##fig8

dev.off()



#




