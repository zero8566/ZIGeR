
##########################################################
###   Section 3.1  
##########################################################
setwd("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/第6次修改/1-code/0-提交纯净版/Code/sim-censor10%/results/")


load("../results/TypeI_1.RData");
library(xtable)

apply(conv,2:3, sum); 

##########calculate non-convergence samples
b1=apply(is.na(pl),2:3, sum);###the LR test
b2=apply(is.na(pw),2:3, sum);  ### the Wald test
b3=apply(is.na(pc),2:3, sum);  ### the score test
b4=apply(is.na(ps),2:3, sum);  ### the new test
# print(ymin[which(is.na(pl))]);
# print(ymax[which(is.na(pl))]);

######## calculate the rejection rate for zero-inflation for all tests with gamma0=0.05
a1=apply((pl>0.90)&(w_hat>0), 2:3,mean,na.rm=T)*100;  ###the LR test
a2=apply(pw>0.95, 2:3,mean,na.rm=T)*100;  ### the Wald test
a3=apply(pc>0.95, 2:3,mean,na.rm=T)*100;  ### the score test
a4=apply(ps>0.95, 2:3,mean,na.rm=T)*100;  ### the new test
pow=matrix(nrow=5,ncol=16)
pow[,(1:4)*4-3]=a1;
pow[,(1:4)*4-2]=a2;
pow[,(1:4)*4-1]=a3;
pow[,(1:4)*4]=a4;


print(xtable(cbind(c(50,100,200,500,1000),pow)))

##### the p value of LR test due to chi-square#####
pl=1-pl;

###########
## Figure 1
pdf("QQ_const_m0.pdf")

op <- par(mfrow = c(4,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
############
# first row
# k for the different means
##############
k=1;kk=5;

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=T,main  = "LR Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Wald Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Score Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=T, main  = "New Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
############
# second row
##############
k=2;kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

############
# third row
##############
k=3; kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));
############
# fourth row
##############
k=4;kk=5

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='s', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='s', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));



title(xlab = "Empirical  p-values",
      ylab = "Theoretical p-values",
      outer = TRUE, line = 3);
par(xpd=NA)
legend(-2.6,-0.47,legend=c("50","100","200","500","1000"), lty = c(2,3,4,5,1), col=c(2,3,4,5,1), horiz = TRUE);

par(op)##fig1

dev.off()







###############################################################################

load("../results/TypeI_2.RData");
apply(conv,2:3, sum); 

##########calculate non-convergence samples
b1=apply(is.na(pl),2:3, sum);###the LR test
b2=apply(is.na(pw),2:3, sum);  ### the Wald test
b3=apply(is.na(pc),2:3, sum);  ### the score test
b4=apply(is.na(ps),2:3, sum);  ### the new test
# print(ymin[which(is.na(pl))]);
# print(ymax[which(is.na(pl))]);

######## calculate the rejection rate for zero-inflation for all tests with gamma0=0.05
a1=apply((pl>0.90)&(w_hat>0), 2:3,mean,na.rm=T)*100;  ###the LR test
a2=apply(pw>0.95, 2:3,mean,na.rm=T)*100;  ### the Wald test
a3=apply(pc>0.95, 2:3,mean,na.rm=T)*100;  ### the score test
a4=apply(ps>0.95, 2:3,mean,na.rm=T)*100;  ### the new test
pow=matrix(nrow=5,ncol=16)
pow[,(1:4)*4-3]=a1;
pow[,(1:4)*4-2]=a2;
pow[,(1:4)*4-1]=a3;
pow[,(1:4)*4]=a4;
print(xtable(cbind(c(50,100,200,500,1000),pow)))

##### the p value of LR test due to chi-square#####
pl=1-pl;
##############
## Figure 2
pdf("QQ_unif_m0.pdf")

op <- par(mfrow = c(4,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
############
# first row
# k for the different means
##############
k=1;kk=5;

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=T,main  = "LR Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Wald Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Score Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=T, main  = "New Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
############
# second row
##############
k=2;kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

############
# third row
##############
k=3; kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));
############
# fourth row
##############
k=4;kk=5

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='s', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
#text( 0.8, .07, paste("a = ", tp[k]));
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='s', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));



title(xlab = "Empirical  p-values",
      ylab = "Theoretical p-values",
      outer = TRUE, line = 3);
par(xpd=NA)
legend(-2.6,-0.47,legend=c("50","100","200","500","1000"), lty = c(2,3,4,5,1), col=c(2,3,4,5,1), horiz = TRUE);

par(op)

dev.off()





################################################################################################################


load("../results/TypeI_3.RData");

apply(conv,2:3, sum); 

##########calculate non-convergence samples
b1=apply(is.na(pl),2:3, sum);###the LR test
b2=apply(is.na(pw),2:3, sum);  ### the Wald test
b3=apply(is.na(pc),2:3, sum);  ### the score test
b4=apply(is.na(ps),2:3, sum);  ### the new test
print(ymin[which(is.na(pl))]);
print(ymax[which(is.na(pl))]);

######## calculate the rejection rate for zero-inflation for all tests with gamma0=0.05
a1=apply((pl>0.90)&(w_hat>0), 2:3,mean,na.rm=T)*100;  ###the LR test
a2=apply(pw>0.95, 2:3,mean,na.rm=T)*100;  ### the Wald test
a3=apply(pc>0.95, 2:3,mean,na.rm=T)*100;  ### the score test
a4=apply(ps>0.95, 2:3,mean,na.rm=T)*100;  ### the new test
pow=matrix(nrow=5,ncol=16)
pow[,(1:4)*4-3]=a1;
pow[,(1:4)*4-2]=a2;
pow[,(1:4)*4-1]=a3;
pow[,(1:4)*4]=a4;
print(xtable(cbind(c(50,100,200,500,1000),pow)))

##### the p value of LR test due to chi-square#####
pl=1-pl;
##############
## Figure 3

pdf("QQ_normal_m1.pdf")

op <- par(mfrow = c(4,4),
          oma = c(7,4,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1)
############
# first row
# k for the different means
##############
k=1;kk=5;

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=T,main  = "LR Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Wald Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=T,main  = "Score Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=T, main  = "New Test");
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}

text(0.8,.07, bquote(gamma == .(tp[k])));
############
# second row
##############
k=2;kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

############
# third row
##############
k=3; kk=5
plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='n', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='n', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='n', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));
############
# fourth row
##############
k=4;kk=5

plot(((1:length(sort(pl[,kk,k])))-.5)/length(sort(pl[,kk,k])),sort(pl[,kk,k]) ,
     xaxt='s', yaxt='s',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pl[,jj,k])))-.5)/length(sort(pl[,jj,k])),sort(pl[,jj,k]),lty=jj+1,col=jj+1)
}
#text( 0.8, .07, paste("a = ", tp[k]));
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pw[,kk,k])))-.5)/length(sort(pw[,kk,k])),sort(pw[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pw[,jj,k])))-.5)/length(sort(pw[,jj,k])),sort(pw[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(pc[,kk,k])))-.5)/length(sort(pc[,kk,k])),sort(pc[,kk,k]) ,
     xaxt='s', yaxt='n',ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(pc[,jj,k])))-.5)/length(sort(pc[,jj,k])),sort(pc[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));

plot(((1:length(sort(ps[,kk,k])))-.5)/length(sort(ps[,kk,k])),sort(ps[,kk,k]) ,
     xaxt='s', yaxt='n', ylim=c(0,1),type="l",ann=F);
for (jj in 1:(kk-1)){
  lines(((1:length(sort(ps[,jj,k])))-.5)/length(sort(ps[,jj,k])),sort(ps[,jj,k]),lty=jj+1,col=jj+1)
}
text(0.8,.07, bquote(gamma == .(tp[k])));



title(xlab = "Empirical  p-values",
      ylab = "Theoretical p-values",
      outer = TRUE, line = 3);
par(xpd=NA)
legend(-2.5,-0.47,legend=c("50","100","200","500","1000"), lty = c(2,3,4,5,1), col=c(2,3,4,5,1), horiz = TRUE);

par(op)

dev.off()

###############################################################
