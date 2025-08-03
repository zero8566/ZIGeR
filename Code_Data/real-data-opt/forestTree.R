#bioChemists
rm(list=ls())
setwd("D:/0Adas/Code_and_Data/")

#source("JABES_prog_Rcode.R")
source("./sim/real-likehood.R")
source("./sim/real-bat-test.R")
#source("D:/0Adas/0收纳论文写作/10删失几何零膨胀检验/codes/sim/real-likehood-yinwen.R")


load("./real-data-opt/BCI.data.RData")


data<-BCI.data
hist(data$y)

cut<-17
n <-nrow(BCI.data)  # Sample size.
y<-data$y



# 计算每个唯一值的频数
y_counts <- table(y)

barplot(y_counts, 
        xlab = "Count data of Hybanthus prunifolius on BCI", 
        ylab = "Frequency",
        col = "#10B981",
        border = "black",           # 设置条形图边框颜色
        las = 1,                    # 让y轴标签水平显示
        cex.lab = 1.2,              # 增大坐标轴标签字体
        cex.axis = 1.1,             # 增大坐标轴刻度字体
        cex.main = 1.3,             # 增大标题字体
        ylim = c(0, max(y_counts) * 1.1)  # 稍微扩展y轴范围
)

# 添加网格线
grid(col = "gray80", lty = "dotted")  # 设置网格线颜色和线型

# 添加边框，确保图形被封闭
box(col = "gray30", lwd = 1.2)        # 设置边框颜色和粗细













cov<-x<-as.matrix(rep(1, length(y))) 

outcome=data$y
head(data, n = 5)
dat=data.frame(x=x,y=outcome)
head(dat, n = 5)
a<-glm(outcome~ 1,family=negative.binomial(1))
beta=as.vector(coef(a));



res0<-try(fit0<-stats4::mle(lge0,start=beta,nobs=NROW(dat)));
loglike0<-logLik(fit0)
length(beta)

zz<-coef(fit0)
beta=as.vector(zz);

# 
# 
# res<-optim(beta,lge0)
# 
# zz<-res$par
# beta=as.vector(zz);
# beta_mle= beta;

# try(loglike1<--res$value);

try(z<-test.new(beta,y,cov,cut));
try(ps<-1-pnorm(z)); ### the p value of  new test for zero-inflation; 
try(c<-test.score(beta,y,cov,cut));
try(pc<-1-pnorm(c));

#try(loglike1<--res$value);


res2<-try(fit2<-stats4::mle(lzige0,start =c(0,beta),nobs =NROW(dat),control=list(maxit=1000)));


if(class(res2) != "try-error"&& (fit2@details$convergence==0))
{ 
  omega_hat=unname(coef(fit2)[1]);  
  omega_se=unname(sqrt(diag(vcov(fit2)))[1]);
  w=omega_hat/omega_se; #Wald statistic
  try(pw<-1-pnorm(w));  ##the p value of  Wald test for zero-inflation
  
  try(loglike2<-logLik(fit2));
  l<-2*(loglike2-loglike0);
  try(pl<-1-pchisq(l,df =1));###the p value of LR test
}

#remove(fit2);

aa2_w=c(omega_hat,omega_se);
aa2=round(as.matrix(rbind(t(c(w,l,c,z)),t(c(pw,pl,pc,ps)))),4); 
colnames(aa2) <- c("Wald","LR","score","new");
rownames(aa2)<- c("Statistic", "P-value")
print(aa2)


# res$par
mu=exp(1.3307)
mu




##########计算系数##aic-bic##########
# fit1<-glm(outcome ~ x,family=negative.binomial(1))
# summary(fit1)
# beta_hat1=coef(fit1);
# beta_hat1



#######1---RCGeR###########的AIC和BIC
res0
# 获取系数
coefficients0 <-coef(res0)
# 获取协方差矩阵
cov_matrix0 <- vcov(res0)
# 计算标准误
standard_errors0 <- sqrt(diag(cov_matrix0))
# 计算t统计量
t_statistics0 <- coefficients0 / standard_errors0
# 计算p值（假设为双侧检验）
p_values10 <- 2 * (1 - pt(abs(t_statistics0), df = length(dat$y)- length(coefficients0)))
print(p_values10)

sign10<- ifelse(p_values10 < 0.001, "***", 
                ifelse(p_values10 < 0.01, "**", 
                       ifelse(p_values10 < 0.05, "*", 
                              ifelse(p_values10 < 0.1, ".", ""))))

sign10<-as.data.frame(sign10)

fit00<-cbind(coef(res0),standard_errors0,t_statistics0,p_values10,sign10)

fit00
#####################fit2###############
log_likelihood0 <- logLik(fit0)
# 获取模型的参数个数
k0 <- length(coef(fit0))
# 获取样本大小
n <- NROW(dat)
# 手动计算AIC
aic_value_manual0 <- 2 * k0 - 2 * log_likelihood0
# 手动计算BIC
bic_value_manual0 <- log(n) * k0 - 2 * log_likelihood0



aic_value_manual0 <- format(round(aic_value_manual0, 4), nsmall = 4)
bic_value_manual0 <- format(round(bic_value_manual0, 4), nsmall = 4)

aic_value_manual0
bic_value_manual0








####2---ZIRCGeR##AICBIC#########
res2
# 获取系数
coefficients <-coef(res2)
# 获取协方差矩阵
cov_matrix <- vcov(res2)
# 计算标准误
standard_errors <- sqrt(diag(cov_matrix))
# 计算t统计量
t_statistics <- coefficients / standard_errors
# 计算p值（假设为双侧检验）
p_values1 <- 2 * (1 - pt(abs(t_statistics), df = length(dat$y)- length(coefficients)))
print(p_values1)

sign1<- ifelse(p_values1 < 0.001, "***", 
               ifelse(p_values1 < 0.01, "**", 
                      ifelse(p_values1 < 0.05, "*", 
                             ifelse(p_values1 < 0.1, ".", ""))))

sign1<-as.data.frame(sign1)

fit222<-cbind(coef(res2),standard_errors,t_statistics,p_values1,sign1)

fit222

#####################fit2###############
log_likelihood <- logLik(fit2)
# 获取模型的参数个数
k <- length(coef(fit2))
# 获取样本大小
n <- NROW(dat)
# 手动计算AIC
aic_value_manual <- 2 * k - 2 * log_likelihood
# 手动计算BIC
bic_value_manual <- log(n) * k - 2 * log_likelihood

# summary(fit1)
print(p_values1)

aic_value_manual <- format(round(aic_value_manual, 4), nsmall = 4)
bic_value_manual <- format(round(bic_value_manual, 4), nsmall = 4)

aic_value_manual
bic_value_manual





