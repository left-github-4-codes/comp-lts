#code for table 1-3 for AA1-LTS and AA2-LTS, ltsReg by Y.Zuo on 03/09/22
#initialization

library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

R=1000   #tuning this parameter for different settings
n=300    #tuning this parameter for different settings
p=9      #tuning this parameter for different settings
alpha=2/4
c=0
cut_off=10^{-3}
N=200
gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p)
t_aa1=t_aa2=t_aa3=0

epsilon=0.30; n1=floor(epsilon*n)  #tuning your epsilon based on different settings

for (i in 1:R)
{ #ss=matrix(0.9, nrow=p, ncol=p); diag(ss)<-1   #choose your m1 based on different settings
  #m1=rmvnorm(n, mean=(rep(0, p)),ss)
  #m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.9, 0.9, 1), ncol=2, byrow=T))
  m1=rmvnorm(n, mean=(rep(0, p)), sigma=diag(rep(1,p)))
  if(n1>=1)
   {
     mu=rep(7,p); mu[p]=-2
     m2=rmvnorm(n1, mu,sigma=diag(rep(0.1, p)))
     m1[sample(1:n,n1),]<-m2 
  }
 z=m1
  
  t1=Sys.time()  
  beta_aa1[i,]=AA1_main_lts(z, alpha, c,N, cut_off )
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=AA2_main_lts(z, alpha, c, gamma, cut_off )
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
   t1=Sys.time()  
   fit1<-ltsReg(z[,1:(p-1)], z[,p])
   beta_aa3[i,]=as.numeric(fit1$coefficients) 
   t2=Sys.time()-t1
   t_aa3=t_aa3+t2
}

EMSE_aa1=sum(beta_aa1*beta_aa1)/R
EMSE_aa2=sum(beta_aa2*beta_aa2)/R
EMSE_aa3=sum(beta_aa3*beta_aa3)/R

print(c(R, n,p,c, N,alpha, gamma,epsilon,cut_off) )
print( c(EMSE_aa1, EMSE_aa2, EMSE_aa3))
print(c(t_aa1, t_aa2, t_aa3))