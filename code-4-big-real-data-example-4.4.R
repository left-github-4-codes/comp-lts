#code for real big data set by Y. Zuo on 03/11/22
#code for the big real data set by Y. Zuo on 01/19/22
#rm(list=ls())
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

# library(AER)
# #library(MASS)
# #data(CPS1985)
# data(CPS1985)
# #attach(CPS1985)
# z=CPS1985
# c1=colnames(z)[2:11]
# ctemp=rep("wage",11)
# ctemp[1:10]<-c1
# z=z[ctemp] # switch the 1st colum of z to the last colum since wage is the reponse variable (y) 

library(gamlss.data)
data(plasma)
z=plasma
colnames(z)<-NULL #remove the column names
#data(Boston, package="MASS")
options(warn=-1)

###########################################################################
##### just entire data set and calulate the total time consumed and sample variance
###########################################################################
#z=Boston
n=dim(z)[1]; p=dim(z)[2]
z=matrix(unlist(z), nrow=n, ncol=p)
zz=z
plasma=z

RepN=1000; R=RepN
alpha=1/2 #h will automaticaly be set to be the default value: h=floor((n+p+1)/2) when p>2
c=0
cut_off=10^{-3}
epsilon=cut_off
N=200
gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by LTS
t_aa1=t_aa2=t_aa3=0
#i=1
for (i in 1:RepN)
{
  t1=Sys.time()  
  beta_aa1[i,]=AA1_main_lts(z, alpha, c,N, cut_off )
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=AA2_main_lts(z,alpha, c,gamma, cut_off )
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
}
beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, nrow=RepN, ncol=p)

EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN

print(c(R, n,p,c, N,alpha, gamma,epsilon,cut_off) )
print(c(t_aa1, t_aa2, t_aa3))
print(c(EMSE_aa1, EMSE_aa2, EMSE_aa3))


############################################################################
############## sampling n1 points from the entire data set
############################################################################
z=zz
#z=Boston
#Bos_minus=z[,c(-3,-7)] # delete indus and age two predictor since they are not significant with p-value <0.001
#z=Bos_minus
#n=dim(z)[1]; p=dim(z)[2]
RepN=10000; R=RepN
n1=315 # tune this to 50, 100, 150, 200, or 300
alpha=1/2 
c=0
cut_off=10^{-3}
N=200
gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by ltsReg
t_aa1=t_aa2=t_aa3=0
for (i in 1:RepN)
{
  index=sample(1:n, n1)
  z=zz[index,]

  t1=Sys.time()  
  beta_aa1[i,]=AA1_main_lts(z, alpha, c,N, cut_off )
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=AA2_main_lts(z,alpha, c,gamma, cut_off )
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
  
  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
}

beta_aa1_mean=colMeans(beta_aa1)
deviat_beta_aa1=beta_aa1-matrix(beta_aa1_mean, nrow=RepN, ncol=p)
beta_aa2_mean=colMeans(beta_aa2)
deviat_beta_aa2=beta_aa2-matrix(beta_aa2_mean, nrow=RepN, ncol=p)
beta_aa3_mean=colMeans(beta_aa3)
deviat_beta_aa3=beta_aa3-matrix(beta_aa3_mean, nrow=RepN, ncol=p)

EMSE_aa1=sum(deviat_beta_aa1*deviat_beta_aa1)/RepN
EMSE_aa2=sum(deviat_beta_aa2*deviat_beta_aa2)/RepN
EMSE_aa3=sum(deviat_beta_aa3*deviat_beta_aa3)/RepN

print(c(R, n,p, n1,c, N,alpha, gamma,epsilon,cut_off) )
print(c(t_aa1, t_aa2, t_aa3))
print(c(EMSE_aa1, EMSE_aa2, EMSE_aa3))


########################################################################################
#### Second part (a) dropping age and indus and (b) one part for estimation and the other part for predition
########################################################################################
zzz=plasma #zz used below and changes every time so reserve zz in zzz
z=zzz 
#z=Boston
#Bos_minus=zz[,c(-3,-7)] # delete indus and age two predictor since they are not significant with p-value <0.001
#z=Bos_minus
n=dim(z)[1]; p=dim(z)[2]
RepN=10000; R=RepN
m=150  # tune this 200, 300, 400
alpha=1/2 # one or three?
c=0
cut_off=10^{-3}
N=100
gamma=10^2
beta_aa1=beta_aa2=beta_aa3=matrix(0, nrow=R, ncol=p) #aa3 is actually replace by ltsReg
rss_aa1=rss_aa2=rss_aa3=0
res_aa1=res_aa2=res_aa3=matrix(0, nrow=n-m, ncol=1)
t_aa1=t_aa2=t_aa3=t1=t2=0
#i=1
for (i in 1:RepN)
{
  index=sample(1:n, m)
  z=zzz[index,]
  zz=zzz[-index,] # could drop this and use all data points

  t1=Sys.time()  
  fit1<-ltsReg(z[,1:(p-1)], z[,p])
  beta_aa3[i,]=as.numeric(fit1$coefficients) 
  res_aa3=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  
  if(sum(beta_aa3[i,]!=0)) #ltsReg often output the initialized zero vector for beta
{  
  rss_temp=sum(res_aa3*res_aa3)
  rss_aa3=rss_aa3+rss_temp
  t2=Sys.time()-t1
  t_aa3=t_aa3+t2
    
  t1=Sys.time()  
  beta_aa1[i,]=AA1_main_lts(z, alpha, c, N, cut_off )
  res_aa1=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa1[i,], nrow=p, ncol=1)
  
  rss_temp=sum(res_aa1*res_aa1)
  rss_aa1=rss_aa1+rss_temp
  t2=Sys.time()-t1
  t_aa1=t_aa1+t2
  
  t1=Sys.time()  
  beta_aa2[i,]=AA2_main_lts(z,alpha, c,gamma, cut_off )
  res_aa2=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa2[i,], nrow=p, ncol=1)
 
  rss_temp=sum(res_aa2*res_aa2)
  rss_aa2=rss_aa2+rss_temp
  t2=Sys.time()-t1
  t_aa2=t_aa2+t2
} # end of if 
  # t1=Sys.time()  
  # fit1<-ltsReg(z[,1:(p-1)], z[,p])
  # beta_aa3[i,]=as.numeric(fit1$coefficients) 
  # res_aa3=zz[,p]-as.matrix(cbind(matrix(1, nrow=n-m), zz[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  # #res_aa3=Bos[,p]-as.matrix(cbind(matrix(1, nrow=n), Bos[, 1:(p-1)]))%*%matrix(beta_aa3[i,], nrow=p, ncol=1)
  # 
  # rss_temp=sum(res_aa3*res_aa3)
  # rss_aa3=rss_aa3+rss_temp
  # t2=Sys.time()-t1
  # t_aa3=t_aa3+t2
} # end of for loop

print(c(R, n,p, m,c, N,alpha, gamma,epsilon,cut_off) )
print(c(t_aa1, t_aa2, t_aa3))
print(c(rss_aa1, rss_aa2, rss_aa3)/RepN)
