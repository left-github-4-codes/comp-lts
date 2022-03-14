#code for figure two by Y. Zuo on 03/09/22

n=100; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
m11=m1 #this steps are to fix a random a wanted(nice) data set so that to add text 
m1=m11
par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main_lts(m1, 1/2, 1, 500, 0.001), col=2, lty=2, lwd=1)
#abline(AA1_main_lts(m1,1,500, 0, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)

text(5.8,7, expression("LTS-ltsReg"))
text(6,6., expression("LS"))
text(6.2, 3, expression("LTS-AA1"))

legend(-2.5, 6.5, legend=c("LTS-ltsReg", "LST-AA1", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}

m21=m1 #this steps are to fix a random a wanted(nice) data set so that to add text 
m1=m21
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main_lts(m1,1/2, 1, 500,0.001), col=2, lty=2, lwd=1)
#abline(AA2_main_lts(m1,1/2, 1, 100,0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)

text(2.5, -1.6, expression ('LTS-ltsReg'))
text(6, 5.7, expression('LTS-AA1'))
text(4, -0.9, expression('LS'))

legend(-2.5, 6.5, legend=c("LTS-ltsReg", "LTS-AA1", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')
