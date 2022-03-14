#for example 1.1 Figure one part (a) by Y. Zuo on 03/07/22
##############Part a#####################################################
x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
m1=cbind(x,y)
par(mfrow=c(1,2))
plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
abline(0,0, lty=1, col=1, lwd=1)
abline(0, 1, lty=2, col=2, lwd=1)
text(x, y-0.5, as.character(seq(1:7)))
text(7.5, -.5, expression('L'[1]))
text(7, 6, expression('L'[2]))

p=dim(m1)[2]
plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
beta_ltsReg=as.numeric(fit1$coefficients) 
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
betalts=AA1_main_lts(m1, 1/2, 1, 500, 0.001)
abline(betalts,lty=2, col=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)
text(x, y-0.5, as.character(seq(1:7)))
#text(8, -1, expression('L'[1]))
#text(7, 7.2, expression('L'[2]))
text(6.7, 1.8, expression ('LTS-ltsReg'))
text(6.7, 6.5, expression('LST-AA1'))
text(7.3, 2.5, expression('LS'))
legend(-2.5, 6.5, legend=c(#expression('L'[1]),expression('L'[2]), 
                           "LTS-ltsReg", "LTS-AA1", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')

#####part b#######################################################
######Part (b)
#alpha=1/2     #cut-off value alpha
#N=500               # the total iteration number allowed
#c=0                #c=0,initial beta is rep(0, p);c=1,it is by LS; c=2, it is by LTS.

n=7; p=2#
epsilon=0.3; n1=floor(epsilon*n)
m1=rmvnorm(n, mean=(rep(0, p)),sigma=matrix(c(1, 0.88, 0.88, 1), ncol=2, byrow=T))
#m11=m1 #this steps are to fix a random a wanted(nice) data set so that to add text 
#m1=m11
par(mfrow=c(1,2))  
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit0<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit0$coefficients[[1]], fit0$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main_lts(m1, 1/2, 1, 500, 0.001), col=2, lty=2, lwd=1)
#abline(AA1_main_lts(m1,1,500, 0, 0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)

text(6,4, expression("LTS-ltsReg"))
text(5,5., expression("LS"))
text(3.4, 6, expression("LTS-AA1"))

legend(-2.5, 6.5, legend=c("LTS-ltsReg", "LST-AA1", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')

if(n1>=1)
{
  m2=rmvnorm(n1, mean=c(7,-2),sigma=diag(rep(0.1, p)))
  m1[sample(1:n,n1),]<-m2 
}

#m21=m1 #this steps are to fix a random a wanted(nice) data set so that to add text 
#m1=m21
plot(m1[,1], m1[,2], xlim=c(-3,8), ylim=c(-2, 7), xlab="x-axis", ylab="y-axis", pch=20)
fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=1, lty=1, lwd=1)
abline(AA1_main_lts(m1,1/2, 1, 500,0.001), col=2, lty=2, lwd=1)
fit2<-lm(m1[,2]~m1[,1]) #p is assume to be 2
abline(fit2, col=3, lty=3, lwd=1)

text(4., -1.6, expression ('LTS-ltsReg'))
text(6, 5.7, expression('LTS-AA1'))
text(4, 0, expression('LS'))

legend(-2.5, 6.5, legend=c("LTS-ltsReg", "LTS-AA1", "LS"),    
       col=1:3, lty=1:3, cex=0.6,
       title="Line types", text.font=3, bg='lightblue')
