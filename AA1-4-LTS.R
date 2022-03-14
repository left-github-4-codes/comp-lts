#AA1 for LTS by Y. Zuo on 03/07/22
library(mvtnorm)
library(robustbase)
library(MASS)

get_initial_beta_lts=function(z,c)
{# input: Z is a n by p matrix, first (p-1) columns for x coordinates, pth for the y coordinates
  # output: a vector beta which might be ls or lts estimator
  # c=0, beta=rep(0,p) or any other constant pre-determined beta;
  # c=1, beta is given by LS; c=2, beta is given by LST 
  
  p=dim(z)[2]
  # print("p"); print(p)
  if (c==1)
  { # get a LTS line 
    beta=AA1_main(z, 1, 0, 500, 0.001) # need to call AA1-1.R in LST folder
  }  
  if (c==2)
  { # get a LS line
    fit2<-lm(z[,p]~z[,1:(p-1)]) 
    beta=as.numeric(fit2$coefficients)
  }
  if (c==0)  
  { # get a pre-determined beta such as rep(0,p) or rep(1,p)
    beta=c(rep(0,p))
  }
  return(beta)
}
#
get_weight_wi_lts=function(z, beta_1, alpha) #alternatively function(z, beta, h)
{# Input: Z (n by p) matrix, beta (p by 1) vector, alpha a constant
 # Output: weight w_i=indicator(|r_i-m(r_i)|/sigma(r_i)\leq alpha)  
  
  n=dim(z)[1]; p=dim(z)[2]
  h=floor(alpha*n) +1 #(default is floor((n+p+1)/2))
  if (p>2) {h= floor((n+p+1)/2)}
  Wn=rep(0, n); rn=rep(0, n); rsqn=rep(1, n);index=rep(0,n)
  
  rn=z[,p]-cbind(rep(1, n), z[, 1:(p-1)])%*%as.matrix(matrix(beta_1, nrow=p, ncol=1)) 
  # n by 1 vector
  rsqn=rn*rn  # squares of n residuals, n by 1
  temp=order(rsqn); id=temp[h] # get the id of the hth squared residual
  threshold=rsqn[id]              # the hth squared residual
  index=which(rsqn<=threshold) # which returns indeice meeting the condition

  diagg=rep(0, n) # 1 by n vector with 0 for all elements
  diagg[index]<-1 # change entries to one based on the which returned indices
  return(diagg)   # a weight vector 1 by n
}
interative_beta=function(z, wi)
{# input: data set z (n by p) matrix, and a 1 by n vector wi given by get_weight_wi 
  # output: beta_new=(Xnn Wnn^k Xnn')^{-1}(Xnn Wnn^k Ynn)
  
  n=dim(z)[1]; p=dim(z)[2] #; xnn=wnn=matrix(0, nrow=n, ncol=p)
  xnn=as.matrix(cbind(matrix(1, nrow=n), z[,1:(p-1)]))  #X matrix n by p with 1st column 1  
  
  ynn=as.matrix(z[, p])                                 # Y vector
  wnn=diag(wi)                               # Weight matrix
  
  m=t(xnn)%*%wnn%*%xnn

  if (det(m)==0) {beta_new=ginv(m)%*%t(xnn)%*%wnn%*%ynn}
  else
  {beta_new=solve(t(xnn)%*%wnn%*%xnn)%*%t(xnn)%*%wnn%*%ynn} # weighted LS solution

  return(c(beta_new))
}
#main part
#Following are parameters that need to be tuned every time one runs the program
AA1_main_lts=function(z, alpha, choice_value_4_initial_beta, replication_number, cut_off)
{ # z is a given data matrix n by p, replication_number is the number user provide
  # for the total iteration allowed for the AA1 procedure (usually 100 is enough)
  # choice_value_4_initial_beta usually have three choices, 0, 1, 2, corresponding to
  # a pre-specified initial beta^0, a LS beta, or a beta from LST
  # cut_off is a value to stop the loop of iteration
  # alpha is the number in the defintion of LTS (or h=floor(n+p+1)/2))
  
  n=dim(z)[1]
  p=dim(z)[2]
  
  N=replication_number      # 200 is more than enough in the most cases
  c=choice_value_4_initial_beta  #c=0, a pre-specified beta, eg rep(0, p), rep(1,p)
                                 #c=1, a LS estimator c=2, the LST estimator of Zuo (2022)
  
  threshold= cut_off                                     # stopping criterion
  beta_old=get_initial_beta_lts(z, c)                    #initial beta^0
  wi=get_weight_wi_lts(z, beta_old, alpha)               # weight wi is a 1 by n row vector

  
  for (i in (1:N))
  {    
    beta_new=interative_beta(z, wi) 
    
    diff=beta_new-beta_old; norm_diff=sqrt(sum((diff)*(diff)))
     
    if(norm_diff<=threshold){break}
    else 
    { wi=get_weight_wi_lts(z, beta_new, alpha) # obtain a new weight vector based on the beta_new
    
    beta_old=beta_new    #update beta_old and resady for the next loop
    }
  }
  return(c(beta_new))
} #entire main end

# testing example
# x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
# y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
# m1=cbind(x,y)
# dev.off()
# par(mfrow=c(1,1))
# plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
# abline(0,0, lty=1, col=1, lwd=1)
# abline(0, 1, lty=2, col=2, lwd=1)
# p=dim(m1)[2]
# fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
# beta_ltsReg=as.numeric(fit1$coefficients) 
# abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=3, lty=3, lwd=1)
# lst_beta=AA1_main_lts(m1, 1/2, 1, 500, 0.001)
# lst_beta_1=AA2_main_lts(m1, 1/2, 1, 100, 0.001)
# abline(lst_beta,lty=4, col=4, lwd=1)
# abline(lst_beta_1,lty=5, col=5, lwd=1)
# text(x, y-0.5, as.character(seq(1:7)))
# text(8, -1, expression('L'[1]))
# text(7, 7.2, expression('L'[2]))
# text(7.5, 1.8, expression ('LTS-ltsReg'))
# text(7.5, 6.5, expression('LST-AA1'))
# legend(-2.5, 6.5, legend=c(expression('L'[1]),expression('L'[2]), "LTS-ltsReg", "LTS-AA1"),    
#        col=1:4, lty=1:4, cex=0.8,
#        title="Line types", text.font=3, bg='lightblue')