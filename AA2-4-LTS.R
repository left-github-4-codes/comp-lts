#AA2-4-lts for LTS by Y. Zuo on 03/09/22
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)
# will call functions: get_initial_beta_lts to get initial beta^0, get_weight_wi_lts for the weights
# they are defined in AA1-4-LTS.R #(AA1-1.R or AA1.R or AA11.R)


get_objective_value1=function(z, wi,beta_k) #calculate quantity Q at beta_k
{ #z is a data matrix n by p, in the ith row, the first (p-1) columns for x_i  and pth column for y_i 
  #w_i is a 1 by n row vector, obtained by get_weight_wi_lts function
  #beta_k is the kth step x value (or x^k)
  
  n=dim(z)[1]; p=dim(z)[2] 
  rn=rep(0, n)
  xnn=as.matrix(cbind(matrix(1, nrow=n), z[,1:(p-1)]))          #X matrix n by p with 1st column 1  
  wnn=diag(wi)                                                  #Weight matrix
  rn=z[,p]-as.matrix(cbind(matrix(1, nrow=n), z[, 1:(p-1)]))%*%matrix(beta_k, nrow=p, ncol=1)
  #residuals n by 1 vector                             
  objective_value=t(rn)%*%wnn%*%rn           #O(beta_k)
  return(objective_value)
}
get_objective_value=cmpfun(get_objective_value1)

get_gradient_at_beta_k1=function(z,beta_k,wi) #return nabla_f(beta_k)
{ #z is a data matrix n by p, in the ith row, the first (p-1) columns for x_i  and pth column for y_i 
  #w_i is a 1 by n row vector, obtained by get_weight_wi_lts function
  #beta_k is the x-value at kth step
  
  n=dim(z)[1]; p=dim(z)[2] #; xnn=wnn=matrix(0, nrow=n, ncol=p)
  rn=rep(0, n)
  xnn=as.matrix(cbind(matrix(1, nrow=n), z[,1:(p-1)]))  #X matrix n by p with 1st column 1  
  wnn=diag(wi)                               #Weight matrix
  rn=z[,p]-as.matrix(cbind(matrix(1, nrow=n), z[, 1:(p-1)]))%*%matrix(beta_k, nrow=p, ncol=1)  
  #residuals n by 1 vector
  nabla_f_beta_k=-2*t(xnn)%*%wnn%*%rn        #p by 1 vector-the gradient at beta_k
  #here we have used the 1st-order derivation of O(beta)
  return(nabla_f_beta_k)  
}
get_gradient_at_beta_k=cmpfun(get_gradient_at_beta_k1)

get_modified_Hessian_matrix1=function(z,gamma, wi) #return H(x)+gamma I
{ #z is a data matrix n by p, in the ith row, the first (p-1) columns for x_i  and pth column for y_i 
  #w_i is a 1 by n row vector, obtained by get_weight_wi_lts function
  #gamma is the factor in the [H(x)+gamma I] in the modified Marquardt-Levenberg method
  
  n=dim(z)[1]; p=dim(z)[2] 
  xnn=as.matrix(cbind(matrix(1, nrow=n), z[,1:(p-1)]))  #X matrix n by p with 1st column 1  
  wnn=diag(wi)                               # Weight matrix
  
  hessian=2*t(xnn)%*%wnn%*%xnn               # the 2nd order derivative of O^n(beta)
  hessian_new=hessian + gamma*diag(rep(1,p)) #first term is the 2nd order derivative
  
  return(list(Hessian=hessian, Hessian_new=hessian_new))
}
get_modified_Hessian_matrix=cmpfun(get_modified_Hessian_matrix1)

get_search_direction1=function(modified_hessian, nabla_f_beta_k) #retuen s^k (or Delta_x_{nt}), 
  #the search direction
{ #modified_hessian is the matrix H(x)+gamma I
  #nabla_f_beta_k is the gradient of O at beta_k
  
  inverse_matrix=solve(modified_hessian)
  
  s_k= -1*inverse_matrix%*%nabla_f_beta_k #corresponding to Newton's step Delta_x_{nt} or s^k
  return(s_k)  
}  
get_search_direction=cmpfun(get_search_direction1)

search_a_k1=function(z, beta_k, s_k, wi, nabla_f_beta_k)#backtracking_line_search
{ #z is a data matrix n by p, in the ith row, the first (p-1) columns for x_i  and pth column for y_i 
  #w_i is a 1 by n row vector, obtained by get_weight_wi_lts function
  #beta_k is the kth step x value (or x^k)
  #s_k is the negtive gradient direction, corresponding to Newton's step Delta x_{nt}
  
  Delta_x=s_k; alpha1=1/4; beta1=1/2; t=1    #0<alpha<.5; 0<beta<1, t=1
  beta_k_new=beta_k+t*Delta_x; first_order_value=alpha1*t*t(nabla_f_beta_k)%*%Delta_x
  Q1=get_objective_value(z,wi,beta_k_new)
  Q2=get_objective_value(z,wi,beta_k)+first_order_value
  while(Q1>Q2) 
  {t=beta1*t          #reduce step size
  beta_k_new=beta_k+t*Delta_x; first_order_value=alpha1*t*t(nabla_f_beta_k)%*%Delta_x
  Q1=get_objective_value(z,wi,beta_k_new)
  Q2=get_objective_value(z,wi,beta_k)+first_order_value
  }
  return(t) # this is  a^k, or the search step size (or legth)
}           # it guarantees f(x+tDelta)<=f(x)+first_order_value # the 2nd term is negative
# it guarantees f(x_k+a_ks_k)<f(x_k)
search_a_k=cmpfun(search_a_k1)

AA2_main_lts=function(z, alpha, choice_value_4_initial_beta, gamma, epsilon) 
{#z is a data matrix n by p, in the ith row, the first (p-1) columns for x_i  and pth column for y_i 
  #epsilon is the cut-off value for stopping procedure
  #alpha is the value used in the LST definition
  #gamma is the initial big value to force the hessian to be positive definite 
  # choice_value_4_initial_beta usually have three choices, 0, 1, 2, corresponding to
  # a pre-specified initial beta^0, a LS beta, or a beta from LTS
  
  k=0
  c=choice_value_4_initial_beta
  beta_k=get_initial_beta_lts(z,c)
  wi=get_weight_wi_lts(z,beta_k,alpha) #; print("t0")
  nabla_f_beta_k=get_gradient_at_beta_k(z, beta_k, wi) #; print("t1")
  
  norm_gradient=norm(nabla_f_beta_k, type="2")
  while(norm_gradient>=epsilon & k<=50)#here should add another condition and k<50 to break the loop
  {                                    #while((norm_gradient>=epsilon)&(k<50))
    
    t=get_modified_Hessian_matrix(z, gamma, wi)
    
    hessian_new=t$Hessian_new; hessian=t$Hessian
    
    if(det(hessian)!=0) {henssian_new=hessian}
    s_k=get_search_direction(hessian_new, nabla_f_beta_k) # corresponding to Newton's step Delta_x_{nt}
    
    norm_s_k=norm(s_k, type="2")
    
    if(norm_s_k<epsilon)
    {break}
    quadratic=t(nabla_f_beta_k)%*%s_k #-nabla_f_beta_k (Delta_f_beta_k)^{-1} nabla_f_beta_k
    
    if(quadratic>=0)
    {gamma=2*gamma; 
    }
    else #if (quadratic<0) #on the right descent direction
    { 
      a_k=search_a_k(z,beta_k, s_k, wi, nabla_f_beta_k)#get search step/length
      gamma=gamma/2                                    #reduce gamma
      beta_k=beta_k+a_k*s_k                            # update beta_k
      k=k+1                                           
      
      wi=get_weight_wi_lts(z,beta_k,alpha)
      nabla_f_beta_k=get_gradient_at_beta_k(z, beta_k, wi)
      norm_gradient=norm(nabla_f_beta_k, type="2")
    }
    if(gamma<=epsilon){break}
    
  }
  return(t(beta_k))  
}

# Example to call AA2_main

# x=c(5, 5.5, 4, 3.5, 3, 2.5, -2)
# y=c(-.5, -.5, 6, 4, 2.4, 2, .5)
# m1=cbind(x,y)
# par(mfrow=c(1,1))
# plot(m1, pch=16, xlim=c(-3,8), ylim=c(-2, 7))
# abline(0,0, lty=1, col=1, lwd=1)
# abline(0, 1, lty=2, col=2, lwd=1)
# p=dim(m1)[2]
# fit1<-ltsReg(m1[,1:(p-1)], m1[,p])
# beta_ltsReg=as.numeric(fit1$coefficients) 
# abline(fit1$coefficients[[1]], fit1$coefficients[[2]], col=3, lty=3, lwd=1)
# #lst_beta=AA1_main(m1,1,0,500, 0.001)
# lst_beta_1=AA2_main(m1, 1, 0,200, 0.001)
# #abline(lst_beta,lty=3, col=4, lwd=1)
# abline(lst_beta_1,lty=4, col=5, lwd=1)
