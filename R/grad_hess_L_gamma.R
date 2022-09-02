grad_hess_L_gamma <-
function(Y,X,eta,gamma,I,J)
{
  T=dim(Y)[2]
  q<-dim(gamma)[2]
  eta_matrice=matrix(eta,nrow=I,ncol=T,byrow = TRUE)
  
  eta_repet=X%*%eta_matrice
  
  mu=matrix(0,nrow=(I*J),ncol=T)
  Z=matrix(0,nrow=(I*J),ncol=T)
  E=matrix(0,nrow=(I*J),ncol=T)
  
  mu[,1]=exp(eta_repet[,1])
  E[,1]=(Y[,1]/mu[,1]-1)
  
  Z[,2]=gamma[1:min(q, 2-1)]%*%t(E[,(2-1):max(1, 2-q)])
  mu[,2]=exp(eta_repet[,2]+Z[,2])
  E[,2]=(Y[,2]/mu[,2]-1)
  
  for (t in 3:T)
  {
    Z[,t]=gamma[1:min(q, t-1)]%*%t(E[,(t-1):max(1, t-q)])
    mu[,t]=exp(eta_repet[,t]+Z[,t])
    E[,t]=(Y[,t]/mu[,t])-1
  }
  
  grad_L_gamma=matrix(0,nrow=1,ncol=q)
  hess_L_gamma=matrix(0,nrow=q,ncol=q)
  
  for (q0 in 1:q)
  {
    grad_W = matrix(0,nrow=(I*J),ncol=T)
    grad_W[,2]=E[,(2-1)]
    
    grad_W_q1 = matrix(0,nrow=(I*J),ncol=T)
    grad_W_q1[,2]=E[,(2-1)]
    
    hess_W=matrix(0,nrow=(I*J),ncol=T)
    
    for (t in 3:T)
    {
      jsup=min(q,(t-1))
      if(t-q <= 0){
        grad_W[,t] = 0
      }
      else{
        grad_W[,t] = E[,(t-q0)]
      }
      if(q0 <= jsup)
      {
        for(k in 1:jsup){
          grad_W[,t] = grad_W[,t] - (1 + E[,(t-k)])%*%as.matrix(gamma[k])*grad_W[,(t-k)]
        }
      }
    }
    grad_L_gamma[q0]=sum((Y-mu)*grad_W)
    
    for(q1 in q0:q){
      for (t in 3:T)
      {
        jsup=min(q,(t-1))
        if(t-q <= 0){
          grad_W_q1[,t] = 0
        }
        else{
          grad_W_q1[,t] = E[,(t-q1)]
        }
        if (q1<=jsup)
        {
          for(k in 1:jsup){
            grad_W_q1[,t] = grad_W_q1[,t] - (1 + E[,(t-k)])%*%as.matrix(gamma[k])*grad_W_q1[,(t-k)]
          }
        }
        if (q0+q1<t) 
        {
          A1=-(1+E[,(t-q0)])*grad_W_q1[,(t-q0)]
          A2=-(1+E[,(t-q1)])*grad_W[,(t-q1)]
        }
        else 
        {
          A1=0
          A2=0
        }
        hess_W[,t]=A1+A2
        for (k in 1:jsup){
          hess_W[,t] = hess_W[,t] - as.matrix(gamma[k])%*%((1 + E[,(t-k)])*grad_W[,(t-k)]*grad_W_q1[,(t-k)]  + (1+E[,(t-k)]) * hess_W[,(t-1)])
        }
      }
      A=(Y-mu)*hess_W
      B=mu*grad_W*grad_W_q1
      C=A-B
      L = sum(C)
      hess_L_gamma[q0, q1]=L
      hess_L_gamma[lower.tri(hess_L_gamma)]=hess_L_gamma[upper.tri(hess_L_gamma)]
    }
  }
  return(list(grad_L_gamma=grad_L_gamma,hess_L_gamma=hess_L_gamma))
}
