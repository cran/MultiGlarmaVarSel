grad_hess_L_eta <-
function(Y,X,eta_vect,gamma,I,J)
{
  T=dim(Y)[2]
  q<-dim(gamma)[2]
  eta_matrice=matrix(eta_vect,nrow=I,ncol=T,byrow = TRUE)
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
  
  grad_L_eta=matrix(0,nrow=1,ncol=(I*T))
  
  for (i0 in 1:I)
  {
    for (t0 in 1:T)
    {
      grad_W=matrix(0,nrow=(I*J),ncol=T)
      der_eta_matrice=matrix(0,nrow=I,ncol=T)
      der_eta_matrice[i0,t0]=1
      der_eta_repet=X%*%der_eta_matrice
      if (t0 == 1) {grad_W[((i0-1)*J+1):(i0*J),1]=1}
      for (t in 2:T)
      { 
        jsup=min(q,(t-1))
        grad_W[,t]=der_eta_repet[,t]
        for(k in 1:jsup){
          grad_W[,t] = grad_W[,t] - (1+E[,(t-k)])%*%as.matrix(gamma[k])*grad_W[,(t-k)]
        }
        
      }
      grad_L_eta[(i0-1)*T+t0]=sum(grad_W*(Y-mu))
    }
  }
  
  hess_L_eta=Matrix(0,nrow=(I*T),ncol=(I*T))
  block_hess_list=list()
  block_hess_inv_list=list()
  
  for (i0 in 1:I)
  {
    block_hess_eta=matrix(0,nrow=T,ncol=T)
    
    for (t0 in 1:T)
    {
      for (t1 in t0:T)
      {
        hess_W=matrix(0,nrow=(I*J),ncol=T)
        grad_W_t0=matrix(0,nrow=(I*J),ncol=T)
        grad_W_t1=matrix(0,nrow=(I*J),ncol=T)
        
        der_eta_matrice_t0=matrix(0,nrow=I,ncol=T)
        der_eta_matrice_t0[i0,t0]=1
        der_eta_repet_t0=X%*%der_eta_matrice_t0
        der_eta_matrice_t1=matrix(0,nrow=I,ncol=T)
        der_eta_matrice_t1[i0,t1]=1
        der_eta_repet_t1=X%*%der_eta_matrice_t1
        if (t1==1) grad_W_t1[((i0-1)*J+1):(i0*J),1]=1
        if (t0==1) grad_W_t0[((i0-1)*J+1):(i0*J),1]=1
        for (t in 2:T)
        {
          jsup=min(q,(t-1))
          grad_W_t0[,t]=der_eta_repet_t0[,t]
          grad_W_t1[,t]=der_eta_repet_t1[,t]
          hess_W[,t] = 0
          for(k in 1:jsup){
            grad_W_t0[,t] = grad_W_t0[,t] - (1+E[,(t-k)])*gamma[k]*grad_W_t0[,(t-k)]
            grad_W_t1[,t] = grad_W_t1[,t] - (1+E[,(t-k)])*gamma[k]*grad_W_t1[,(t-k)]
            term1 = (1+E[,(t-k)])*gamma[k]*grad_W_t0[,(t-k)]*grad_W_t1[,(t-k)]
            term2 = -(1+E[,(t-k)])*gamma[k]*hess_W[,(t-k)]
            hess_W[,t] = hess_W[,t] + term1 + term2
          }
        }
        A=(Y-mu)*hess_W
        B=mu*grad_W_t0*grad_W_t1
        C=A-B
        block_hess_eta[t0,t1]=sum(C[((i0-1)*J+1):(i0*J),])
      }
    }
    
    block_hess_list[[i0]]=(block_hess_eta+t(block_hess_eta))/2
    res_svd=svd(block_hess_list[[i0]])
    ind_vp_not_null=which(round(res_svd$d,digits=6)!=0)
    #block_hess_inv_list[[i0]]=res_svd$v[,ind_vp_not_null]%*%diag(1/(res_svd$d[ind_vp_not_null]), length(ind_vp_not_null))%*%(t(res_svd$u[,ind_vp_not_null]))
  }
  hess_L_eta=bdiag(block_hess_list)
  #hess_L_eta_inv=bdiag(block_hess_inv_list)
  
  return(list(grad_L_eta=grad_L_eta,hess_L_eta=hess_L_eta))
}
