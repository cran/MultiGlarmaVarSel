NR_gamma <-
function(Y,X,eta,gamma,I,J,n_iter=100)
{
  q=dim(gamma)[2]
  gamma_current = gamma
  eps=10^-6
  error=1+eps
  iter=1
  
  while((error>eps)&&(iter<=n_iter)){
    res=grad_hess_L_gamma(Y,X,eta,gamma_current,I,J)
    Hessienne = res$hess_L_gamma
    if (is.na(sum(Hessienne)) || is.infinite(sum(Hessienne))){
      gamma_new=gamma_current
      break}
    
    res_svd=svd(res$hess_L_gamma)
    U=res_svd$u
    V=res_svd$v
    ind_vp_not_null=which(round(res_svd$d,digits=6)!=0)
    Lambda=diag(1/res_svd$d[ind_vp_not_null],length(ind_vp_not_null))
    hess_inv=V[,ind_vp_not_null]%*%Lambda%*%t(U[,ind_vp_not_null])
    gamma_new=gamma_current-as.numeric(hess_inv%*%t(res$grad_L_gamma))
    
    error=max(abs(gamma_new-gamma_current))
    gamma_current = gamma_new
    iter=iter+1
  }
  return(matrix(gamma_new, nrow = 1, ncol = q))
}
