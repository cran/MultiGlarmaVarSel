variable_selection <-
function(Y, X, gamma, k_max=1, n_iter=100, method="min", nb_rep_ss=1000,threshold=0.6)
{
    gamma0=gamma
    I=dim(X)[2]
    J=dim(X)[1]/I
    T=dim(Y)[2]
    eta_glm_mat_0 = matrix(0,ncol=T,nrow=I)
    for (t in 1:T)
    {
      result_glm_0 = glm(Y[,t]~X-1,family=poisson(link='log'))
      eta_glm_mat_0[,t]=as.numeric(result_glm_0$coefficients)
    }
    eta0 = round(as.numeric(t(eta_glm_mat_0)),digits=6)
    
    for(k in 1:k_max){
      gamma_est = NR_gamma(Y,X,eta0,gamma0,I,J,n_iter=100)
      grad_hess_res_est = grad_hess_L_eta(Y,X,eta0,gamma_est,I,J)
      Gradient=grad_hess_res_est$grad_L_eta
      Hessienne=grad_hess_res_est$hess_L_eta
      res_svd=svd(-Hessienne)
      U=res_svd$u
      ind_vp_not_null=which(round(res_svd$d,digits=10)!=0)
      Lambda_rac=diag(sqrt(res_svd$d[ind_vp_not_null]))
      Lambda_rac_inv=diag(1/sqrt(res_svd$d[ind_vp_not_null]))
      Y_eta=Lambda_rac_inv%*%t(U[,ind_vp_not_null])%*%t(Gradient)+Lambda_rac%*%t(U[,ind_vp_not_null])%*%eta0
      X_eta=Lambda_rac%*%t(U[,ind_vp_not_null])  
      
      if(method=='min'){
        all_lambda=glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=TRUE)$lambda
        lambda_min = c(min(all_lambda))
      }else if(method=='cv'){
        lambda_min=cv.glmnet(X_eta,Y_eta,family="gaussian",alpha=1,parallel=TRUE)$lambda.min
      }
      
      p = length(eta0)
      b_sort_matrix = matrix(NA, nrow = nb_rep_ss, ncol = floor(p/2))
      for(j in 1:nb_rep_ss){
        b_sort_matrix[j,] <- sort(sample(1:p,floor(p/2)))
      }
      
      res.cum = rep(0,(p+1))
      for(j in 1:nb_rep_ss){
        b_sort = b_sort_matrix[j,]
        resultat_glmnet=glmnet(X_eta[b_sort,],Y_eta[b_sort],family="gaussian",alpha=1,lambda=lambda_min)
        ind_glmnet=which(as.numeric(resultat_glmnet$beta)!=0)
        res.cum = res.cum + tabulate(ind_glmnet,(p+1))
      }
      freq=res.cum/nb_rep_ss
      estim_active=which(freq>=threshold)
      
      eta_est_mat = matrix(0,ncol=T,nrow=I)
      eta_vec_temporary=rep(0,(I*T))
      eta_vec_temporary[estim_active]=1
      eta_mat_temporary=matrix(eta_vec_temporary,nrow=I,ncol=T,byrow=TRUE)
      
      for (t in 1:T)
      {
        active_coeffs = which(eta_mat_temporary[,t]!=0,arr.ind = T)
        if (length(active_coeffs) > 0){
          
          result_glm = glm(Y[,t]~X[,active_coeffs]-1,family=poisson(link='log'))
          eta_est_mat[active_coeffs,t]=as.numeric(result_glm$coefficients)
        }
      }
      eta_est = round(as.numeric(t(eta_est_mat)),digits=6)
      eta0 = eta_est
      gamma0 = gamma_est
    }
    
    return(list(estim_active=estim_active,eta_est=eta_est,gamma_est=gamma_est))
  }
