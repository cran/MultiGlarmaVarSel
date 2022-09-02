\name{MultiGlarmaVarSel-package}
\alias{MultiGlarmaVarSel-package}
\alias{MultiGlarmaVarSel}
\docType{package}
\title{
\packageTitle{MultiGlarmaVarSel}
}
\description{
MultiGlarmaVarSel consists of four functions: "variable_selection.R", 
"grad_hess_L_gamma.R", "grad_hess_L_eta.R", and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\details{
This package consists of four functions: "variable_selection.R", 
"grad_hess_L_gamma.R", "grad_hess_L_eta.R" and "NR_gamma.R"
For further information on how to use these functions, 
we refer the reader to the vignette of the package.
}
\author{
Marina Gomtsyan

Maintainer: Marina Gomtsyan <marina.gomtsyan@agroparistech.fr>
}
\references{
M. Gomtsyan et al. "Variable selection in sparse multivariate GLARMA models: Application to germination control by environment", arXiv:2208.14721
}
\keyword{package}
\examples{
data(Y)
I=3
J=100
T=dim(Y)[2]
q=1
X=matrix(0,nrow=(I*J),ncol=I)
for (i in 1:I)
{
  X[((i-1)*J+1):(i*J),i]=rep(1,J)
}
gamma_0 = matrix(0, nrow = 1, ncol = q)
result=variable_selection(Y, X, gamma_0, k_max=1, 
n_iter=100, method="min", nb_rep_ss=1000, threshold=0.6)
estim_active = result$estim_active
eta_est = result$eta_est
gamma_est = result$gamma_est
}
