\name{grad_hess_L_gamma}
\alias{grad_hess_L_gamma}
\title{
Gradient and Hessian of the log-likelihood with respect to gamma
}
\description{
This function calculates the gradient and Hessian of the log-likelihood with
respect to gamma
}
\usage{
grad_hess_L_gamma(Y, X, eta, gamma, I, J)
}
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
}
  \item{eta}{
Initial eta vector
}
  \item{gamma}{
Initial gamma vector
}
  \item{I}{
Number of conditions
}
  \item{J}{
Number of replications
}
}

\value{
\item{grad_L_gamma}{
Vector of the gradient of L with respect to gamma
}
\item{hess_L_gamma}{
Matrix of the Hessian of L with respect to gamma
}
}
\references{
M. Gomtsyan et al. "Variable selection in sparse multivariate GLARMA models: Application to germination control by environment", arXiv:2208.14721
}
\author{
Marina Gomtsyan

Maintainer: Marina Gomtsyan <marina.gomtsyan@agroparistech.fr>
}
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
eta_glm_mat_0 = matrix(0,ncol=T,nrow=I)
for (t in 1:T)
{
  result_glm_0 = glm(Y[,t]~X-1,family=poisson(link='log'))
  eta_glm_mat_0[,t]=as.numeric(result_glm_0$coefficients)
}
eta_0 = round(as.numeric(t(eta_glm_mat_0)),digits=6)
result = grad_hess_L_gamma(Y, X, eta_0, gamma_0, I, J)
grad = result$grad_L_gamma
Hessian = result$hess_L_gamma
}
