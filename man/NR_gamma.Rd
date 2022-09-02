\name{NR_gamma}
\alias{NR_gamma}
\title{
Newton-Raphson method for estimation of gamma
}
\description{
This function estimates gamma with Newton-Raphson method
}
\usage{
NR_gamma(Y, X, eta, gamma, I, J, n_iter = 100)
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
  \item{n_iter}{
Number of iterations of the algorithm. Default=100
}
}
\value{
Estimated gamma vector
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
gamma_est=NR_gamma(Y, X, eta_0, gamma_0, I, J, n_iter = 100)
}

