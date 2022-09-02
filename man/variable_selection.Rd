\name{variable_selection}
\alias{variable_selection}
\title{
Variable selection
}
\description{
This function performs variable selection, estimates a new vector eta and a new vector gamma
}
\usage{
variable_selection(Y, X, gamma, k_max = 1, n_iter = 100, 
method = "min", nb_rep_ss = 1000, threshold = 0.6)
}
\arguments{
  \item{Y}{
Observation matrix
}
  \item{X}{
Design matrix
}
  \item{gamma}{
Initial gamma vector
}
  \item{k_max}{
Number of iteration to repeat the whole algorithm
}
  \item{n_iter}{
Number of iteration for Newton-Raphson algorithm
}
  \item{method}{
Stability selection method: "min" or "cv". In "min" the smallest lambda is chosen, in "cv" cross-validation lambda is chosen for stability selection. The default is "min"
}
  \item{nb_rep_ss}{
Number of replications in stability selection step. The default is 1000
}
  \item{threshold}{
Threshold for stability selection. The default is 0.9
}
}
\value{
\item{estim_active}{Vector of stimated active coefficients}
\item{eta_est}{Vector of estimated eta values}
\item{gamma_est}{Vector of estimated gamma values}
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
result=variable_selection(Y, X, gamma_0, k_max=1, 
n_iter=100, method="min", nb_rep_ss=1000, threshold=0.6)
estim_active = result$estim_active
eta_est = result$eta_est
gamma_est = result$gamma_est
}
