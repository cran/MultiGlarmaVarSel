## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)
library(MultiGlarmaVarSel)
library(ggplot2)
library(formatR)
set.seed(12345)

## ----Y------------------------------------------------------------------------
data(Y)

## ----dimensions---------------------------------------------------------------
I=3
J=100
T=dim(Y)[2]
q=1

## ----gamma--------------------------------------------------------------------
gamma = matrix(c(0.5), nrow = 1, ncol = q)

## ----eta----------------------------------------------------------------------
active=c(2, 24, 33, 34, 37, 41)
non_active = setdiff(1:(I*T),active)
eta_true=rep(0,(I*T))
eta_true[active]=c(0.63, 2.62, 1.69, 2.27, 0.72, 0.41)

## ----X------------------------------------------------------------------------
X=matrix(0,nrow=(I*J),ncol=I)
for (i in 1:I)
{
  X[((i-1)*J+1):(i*J),i]=rep(1,J)
}

## ----initial------------------------------------------------------------------
gamma_0 = matrix(0, nrow = 1, ncol = q)
eta_glm_mat_0 = matrix(0,ncol=T,nrow=I)
for (t in 1:T)
{
  result_glm_0 = glm(Y[,t]~X-1,family=poisson(link='log'))
  eta_glm_mat_0[,t]=as.numeric(result_glm_0$coefficients)
}
eta_0 = round(as.numeric(t(eta_glm_mat_0)),digits=6)

## ----gamma_est----------------------------------------------------------------
gamma_est=NR_gamma(Y, X, eta_0, gamma_0, I, J, n_iter = 100)
cat("Estimated gamma: ", gamma_est, "\n")

## ----variableSelection, tidy=TRUE, tidy.opts=list(width.cutoff=60)------------
result=variable_selection(Y, X, gamma_est, k_max=1, n_iter=100, method="min", nb_rep_ss=1000, threshold=0.67)
estim_active = result$estim_active
eta_est = result$eta_est
gamma_est = result$gamma_est

## ----output, echo=FALSE, eval=TRUE--------------------------------------------
eta_data = data.frame(eta_est)
eta_data$Condition = c(rep(1,T), rep(2,T), rep(3,T))
eta_data$t = c(rep(seq(1,T,1),I))
colnames(eta_data)[1] <- "eta"
eta_data = eta_data[eta_data$eta!=0,]
eta_true_data = data.frame(eta_true)
colnames(eta_true_data)[1] <- "eta"
eta_true_data$Condition = c(rep(1,T), rep(2,T), rep(3,T))
eta_true_data$t = c(rep(seq(1,T,1),I))
eta_true_data = eta_true_data[eta_true_data$eta !=0,]

eta_data = eta_data[,2:3]
eta_data_coef = noquote(apply(eta_data, 1, paste, collapse=","))
eta_data_coef_print = noquote(paste0("(",eta_data_coef,")"))
eta_true_data = eta_true_data[,2:3]
eta_true_coef = noquote(apply(eta_true_data, 1, paste, collapse=","))
eta_true_coef_print = noquote(paste0("(",eta_true_coef,")"))


cat("True active coefficient pairs: ", eta_true_coef_print, "\n")
cat("Estimated active coefficient pairs: ", eta_data_coef_print, "\n")
cat("True values of elements : ", eta_true[active], "\n")
cat("Estimated values of elements: ", round(eta_est[active], digits = 2), "\n")

## ----plot, fig.align="center", fig.width=5, fig.height=3.5, tidy=FALSE, tidy.opts=list(width.cutoff=52)----
#First, we make a dataset of estimated etas
eta_df = data.frame(eta_est)
eta_df$t = c(rep(seq(1,T,1),I))
eta_df$I = c(rep(1,T), rep(2,T), rep(3,T))
colnames(eta_df)[1] <- "eta"
eta_df = eta_df[eta_df$eta!=0,]
#Next, we make a dataset of true etas
eta_t_df = data.frame(eta_true)
colnames(eta_t_df)[1] <- "eta"
eta_t_df$I = c(rep(1,T), rep(2,T), rep(3,T))
eta_t_df$t = c(rep(seq(1,T,1),I))
eta_t_df = eta_t_df[eta_t_df$eta !=0,]
#Finally, we plot the result
plot = ggplot()+
  geom_point(data = eta_df, aes(x=t, y=I, color=eta), pch=20, size=3, stroke = 1)+
  geom_point(data= eta_t_df, aes(x=t, y=I, color=eta), pch=4, size=4.5)+
  scale_color_gradient2(name=expression(hat(eta)), midpoint=0, 
                        low="steelblue", mid = "white", high ="red")+
  theme_bw()+ylab('I')+xlab('T')+
  theme(legend.position = "bottom")+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.title = element_text(size = 15, face="bold"))+
  theme(legend.text = element_text(size = 7, color="black"))+
  scale_y_continuous(breaks=seq(1, I, 1), limits=c(0, I))+
  scale_x_continuous(breaks=c(1, seq(10, T, 10)), limits=c(0, T))+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.text=element_text(size=10, color="black"))+
  theme(axis.title=element_text(size=15,face="bold"))
plot

