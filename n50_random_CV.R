rm(list=ls(all=TRUE))
source("functions.R")
set.seed(1)

#######################################################
#simulation-dependent parameters (only this changes accross different simulations)
#######################################################
n=50
sX="random"
sigma_0=1
theta_0=3
estimator="CV"

#######################################################
#Simulation-independent parameters (this does not change accross different simulations) 
#######################################################
lim_theta = c(0.1,10)
lim_sigma = c(0.3,30)
n_descent=5
n_iter=1000
nmc=2000
period_message=100

#######################################################
#Simulation (this does not change accross different simulations) 
#######################################################
if (sX == "regular") {
  X = (1:n)/n
}
if (sX == "maximal") {
  X = maximal_variance_design(n)
}
if (sX == "random") {
  X = sort(runif(n))
}
m_hat_par = matrix(nrow=nmc,ncol=2,data=0)
R = sigma_0^2 * exp( -  theta_0*abs(outer(X,X,'-')) )
cR = chol(R)
for (i in 1:nmc) {
  if (  floor(i/period_message) == i/period_message  ) {
    cat("Estimation",i,"\n")
  }
  y = t(cR)%*%rnorm(n)
  m_hat_par[i,] = estim_CV(X=X,y=y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_descent=n_descent,
                              n_iter=n_iter)
}

######################################################
#Save of the results
#####################################################
name_ = paste0("n",n,"_X",sX,"_",estimator)
name_file = paste0(name_,".Rdata")
save(m_hat_par,n,theta_0,sigma_0,X,file = name_file)
