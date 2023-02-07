rm(list=ls(all=TRUE))
source("functions.R")
set.seed(1)

#######################################################
#simulation-dependent parameters (only this changes accross different simulations)
#######################################################
n=12
sX="maximal"
sigma_0=1
theta_0=3
estimator="CV"
remark = ""

#######################################################
#Simulation-independent parameters (this does not change accross different simulations) 
#######################################################
lim_theta = c(0.1,10)
lim_sigma = c(0.3,30)
if (remark == "theta_known") {
  lim_theta = c(theta_0,theta_0)
}
if (remark == "sigma_known") {
  lim_sigma = c(sigma_0,sigma_0)
}
n_descent=5
n_iter=1000
nmc=2000
period_message=100
ratio=3
alpha=0.5

#######################################################
#Simulation (this does not change accross different simulations) 
#######################################################
if (sX == "regular") {
  X = seq(from=0,to=1,length=n)
}
if (sX == "maximal") {
  X = maximal_variance_design(n)
}
if (sX == "random") {
  X = sort(runif(n))
}
if (sX == "minimal") {
  X = minimal_variance_design(n,alpha)
}
m_hat_par = matrix(nrow=nmc,ncol=2,data=0)
R = sigma_0^2 * exp( -  theta_0*abs(outer(X,X,'-')) )
cR = chol(R)
for (i in 1:nmc) {
  if (  floor(i/period_message) == i/period_message  ) {
    cat("Estimation",i,"\n")
  }
  y = t(cR)%*%rnorm(n)
  if (remark == "") {
    m_hat_par[i,] = estim_CV(X=X,y=y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_descent=n_descent,
                             n_iter=n_iter)
  }
  if (remark == "theta_known")  {
    m_hat_par[i,] = estim_CV_theta_known(X=X,y=y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_descent=n_descent,
                             n_iter=n_iter)
  }
  if (remark == "sigma_known")  {
    m_hat_par[i,] = estim_CV_sigma_known(X=X,y=y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_descent=n_descent,
                                         n_iter=n_iter)
  }
}

######################################################
#Save of the results
#####################################################
name_ = paste0("n",n,"_X",sX,"_",estimator)
name_file = paste0(name_,".Rdata")
save(m_hat_par,n,theta_0,sigma_0,X,file = name_file)
