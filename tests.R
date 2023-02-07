rm(list=ls(all=TRUE))
source("functions.R")
set.seed(1)

#######################################################
#Test of the explicit formula with Markovian properties
#######################################################

X = c(0.12,0.3,-0.4,12,-5,0.76,0.18,0,-0.92,64)
y = c(1,2.3,-1,0.33,5.23,-3.33,-34,32,-3,-3.5435)
sigma = 5.56
theta=0.584

cat("evaluation loo criterion with Markov: ", crit_CV(theta,sigma,X,y) , "\n")
cat("evaluation loo criterion with matrix: ", crit_CV_matrix(theta,sigma,X,y) , "\n")

#######################################################
#Test of the gradient descent
#######################################################

X = runif(n=100)
y = cos(3*exp(2*X) - 2*X^2 + 1) - X
lim_theta = c(0.01,1)
lim_sigma = c(0.01,10)
n_descent=30
n_iter=10000
n_sampling=10000

hat_par_gradient = estim_CV(X,y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_descent=n_descent,
                            n_iter=n_iter)
opt_crit_gradient = crit_CV(theta=hat_par_gradient[1],sigma=hat_par_gradient[2],X=X,y=y)
hat_par_sampling = estim_CV_random(X,y,lim_theta=lim_theta,lim_sigma=lim_sigma,n_sampling=n_sampling)
opt_crit_sampling = crit_CV(theta=hat_par_sampling[1],sigma=hat_par_sampling[2],X=X,y=y)


cat("estimated theta sigma with gradient-based minimization",hat_par_gradient,"\n")
cat("optimum criterion with gradient-based minimization",opt_crit_gradient,"\n")
cat("estimated theta sigma with sampling-based minimization",hat_par_sampling,"\n")
cat("optimum criterion with sampling-based minimization",opt_crit_sampling,"\n")

