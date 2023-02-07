
crit_CV = function(theta,sigma,X,y) {
  #Compute the log predictive probability LOO criterion, see Rasmussen 06
  #for the exponential covariance model in dimension 1
  #makes use of the explicit formula using Markovian property
  #
  #INPUTS
  #theta, sigma: covariance parameters, for covariance function sigma^2 exp(- theta t )
  #X,y: vectors of observation points and observations
  #
  #RETURNS
  #log predictive probability LOO criterion
  #
  #Preparations
  n = length(y)
  iX = sort(X,index.return=TRUE,decreasing=FALSE)$ix
  X = X[iX]
  y = y[iX]
  Delta2 = X[2] - X[1]
  Deltan = X[n] - X[n-1]
  Delta2_nm1 =  X[2:(n-1)] - X[ 1:(n-2) ]
  Delta3_n =  X[3:n] - X[ 2:(n-1) ]
  #
  #The log part of the criterion
  crit_log = log(1-exp(-2*theta*Delta2)) +log(1-exp(-2*theta*Deltan))+n*log(sigma^2)- sum( log( 1/(1-exp(-2*theta*Delta2_nm1))      + exp(-2*theta*Delta3_n)/(1-exp(-2*theta*Delta3_n))) )
  #
  #The data-based part of the criterion
  crit_err = ((y[1] - exp(-theta*Delta2)*y[2])^2) / (sigma^2*(1-exp(-2*theta*Delta2)))
  crit_err=crit_err +  ((y[n] - exp(-theta*Deltan)*y[n-1])^2)/(sigma^2*(1-exp(-2*theta*Deltan)))
  crit_err=crit_err + (1/sigma^2) * sum(  (1/(1-exp(-2*theta*Delta2_nm1)) + exp(-2*theta*Delta3_n)/(1 - exp(-2*theta*Delta3_n)) ) *( y[2:(n-1)] - ( (exp(-theta*Delta2_nm1)/(1-exp(-2*theta*Delta2_nm1))*y[1:(n-2)] + exp(-theta*Delta3_n)/(1-exp(-2*theta*Delta3_n))*y[3:n]) ) /(1/(1-exp(-2*theta*Delta2_nm1)) + exp(-2*theta*Delta3_n)/(1 - exp(-2*theta*Delta3_n)) ) )^2)
  #
  #Returns of the result
  crit_log + crit_err
}

crit_CV_matrix = function(theta,sigma,X,y) {
  #Compute the log predictive probability LOO criterion, see Rasmussen 06
  #for the exponential covariance model in dimension 1
  #uses the general matrix-based formula
  #
  #INPUTS
  #theta, sigma: covariance parameters, for covariance function sigma^2 exp(- theta t )
  #X,y: vectors of observation points and observations
  #
  #RETURNS
  #log predictive probability LOO criterion
  #
  n = length(y)
  R = sigma^2 * exp( -  theta*abs(outer(X,X,'-')) )
  errloo = solve(diag(diag(solve(R))))%*% solve(R) %*% matrix(nrow = n,ncol=1,data=y)
  sd2loo = 1/diag(solve(R))
  sum(log(sd2loo)) + sum( (errloo^2)/sd2loo  )
}


estim_CV = function(X,y,lim_theta,lim_sigma,n_descent,n_iter) {
  #return the optimizers of the loo log predictive probability
  #optimization by best among gradient descents
  #
  #INPUTS
  #X,y: vectors of observation points and observations
  #lim_theta: vector of minimum and maximum values for optimization w.r.t. theta
  #lim_sigma: same for sigma
  #n_descent: number of gradient descent (the one with best final score is kept)
  #n_iter: maximal number of iterations for the gradient descents
  #
  #RETURNS
  #a list with
  #hat_theta: the theta corresponding to the minimizer of the loo log predictive probability
              #w.r.t (theta,sigma)
  #hat_sigma: same for sigma
  #
  m_hat_par_descents = matrix(nrow=n_descent,ncol=2,data=0)
  v_crit_descents = 0*(1:n_descent)
  #
  to_min <-function( par, X, y )  {   # the function to minimize; par = (theta,sigma)
    crit_CV(par[1],par[2],X,y)
  }
  #
  for (i in 1:n_descent)  {  #Loop on the gradient descents
    par_init = c( runif(n=1,min=lim_theta[1],max=lim_theta[2]) , runif(n=1,min=lim_sigma[1],max=lim_sigma[2]))
    resOpt = optim( par=par_init , f=to_min , method="L-BFGS-B" , 
                    lower = c(lim_theta[1],lim_sigma[1]) , upper = c(lim_theta[2],lim_sigma[2]),
                    control = list(maxit=n_iter),X=X,y=y)
    m_hat_par_descents[i,] = resOpt$par
    v_crit_descents[i] = resOpt$value
  }
  ind_opt = which.min(v_crit_descents)
  m_hat_par_descents[ind_opt,]
}

estim_CV_random = function(X,y,lim_theta,lim_sigma,n_sampling) {
  #return the optimizers of the loo log predictive probability
  #optimization by best among random sampling of parameters
  #
  #INPUTS
  #X,y: vectors of observation points and observations
  #lim_theta: vector of minimum and maximum values for optimization w.r.t. theta
  #lim_sigma: same for sigma
  #n_sampling: number of randomly sampled points (the one with best final score is kept)
  #
  #RETURNS
  #a list with
  #hat_theta: the theta corresponding to the minimizer of the loo log predictive probability
              #w.r.t (theta,sigma)
  #hat_sigma: same for sigma
  #
  v_sampling_theta = runif(n=n_sampling,min = lim_theta[1],max=lim_theta[2])
  v_sampling_sigma = runif(n=n_sampling,min = lim_sigma[1],max=lim_sigma[2])
  m_par_sampling = cbind(v_sampling_theta,v_sampling_sigma)
  v_crit_sampling = 0*(1:n_descent)
  #
  for (i in 1:n_sampling)  {  #Loop on the criterion evaluations
    v_crit_sampling[i] = crit_CV(m_par_sampling[i,1],m_par_sampling[i,2],X,y)
  }
  ind_opt = which.min(v_crit_sampling)
  m_par_sampling[ind_opt,]
}

tau_n = function(X) {
  #Compute the tau_n (cf our paper) as a function of the vector X
  #     of observation points
  X = sort(X)
  n=length(X)
  delta_im1 = X[2:(n-2)] - X[1:(n-3)]
  delta_i = X[3:(n-1)] - X[2:(n-2)]
  delta_ip1 = X[4:n] - X[3:(n-1)]
  sum_of_square = sum( ( (delta_ip1/(delta_i+delta_ip1)) +
                          (delta_im1/(delta_i+delta_im1)) )^2 ) 
  sum_2 = sum( (delta_i*delta_ip1)/((delta_i+delta_ip1)^2) )
  sqrt((2/n)*(sum_of_square+2*sum_2))
}

maximal_variance_design = function(n) {
  #Compute the vector X of size n giving the maximum possible
  #variance for the CV estimator (cf our paper)
  X = 0*(1:(2*n+2))
  gamma_n = 1/n
  X_odd = seq(from=0,by=2/n,length=n)
  X_even = (2*gamma_n)/n + seq(from=0,by=2/n,length=n)
  X[seq(from=1,by=2,to=n)] = X_odd[1:floor((n+1)/2)]
  X[seq(from=2,by=2,to=n)] = X_even[1:floor(n/2)]
  X[n] = 1
  X[1:n]
}

estim_CV_theta_known = function(X,y,lim_theta,lim_sigma,n_descent,n_iter) {
  #return the optimizers of the loo log predictive probability
  #optimization by best among gradient descents
  #
  #INPUTS
  #X,y: vectors of observation points and observations
  #lim_theta: vector of minimum and maximum values for optimization w.r.t. theta
  #        minimal and maximal must be equal here
  #lim_sigma: same for sigma
  #n_descent: number of gradient descent (the one with best final score is kept)
  #n_iter: maximal number of iterations for the gradient descents
  #
  #RETURNS
  #a list with
  #hat_theta: the theta corresponding to the minimizer of the loo log predictive probability
  #w.r.t (theta,sigma)
  #hat_sigma: same for sigma
  #
  m_hat_par_descents = matrix(nrow=n_descent,ncol=2,data=0)
  v_crit_descents = 0*(1:n_descent)
  #
  to_min <-function( par,theta, X, y )  {   # the function to minimize for theta known; par = (sigma)
    crit_CV(theta,par[1],X,y)
  }
  #
  for (i in 1:n_descent)  {  #Loop on the gradient descents
    par_init =  runif(n=1,min=lim_sigma[1],max=lim_sigma[2])
    resOpt = optim( par=par_init , f=to_min , method="L-BFGS-B" , 
                    lower = lim_sigma[1] , upper = lim_sigma[2],
                    control = list(maxit=n_iter),theta=lim_theta[1],
                    X=X,y=y)
    m_hat_par_descents[i,] = c(lim_theta[1],resOpt$par)
    v_crit_descents[i] = resOpt$value
  }
  ind_opt = which.min(v_crit_descents)
  m_hat_par_descents[ind_opt,]
}

estim_CV_sigma_known = function(X,y,lim_theta,lim_sigma,n_descent,n_iter) {
  #return the optimizers of the loo log predictive probability
  #optimization by best among gradient descents
  #
  #INPUTS
  #X,y: vectors of observation points and observations
  #lim_theta: vector of minimum and maximum values for optimization w.r.t. theta
  #lim_sigma: same for sigma. Minimal and maximal must be equal here
  #n_descent: number of gradient descent (the one with best final score is kept)
  #n_iter: maximal number of iterations for the gradient descents
  #
  #RETURNS
  #a list with
  #hat_theta: the theta corresponding to the minimizer of the loo log predictive probability
  #w.r.t (theta,sigma)
  #hat_sigma: same for sigma
  #
  m_hat_par_descents = matrix(nrow=n_descent,ncol=2,data=0)
  v_crit_descents = 0*(1:n_descent)
  #
  to_min <-function( par,sigma, X, y )  {   # the function to minimize for theta known; par = (sigma)
    crit_CV(par[1],sigma,X,y)
  }
  #
  for (i in 1:n_descent)  {  #Loop on the gradient descents
    par_init =  runif(n=1,min=lim_theta[1],max=lim_theta[2])
    resOpt = optim( par=par_init , f=to_min , method="L-BFGS-B" , 
                    lower = lim_theta[1] , upper = lim_theta[2],
                    control = list(maxit=n_iter),sigma=lim_sigma[1],
                    X=X,y=y)
    m_hat_par_descents[i,] = c(resOpt$par,lim_sigma[1])
    v_crit_descents[i] = resOpt$value
  }
  ind_opt = which.min(v_crit_descents)
  m_hat_par_descents[ind_opt,]
}


minimal_variance_design = function(n,alpha) {
  #Compute the vector X of size n giving the minimum possible
  #variance for the CV estimator (cf our paper)
  vdelta = 0*(1:(n))
  for (i in (floor(n^alpha)+1):(n)) {
    vdelta[i] = 1/factorial(i)
  }
  rn = sum(vdelta)
  for (i in 2:floor(n^alpha)) {
    vdelta[i] = (1-rn) / (floor(n^alpha) -1)
  }
  X = 0*(1:n)
  for (i in 2:n) {
    X[i] = X[i-1] + vdelta[i]
  }
  X
}
