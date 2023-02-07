rm(list=ls(all=TRUE))
set.seed(1)



X = c(0.2,0.33,0.43,0.65,0.77,0.99,1)
y = c(1,2.3,-1,0.33,5.23,-3.33,-34)
sigma = 1
theta=0.34

n = length(y)

R = sigma^2 * exp( -  theta*abs(outer(X,X,'-')) )

errloo = solve(diag(diag(solve(R))))%*% solve(R) %*% matrix(nrow = n,ncol=1,data=y)
sd2loo = 1/diag(solve(R))

crit1_log = sum(log(sd2loo)) 
crit1_err =  sum( (errloo^2)/sd2loo  )

cat("first evaluation loo log criterion: ",crit1_log, "\n")
cat("first evaluation loo err criterion: ",crit1_err, "\n")

Delta2 = X[2] - X[1]
Deltan = X[n] - X[n-1]
Delta2_nm1 =  X[2:(n-1)] - X[ 1:(n-2) ]
Delta3_n =  X[3:n] - X[ 2:(n-1) ]

crit2_log = log(1-exp(-2*theta*Delta2)) +log(1-exp(-2*theta*Deltan))+n*log(sigma^2)- sum( log( 1/(1-exp(-2*theta*Delta2_nm1))      + exp(-2*theta*Delta3_n)/(1-exp(-2*theta*Delta3_n))) )

crit2_err = ((y[1] - exp(-theta*Delta2)*y[2])^2) / (sigma^2*(1-exp(-2*theta*Delta2)))
crit2_err=crit2_err +  ((y[n] - exp(-theta*Deltan)*y[n-1])^2)/(sigma^2*(1-exp(-2*theta*Deltan)))
crit2_err=crit2_err + (1/sigma^2) * sum(  (1/(1-exp(-2*theta*Delta2_nm1)) + exp(-2*theta*Delta3_n)/(1 - exp(-2*theta*Delta3_n)) ) *( y[2:(n-1)] - ( (exp(-theta*Delta2_nm1)/(1-exp(-2*theta*Delta2_nm1))*y[1:(n-2)] + exp(-theta*Delta3_n)/(1-exp(-2*theta*Delta3_n))*y[3:n]) ) /(1/(1-exp(-2*theta*Delta2_nm1)) + exp(-2*theta*Delta3_n)/(1 - exp(-2*theta*Delta3_n)) ) )^2)

cat("second evaluation loo log criterion: ",crit2_log, "\n")
cat("second evaluation loo err criterion: ",crit2_err, "\n")

