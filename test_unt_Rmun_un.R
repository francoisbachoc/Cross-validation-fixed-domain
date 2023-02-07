rm(list=ls(all=TRUE))
set.seed(1)


n=1000
X = runif(n=n)
sigma = 1
theta=3.4

R = sigma^2 * exp( -  theta*abs(outer(X,X,'-')) )
un = seq(from=1,to=1,length=n)

info_exact = t(un)%*%solve(R)%*%un
info_assympt = 1+0.5*theta

cat("info exacte: ",info_exact, "\n")
cat("info assymptotique: ",info_assympt, "\n")


