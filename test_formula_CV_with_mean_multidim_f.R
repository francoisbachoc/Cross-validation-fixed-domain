rm(list=ls(all=TRUE))
set.seed(1)


i = 3
n=7
y = runif(n=n)
X = runif(n=n)
p=4
mF = matrix(nrow=n,ncol=p,data=runif(n*p))
beta_0 = runif(n=p)
z = y+ mF%*%beta_0
sigma = 1
theta=0.14

R = sigma^2 * exp( -  theta*abs(outer(X,X,'-')) )
iR = solve(R)
Ri = R[-i,-i]
iRi = solve(Ri)
ri = R[-i,i]
Q = iR - iR%*%mF%*%solve(t(mF)%*%iR%*%mF)%*%t(mF)%*%iR

hatbeta_i = solve(t(mF[-i,])%*%iRi%*%mF[-i,])%*%t(mF[-i,])%*%iRi%*%z[-i]

predi_1 = t(mF[i,])%*%hatbeta_i + t(ri)%*%iRi%*%( z[-i] - mF[-i,] %*%hatbeta_i)
sd2_i = 1/Q[i,i]

critloo_log1 = log(sd2_i) 
critloo_err1 = ((z[i] - predi_1)^2)/sd2_i

critloo_log_y = log(1/(iR[i,i])) 
critloo_err_y = ( ((1/(iR[i,i]))*(iR%*%y)[i])^2 )/(1/(iR[i,i]))
epsi =  (t(mF[i,])-t(ri)%*%iRi%*%mF[-i,]  )%*%(beta_0 - hatbeta_i )
epsi_bar = (iR%*%mF%*%solve(t(mF)%*%iR%*%mF)%*%t(mF)%*%iR)[i,i]

rest_log = -log( iR[i,i] - epsi_bar ) + log(iR[i,i]) 
rest_err = iR[i,i]*epsi^2 + 2 * iR[i,i] * epsi * (y[i] - t(ri)%*%iRi%*%y[-i]  )  - epsi_bar * (y[i] - t(ri)%*%iRi%*%y[-i]  + epsi)^2 
  


cat("first evaluation loo log criterion: ",critloo_log1, "\n")
cat("second evaluation loo log criterion: ",critloo_log_y + rest_log, "\n")
cat("first evaluation loo err criterion: ",critloo_err1, "\n")
cat("second evaluation loo err criterion: ",critloo_err_y + rest_err, "\n")

