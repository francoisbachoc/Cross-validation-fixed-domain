rm(list=ls(all=TRUE))
set.seed(1)


i = 3
X = c(0.2,0.33,0.43,0.65,0.77,0.99,1)
y = c(1,2.3,-1,0.33,5.23,-3.33,-34)
beta_0 = 2
z = y+ beta_0
sigma = 1
theta=0.14
one = seq(from=1,to=1,length=length(X))

R = sigma^2 * exp( -  theta*abs(outer(X,X,'-')) )
iR = solve(R)
Ri = R[-i,-i]
iRi = solve(Ri)
ri = R[-i,i]
Q = iR - iR%*%one%*%solve(t(one)%*%iR%*%one)%*%t(one)%*%iR

hatbeta_i = solve(t(one[-i])%*%iRi%*%one[-i])%*%t(one[-i])%*%iRi%*%z[-i]

predi_1 = hatbeta_i + t(ri)%*%iRi%*%( z[-i] - hatbeta_i )
sd2_i = 1/Q[i,i]

critloo_log1 = log(sd2_i) 
critloo_err1 = ((z[i] - predi_1)^2)/sd2_i

critloo_log_y = log(1/(iR[i,i])) 
critloo_err_y = ( ((1/(iR[i,i]))*(iR%*%y)[i])^2 )/(1/(iR[i,i]))
epsi = (beta_0 - hatbeta_i ) * (1-t(ri)%*%iRi%*%one[-i]  )
epsi_bar = (iR%*%one%*%solve(t(one)%*%iR%*%one)%*%t(one)%*%iR)[i,i]

rest_log = -log( iR[i,i] - epsi_bar ) + log(iR[i,i]) 
rest_err = iR[i,i]*epsi^2 + 2 * iR[i,i] * epsi * (y[i] - t(ri)%*%iRi%*%y[-i]  )  - epsi_bar * (y[i] - t(ri)%*%iRi%*%y[-i]  + epsi)^2 
  


cat("first evaluation loo log criterion: ",critloo_log1, "\n")
cat("second evaluation loo log criterion: ",critloo_log_y + rest_log, "\n")
cat("first evaluation loo err criterion: ",critloo_err1, "\n")
cat("second evaluation loo err criterion: ",critloo_err_y + rest_err, "\n")

