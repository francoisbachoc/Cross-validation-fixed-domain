rm(list=ls(all=TRUE))
source("functions.R")
set.seed(1)

#######################################################
#simulation-dependent parameters (only this changes accross different simulations)
#######################################################
n=12
alpha=0.5
X_min = minimal_variance_design(n,alpha)
pdf(file="minimal_design.pdf")
plot(X_min,0*X_min+1,type="p",main="",xlab="",ylab="",,pch=21,cex=1.2,bg = "blue",yaxt='n',cex.axis=2.2)
dev.off()

X_eq = seq(from=0,to=1,length=n)
pdf(file="regular_design.pdf")
plot(X_eq,0*X_min+1,type="p",main="",xlab="",ylab="",,pch=21,cex=1.2,bg = "blue",yaxt='n',cex.axis=2.2)
dev.off()

gamma=0.5
X_max = maximal_variance_design(n)
pdf(file="maximal_design.pdf")
plot(X_max,0*X_min+1,type="p",main="",xlab="",ylab="",,pch=21,cex=1.2,bg = "blue",yaxt='n',cex.axis=2.2)
dev.off()



