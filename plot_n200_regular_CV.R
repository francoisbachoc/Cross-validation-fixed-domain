rm(list=ls(all=TRUE))
source("functions.R")
set.seed(1)


#######################################################
#simulation-dependent parameters (only this changes accross different simulations)
#######################################################
name_ = "n200_Xregular_CV"


#######################################################
#Figure (this does not change accross different simulations) 
#######################################################
load(paste0(name_,".Rdata"))
to_plot = sqrt(n)*(1/(theta_0*sigma_0^2))*
  (m_hat_par[,1]*m_hat_par[,2]^2 - theta_0*sigma_0^2)

pdf(paste0(name_,".pdf"), useDingbats=F )
opar <- par(lwd=4)
xlim=c(min(to_plot),max(to_plot))
hist(to_plot,breaks=30,xlab="",ylab="",main="n=200, regular",
     freq=FALSE,font.axis=2,cex.axis=2,cex=3.5,cex.lab=2,font.lab=2,cex.main=2,
     font.main=2,xlim=c(-6,10),ylim=c(0,0.4))
x_plot = seq(from=xlim[1],to=xlim[2],length=1000)
points(x=x_plot,dnorm(x_plot,sd=tau_n(X)),
       col="red",lwd=4,type="l")
par(opar)
dev.off()  





