rm(list=ls(all=TRUE))
set.seed(12)
library('Rcpp')
# library('MCMCpack')

setwd('U:\\GIT_models\\git_LDA_fungi') 
source('gibbs LDA ordinal aux.R')
source('gibbs LDA ordinal main function.R')
sourceCpp('LDA_ordinal_rcpp.cpp')
dat=data.matrix(read.csv('fake data.csv',as.is=T))
ngibbs=10000
ncomm=10
prop.burn=0.5
gamma1=0.1#0.01
res=LDA_ordinal(dat=dat,ncomm=ncomm,ngibbs=ngibbs,prop.burn=prop.burn,gamma1=gamma1)
#---------------------------------------------------
nloc=nrow(dat)
nspp=ncol(dat)

#look at logl
seq1=1:ngibbs
seq1=(ngibbs*0.7):ngibbs
plot(res$logl[seq1],type='l')

#look at theta
theta.estim=matrix(res$theta[nrow(res$theta),],nloc,ncomm)
boxplot(theta.estim)

#black, red, green, blue, cyan
# seq.comm=1:5
seq.comm=c(3,4,5,1,2)
theta.estim1=theta.estim[,seq.comm]
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:length(seq.comm)){
  lines(1:nloc,theta.estim1[,i],col=i)
  lines(1:nloc,theta.true[,i],col=i,lty=3)
}

#look at phi
phi.estim=matrix(colMeans(res$phi),ncomm,nspp)[seq.comm,]
rango=range(phi.estim,phi.true)
plot(phi.true,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango)

#look at breaks
break2=break1.true[-c(1,length(break1.true))]
break.estim=colMeans(res$break1)
rango=range(break.estim,break2)
plot(break2,break.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')

#look to see if breaks are moving
plot(res$break1[,2],type='l')
plot(res$break1[,4],type='l')