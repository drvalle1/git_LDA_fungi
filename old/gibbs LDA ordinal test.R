rm(list=ls(all=TRUE))
set.seed(4)
library('Rcpp')

setwd('U:\\GIT_models\\git_LDA_fungi') 
source('gibbs LDA ordinal aux.R')
source('gibbs LDA ordinal function.R')
sourceCpp('LDA_ordinal_rcpp.cpp')

dat=data.matrix(read.csv('fake data5.csv',as.is=T))
ngibbs=10000
ncommun=20
res=LDA_ordinal(dat=dat,ncomm.max=ncommun,ngibbs=ngibbs,prop.burn=0.9,gamma=0.1)

#---------------------------------------------------
nloc=nrow(dat)
nspp=ncol(dat)

#look at logl
seq1=(ngibbs*0.2):ngibbs
plot(res$logl[seq1],type='l')

#look at theta
theta.estim=matrix(colMeans(res$theta),nloc,ncommun)
boxplot(theta.estim)

seq.comm=1:5
seq.comm=c(1,2,5,3,4)
theta.estim1=theta.estim[,seq.comm]
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:5){
  lines(1:nloc,theta.estim1[,i],col=i)
}

#look at phi
phi.estim=matrix(colMeans(res$phi),ncommun,nspp)[seq.comm,]
rango=range(phi.estim,phi.true)
plot(phi.true,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango)

#look at breaks
break2=break1.true[-c(1,length(break1.true))]
break.estim=colMeans(res$break1)
rango=range(break.estim,break2)
plot(break2,break.estim,xlim=rango,ylim=rango)
lines(rango,rango)
