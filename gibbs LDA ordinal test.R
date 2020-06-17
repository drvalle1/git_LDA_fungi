rm(list=ls(all=TRUE))
set.seed(6)
library('Rcpp')
# library('MCMCpack')

setwd('U:\\GIT_models\\git_LDA_fungi') 
source('gibbs LDA ordinal aux.R')
source('gibbs LDA ordinal function.R')
sourceCpp('LDA_ordinal_rcpp.cpp')

thetas.p.lprob=read.csv('potential thetas lprob.csv',as.is=T)$x
thetas.pot=data.matrix(read.csv('potential thetas.csv',as.is=T))
dat=data.matrix(read.csv('fake data3.csv',as.is=T))
ngibbs=1000
ncomm.max=ncol(thetas.pot); ncomm.max
prop.burn=0.9
res=LDA_ordinal(dat=dat,ncomm.max=ncomm.max,ngibbs=ngibbs,prop.burn=prop.burn,
                thetas.p.lprob=thetas.p.lprob,thetas.pot=thetas.pot)

#---------------------------------------------------
nloc=nrow(dat)
nspp=ncol(dat)

#look at logl
seq1=1:ngibbs
seq1=(ngibbs*0.8):ngibbs
plot(res$logl[seq1],type='l')

#look at theta
theta.estim=matrix(colMeans(res$theta),nloc,ncomm.max)
boxplot(theta.estim)

#black, red, green, blue, cyan
seq.comm=1:5
# seq.comm=c(4,2,1)
theta.estim1=theta.estim[,seq.comm]
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:length(seq.comm)){
  lines(1:nloc,theta.estim1[,i],col=i)
}

#look at phi
phi.estim=matrix(colMeans(res$phi),ncomm.max,nspp)[seq.comm,]
rango=range(phi.estim,phi.true)
plot(phi.true,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango)

#look at breaks
break2=break1.true[-c(1,length(break1.true))]
break.estim=colMeans(res$break1)
rango=range(break.estim,break2)
plot(break2,break.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')
lines(rango,rango)
