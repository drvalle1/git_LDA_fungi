# rm(list=ls(all=TRUE))
set.seed(4)
library('Rcpp')

setwd('U:\\GIT_models\\git_LDA_fungi') 
source('gibbs LDA fungi functions.R')
sourceCpp('aux1.cpp')
dat=data.matrix(read.csv('fake data5.csv',as.is=T))

nloc=nrow(dat)
nspp=ncol(dat)
ncommun=10
gamma=0.1

#set initial values
vlk=cbind(matrix(0.5,nloc,ncommun-1),1)
phi=matrix(0,ncommun,nspp)

ones.nloc=rep(1,nloc)
theta=convertVtoTheta(vlk,ones.nloc)
ones.nspp=rep(1,nspp)

#things for MH algorithm
n.vlk=(ncommun-1)*nloc
lo.vlk=rep(0,n.vlk)
hi.vlk=rep(1,n.vlk)
n.phi=ncommun*nspp

#gibbs stuff
ngibbs=10000
accept.output=50
vec.theta=matrix(NA,ngibbs,nloc*ncommun)
vec.phi=matrix(NA,ngibbs,ncommun*nspp)

jump1=list(vlk=matrix(1,nloc,ncommun-1),
           phi=matrix(1,ncommun,nspp))
accept1=list(vlk=matrix(0,nloc,ncommun-1),
             phi=matrix(0,ncommun,nspp))
param=list(theta=theta,phi=phi,vlk=vlk)

#core MCMC algorithm
for (i in 1:ngibbs){
  print(i)
  tmp=sample.vlk(param,jump1$vlk)
  accept1$vlk=accept1$vlk+tmp$accept
  param$vlk=tmp$vlk
  param$theta=convertVtoTheta(param$vlk,ones.nloc)
  # param$theta=theta.true
    
  tmp=sample.phi(param,jump1$phi)
  accept1$phi=accept1$phi+tmp$accept
  param$phi=tmp$phi
  # param$phi=phi.true
  
  #adaptive piece
  if (i%%accept.output==0 & i<1000){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
  }
  
  #store results
  vec.theta[i,]=param$theta
  vec.phi[i,]=param$phi
}

#---------------------------------------------------
seq1=(ngibbs*0.9):ngibbs
theta.estim=matrix(colMeans(vec.theta[seq1,]),nloc,ncommun)
boxplot(theta.estim)

seq.comm=1:5
seq.comm=c(2,5,3,1,4)
theta.estim1=theta.estim[,seq.comm]
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:5){
  lines(1:nloc,theta.estim1[,i],col=i)
}

phi.estim=matrix(colMeans(vec.phi[seq1,]),ncommun,nspp)[seq.comm,]
rango=range(phi.estim,phi.true)
plot(phi.true,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango)
