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
uni=sort(unique(as.numeric(dat)))
nuni=length(uni)

#set initial values
vlk=cbind(matrix(0.5,nloc,ncommun-1),1)
phi=matrix(0,ncommun,nspp)
break1=seq(from=0,to=0.5,length.out=nuni-1)
ones.nloc=rep(1,nloc)
theta=convertVtoTheta(vlk,ones.nloc)
ones.nspp=rep(1,nspp)

#things for MH algorithm
n.vlk=(ncommun-1)*nloc
lo.vlk=rep(0,n.vlk)
hi.vlk=rep(1,n.vlk)
n.phi=ncommun*nspp

#useful stuff
indicator=list()
n.indicator=rep(NA,nuni)
z=matrix(NA,nloc,nspp)
for (i in 1:nuni){
  cond=dat==uni[i]
  ind=which(cond)
  indicator[[i]]=ind
  z[ind]=i/100
  n.indicator[i]=length(ind)
}

#gibbs stuff
ngibbs=10000
accept.output=50
vec.theta=matrix(NA,ngibbs,nloc*ncommun)
vec.phi=matrix(NA,ngibbs,ncommun*nspp)
vec.break1=matrix(NA,ngibbs,nuni-1)
vec.logl=matrix(NA,ngibbs,1)

jump1=list(vlk=matrix(1,nloc,ncommun-1),break1=rep(1,nuni-1))
accept1=list(vlk=matrix(0,nloc,ncommun-1),break1=rep(0,nuni-1))
param=list(theta=theta,phi=phi,vlk=vlk,break1=break1,z=z)

#core MCMC algorithm
options(warn=2)
for (i in 1:ngibbs){
  print(i)
  tmp=sample.vlk(param,jump1$vlk)
  accept1$vlk=accept1$vlk+tmp$accept
  param$vlk=tmp$vlk
  param$theta=convertVtoTheta(param$vlk,ones.nloc)
  # param$theta=theta.true
    
  param$phi=sample.phi(param)
  # param$phi=phi.true
  
  tmp=sample.break(param,jump1$break1)
  accept1$break1=accept1$break1+tmp$accept
  param$break1=tmp$break1
  # param$break1=break1.true[-c(1,length(break1.true))]
  
  param$z=sample.z(param)
  # param$z=z.true
  
  #adaptive piece
  if (i%%accept.output==0 & i<1000){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
  }
  
  #store results
  vec.theta[i,]=param$theta
  vec.phi[i,]=param$phi
  vec.break1[i,]=param$break1
  vec.logl[i]=get.marg.logl(param$theta,param$phi,param$break1)
}

#---------------------------------------------------
#look at theta
seq1=(ngibbs*0.8):ngibbs
theta.estim=matrix(colMeans(vec.theta[seq1,]),nloc,ncommun)
boxplot(theta.estim)

# seq.comm=1:5
seq.comm=c(4,1,2,3,5)
theta.estim1=theta.estim[,seq.comm]
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:5){
  lines(1:nloc,theta.estim1[,i],col=i)
}

#look at phi
phi.estim=matrix(colMeans(vec.phi[seq1,]),ncommun,nspp)[seq.comm,]
rango=range(phi.estim,phi.true)
plot(phi.true,phi.estim,ylim=rango,xlim=rango)
lines(rango,rango)

#look at z
rango=range(param$z,z.true)
plot(z.true,param$z,xlim=rango,ylim=rango)
lines(rango,rango)

#look at breaks
break2=break1.true[-c(1,length(break1.true))]
rango=range(param$break1,break2)
plot(break2,param$break1,xlim=rango,ylim=rango)
lines(rango,rango)
