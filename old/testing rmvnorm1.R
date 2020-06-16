rm(list=ls(all=TRUE))
set.seed(1)
setwd('U:\\GIT_models\\git_LDA_fungi') 
source('gibbs LDA fungi functions.R')

seq1=0:10
tmp=outer(seq1,seq1,'-')
rho=0.9
sig2=1
sigma=sig2*(rho^abs(tmp))

n=100000
z1=rmvnorm1(n,sigma)
sigma.estim=var(z1)

rango=range(c(sigma,sigma.estim))
plot(sigma,sigma.estim,xlim=rango,ylim=rango)
lines(rango,rango)