rm(list=ls(all=TRUE))
setwd('U:\\GIT_models\\git_LDA_fungi')

#generate theta for 5 groups
ncomm=5
lo=0.00001
seq1=seq(from=lo,to=1-lo,length.out=60)
combo=expand.grid(x1=seq1,x2=seq1,x3=seq1,x4=seq1)
soma=rowSums(combo)
combo$fim=1-soma
cond=combo$fim >= 0
combo1=combo[cond,]; dim(combo); dim(combo1)
unique(rowSums(combo1))

#output results
colnames(combo1)=paste0('theta',1:ncomm)
write.csv(combo1,'potential thetas.csv',row.names=F)

#------------------------------------------
#calculate corresponding probabilities based on TSP
gamma1=0.1
unico=unique(unlist(combo1)); n.unico=length(unico); n.unico
lprobs=dbeta(unico,1,gamma1,log=T)
plot(x=unico,y=lprobs,type='h')

#convert thetas to vs
theta=combo1
n.theta=nrow(theta)
v=matrix(NA,n.theta,ncomm-1)
v[,1]=theta[,1]
tmp=rep(1,n.theta)
for (i in 2:(ncomm-1)){
  tmp=tmp*(1-v[,i-1])
  v[,i]=theta[,i]/tmp
}
range(v)

#checking
theta.estim=matrix(NA,n.theta,ncomm)
tmp=rep(1,n.theta)
for(i in 1:(ncomm-1)){
  theta.estim[,i]=v[,i]*tmp  
  tmp=tmp*(1-v[,i])
}
theta.estim[,ncomm]=tmp
hist(data.matrix(theta)-theta.estim)

#calculate probabilities
lprobs.mat=dbeta(v,1,gamma1,log=T)
lprobs.tot=rowSums(lprobs.mat)
write.csv(lprobs.tot,'potential thetas lprob.csv',row.names=F)
