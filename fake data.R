rm(list=ls(all=TRUE))
set.seed(291)

nloc=100
nspp=200
ncommun=5
base=floor(nloc/(ncommun-1))

#generate thetas
x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)
  
init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)
  
theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)

#eliminate zeroes
cond=theta<0.000001
theta[cond]=0.000001
theta=theta/rowSums(theta)

theta.true=theta
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate phi  
phi=matrix(NA,ncommun,nspp)
mu.large=1.5
mu.small=-1
ind=matrix(rbinom(ncommun*nspp,size=1,prob=0.2),ncommun,nspp)
for (i in 1:ncommun){
  ind1=ind[i,]
  cond=ind1==1
  tmp=rep(NA,nspp)
  tmp[cond ]=rnorm(sum( cond),mean=mu.large,sd=0.6)
  tmp[!cond]=rnorm(sum(!cond),mean=mu.small,sd=0.4)
  phi[i,]=tmp
}
image(phi)
phi.true=phi
phi[1:5,1:5]
hist(phi)
ind[1:5,1:5]

#how many unique species does each community have?
teste=rep(0,ncommun)
for (i in 1:nspp){
  tmp=ind[,i]
  if (sum(tmp)==1){
    ind1=which(tmp==1)
    teste[ind1]=teste[ind1]+1
  }
}
teste

#calculate probabilities
medias=theta%*%phi; dim(medias)
z.true=z=matrix(rnorm(nloc*nspp,mean=medias,sd=1),nloc,nspp)

#generate actual observations y
range(z)
max1=max(phi.true); max1
break1.true=break1=c(-Inf,seq(from=0,to=3.1,by=0.1),Inf)
y=matrix(NA,nloc,nspp)
for (i in 2:length(break1)){
  cond=z>break1[i-1] & z<break1[i]  
  y[cond]=i-2
}
tmp=table(y); tmp; 
length(tmp); length(break1)

setwd('U:\\GIT_models\\git_LDA_fungi') 
write.csv(y,'fake data.csv',row.names=F)    
