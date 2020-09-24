rmvnorm1=function(n,sigma){
  ev <- eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = T) %*% R
}
#------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-----------------------------------------------------------------------------------------------
acceptMH.indicator <- function(p0,p1){  
  nz           <- length(p0)  #no. to accept
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  which(z < a)
  
}
#-----------------------------------------------------------------------------------------------
rdirichlet1=function(alpha){
  nrow1=nrow(alpha)
  ncol1=ncol(alpha)
  x=matrix(rgamma(nrow1*ncol1,alpha),nrow1,ncol1)
  soma=rowSums(x)
  x/soma
}
#-----------------------------------------------------------------------------------------------
ldirichlet=function(alpha,x){
  #calculate dirichlet
  p1=lgamma(rowSums(alpha))
  p2=-rowSums(lgamma(alpha))
  p3=rowSums((alpha-1)*log(x))
  p1+p2+p3
}
#-----------------------------------------------------------------------------------------------
v.to.theta=function(vmat,ncommun,nloc){
  theta=matrix(NA,nloc,ncommun)
  tmp=rep(1,nloc)
  for (j in 1:(ncommun-1)){
    theta[,j]=vmat[,j]*tmp
    tmp=tmp*(1-vmat[,j])
  }
  theta[,ncommun]=tmp
  # rowSums(theta)
  theta
}

#------------------------------------------------------------------------------
theta.to.v=function(theta,nloc,ncommun){
  vmat=matrix(NA,nloc,ncommun-1)
  tmp=rep(1,nloc)  
  for (i in 1:(ncommun-1)){
    vmat[,i]=theta[,i]/tmp
    tmp=tmp*(1-vmat[,i])
  }
  vmat
}
#------------------------------------------------------------------------------
sample.theta=function(nloc,vmat,ncommun,phi,nspp,jump1,break1,nuni,indicator,gamma1){
  vorig=vold=vmat
  tmp=tnorm(nloc*(ncommun-1),lo=0,hi=1,mu=vold,sig=jump1)
  v.prop=matrix(tmp,nloc,ncommun-1)
  tmp=(pnorm(1,mean=vold,sd=jump1)  -pnorm(0,mean=vold,sd=jump1))/
      (pnorm(1,mean=v.prop,sd=jump1)-pnorm(0,mean=v.prop,sd=jump1))
  lprob.propose=matrix(log(tmp),nloc,ncommun-1)
  
  #get priors
  lprob.prior.old=dbeta(vold,1,gamma1,log=T)
  lprob.prior.new=dbeta(v.prop,1,gamma1,log=T)

  for (j in 1:(ncommun-1)){
    vnew=vold
    vnew[,j]=v.prop[,j]
    theta.old=v.to.theta(vmat=vold,nloc=nloc,ncommun=ncommun)
    theta.new=v.to.theta(vmat=vnew,nloc=nloc,ncommun=ncommun)
    llk.old=rowSums(get.marg.logl(media=theta.old%*%phi,break1=break1,nuni=nuni,indicator=indicator,
                                  nloc=nloc,nspp=nspp))
    llk.new=rowSums(get.marg.logl(media=theta.new%*%phi,break1=break1,nuni=nuni,indicator=indicator,
                                  nloc=nloc,nspp=nspp))
    pold=llk.old+lprob.prior.old[,j]
    pnew=llk.new+lprob.prior.new[,j]+lprob.propose[,j]
    
    k=acceptMH(p0=pold,p1=pnew,x0=vold[,j],x1=vnew[,j],BLOCK=F)
    vold[,j]=k$x
  }
  theta=v.to.theta(vmat=vold,nloc=nloc,ncommun=ncommun)
  list(vmat=vold,accept=vold!=vorig,theta=theta)
}
#-----------------------------------------------------------------------------------------------
sample.phi=function(param,ncommun,nspp){
  theta.t.theta=t(param$theta)%*%param$theta
  prec=theta.t.theta+diag(1,ncommun)
  var1=solve(prec)
  w=t(rmvnorm1(nspp,sigma=var1))
  media=var1%*%t(param$theta)%*%param$z
  media+w
}
#-----------------------------------------------------------------------------------------------
sample.break=function(param,jump,nuni,indicator,nloc,nspp){
  media=param$theta%*%param$phi
  diffs.old=diff(param$break1)
  
  for (i in 1:length(diffs.old)){ #the first break is set to 0 for identifiability purposes
    diffs.new=diffs.old
    diffs.new[i]=abs(rnorm(1,mean=diffs.old[i],sd=jump[i])) #reflect off of zero
    
    #to avoid numerical issues
    if (diffs.new[i]<0.0000001) diffs.new[i]=0.0000001
    
    #get implied breakpoints
    break1.old=c(0,cumsum(diffs.old))
    break1.new=c(0,cumsum(diffs.new))
    
    #get loglikel
    llikel.old=sum(get.marg.logl(media=media,break1=break1.old,nuni=nuni,indicator=indicator,
                                 nloc=nloc,nspp=nspp))
    llikel.new=sum(get.marg.logl(media=media,break1=break1.new,nuni=nuni,indicator=indicator,
                                 nloc=nloc,nspp=nspp))
    
    #get priors
    prior.old=0 #uniform prior
    prior.new=0 #uniform prior
    
    #accept or reject
    k=acceptMH(p0=llikel.old+prior.old,
               p1=llikel.new+prior.new,
               x0=diffs.old[i],x1=diffs.new[i],F)
    diffs.old[i]=k$x
  }
  break1=c(0,cumsum(diffs.old))
  list(break1=break1,accept=break1!=param$break1)
}
#-----------------------------------------------------------------------------------------------
sample.z=function(param,nuni,indicator,n.indicator){
  media=param$theta%*%param$phi
  break1=c(-100,param$break1,100) #to avoid numerical issues
  z=param$z
  for (i in 1:nuni){
    ind=indicator[[i]]
    z[ind]=tnorm(n.indicator[i],lo=break1[i],hi=break1[i+1],mu=media[ind],sig=1)
  }
  z
}
#-----------------------------------------------------------------------------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.0001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#-----------------------------------------------------------------------------------------------
get.marg.logl=function(media,break1,nuni,indicator,nloc,nspp){
  break2=c(-Inf,break1,Inf)
  prob=matrix(NA,nloc,nspp)
  for (i in 1:nuni){
    #get probabilities
    ind1=indicator[[i]]
    tmp=pnorm(break2[i+1]-media[ind1])-pnorm(break2[i]-media[ind1])
    prob[ind1]=log(tmp)
  }
  prob
}
#----------------------------------------------------------------------------
# sample.break.sum=function(param,jump,nuni,indicator,nloc,nspp){
#   media=param$theta%*%param$phi
#   param$break1[1]=0
#   break1.old=param$break1
#   
#   #lower bound avoids choosing a number that would result in param$break1[2]<0
#   lo1=-param$break1[2]
#   tmp=tnorm(1,lo=lo1,hi=Inf,mu=0,sig=jump) 
#   break1.new=break1.old+tmp
#   break1.new[1]=0
#   
#   #get probabilities
#   pold=sum(get.marg.logl(media=media,break1=break1.old,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp))+
#     pnorm(break1.new[2],mean=0,sd=jump,log=T)
#   pnew=sum(get.marg.logl(media=media,break1=break1.new,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp))+
#     pnorm(break1.old[2],mean=0,sd=jump,log=T)
#   
#   #accept or reject
#   k=acceptMH(pold,pnew,0,1,F)
#   if (k$x==1) break1.old=break1.new
#   list(break1=break1.old,accept=k$x)
# }
#-----------------------------------------------------------------------------------------------
sample.break.mult=function(param,jump,nuni,indicator,nloc,nspp){
  media=param$theta%*%param$phi
  param$break1[1]=0
  break1.old=param$break1
  
  #lower bound avoids choosing a number that would result in param$break1[2]<0
  tmp=tnorm(1,lo=0,hi=Inf,mu=1,sig=jump) 
  break1.new=break1.old*tmp
  
  #adjustment for truncated proposal
  fix1=dnorm(1/tmp,mean=1,sd=jump,log=T)-dnorm(tmp,mean=1,sd=jump,log=T)
  
  #get probabilities
  pold=sum(get.marg.logl(media=media,break1=break1.old,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp))
  pnew=sum(get.marg.logl(media=media,break1=break1.new,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp))

  #accept or reject
  k=acceptMH(pold,pnew+fix1,0,1,F)
  if (k$x==1) break1.old=break1.new
  list(break1=break1.old,accept=k$x)
}
