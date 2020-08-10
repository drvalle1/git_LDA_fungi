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
sample.theta=function(nloc,theta,ncommun,phi,nspp,jump1,break1,nuni,indicator){
  jump2=matrix(1/jump1,nloc,ncommun) #bigger values of jump2==smaller variance
  theta.old=theta
  
  #create theta.new
  alpha.old.to.new=theta.old*jump2
  theta.new=rdirichlet1(alpha=alpha.old.to.new)
  alpha.new.to.old=theta.new*jump2
  
  #for numerical stability
  cond=theta.new<0.0000001
  theta.new[cond]=0.0000001
  theta.new=theta.new/rowSums(theta.new)
  
  #get jump probabilities
  jump.old.to.new=ldirichlet(alpha=alpha.old.to.new,x=theta.new)
  jump.new.to.old=ldirichlet(alpha=alpha.new.to.old,x=theta.old)
  
  #MH algorithm
  media.old=theta.old%*%phi
  media.new=theta.new%*%phi
  
  #get loglikel
  llk.old=rowSums(get.marg.logl(media=media.old,break1=break1,nuni=nuni,indicator=indicator,
                                nloc=nloc,nspp=nspp))
  llk.new=rowSums(get.marg.logl(media=media.new,break1=break1,nuni=nuni,indicator=indicator,
                                nloc=nloc,nspp=nspp))

  #get prior probabilities
  # alpha.mat=matrix(theta.prior,nloc,ncommun)
  prior.old=0#ldirichlet(alpha=alpha.mat,x=theta.old)
  prior.new=0#ldirichlet(alpha=alpha.mat,x=theta.new)
  
  ind=acceptMH.indicator(p0=llk.old+jump.old.to.new+prior.old,
                         p1=llk.new+jump.new.to.old+prior.new)
  accept=rep(0,nloc)
  if (length(ind)!=0) {
    theta.old[ind,]=theta.new[ind,]
    accept[ind]=1
  }
  list(theta=theta.old,accept=accept)

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
