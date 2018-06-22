rmvnorm1=function(n,sigma){
  ev <- eigen(sigma, symmetric = TRUE)
  R = t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values,0))))
  matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = T) %*% R
}
#------------------------------------
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
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
get.logl=function(theta,phi,z){
  media=theta%*%phi
  dnorm(z,media,sd=1,log=T)
}
#-----------------------------------------------------------------------------------------------
sample.vlk=function(param,jump){
  vlk.orig=vlk.old=param$vlk
  
  #things for MH
  tmp=tnorm(n.vlk,lo=lo.vlk,hi=hi.vlk,mu=vlk.old[,-ncommun],sig=jump)
  vlk.prop=matrix(tmp,nloc,ncommun-1)
  tmp=fix.MH(lo=lo.vlk,hi=hi.vlk,old1=vlk.old[,-ncommun],new1=vlk.prop,jump)
  fix1=matrix(tmp,nloc,ncommun-1)
  
  #priors
  prior.old=dbeta(vlk.old[,-ncommun],1,gamma,log=T)
  prior.new=dbeta(vlk.prop,1,gamma,log=T)
  
  for (i in 1:(ncommun-1)){
    vlk.new=vlk.old
    vlk.new[,i]=vlk.prop[,i]
    
    theta.old=convertVtoTheta(vlk.old,ones.nloc)
    theta.new=convertVtoTheta(vlk.new,ones.nloc)
    pold=get.logl(theta.old,param$phi,param$z)
    pnew=get.logl(theta.new,param$phi,param$z)
    pold1=(pold%*%ones.nspp)+prior.old[,i]
    pnew1=(pnew%*%ones.nspp)+prior.new[,i]
    
    k=acceptMH(pold1,pnew1+fix1[,i],vlk.old[,i],vlk.new[,i],F)
    vlk.old[,i]=k$x
  }
  
  list(vlk=vlk.old,accept=vlk.old[,-ncommun]!=vlk.orig[,-ncommun])
}
#-----------------------------------------------------------------------------------------------
sample.phi=function(param){
  theta.t.theta=t(param$theta)%*%param$theta
  prec=theta.t.theta+diag(1,ncommun)
  var1=solve(prec)
  w=t(rmvnorm1(nspp,sigma=var1))
  media=var1%*%t(param$theta)%*%param$z
  media+w
}
#-----------------------------------------------------------------------------------------------
sample.break=function(param){
  util=matrix(NA,nuni,2)  
  for (i in 1:nuni){
    util[i,]=range(param$z[indicator[[i]]])
  }
  colnames(util)=c('min1','max1')
  lo1=util[-nuni,'max1']
  hi1=util[-1,'min1']
  runif(nuni-1,min=lo1,max=hi1)
}
#-----------------------------------------------------------------------------------------------
sample.z=function(param){
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
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
