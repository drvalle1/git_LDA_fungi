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
get.logl=function(theta,phi){
  media=theta%*%phi
  dnorm(dat,mean=media,sd=1,log=T)
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
    pold=get.logl(theta.old,param$phi)
    pnew=get.logl(theta.new,param$phi)
    pold1=(pold%*%ones.nspp)+prior.old[,i]
    pnew1=(pnew%*%ones.nspp)+prior.new[,i]
    
    k=acceptMH(pold1,pnew1+fix1[,i],vlk.old[,i],vlk.new[,i],F)
    vlk.old[,i]=k$x
  }
  
  list(vlk=vlk.old,accept=vlk.old[,-ncommun]!=vlk.orig[,-ncommun])
}
#-----------------------------------------------------------------------------------------------
sample.phi=function(param,jump){
  phi.orig=phi.old=param$phi
  
  #things for MH
  tmp=rnorm(n.phi,mean=phi.old,sd=jump)
  phi.prop=matrix(tmp,ncommun,nspp)
  
  #priors
  prior.old=dnorm(phi.old,mean=0,sd=1,log=T)
  prior.new=dnorm(phi.prop,mean=0,sd=1,log=T)
  
  for (i in 1:ncommun){
    phi.new=phi.old
    phi.new[i,]=phi.prop[i,]
    
    pold=get.logl(param$theta,phi.old)
    pnew=get.logl(param$theta,phi.new)
    pold1=(ones.nloc%*%pold)+prior.old[i,]
    pnew1=(ones.nloc%*%pnew)+prior.new[i,]
    
    k=acceptMH(pold1,pnew1,phi.old[i,],phi.new[i,],F)
    phi.old[i,]=k$x
  }
  
  list(phi=phi.old,accept=phi.old!=phi.orig)
}
#----------------------------
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
