LDA_ordinal=function(dat,ncomm,ngibbs,prop.burn,gamma1){
  nloc=nrow(dat)
  nspp=ncol(dat)
  ncommun=ncomm
  uni=sort(unique(as.numeric(dat)))
  nuni=length(uni)
  
  #set initial values
  phi=matrix(0,ncommun,nspp)
  break1=seq(from=0,to=1,length.out=nuni-1)
  ones.nloc=rep(1,nloc)
  ones.nspp=rep(1,nspp)
  theta=matrix(1/ncommun,nloc,ncommun)
  vmat=theta.to.v(theta=theta,nloc=nloc,ncommun=ncommun)
  # theta.tmp=v.to.theta(vmat=vmat,nloc=nloc,ncommun=ncommun)
  
  #things for MH algorithm
  n.phi=ncommun*nspp
  
  #useful stuff
  lo.prob=0.000001
  hi.prob=1-lo.prob
  indicator=list()
  n.indicator=rep(NA,nuni)
  break2=c(-0.1,break1,1.1)
  z=matrix(NA,nloc,nspp)
  for (i in 1:nuni){
    cond=dat==uni[i]
    ind=which(cond)
    indicator[[i]]=ind
    z[ind]=mean(break2[i:(i+1)])
    n.indicator[i]=length(ind)
  }

  #gibbs stuff
  accept.output=50
  nkeep=ceiling(ngibbs*(1-prop.burn))
  vec.theta=matrix(NA,nkeep,nloc*ncommun)
  vec.phi=matrix(NA,nkeep,ncommun*nspp)
  vec.break1=matrix(NA,nkeep,nuni-1)
  vec.logl=matrix(NA,ngibbs,1)
  
  jump1=list(break1=rep(0.01,nuni-1),vmat=matrix(0.2,nloc,ncommun-1),break1.mult=0.01)
  accept1=list(break1=rep(0,nuni-1),vmat=matrix(0,nloc,ncommun-1),break1.mult=0)
  param=list(theta=theta,vmat=vmat,phi=phi,break1=break1,z=z)
  
  #core MCMC algorithm
  options(warn=2)
  oo=1
  for (i in 1:ngibbs){
    print(i)
    
    #sample breaks
    tmp=sample.break(param=param,jump=jump1$break1,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp)
    accept1$break1=accept1$break1+tmp$accept
    param$break1=tmp$break1
    # param$break1=break1.true[-c(1,length(break1.true))]    
    
    tmp=sample.break.mult(param=param,jump=jump1$break1.mult,nuni=nuni,indicator=indicator,nloc=nloc,nspp=nspp)
    accept1$break1.mult=accept1$break1.mult+tmp$accept
    param$break1=tmp$break1
    
    tmp=sample.theta(nloc=nloc,ncommun=ncommun,nspp=nspp,vmat=param$vmat,
                     phi=param$phi,jump1=jump1$vmat,break1=param$break1,nuni=nuni,
                     indicator=indicator,gamma1=gamma1)
    param$theta=tmp$theta
    param$vmat=param$vmat
    accept1$vmat=accept1$vmat+tmp$accept
    # param$theta=theta.true
    
    param$phi=sample.phi(param,ncommun,nspp)
    # param$phi=rbind(phi.true,matrix(0.001,ncommun-3,nspp))
    
    param$z=sample.z(param,nuni,indicator,n.indicator)
    # param$z=z.true
    
    if (i%%accept.output==0 & i<ngibbs*prop.burn){
      #adaptive MH
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #re-shuffle groups
    if (i%%200==0 & i<ngibbs*prop.burn){
      med=apply(param$theta,2,median)
      ind=order(med,decreasing=T)
      param$phi=param$phi[ind,]
      param$theta=param$theta[,ind]
      param$vmat=theta.to.v(theta=param$theta,nloc=nloc,ncommun=ncommun)
      boxplot(param$theta,main=i)
      jump1$vmat[]=0.05
    }
    
    #store results
    media=param$theta%*%param$phi
    vec.logl[i]=sum(get.marg.logl(media=media,break1=param$break1,nuni=nuni,indicator=indicator,
                                  nloc=nloc,nspp=nspp))
    
    if (i>(ngibbs*prop.burn)){
      vec.theta[oo,]=param$theta
      vec.phi[oo,]=param$phi
      vec.break1[oo,]=param$break1
      oo=oo+1
    }
  }
  list(theta=vec.theta,phi=vec.phi,break1=vec.break1,logl=vec.logl)
}
