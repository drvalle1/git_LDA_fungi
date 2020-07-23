LDA_ordinal=function(dat,ncomm,ngibbs,prop.burn){
  nloc=nrow(dat)
  nspp=ncol(dat)
  ncommun=ncomm
  uni=sort(unique(as.numeric(dat)))
  nuni=length(uni)
  
  #set initial values
  phi=matrix(0,ncommun,nspp)
  break1=seq(from=0,to=10,length.out=nuni-1)
  ones.nloc=rep(1,nloc)
  ones.nspp=rep(1,nspp)
  theta=matrix(1/ncommun,nloc,ncommun)
    
  #things for MH algorithm
  n.phi=ncommun*nspp
  
  #useful stuff
  lo.prob=0.000001
  hi.prob=1-lo.prob
  indicator=list()
  n.indicator=rep(NA,nuni)
  break2=c(-1,break1,11)
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
  
  jump1=list(break1=rep(1,nuni-1),
             theta=rep(0.1,nloc))
  accept1=list(break1=rep(0,nuni-1),
               theta=rep(0,nloc))
  param=list(theta=theta,phi=phi,break1=break1,z=z)
  
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
    
    tmp=sample.theta(nloc=nloc,ncommun=ncommun,nspp=nspp,theta=param$theta,
                     phi=param$phi,jump1=jump1$theta,break1=param$break1,nuni=nuni,
                     indicator=indicator)
    param$theta=tmp$theta
    accept1$theta=accept1$theta+tmp$accept
    # param$theta=theta.true
    
    param$phi=sample.phi(param,ncommun,nspp)
    # param$phi=rbind(phi.true,matrix(0.001,ncommun-3,nspp))
    
    param$z=sample.z(param,nuni,indicator,n.indicator)
    # param$z=z.true
    
    if (i%%accept.output==0 & i<(ngibbs/2)){
      #adaptive MH
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
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
