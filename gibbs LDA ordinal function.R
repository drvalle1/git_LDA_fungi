LDA_ordinal=function(dat,ncomm.max,ngibbs,prop.burn,thetas.p.lprob,thetas.pot){
  nloc=nrow(dat)
  nspp=ncol(dat)
  ncommun=ncomm.max
  uni=sort(unique(as.numeric(dat)))
  nuni=length(uni)
  
  #set initial values
  phi=matrix(0,ncommun,nspp)
  break1=seq(from=0,to=0.1,length.out=nuni-1)
  ones.nloc=rep(1,nloc)
  ones.nspp=rep(1,nspp)
  
  #things for MH algorithm
  n.phi=ncommun*nspp
  
  #useful stuff
  lo.prob=0.000001
  hi.prob=1-lo.prob
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
  accept.output=50
  nkeep=ceiling(ngibbs*(1-prop.burn))
  vec.theta=matrix(NA,nkeep,nloc*ncommun)
  vec.phi=matrix(NA,nkeep,ncommun*nspp)
  vec.break1=matrix(NA,nkeep,nuni-1)
  vec.logl=matrix(NA,ngibbs,1)
  
  jump1=list(break1=rep(1,nuni-1),
             break1.sum=0.2,
             break1.mult=0.05)
  accept1=list(break1=rep(0,nuni-1),
               break1.sum=0,
               break1.mult=0)
  param=list(theta=NA,phi=phi,break1=break1,z=z)
  
  #core MCMC algorithm
  options(warn=2)
  oo=1
  for (i in 1:ngibbs){
    print(i)
    param$theta=sample.theta(nloc=nloc,ncommun=ncommun,nspp=nspp,thetas.pot=thetas.pot,
                             phi=param$phi,z=param$z,thetas.p.lprob=thetas.p.lprob)
    # param$theta=theta.true
    
    param$phi=sample.phi(param,ncommun,nspp)
    # param$phi=rbind(phi.true,matrix(0.001,ncommun-3,nspp))
    
    #sample breaks
    tmp=sample.break(param,jump1$break1,nuni,indicator)
    accept1$break1=accept1$break1+tmp$accept
    param$break1=tmp$break1

    tmp=sample.break.sum(param,jump1$break1.sum,nuni,indicator)
    accept1$break1.sum=accept1$break1.sum+tmp$accept
    param$break1=tmp$break1

    tmp=sample.break.mult(param,jump1$break1.mult,nuni,indicator)
    accept1$break1.mult=accept1$break1.mult+tmp$accept
    param$break1=tmp$break1
    # param$break1=break1.true[-c(1,length(break1.true))]
    
    param$z=sample.z(param,nuni,indicator,n.indicator)
    # param$z=z.true
    
    if (i%%accept.output==0 & i<(ngibbs/2)){
      #adaptive MH
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
      
      #re-order results
      med=apply(param$theta,2,median)
      ind=order(med,decreasing=T)
      param$theta=param$theta[,ind] 
      param$phi=param$phi[ind,]
    }
    
    #store results
    vec.logl[i]=get.marg.logl(param$theta,param$phi,param$break1,nuni,indicator)
    
    if (i>(ngibbs*prop.burn)){
      vec.theta[oo,]=param$theta
      vec.phi[oo,]=param$phi
      vec.break1[oo,]=param$break1
      oo=oo+1
    }
  }
  list(theta=vec.theta,phi=vec.phi,break1=vec.break1,logl=vec.logl)
}
