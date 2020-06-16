library('microbenchmark')

func=function(){
  pp=matrix(z[i,],n.thetas.pot,nspp,byrow=T)
  p1=-(1/2)*rowSums((ztmp-media)^2)  
}

microbenchmark(
  func(),
  -(1/2)*CalcSqDiff(z=z[i,], media=media),
  times=10
)
 
microbenchmark(
  SampleIndTheta(z=z, media=media, lprior=thetas.p.lprob, runif1=runif(nloc)),
  sample.theta(nloc=nloc,ncommun=ncommun,nspp=nspp,thetas.pot=thetas.pot,
               phi=phi,z=z,thetas.p.lprob=thetas.p.lprob),
  times=10
)

