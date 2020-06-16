theta=matrix(NA,nloc,ncomm)
for (i in 1:nloc){
  media=theta.pot%*%phi
  ztmp=matrix(z[i,],nloc,nspp,byrow=T)
  prob=-(1/2)*rowsSums((ztmp-media)^2)
  prob1=exp(prob)
  prob2=prob1/rowSums(prob1)
  ind=rmultinom(1,size=1,prob=prob2)
  theta[i,]=theta.pot[ind,]
}
theta