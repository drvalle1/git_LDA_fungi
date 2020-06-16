res=matrix(NA,100,2)
for (i in 1:100){
  theta.old=rdirichlet(1,alpha=rep(1,ncommun))
  theta.new=rdirichlet(1,alpha=rep(1,ncommun))
  
  res[i,1]=log.ddirichlet(alpha=theta.old,theta=theta.new)
  res[i,2]=log(ddirichlet(x=theta.new,alpha=theta.old))
}
hist(res[,1]-res[,2])