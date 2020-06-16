#to test function "convertThetatoV"

set.seed(1)
theta.new=rdirichlet(nloc,alpha=rep(1,ncommun))

#convert from theta to v to get priors
vlk.new=convertThetatoV(theta.new)
vlk.new[,ncommun]=1

#compare to R code
v=matrix(NA,nloc,ncommun)
v[,1]=theta.new[,1]
v[,ncommun]=1
pred=1-v[,1]
for (i in 2:ncommun){
  v[,i]=theta.new[,i]/pred
  pred=pred*(1-v[,i])
}

plot(v,vlk.new,xlim=rango,ylim=rango)

#compare the generated theta
teste=convertVtoTheta(vlk.new,ones.nloc)

rango=c(0,1)
plot(theta.new,teste,xlim=rango,ylim=rango)
lines(rango,rango)

