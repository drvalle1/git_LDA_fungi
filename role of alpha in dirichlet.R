#smaller values of alpha = greater variability
z=rdirichlet(1000,rep(100,ncommun))
apply(z,2,var)

z=rdirichlet(1000,rep(1,ncommun))
apply(z,2,var)

z=rdirichlet(1000,rep(.01,ncommun))
apply(z,2,var)