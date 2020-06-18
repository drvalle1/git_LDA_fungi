rm(list=ls(all=TRUE))
setwd('U:\\GIT_models\\git_LDA_fungi')

seq1=seq(from=0,to=1,length.out=21); seq1

#generate theta for 5 groups
ncomm=5

combo=expand.grid(x1=seq1,x2=seq1,x3=seq1,x4=seq1)
soma=rowSums(combo)
combo$fim=1-soma
cond=combo$fim >= 0
combo1=combo[cond,]; dim(combo); dim(combo1)

#cannot never be 0 or 1
lo=0.000001
combo1[combo1<lo]=lo
combo1[combo1>1-lo]=1-lo
combo1=combo1/rowSums(combo1) #re-normalize
unique(rowSums(combo1))
apply(combo1,2,range)
boxplot(combo1)

#output results
colnames(combo1)=paste0('theta',1:ncomm)
write.csv(combo1,'potential thetas.csv',row.names=F)

#------------------------------------------
#calculate corresponding probabilities based on dirichlet 0.1
gamma1=0.1
lprobs.tot=rowSums((0.1-1)*log(combo1))
range(lprobs.tot)

write.csv(lprobs.tot,'potential thetas lprob.csv',row.names=F)
