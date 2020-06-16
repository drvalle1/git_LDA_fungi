rm(list=ls(all=TRUE))
setwd('U:\\GIT_models\\git_LDA_fungi')

seq1=seq(from=0,to=1,length.out=21); seq1

#generate theta for 5 groups
ncomm=3
lo=0.00001
seq1=seq(from=lo,to=1-lo,length.out=61)
combo=expand.grid(x1=seq1,x2=seq1)
soma=rowSums(combo)
combo$fim=1-soma
cond=combo$fim >= 0
combo1=combo[cond,]; dim(combo); dim(combo1)
unique(rowSums(combo1))

#output results
colnames(combo1)=paste0('theta',1:ncomm)
write.csv(combo1,'potential thetas.csv',row.names=F)
