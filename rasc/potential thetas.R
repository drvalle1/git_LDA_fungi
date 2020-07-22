rm(list=ls(all=TRUE))

seq1=seq(from=0,to=1,length.out=41); seq1

#generate thetas for 3 groups
for (ncomm in 2:5){
  if (ncomm==2) combo=data.frame(x1=seq1)
  if (ncomm==3) combo=expand.grid(x1=seq1,x2=seq1)
  if (ncomm==4) combo=expand.grid(x1=seq1,x2=seq1,x3=seq1)
  if (ncomm==5) combo=expand.grid(x1=seq1,x2=seq1,x3=seq1,x4=seq1)
  
  soma=rowSums(combo)
  combo$fim=1-soma
  cond=combo$fim >= 0
  combo1=combo[cond,]; dim(combo); 
  print(c(ncomm,dim(combo1)))
  
  #output results
  colnames(combo1)=paste0('theta',1:ncomm)
  nome=paste0('potential thetas',ncomm,'.csv')
  
  setwd('U:\\GIT_models\\git_LDA_fungi')
  write.csv(combo1,nome,row.names=F)
}
