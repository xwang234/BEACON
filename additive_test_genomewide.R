#!/usr/bin/env Rscirpt
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Dong23SNPs.RData")
# salloc -t 6-1 -n 51 mpirun -n 1 R --interactive
#salloc -t 6-1 --mem-per-cpu 32G -n 51  mpirun -n 1 R --interactive
library(CGEN)
library(readxl)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")
#genotype data needs to be 0,1,2
for (i in 1:ncol(allgenotype))
{
  idx=which(!allgenotype[,i] %in% c(0,1,2))
  if (length(idx)>0)
  {
    if (length(idx)>nrow(allgenotype)/100) 
    { 
      print(i)
      print(length(idx)/nrow(allgenotype))
    }
    idx1=which(allgenotype[idx,i]>0 & allgenotype[idx,i]<0.5)
    allgenotype[idx[idx1],i]=0
    idx1=which(allgenotype[idx,i]>=0.5 & allgenotype[idx,i]<1.5)
    allgenotype[idx[idx1],i]=1
    idx1=which(allgenotype[idx,i]>=1.5)
    allgenotype[idx[idx1],i]=2
  }
}
nprobes <- nrow(allgenotype)

sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

#fam is for the genotype data
allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid
all(allsamples==fam$V2) #T

sampletable$BEEA=NA

sampletable$BEEA[which(sampletable$phenoBE_bc==1 | sampletable$phenoEA_bc==1)]=0

sampletable$BEEA[which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2)]=1

table(sampletable$BEEA,useNA="ifany")

sampletable$BE=I(sampletable$phenoBE_bc==2)

sampletable$EA=I(sampletable$phenoEA_bc==2)

sampletable$bmi_cat <- ifelse(sampletable$bmi_recent_healthy>=30,1,0)

table(sampletable$bmi_cat,useNA="ifany")

sampletable$sex=as.factor(sampletable$sex)

exposures=c("bmi_cat","recurrent_HB_RF","cig_smk_ever")

sampletable1=sampletable[,which(colnames(sampletable) %in% c("age","sex","BE","EA","BEEA",exposures))]
sampletable1$age[is.na(sampletable1$age)]=mean(sampletable1$age,na.rm=T)
compute_pvalue=function(i=1,response="BE",genodat)
{
  res=data.frame(matrix(NA,1,3),stringsAsFactors = F)
  colnames(res)=exposures
  dat1=cbind(sampletable1,snp=as.numeric(genodat[i,]))
  idx=complete.cases(dat1[,c("snp","age","sex",response,exposures[1])])
  dat2=dat1[idx,]
  if (length(unique(dat2$snp))>1) #snp needs to have at least 2 levels
  {
    tmp1=additive.test(data=dat1,response.var = response,snp.var ="snp",main.vars = c("age","sex"),exposure.var = exposures[1],op=list(genetic.model=1))
    res[,1]=tmp1$pval.add
  }
  idx=complete.cases(dat1[,c("snp","age","sex",response,exposures[2])])
  dat2=dat1[idx,]
  if (length(unique(dat2$snp))>1) 
  {
    tmp1=additive.test(data=dat1,response.var = response,snp.var ="snp",main.vars = c("age","sex"),exposure.var = exposures[2],op=list(genetic.model=1))
    res[,2]=tmp1$pval.add
  }
  idx=complete.cases(dat1[,c("snp","age","sex",response,exposures[3])])
  dat2=dat1[idx,]
  if (length(unique(dat2$snp))>1)
  {
    tmp1=additive.test(data=dat1,response.var = response,snp.var ="snp",main.vars = c("age","sex"),exposure.var = exposures[3],op=list(genetic.model=1))
    res[,3]=tmp1$pval.add
  }
  return(res)
}


mpi_compute_pvalue <- function(nsnps=2*njob,response="BE"){
  
  nchunks=ceiling(nprobes/nsnps)
  print(paste0("number of total:",nchunks))
  rows=1:nprobes
  out=NULL
  for (ii in 1:nchunks)
  {
    if (ii%% 500==0)
    {
      mytime=as.character(Sys.time())
      cat(round(ii/nchunks,2),mytime,"..")
    }  
      
    if (ii<nchunks)
    {
      seqs=rows[((ii-1)*nsnps+1):(ii*nsnps)]
    }else
    {
      seqs=rows[((ii-1)*nsnps+1):length(rows)]
    }
    out1 <- matrix(NA,length(seqs),3)
    rownames(out1)=rownames(allgenotype)[seqs]
    genodat=allgenotype[seqs,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genodat=genodat,response=response,job.num=njob)
    out1[,1:3]=matrix(unlist(tmp),ncol=3,byrow = T)
    out=rbind(out,out1)
  }
  colnames(out)=exposures
  return(out)
}

library('Rmpi')
njob <- mpi.universe.size() - 1L
mpi.spawn.Rslaves(nslaves=njob,needlog=F)
mpi.bcast.Robj2slave(sampletable1)
mpi.bcast.Robj2slave(exposures)
mpi.bcast.Robj2slave(compute_pvalue)
mpi.remote.exec(library(CGEN))


print(Sys.time())
BEres=mpi_compute_pvalue()
save(BEres,file="../result/res_additive_test_genomewide.RData")
print(Sys.time())
EAres=mpi_compute_pvalue(response="EA")
save(BEres,EAres,file="../result/res_additive_test_genomewide.RData")
print(Sys.time())
BEEAres=mpi_compute_pvalue(response="BEEA")
save(BEres,EAres,BEEAres,file="../result/res_additive_test_genomewide.RData")
print(Sys.time())
print("done")