#!/usr/bin/env Rscirpt

#request large node 
#salloc -t 2-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript BCA_SKAT_MPI.R

rm(list=ls())

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/skat_gene_snp_bca_filteredMAF_19Feb2015.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid

Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)





nprobes <- nrow(genesnptable)

sink("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/MPI_output.txt") #changed-

### Prepare Rmpi environment, DO NOT CHANGE ###
Sys.getenv(c("SLURM_SUBMIT_DIR"))
Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST","SLURM_NNODES", "SLURM_NTASKS", 
             "SLUR M_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","SLURM_NTASKS_PER_NODE", 
             "SLURM_TASK_PID", "SLURM_ PARTITION"))
#
library('Rmpi')
njob <- mpi.universe.size() - 1L
mpi.spawn.Rslaves(nslaves=njob,needlog = F)
.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}
#
starttime<-Sys.time();starttime

compute_pvalue=function(i=1,genetable,pheno,genodat,covdat)
{
  result=data.frame(p=NA)
  selectedsnps <- unlist(strsplit(genetable[i,2],",",fixed=T)) 
  
  idx=match(selectedsnps,rownames(genodat))
  
  if (sum(is.na(idx))!=length(idx)) {
    library(SKAT)
    idx <- idx[!is.na(idx)]
    Z=t(genodat[idx,,drop=F])
    colnames(Z) <- row.names(genodat)[idx]
    obj.s<-SKAT_Null_Model(pheno ~ as.matrix(covdat),out_type="D")
    out=SKATBinary(Z, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
    result$p <- out$p.value
  }  
  
  return(result)
}

mpi.bcast.Robj2slave(sampletable)
#mpi.bcast.Robj2slave(genotype)
mpi.bcast.Robj2slave(Covariate)
mpi.bcast.Robj2slave(compute_pvalue)

mpi_compute_pvalue <- function(ngenes=40){
  
  nchunks=ceiling(nprobes/ngenes)
  print(paste0("number of total:",nchunks))
  rows=1:nprobes
  out=NULL
  for (ii in 1:nchunks)
  {
    if (ii%% 50==0)cat(ii,"..")
    if (ii<nchunks)
    {
      seqs=rows[((ii-1)*ngenes+1):(ii*ngenes)]
    }else
    {
      seqs=rows[((ii-1)*ngenes+1):length(rows)]
    }
    
    out1 <- matrix(NA,length(seqs),4)
    rownames(out1)=genesnptable$gene[seqs]
    cgenes <- genesnptable$gene[seqs]
    csnps <-  unique(unlist(strsplit(genesnptable[genesnptable$gene%in%cgenes,2],",",fixed=T))) 
    
    genodat=allgenotype[rownames(allgenotype)%in% csnps,]
    genetable <- genesnptable[seqs,]
    
    pheno=1*(sampletable$phenoBE_bca==2)
    
    idx=!is.na(pheno) & sampletable$recurrent_HB_RF==1
    pheno=pheno[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,genodat=genodat1,covdat=covdat,job.num=njob)
    out1[,1]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=1*(sampletable$phenoEA_bca==2)
    
    idx=!is.na(pheno)& sampletable$recurrent_HB_RF==1
    pheno=pheno[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,genodat=genodat1,covdat=covdat,job.num=njob)
    out1[,2]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=1*(sampletable$phenoEABE_bca==2)
    
    idx=!is.na(pheno)& sampletable$recurrent_HB_RF==1
    pheno=pheno[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,genodat=genodat1,covdat=covdat,job.num=njob)
    out1[,3]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=sampletable$phenoEA_bca
    pheno[is.na(pheno)] = 3
    pheno[pheno==1] = NA
    pheno=1*(pheno==2)
    
    idx=!is.na(pheno)& sampletable$recurrent_HB_RF==1
    pheno=pheno[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,genodat=genodat1,covdat=covdat,job.num=njob)
    out1[,4]=matrix(unlist(tmp),ncol=1,byrow = T)
    out=rbind(out,out1)
  }
  return(out)
}



Sys.time()
SKATout=mpi_compute_pvalue()

sink()

save(SKATout,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_SKAT_output.RData") 

