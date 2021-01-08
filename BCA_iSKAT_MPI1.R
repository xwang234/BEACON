#!/usr/bin/env Rscirpt
# library(devtools)
# install_github("lin-lab/iSKAT-GESAT")
#request large node 
#salloc -t 2-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript BCA_iSKAT_MPI1.R

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
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)
Covariate$ev1_bca[is.na(Covariate$ev1_bca)]=mean(Covariate$ev1_bca,na.rm = T)
Covariate$ev2_bca[is.na(Covariate$ev2_bca)]=mean(Covariate$ev2_bca,na.rm = T)
Covariate$ev3_bca[is.na(Covariate$ev3_bca)]=mean(Covariate$ev3_bca,na.rm = T)
Covariate$ev4_bca[is.na(Covariate$ev4_bca)]=mean(Covariate$ev4_bca,na.rm = T)





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

compute_pvalue=function(i=1,genetable,pheno,E1,genodat1,covdat)
{
  result=data.frame(p=NA)
  
  selectedsnps <- unlist(strsplit(genetable[i,2],",",fixed=T)) 
  
  idx=match(selectedsnps,rownames(genodat1))
  
  if (sum(is.na(idx))!=length(idx)) {
    library(iSKAT)
    idx <- idx[!is.na(idx)]
    Z=t(genodat1[idx,,drop=F])
    colnames(Z) <- row.names(genodat1)[idx]
    pheno <- as.matrix(pheno)
    E1 <- as.matrix(E1)
    #catch erros
    out=tryCatch(
      {
        GESAT(Z, pheno, E1, X=as.matrix(covdat),out_type="D")
      },
      error=function(e)
      {
        return(F)
      }
    )
    if (class(out)=="logical")
    {
      result$p=NA
    }else
    {
      result$p <- out$pvalue
    }
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
    
    E=sampletable$recurrent_HB_RF
    
    idx=!is.na(pheno) & !is.na(E)
    pheno=pheno[idx]
    E1 <- E[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,E1=E1,genodat1=genodat1,covdat=covdat,job.num=njob)
    out1[,1]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=1*(sampletable$phenoEA_bca==2)
    idx=!is.na(pheno)& !is.na(E)
    pheno=pheno[idx]
    E1 <- E[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,E1=E1,genodat1=genodat1,covdat=covdat,job.num=njob)
    out1[,2]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=1*(sampletable$phenoEABE_bca==2)
    idx=!is.na(pheno)& !is.na(E)
    pheno=pheno[idx]
    E1 <- E[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,E1=E1,genodat1=genodat1,covdat=covdat,job.num=njob)
    out1[,3]=matrix(unlist(tmp),ncol=1,byrow = T)
    
    pheno=sampletable$phenoEA_bca
    pheno[is.na(pheno)] = 3
    pheno[pheno==1] = NA
    pheno=1*(pheno==2)
    
    idx=!is.na(pheno)& !is.na(E)
    pheno=pheno[idx]
    E1 <- E[idx]
    genodat1=genodat[,idx]
    covdat=Covariate[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,genetable=genetable,pheno=pheno,E1=E1,genodat1=genodat1,covdat=covdat,job.num=njob)
    out1[,4]=matrix(unlist(tmp),ncol=1,byrow = T)
    out=rbind(out,out1)
  }
  return(out)
}



Sys.time()
iSKATout=mpi_compute_pvalue()
Sys.time()
sink()

save(iSKATout,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_iSKAT_output.RData") 

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_iSKAT_output.RData")

#iSKATout1 <- iSKATout[!is.na(iSKATout[,1]),]
#plot(-log((nrow(iSKATout1):1)/nrow(iSKATout1),base=10),-log(iSKATout1[order(iSKATout1[,1],decreasing=T),1],base=10),xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#iSKATout1 <- iSKATout[!is.na(iSKATout[,2]),]
#plot(-log((nrow(iSKATout1):1)/nrow(iSKATout1),base=10),-log(iSKATout1[order(iSKATout1[,2],decreasing=T),2],base=10),xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#iSKATout1 <- iSKATout[!is.na(iSKATout[,3]),]
#plot(-log((nrow(iSKATout1):1)/nrow(iSKATout1),base=10),-log(iSKATout1[order(iSKATout1[,3],decreasing=T),3],base=10),xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#iSKATout1 <- iSKATout[!is.na(iSKATout[,4]),]
#plot(-log((nrow(iSKATout1):1)/nrow(iSKATout1),base=10),-log(iSKATout1[order(iSKATout1[,4],decreasing=T),4],base=10),xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)