
############################################
#### genetic risk prediction for BE/EA #####
#### grant application 2020            #####
############################################

#request large node 
#salloc -t 2-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript BCA_GWAS_bmi_interaction_MPI.R

rm(list=ls())
#install.packages("WGCNA") #bigtranspose
#install.packages("/fh/fast/dai_j/CancerGenomics/Tools/Rpackages/plink2R-master/plink2R_1.1.tar.gz", repos = NULL, type = "source")
# library(plink2R)
# dat_bca=read_plink("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filteredMAF_20Feb2015",impute="avg")
# bim=dat_bca$bim
# genotype=dat_bca$bed 
# fam=dat_bca$fam
# #rm(dat_bca)
# idx=bim$V1 %in% 1:23 #23:X
# bim=bim[idx,]
# genotype=genotype[,idx]
# allgenotype=genotype
# save(allgenotype,fam,bim,file="../result/bca_filteredMAF_20Feb2015.RData",version=2)
# library("WGCNA")
# tmp=transposeBigData(allgenotype, blocksize = 20000)
# all(as.numeric(tmp[2500,])==allgenotype[,2500])
# allgenotype=as.data.frame(tmp)
# save(allgenotype,fam,bim,file="../result/bca_filteredMAF_20Feb2015_t.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid

Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca","bmi_recent_healthy")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)



nprobes <- nrow(allgenotype)


sink("../result/MPI_output.txt") #changed-

### Prepare Rmpi environment, DO NOT CHANGE ###
Sys.getenv(c("SLURM_SUBMIT_DIR"))
Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST","SLURM_NNODES", "SLURM_NTASKS", 
             "SLUR M_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","SLURM_NTASKS_PER_NODE", 
             "SLURM_TASK_PID", "SLURM_ PARTITION"))
#
library('Rmpi')
njob <- mpi.universe.size() - 1L
mpi.spawn.Rslaves(nslaves=njob)
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

compute_pvalue=function(i=1,pheno,genodat,covdat)
{
  res=data.frame(p1=NA,beta1=NA,p2=NA,beta2=NA,p3=NA)
  genotype=as.numeric(genodat[i,])
  #covdat$cig_smk_ever <- 1-covdat$cig_smk_ever
  GE <- genotype*covdat$bmi_recent_healthy
  fit1 <- glm(I(pheno==2)~GE+genotype+bmi_recent_healthy+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=covdat)
  fit2 <- glm(I(pheno==2)~bmi_recent_healthy+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=covdat)
  res$p3 <- pchisq(fit2$deviance-fit1$deviance,df=2,lower.tail = F)
  if ("GE" %in% rownames(summary(fit1)$coef))
  {
    res$p1 <- summary(fit1)$coef[2,4]
    res$beta1 <- summary(fit1)$coef[2,1]
    res$p2 <- summary(fit1)$coef[3,4]
    res$beta2 <- summary(fit1)$coef[3,1]
  }
  return(res)
}

mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(Covariate)
mpi.bcast.Robj2slave(compute_pvalue)

mpi_compute_pvalue <- function(nsnps=1000){
  
  nchunks=ceiling(nprobes/nsnps)
  print(paste0("number of total:",nchunks))
  rows=1:nprobes
  out=NULL
  for (ii in 1:nchunks)
  {
    if (ii%% 50==0)cat(ii,"..")
    if (ii<nchunks)
    {
      seqs=rows[((ii-1)*nsnps+1):(ii*nsnps)]
    }else
    {
      seqs=rows[((ii-1)*nsnps+1):length(rows)]
    }
    out1 <- matrix(NA,length(seqs),20)
    rownames(out1)=rownames(allgenotype)[seqs]
    Covariate1 <- Covariate
    
    pheno=sampletable$phenoBE_bca
    genodat=allgenotype[seqs,]
    idx=!is.na(pheno) & !is.na(Covariate1$bmi_recent_healthy)
    pheno=pheno[idx]
    genodat=genodat[,idx]
    covdat=Covariate1[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,pheno=pheno,genodat=genodat,covdat=covdat,job.num=njob)
    out1[,1:5]=matrix(unlist(tmp),ncol=5,byrow = T)
    
    pheno=sampletable$phenoEA_bca
    genodat=allgenotype[seqs,]
    idx=!is.na(pheno)& !is.na(Covariate1$bmi_recent_healthy)
    pheno=pheno[idx]
    genodat=genodat[,idx]
    covdat=Covariate1[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,pheno=pheno,genodat=genodat,covdat=covdat,job.num=njob)
    out1[,6:10]=matrix(unlist(tmp),ncol=5,byrow = T)
    
    pheno=sampletable$phenoEABE_bca
    genodat=allgenotype[seqs,]
    idx=!is.na(pheno)& !is.na(Covariate1$bmi_recent_healthy)
    pheno=pheno[idx]
    genodat=genodat[,idx]
    covdat=Covariate1[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,pheno=pheno,genodat=genodat,covdat=covdat,job.num=njob)
    out1[,11:15]=matrix(unlist(tmp),ncol=5,byrow = T)
    
    pheno=sampletable$phenoEA_bca
    pheno[is.na(pheno)] = 3
    pheno[pheno==1] = NA
    genodat=allgenotype[seqs,]
    idx=!is.na(pheno)& !is.na(Covariate1$bmi_recent_healthy)
    pheno=pheno[idx]
    genodat=genodat[,idx]
    covdat=Covariate1[idx,]
    tmp=mpi.parSapply(X=1:length(seqs),FUN=compute_pvalue,pheno=pheno,genodat=genodat,covdat=covdat,job.num=njob)
    out1[,16:20]=matrix(unlist(tmp),ncol=5,byrow = T)
    
    out=rbind(out,out1)
  }
  return(out)
}



Sys.time()
BCAout=mpi_compute_pvalue()

sink()

save(BCAout,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_bmi_interaction_output.RData") 


## quit program
mpi.close.Rslaves()
mpi.quit()
quit()



#rm(list=ls())
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_bmi_interaction_output.RData")

#mpoint <- ceiling(quantile(1:nrow(BCAout),0.95)/200)
#pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(BCAout))

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,5],decreasing=T),5],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,6],decreasing=T),6],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,8],decreasing=T),8],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,10],decreasing=T),10],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,11],decreasing=T),11],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,13],decreasing=T),13],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,15],decreasing=T),15],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,20],decreasing=T),20],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

#plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,7],decreasing=T),7],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)
