
############################################
#### genetic risk prediction for BE/EA #####
#### grant application 2020            #####
############################################

rm(list=ls())
#install.packages("/fh/fast/dai_j/CancerGenomics/Tools/Rpackages/plink2R-master/plink2R_1.1.tar.gz", repos = NULL, type = "source")
library(plink2R)
dat_bca=read_plink("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filteredMAF_20Feb2015",impute="avg")
bim=dat_bca$bim
genotype=dat_bca$bed 
fam=dat_bca$fam
rm(dat_bca)

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)

sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
#eigenstratmatrix=readeigenstrat()
allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid

# idx=match(allsamples,colnames(eigenstratmatrix))
# eigenstratmatrix=eigenstratmatrix[,idx]
# sampletable=cbind(sampletable,t(eigenstratmatrix[1:4,]))
Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3","pc4")]
Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)



nprobes <- ncol(genotype)
njob=90
library(Rmpi)



sink("~/SCHARPHOME/ScharpFile/Article/Grant/BEprogression/MPI_output.txt")

### Prepare Rmpi environment, DO NOT CHANGE ###
Sys.getenv(c("SLURM_SUBMIT_DIR"))
Sys.getenv(c("HOST", "SLURM_JOB_ID", "SLURM_NODELIST","SLURM_NNODES", "SLURM_NTASKS", 
             "SLUR M_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","SLURM_NTASKS_PER_NODE", 
             "SLURM_TASK_PID", "SLURM_ PARTITION"))
#
library('Rmpi')
nWorkers <- mpi.universe.size() - 1L
mpi.spawn.Rslaves(nslaves=nWorkers)
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









mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(genotype)
mpi.bcast.Robj2slave(Covariate)



compute_pvalue <- function(r,n.chunk,nprobes){
  
 
  nblock <- floor(nprobes/(n.chunk))
  ncount <- rep(nblock,n.chunk)
  if (nprobes>sum(ncount))  ncount[1:(nprobes-sum(ncount))] <- ncount[1:(nprobes-sum(ncount))]+ rep(1,nprobes-sum(ncount))
  chunks <- c(0,cumsum(ncount))
  chunkr <- seq(chunks[r]+1,chunks[r+1])
  nr <- length(chunkr)
  
  out <- matrix(NA,nr,8)
  for (i in 1:nr) {
  #if (i%%1000==0) cat(i,"..")
  dat1 <- cbind(sampletable$phenoBE_bca,genotype[,i],Covariate)
  dat1 <- data.frame(dat1[!is.na(dat1[,1]),])
  names(dat1)[1:2] <- c("pheno","genotype")
  fit1 <- glm(I(pheno==2)~genotype+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=dat1)
  out[i,1] <- summary(fit1)$coef[2,4]
  out[i,2] <- summary(fit1)$coef[2,1]
  
  dat2 <- cbind(sampletable$phenoEA_bca,genotype[,i],Covariate)
  dat2 <- data.frame(dat2[!is.na(dat2[,1]),])
  names(dat2)[1:2] <- c("pheno","genotype")
  fit2 <- glm(I(pheno==2)~genotype+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=dat2)
  out[i,3] <- summary(fit2)$coef[2,4]
  out[i,4] <- summary(fit2)$coef[2,1]
  
  dat3 <- cbind(sampletable$phenoEABE_bca,genotype[,i],Covariate)
  dat3 <- data.frame(dat3[!is.na(dat3[,1]),])
  names(dat3)[1:2] <- c("pheno","genotype")
  fit3 <- glm(I(pheno==2)~genotype+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=dat3)
  out[i,5] <- summary(fit3)$coef[2,4]
  out[i,6] <- summary(fit3)$coef[2,1]
  
  
  dat4 <- cbind(1-(is.na(sampletable$phenoEA_bca)| sampletable$phenoEA_bca==1),genotype[,i],Covariate)
  dat4 <- data.frame(dat4[!is.na(dat4[,1]),])
  names(dat4)[1:2] <- c("pheno","genotype")
  fit4 <- glm(pheno~genotype+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=dat4)
  out[i,7] <- summary(fit4)$coef[2,4]
  out[i,8] <- summary(fit4)$coef[2,1]
  }  
list(out)
}

mpi.bcast.Robj2slave(compute_pvalue)

#each run compute n rows

Sys.time()
outlist<-mpi.parSapply(1:njob,FUN=compute_pvalue,n.chunk=njob,nprobes,job.num=njob)
BCAout <- do.call(rbind,outlist)  


sink()

save(BCAout,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output.RData") 
mpi.close.Rslaves()
mpi.quit()
quit()


