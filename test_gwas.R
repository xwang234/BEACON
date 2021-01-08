
library(data.table)
if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  for (i in 1:ncol(sampletable))
  {
    idx=which(sampletable[,i]==-9)
    sampletable[idx,i]=NA
  }
}

fam=read.table("../result/bca_filtered_04Apr2019.fam",stringsAsFactors = F)
ped=fread("../result/bca_filtered_04Apr2019.traw",header = T)
ped=as.data.frame(ped)
all(colnames(ped)[7:ncol(ped)]==paste0(fam$V1,"_",fam$V2)) #T

COsamples=intersect(sampletable$localid[sampletable$phenoBE_bc==1 | sampletable$phenoEA_bc==1],fam$V2)
BEsamples=intersect(sampletable$localid[sampletable$phenoBE_bc==2],fam$V2)
EAsamples=intersect(sampletable$localid[sampletable$phenoEA_bc==2],fam$V2)

COsamples=intersect(COsamples,sampletable$localid[sampletable$recurrent_HB_RF==1])
BEsamples=intersect(BEsamples,sampletable$localid[sampletable$recurrent_HB_RF==1])
EAsamples=intersect(EAsamples,sampletable$localid[sampletable$recurrent_HB_RF==1])

allsamples=c(COsamples,BEsamples,EAsamples)
idx=match(allsamples,fam$V2)
peddat=ped[,6+idx]
idx=match(allsamples,sampletable$localid)
clinicaldat=sampletable[idx,]
readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  eigsamples=read.table(eigsampfile,stringsAsFactors = F)
  eigsamples=eigsamples$V2
  if (opt==1)
  {
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=eigsamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}

eigenstratmatrix=readeigenstrat(opt=2)
eigenstratmatrix=t(eigenstratmatrix)
idx=match(allsamples,rownames(eigenstratmatrix))
clinicaldat$pc1=eigenstratmatrix[idx,1]
clinicaldat$pc2=eigenstratmatrix[idx,2]
clinicaldat$pc3=eigenstratmatrix[idx,3]
clinicaldat$pc4=eigenstratmatrix[idx,4]
clinicaldat$sex=as.factor(clinicaldat$sex)

table(clinicaldat$phenoEABE_bc,useNA="ifany")
# 1    2 <NA> 
#   348 2038    3 
which(is.na(clinicaldat$phenoEABE_bc))
# [1]  141  496 2348
clinicaldat[c(141,496,2348),c(3,4,5)]
# phenoBE_bc phenoEA_bc phenoEABE_bc
# 4031         NA          1           NA
# 667           2         NA           NA
# 7958         NA          2           NA
clinicaldat$phenoEABE_bc[141]=1
clinicaldat$phenoEABE_bc[c(496,2348)]=2
table(clinicaldat$phenoEABE_bc,useNA="ifany")
# 1    2 
# 349 2040

save(peddat,clinicaldat,file="../result/reflux_gwas_dat.RData")
if (!exists("peddat"))
  load("../result/reflux_gwas_dat.RData")
pvalue=function(i=1)
{
  x=as.numeric(peddat[i,])
  fit1 <- glm(I(phenoEABE_bc==2)~x+ age + factor(sex) + pc1 + pc2 + pc3 + pc4,family=binomial,data=clinicaldat)
  res=NA
  if ("x" %in% rownames(summary(fit1)$coefficient))
    res=summary(fit1)$coefficient[2,4]
  return(res)
}

#salloc -t 1-1 -n 24 mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
library(Rmpi)
mpi_pvalue=function()
{
  njobs=mpi.universe.size() - 1
  print(njobs)
  if (mpi.comm.size()==0)
  {
    mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
  }
  mpi.remote.exec(load("../result/reflux_gwas_dat.RData"))
  mpi.bcast.Robj2slave(pvalue)
  len=10000
  n=nrow(peddat)
  nchunks=ceiling(n/len)
  print(paste0("total number of iterations: ",nchunks))
  allres=NULL
  for (i in 1:nchunks)
  {
      cat(i,'..')
      start=(i-1)*len+1
      end=min(start+len-1,n)
      res=mpi.parSapply(X=start:end,FUN=pvalue,job.num=njobs)
      allres=c(allres,res)
  }
  save(allres,file="../result/test_gwas_res.RData")
}