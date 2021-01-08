#!/usr/bin/env Rscript
#to apply SKAT-TWAS on GERD=1
#salloc -t 1-1 -n 51 -n 1 Rscript ./tryscat_GERD.R &
#salloc -t 6-1 -n 51 mpirun -n 1 R --interactive
  

#un-comment each datasets to run
#source("tryscat_GERD.R")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11" 
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11"
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11"

load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/bca_extractgenotype.RData"))
tmp=colnames(bcagenotype)
if (grepl("_",tmp[1]))
{
  tmp=strsplit(tmp,"_")
  tmp1=sapply(1:length(tmp),function(x){
    tmp1=tmp[[x]]
    paste0(tmp1[2:length(tmp1)],collapse = "_")
  })
  colnames(bcagenotype)=tmp1 #use localid
}

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
table(sampletable$recurrent_HB_RF)
# 0    1 
# 3528 2432 
table(sampletable$recurrent_HB_RF,sampletable$phenoBE_bca)
# 1    2
# 0 1446 1058
# 1  348 1186
table(sampletable$recurrent_HB_RF,sampletable$phenoEA_bca)
# 1    2
# 0 1448  965
# 1  349  854
table(sampletable$recurrent_HB_RF,sampletable$phenoEABE_bca)
# 1    2
# 0 1447 2022
# 1  348 2038

sampletable=sampletable[which(sampletable$recurrent_HB_RF==1),]

allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,colnames(bcagenotype))
bcagenotype=bcagenotype[,idx]
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid
idx=match(allsamples,sampletable$localid)

Covariate=sampletable[idx,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate)[which(colnames(Covariate)=="ev1_bca")]="pc1"
colnames(Covariate)[which(colnames(Covariate)=="ev2_bca")]="pc2"
colnames(Covariate)[which(colnames(Covariate)=="ev3_bca")]="pc3"
colnames(Covariate)[which(colnames(Covariate)=="ev4_bca")]="pc4"
rownames(Covariate)=sampletable$localid[idx]
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)
#Covariate$sex=factor(Covariate$sex)
library(SKAT)


skat_p=function(Z,Covariate,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1,na.rm=T)
  if (sum(tmp!=0) >1)
  {
    obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
    out2<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #0.0004209424
    res=out2$p.value
    return(res)
  }else
  {
    return(NA)
  }
}

compute_p_arow=function(i)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  gene=genes[i]
  res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA,stringsAsFactors = F)
  rownames(res)=gene
  selectedsnps=unlist(strsplit(genemodel$selectedsnps[which(rownames(genemodel)==gene)],"|",fixed=T))
  idx=match(selectedsnps,rownames(bcagenotype))
  idxtocorrect=which(is.na(idx))
  
  if (length(idxtocorrect)>0)
  {
    for (j in idxtocorrect)
    {
      tmp=unlist(strsplit(selectedsnps[j],"_"))
      selectedsnps[j]=paste0(tmp[c(1,3,2)],collapse = "_")
    }
    idx=match(selectedsnps,rownames(bcagenotype))
  }
  
  Z=t(bcagenotype[idx,,drop=F])
  if (length(idxtocorrect)>0)
  {
    Z[idxtocorrect,]=2-Z[idxtocorrect,]
  }
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bca==2) #case
  idx2=which(sampletable$phenoBE_bca==1)
  res$BE_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoEA_bca==1)
  res$EA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoEA_bca==2)
  idx2=which(sampletable$phenoBE_bca==2)
  res$BEA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoBE_bca==2 | sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoBE_bca==1)
  res$BEEA_p=skat_p(Z,Covariate,idx1,idx2)
  return(res)
}

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
if (mpi.comm.size()==0) mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(outfolder)
mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(allsamples)
mpi.bcast.Robj2slave(Covariate)
mpi.bcast.Robj2slave(skat_p)
mpi.bcast.Robj2slave(compute_p_arow)
mpi.remote.exec(library(SKAT))

#add chr
get_chr_model=function(genemodel=res_min)
{
  genemodel=genemodel[genemodel$glmflag==1,]
  genemodel$chr=genemodel$pos=NA
  for (i in 1:nrow(genemodel))
  {
    tmp=unlist(strsplit(genemodel$selectedsnps[i],"|",fixed=T))
    tmp=unlist(strsplit(tmp,":"))
    genemodel$chr[i]=tmp[1]
    tmp=unlist(strsplit(tmp[2],"_"))
    genemodel$pos[i]=as.numeric(tmp[1])
  }
  return(genemodel)
}
genemodel=get_chr_model()

library(GenomicRanges)
bcabim$V1=as.character(bcabim$V1)
gr_allbcagenotype=GRanges(seqnames = bcabim$V1,ranges = IRanges(start=bcabim$V4,width=1))

mpi_compute_p=function(genemodel=genemodel)
{
  genemodel=genemodel[genemodel$glmflag==1,]
  bcabim$V1=as.character(bcabim$V1)
  allgenemodel=genemodel
  allbcagenotype=bcagenotype
  resall=NULL
  chrs=unique(allgenemodel$chr)
  for (chr in chrs)
  {
    res=NULL
    print(paste0(chr,"---"))
    genemodel=allgenemodel[allgenemodel$chr==chr,]
    idx=order(genemodel$chr,genemodel$pos)
    genemodel=genemodel[idx,]
    gr_genemodel=GRanges(seqnames = genemodel$chr,ranges=IRanges(start=genemodel$pos,width = 1))
    mpi.bcast.Robj2slave(genemodel)
    genes=rownames(genemodel)
    rows=1:nrow(genemodel)
    #
    n=100 #work on a subset of genes each time, to save memory
    nchunks=ceiling(length(rows)/n)
    print(paste0("number of total:",nchunks))
    for (i in 1:nchunks)
    {
      #if (i %% 5==0) cat(i,"..")
      if (i<nchunks)
      {
        seq=rows[((i-1)*n+1):(i*n)]
      }else
      {
        seq=rows[((i-1)*n+1):length(rows)]
      }
      gr_chunk=GRanges(seqnames = chr,ranges=IRanges(start=max(c(genemodel$pos[seq[1]]-1e6,1)),end=genemodel$pos[seq[length(seq)]]+1e6)) #
      tmp=distance(gr_allbcagenotype,gr_chunk)
      idx=which(tmp==0)
      bcagenotype=allbcagenotype[idx,]
      mpi.bcast.Robj2slave(bcagenotype)
      tmp=mpi.parSapply(X=seq,FUN=compute_p_arow,job.num=njobs)
      if (length(tmp) %% 4 !=0) print(i)
      res1=as.data.frame(matrix(unlist(tmp),ncol=4,byrow = T))
      #rownames(res1)=rownames(genemodel)[seq]
      res=rbind(res,res1)
    }
    rownames(res)=rownames(genemodel)
    resall=rbind(resall,res)
  }
  
  return(resall)
}

skat_min2=mpi_compute_p(genemodel=genemodel)

#save skat result
save(skat_min2,file=paste0(outfolder,"/skat_GERD1_res.RData"))

# mpi.close.Rslaves()
# mpi.quit()
# quit()