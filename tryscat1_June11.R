#!/usr/bin/env Rscript
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript ./tryscat1_June11.R &
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 R --interactive

#un-comment each datasets to run
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11" 
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11"

library(SKAT)
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
#load(paste0(outfolder,"/bca_predict_geneexp.RData"))

#SKAT function
#SKATBinary called SKAT, it also includes quality checks on genotypes (impute and flip)
# skat_p=function(Z,Covariate,idx1,idx2)
# {
#   y=c(rep(1,length(idx1)),rep(0,length(idx2)))
#   
#   Z1=Z[c(idx1,idx2),,drop=F]
#   Covariate1=as.matrix(Covariate[c(idx1,idx2),])
#   tmp=rowSums(Z1)
#   if (sum(tmp!=0) >1)
#   {
#     obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
#     #out<-SKAT(Z1, obj.s,weights.beta=c(1,25),r.corr=0) #0.04183234
#     #out<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #1.207545e-05
#     #out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
#     #out=SKATBinary(Z1, obj.s, weights.beta=c(1,25),kernel = "linear.weighted",r.corr=0) #0.04183234
#     out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
#     #out=SKATBinary(Z1, obj.s, method="SKATO", weights.beta=c(1,1),kernel = "linear.weighted") #4.244548e-05
#     #out<-SKATBinary(Z1, obj.s, weights.beta=c(1,25),kernel = "IBS.weighted",r.corr=0) #0.06059779
#     #out<-SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "2wayIX",r.corr=0) #5.994431e-05
#     #out<-SKAT_CommonRare(Z1, obj.s) #0.000466221
#     return(out$p.value)
#   }else
#   {
#     return(NA)
#   }
# }

#check multipe options:
#weights.beta:We suggest setting a1 = 1 and a2 = 25 because it increases the weight of rare variants while still putting decent nonzero weights for variants with MAF 1%–5%.
#             a1 = a2 = 1 corresponds to wj = 1, that is all variants are weighted equally
#r.corr: ρ=0 results in the original weighted linear kernel SKAT, and ρ=1 results in the weighted burden test.
#method: "SKATO" and "optimal.adj" represent a SKAT-O based on an unified approach, and "optimal" is an old version of the implementation of SKAT-O.
# skat_ps=function(Z,Covariate,idx1,idx2)
# {
#   y=c(rep(1,length(idx1)),rep(0,length(idx2)))
#   
#   Z1=Z[c(idx1,idx2),,drop=F]
#   Covariate1=as.matrix(Covariate[c(idx1,idx2),])
#   tmp=rowSums(Z1,na.rm=T)
#   if (sum(tmp!=0) >1)
#   {
#     obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
#     out1<-SKAT(Z1, obj.s,weights.beta=c(1,25),r.corr=0,is_dosage=T) #0.0007251699
#     out2<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #0.0004209424
#     out3<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") #0.001014649 
#     out4<-SKAT_CommonRare(Z1, obj.s) #0.001888707
#     res=data.frame(p1=out1$p.value,p2=out2$p.value,p3=out3$p.value,p4=out4$p.value,stringsAsFactors = F)
#     return(res)
#   }else
#   {
#     return(data.frame(p1=NA,p2=NA,p3=NA,p4=NA,stringsAsFactors = F))
#   }
# }

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
#genomewide:
# sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
# sampletable$site[sampletable$site=="NA"]=NA
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
#eigenstratmatrix=readeigenstrat()
allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,colnames(bcagenotype))
bcagenotype=bcagenotype[,idx]
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid
idx=match(allsamples,sampletable$localid)

# idx=match(allsamples,colnames(eigenstratmatrix))
# eigenstratmatrix=eigenstratmatrix[,idx]
# sampletable=cbind(sampletable,t(eigenstratmatrix[1:4,]))
#Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3","pc4")]
Covariate=sampletable[idx,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate)[which(colnames(Covariate)=="ev1_bca")]="pc1"
colnames(Covariate)[which(colnames(Covariate)=="ev2_bca")]="pc2"
colnames(Covariate)[which(colnames(Covariate)=="ev3_bca")]="pc3"
colnames(Covariate)[which(colnames(Covariate)=="ev4_bca")]="pc4"
rownames(Covariate)=sampletable$localid[idx]
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)


#genemodel,bcagenotype,Covariate were loaded,used for mpi
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
# #add weights
# compute_p_withcoeff_arow=function(i,genemodel=res_min,bcagenotype,Covariate)
# {
#   genes=rownames(genemodel)[genemodel$glmflag==1]
#   gene=genes[i]
#   res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA,stringsAsFactors = F)
#   rownames(res)=gene
#   selectedsnps=unlist(strsplit(genemodel$selectedsnps[which(rownames(genemodel)==gene)],"|",fixed=T))
#   selectedcoefs=as.numeric(unlist(strsplit(genemodel$selectedsnps_coeff[which(rownames(genemodel)==gene)],"|",fixed=T)))
#   idx=match(selectedsnps,rownames(bcagenotype))
#   idxtocorrect=which(is.na(idx))
#   
#   if (length(idxtocorrect)>0)
#   {
#     for (j in idxtocorrect)
#     {
#       tmp=unlist(strsplit(selectedsnps[j],"_"))
#       selectedsnps[j]=paste0(tmp[c(1,3,2)],collapse = "_")
#     }
#     idx=match(selectedsnps,rownames(bcagenotype))
#   }
#   
#   Z=t(bcagenotype[idx,,drop=F])
#   if (length(idxtocorrect)>0)
#   {
#     Z[idxtocorrect,]=2-Z[idxtocorrect,]
#   }
#   rownames(Z)=colnames(bcagenotype)
#   idx=match(allsamples,rownames(Z))
#   Z=Z[idx,,drop=F]
#   Z=t(t(Z)*abs(selectedcoefs)) #add the coeff here
#   #print(rankMatrix(Z)[[1]])
#   idx1=which(sampletable$phenoBE_bc==2) #case
#   idx2=which(sampletable$phenoBE_bc==1)
#   res$BE_p=skat_p(Z,Covariate,idx1,idx2)
#   idx1=which(sampletable$phenoEA_bc==2) #case
#   idx2=which(sampletable$phenoEA_bc==1)
#   res$EA_p=skat_p(Z,Covariate,idx1,idx2)
#   idx1=which(sampletable$phenoBE_bc==2)
#   idx2=which(sampletable$phenoEA_bc==2)
#   res$BEA_p=skat_p(Z,Covariate,idx1,idx2)
#   idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
#   idx2=which(sampletable$phenoBE_bc==1)
#   res$BEEA_p=skat_p(Z,Covariate,idx1,idx2)
#   return(res)
# }


# #run sequentially
# compute_p_allrows=function(genemodel=res_min,bcagenotype,Covariate,opt=1)
# {
#   genes=rownames(genemodel)[genemodel$glmflag==1]
#   res=data.frame(BE_p=rep(NA,length(genes)),EA_p=rep(NA,length(genes)),BEA_p=rep(NA,length(genes)),BEEA_p=rep(NA,length(genes)),stringsAsFactors = F)
#   print(paste0("total genes: ",length(genes)))
#   for (i in 1:length(genes))
#   {
#     if (i %%100 ==0) cat(i,"..")
#     res[i,]=compute_p_arow(i,genemodel,bcagenotype,Covariate,opt=opt)
#   }
#   rownames(res)=genes
#   return(res)
# }
# skat_min=compute_p_allrows(genemodel=res_min,bcagenotype,Covariate)
# skat_1se=compute_p_allrows(genemodel = res_1se,bcagenotype,Covariate)

#parallel version
#salloc -t 1-1 --mem-per-cpu 48G -n 15 --partition=largenode mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
#salloc -t 3-1 -n 75 mpirun -n 1 R --interactive

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(outfolder)
#mpi.remote.exec(load(paste0(outfolder,"/preidiction_michigan_model.RData")))
#mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(allsamples)
mpi.bcast.Robj2slave(Covariate)
mpi.bcast.Robj2slave(skat_p)
mpi.bcast.Robj2slave(compute_p_arow)
#mpi.bcast.Robj2slave(compute_p_withcoeff_arow)
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
    n=200 #work on a subset of genes each time, to save memory
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
save(skat_min2,file=paste0(outfolder,"/skat_res.RData"))

mpi.close.Rslaves()
mpi.quit()
quit()



