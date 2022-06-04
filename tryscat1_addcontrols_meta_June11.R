#!/usr/bin/env Rscript
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript ./tryscat1_June11.R &
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 R --interactive

#un-comment each datasets to run
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11" 
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
# #outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19"
library(SKAT)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/bca_extractgenotype.RData"))
# tmp=colnames(bcagenotype)
# if (grepl("_",tmp[1]))
# {
#   tmp=strsplit(tmp,"_")
#   tmp1=sapply(1:length(tmp),function(x){
#     tmp1=tmp[[x]]
#     paste0(tmp1[2:length(tmp1)],collapse = "_")
#   })
#   colnames(bcagenotype)=tmp1 #use localid
# }
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
  Covariate1=Covariate[c(idx1,idx2),]
  tmp=rowSums(Z1,na.rm=T)
  if (sum(tmp!=0) >1)
  {
    obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
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

#add covariate table (pc1-pc4)
readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype.pedind",
                        nskip=16)
{
  eigsamples=read.table(eigsampfile,stringsAsFactors = F)
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  colnames(tmp)=paste0("pc",1:ncol(tmp))
  rownames(tmp)=eigsamples$V2
  tmp$sex="M"
  tmp$sex[eigsamples$V5==2]="F"
  return(tmp)
}

#Cambridge
covariatetable=readeigenstrat(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC.pca",
                              eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC.pedind")
#BEACON
covariatetable=readeigenstrat(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.pca",
                              eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.pedind")
rownames(covariatetable)=gsub("SEP","",rownames(covariatetable))
#all(colnames(predict_min)[3:ncol(predict_min)] %in% rownames(covariatetable))
#add case/control
tmp=covariatetable
tmp$phenoBE_bca=tmp$phenoEA_bca=tmp$phenoEABE_bca=1
comsamples=intersect(sampletable$localid,rownames(tmp))
idx1=match(comsamples,rownames(tmp))
idx2=match(comsamples,sampletable$localid)
tmp$phenoBE_bca[idx1]=sampletable$phenoBE_bca[idx2]
tmp$phenoEA_bca[idx1]=sampletable$phenoEA_bca[idx2]
tmp$phenoEABE_bca[idx1]=sampletable$phenoEABE_bca[idx2]
tmp$phenoBE_bca[tmp$phenoBE_bca==-9]=NA
tmp$phenoEA_bca[tmp$phenoEA_bca==-9]=NA
tmp$phenoEABE_bca[tmp$phenoEABE_bca==-9]=NA
covariatetable=tmp
covariatetable$sex=factor(covariatetable$sex)

allsamples=intersect(colnames(bcagenotype),rownames(covariatetable))
idx=match(allsamples,colnames(bcagenotype))
bcagenotype=bcagenotype[,idx]
idx=match(allsamples,rownames(covariatetable))
covariatetable=covariatetable[idx,]


# idx=match(allsamples,colnames(eigenstratmatrix))
# eigenstratmatrix=eigenstratmatrix[,idx]
# sampletable=cbind(sampletable,t(eigenstratmatrix[1:4,]))
#Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3","pc4")]
#try different number of PCs
Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","sex")]
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","sex")]
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","sex")]

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
    Z[,idxtocorrect]=2-Z[,idxtocorrect]
  }
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  res$BE_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoEA_bca==1)
  res$EA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoEA_bca==2)
  idx2=which(covariatetable$phenoBE_bca==2)
  res$BEA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoBE_bca==2 | covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
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
#ml OpenMPI/3.1.4-GCC-8.3.0
#salloc -t 3-1 --constraint=gizmok -n 15 mpirun -n 1 R --interactive

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(outfolder)
#mpi.remote.exec(load(paste0(outfolder,"/preidiction_michigan_model.RData")))
#mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(covariatetable)
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

all(rownames(Covariate)==colnames(bcagenotype))
all(rownames(covariatetable)==colnames(bcagenotype))
skat_min2_cambridge=mpi_compute_p(genemodel=genemodel)
save(skat_min2_cambridge,file=paste0(outfolder,"/skat_meta_res.RData"))
skat_min2_beacon=mpi_compute_p(genemodel=genemodel)

#save skat result
save(skat_min2_cambridge,skat_min2_beacon,file=paste0(outfolder,"/skat_meta_res.RData"))
quantile(skat_min2[,1])
quantile(skat_min2[,2])
quantile(skat_min2[,3])
quantile(skat_min2[,4])

mpi.close.Rslaves()
mpi.quit()
quit()

meta_skatp=function(dat1=skat_min2_cambridge,dat2=skat_min2_beacon)
{
  comgenes=intersect(rownames(dat1),rownames(dat2))
  if (length(comgenes) <nrow(dat1)) warning("something could be wrong")
  dat1=dat1[match(comgenes,rownames(dat1)),]
  dat2=dat2[match(comgenes,rownames(dat2)),]
  resp=respadj=data.frame(matrix(NA,nrow=nrow(dat1),ncol=ncol(dat1)))
  rownames(resp)=rownames(respadj)=rownames(dat1)
  for (i in 1:ncol(dat1))
  {
    tmp=cbind(dat1[,i],dat2[,i])
    tmp1=fisher.method(tmp, na.rm=TRUE)
    resp[,i]=tmp1$p.value
    respadj[,i]=tmp1$p.adj
  }
  return(list(resp=resp,respadj=respadj))
  
}
meta_p=meta_skatp()
save(skat_min2_cambridge,skat_min2_beacon,meta_p,file=paste0(outfolder,"/skat_meta_res.RData"))
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
qqplot_twas=function(dat=skat_min2,r2cutoff=0)
{
  skat_min2_code=dat[rownames(dat) %in% proteingenes,]
  idx=match(rownames(skat_min2_code),rownames(res_min))
  skat_min2_code=skat_min2_code[res_min$r2[idx]>=r2cutoff,]
  
  par(mar=c(6,6,3,1))
  par(mfrow=c(2,2))
  qqplot(skat_min2_code[,1])
  qqplot(skat_min2_code[,2])
  qqplot(skat_min2_code[,3])
  qqplot(skat_min2_code[,4])
}
qqplot_twas(dat=meta_p$resp)
qqplot_twas(dat=meta_p$resp,r2cutoff = 0.05)
qqplot_twas(dat=meta_p$resp,r2cutoff = 0.1)
qqplot_twas(dat=skat_min2_beacon)
qqplot_twas(dat=skat_min2_cambridge)