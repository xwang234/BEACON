
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_April18" #models stored for TCGA, came from prediction_michigan_models7_TCGA.R
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18" #models stored for GTEx, came from prediction_michigan_models6_GTEx.R
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_April18" #models stored for GTEx mucosa


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

# if (!exists("sampletable"))
# {
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  # sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  # sampletable$site[sampletable$site=="NA"]=NA
  library(readxl)
  sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
  # geneexpsamplenames=strsplit(colnames(predict_min)[3:ncol(predict_min)],"_")
  # geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  #   tmp=geneexpsamplenames[[x]]
  #   paste0(tmp[2:length(tmp)],collapse = "_")
  # })
# }

# readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
#                         eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
#                         thesamples=geneexpsamplenames,nskip=16,opt=1)
# {
#   tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
#   if (opt==1)
#   {
#     eigsamples=read.table(eigsampfile,stringsAsFactors = F)
#     eigsamples=eigsamples$V2
#     idx=match(thesamples,eigsamples)
#     tmp1=tmp[idx,]
#     tmp1=as.data.frame(t(tmp1))
#     colnames(tmp1)=thesamples
#   }else #don't need to change sample names
#   {
#     tmp1=as.data.frame(t(tmp))
#     colnames(tmp1)=thesamples
#   }
#   rownames(tmp1)=paste0("pc",1:nrow(tmp1))
#   return(tmp1)
# }
# 
# eigenstratmatrix=readeigenstrat()

removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  dat=dat[!idxconst,]
}

#SKAT function
#SKATBinary called SKAT, it also includes quality checks on genotypes (impute and flip)
skat_p=function(Z,Covariate,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  if (sum(tmp!=0) >1)
  {
    obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
    #out<-SKAT(Z1, obj.s,weights.beta=c(1,25),r.corr=0) #0.04183234
    #out<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #1.207545e-05
    #out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
    #out=SKATBinary(Z1, obj.s, weights.beta=c(1,25),kernel = "linear.weighted",r.corr=0) #0.04183234
    out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
    #out=SKATBinary(Z1, obj.s, method="SKATO", weights.beta=c(1,1),kernel = "linear.weighted") #4.244548e-05
    #out<-SKATBinary(Z1, obj.s, weights.beta=c(1,25),kernel = "IBS.weighted",r.corr=0) #0.06059779
    #out<-SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "2wayIX",r.corr=0) #5.994431e-05
    #out<-SKAT_CommonRare(Z1, obj.s) #0.000466221
    return(out$p.value)
  }else
  {
    return(NA)
  }
}

#check multipe options:
#weights.beta:We suggest setting a1 = 1 and a2 = 25 because it increases the weight of rare variants while still putting decent nonzero weights for variants with MAF 1%–5%.
#             a1 = a2 = 1 corresponds to wj = 1, that is all variants are weighted equally
#r.corr: ρ=0 results in the original weighted linear kernel SKAT, and ρ=1 results in the weighted burden test.
#method: "SKATO" and "optimal.adj" represent a SKAT-O based on an unified approach, and "optimal" is an old version of the implementation of SKAT-O.
skat_ps=function(Z,Covariate,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1,na.rm=T)
  if (sum(tmp!=0) >1)
  {
    obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
    out1<-SKAT(Z1, obj.s,weights.beta=c(1,25),r.corr=0,is_dosage=T) #0.0007251699
    out2<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #0.0004209424
    out3<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") #0.001014649 
    out4<-SKAT_CommonRare(Z1, obj.s) #0.001888707
    res=data.frame(p1=out1$p.value,p2=out2$p.value,p3=out3$p.value,p4=out4$p.value,stringsAsFactors = F)
    return(res)
  }else
  {
    return(data.frame(p1=NA,p2=NA,p3=NA,p4=NA,stringsAsFactors = F))
  }
}

#genomewide:
# sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
# sampletable$site[sampletable$site=="NA"]=NA
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
#eigenstratmatrix=readeigenstrat()
allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid
idx=match(allsamples,colnames(bcagenotype))
bcagenotype=bcagenotype[,idx]
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


# #opt: which option of skat to use (1-4)
# compute_p_arow_old=function(i,genemodel=res_min,bcagenotype,Covariate,opt=1)
# {
#   genes=rownames(genemodel)[genemodel$glmflag==1]
#   gene=genes[i]
#   res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA,stringsAsFactors = F)
#   rownames(res)=gene
#   selectedsnps=unlist(strsplit(genemodel$selectedsnps[which(rownames(genemodel)==gene)],"|",fixed=T))
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
#    Z=t(bcagenotype[idx,,drop=F])
#   if (length(idxtocorrect)>0)
#   {
#     Z[idxtocorrect,]=2-Z[idxtocorrect,]
#   }
#   rownames(Z)=colnames(bcagenotype)
#   idx=match(allsamples,rownames(Z))
#   Z=Z[idx,,drop=F]
#   #print(rankMatrix(Z)[[1]])
#   idx1=which(sampletable$phenoBE_bc==2) #case
#   idx2=which(sampletable$phenoBE_bc==1)
#   tmp=skat_ps(Z,Covariate,idx1,idx2)
#   res$BE_p=tmp[,opt]
#   idx1=which(sampletable$phenoEA_bc==2) #case
#   idx2=which(sampletable$phenoEA_bc==1)
#   tmp=skat_ps(Z,Covariate,idx1,idx2)
#   res$EA_p=tmp[,opt]
#   idx1=which(sampletable$phenoBE_bc==2)
#   idx2=which(sampletable$phenoEA_bc==2)
#   tmp=skat_p(Z,Covariate,idx1,idx2)
#   res$BEA_ps=tmp[,opt]
#   idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
#   idx2=which(sampletable$phenoBE_bc==1)
#   tmp=skat_ps(Z,Covariate,idx1,idx2)
#   res$BEEA_p=tmp[,opt]
#   return(res)
# }

#genemodel,bcagenotype,Covariate were loaded,used for mpi
compute_p_arow=function(i,opt=1)
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
  tmp=skat_ps(Z,Covariate,idx1,idx2)
  res$BE_p=tmp[,opt]
  idx1=which(sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoEA_bca==1)
  tmp=skat_ps(Z,Covariate,idx1,idx2)
  res$EA_p=tmp[,opt]
  idx1=which(sampletable$phenoBE_bca==2)
  idx2=which(sampletable$phenoEA_bca==2)
  tmp=skat_ps(Z,Covariate,idx1,idx2)
  res$BEA_p=tmp[,opt]
  idx1=which(sampletable$phenoBE_bca==2 | sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoBE_bca==1)
  tmp=skat_ps(Z,Covariate,idx1,idx2)
  res$BEEA_p=tmp[,opt]
  return(res)
}
#add weights
compute_p_withcoeff_arow=function(i,genemodel=res_min,bcagenotype,Covariate)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  gene=genes[i]
  res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA,stringsAsFactors = F)
  rownames(res)=gene
  selectedsnps=unlist(strsplit(genemodel$selectedsnps[which(rownames(genemodel)==gene)],"|",fixed=T))
  selectedcoefs=as.numeric(unlist(strsplit(genemodel$selectedsnps_coeff[which(rownames(genemodel)==gene)],"|",fixed=T)))
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
  Z=t(t(Z)*abs(selectedcoefs)) #add the coeff here
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  res$BE_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoEA_bc==2) #case
  idx2=which(sampletable$phenoEA_bc==1)
  res$EA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoBE_bc==2)
  idx2=which(sampletable$phenoEA_bc==2)
  res$BEA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  res$BEEA_p=skat_p(Z,Covariate,idx1,idx2)
  return(res)
}


#run sequentially
compute_p_allrows=function(genemodel=res_min,bcagenotype,Covariate,opt=1)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  res=data.frame(BE_p=rep(NA,length(genes)),EA_p=rep(NA,length(genes)),BEA_p=rep(NA,length(genes)),BEEA_p=rep(NA,length(genes)),stringsAsFactors = F)
  print(paste0("total genes: ",length(genes)))
  for (i in 1:length(genes))
  {
    if (i %%100 ==0) cat(i,"..")
    res[i,]=compute_p_arow(i,genemodel,bcagenotype,Covariate,opt=opt)
  }
  rownames(res)=genes
  return(res)
}
skat_min=compute_p_allrows(genemodel=res_min,bcagenotype,Covariate)
skat_1se=compute_p_allrows(genemodel = res_1se,bcagenotype,Covariate)

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
mpi.bcast.Robj2slave(removeconstrows)
mpi.bcast.Robj2slave(skat_ps)
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

mpi_compute_p=function(genemodel=genemodel,opt=1)
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
      tmp=mpi.parSapply(X=seq,FUN=compute_p_arow,opt=opt,job.num=njobs)
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

skat_min1=mpi_compute_p(genemodel=genemodel,opt=1)
skat_min2=mpi_compute_p(genemodel=genemodel,opt=2)
skat_min3=mpi_compute_p(genemodel=genemodel,opt=3)
skat_min4=mpi_compute_p(genemodel=genemodel,opt=4)
#save skat result
save(skat_min1,skat_min2,skat_min3,skat_min4,file=paste0(outfolder,"/skat_res.RData"))
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
idx=skat_min$gene %in% proteingenes
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
# #multiply coeff
# mpi_compute_p_coeff=function(genemodel=res_min)
# {
#   genemodel=genemodel[genemodel$glmflag==1,]
#   genes=rownames(genemodel)
#   rows=1:nrow(genemodel)
#   res=NULL
#   n=njobs
#   nchunks=ceiling(length(rows)/n)
#   print(paste0("number of total:",nchunks))
#   for (i in 1:nchunks)
#   {
#     if (i %% 10==0) cat(i,"..")
#     if (i<nchunks)
#     {
#       seq=rows[((i-1)*n+1):(i*n)]
#     }else
#     {
#       seq=rows[((i-1)*n+1):length(rows)]
#     }
#     tmp=mpi.parSapply(X=seq,FUN=compute_p_withcoeff_arow,genemodel=genemodel,bcagenotype=bcagenotype,X=X,job.num=njobs)
#     if (length(tmp) %% 4 !=0) print(i)
#     res1=matrix(unlist(tmp),ncol=4,byrow = T)
#     #rownames(res1)=rownames(genemodel)[seq]
#     res=rbind(res,res1)
#   }
#   res=as.data.frame(res)
#   rownames(res)=rownames(genemodel)
#   return(res)
# }
# skat_1se=mpi_compute_p_coeff(genemodel = res_1se)
# skat_min=mpi_compute_p_coeff(genemodel = res_min)
# save(skat_1se,skat_min,file=paste0(outfolder,"/skat_coef_abs_res.RData"))

# numsnps=res_min$numselectedsnp[match(skat_min_code$gene,rownames(res_min))]
# idx1=skat_min$BEEA_p[idx]<1e-2
# t.test(numsnps[idx1],numsnps[!idx1])
# quantile(numsnps[!idx1],na.rm=T)
# quantile(numsnps[idx1],na.rm=T)
# 
# numsnps=res_min$numselectedsnp[match(skat_min$gene,rownames(res_min))]
# idx1=skat_min$BEEA_p<5e-2
# t.test(numsnps[idx1],numsnps[!idx1])
# quantile(numsnps[!idx1],na.rm=T)
# quantile(numsnps[idx1],na.rm=T)
# rownames(skat_min_code)=skat_min_code$gene #for GTEx
# rownames(skat_min)=skat_min$gene
compute_fwer_fdr=function(dat=skat_min_code,cutoff=0.05)
{
  compute_fdr=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method="fdr")
    if (sum(tmp<cutoff1,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff)]
      res=res[order(res)]
    }
    return(res)
  }
  compute_fwer=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method = "bonferroni")
    if (sum(tmp<cutoff,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff1)]
      res=res[order(res)]
    }
    return(res)
  }
  BE_fdr=compute_fdr(dat,namecol="BE_p")
  BE_fwer=compute_fwer(dat,namecol="BE_p")
  EA_fdr=compute_fdr(dat,namecol="EA_p")
  EA_fwer=compute_fwer(dat,namecol="EA_p")
  BEA_fdr=compute_fdr(dat,namecol="BEA_p")
  BEA_fwer=compute_fwer(dat,namecol="BEA_p")
  BEEA_fdr=compute_fdr(dat,namecol="BEEA_p")
  BEEA_fwer=compute_fwer(dat,namecol="BEEA_p")
  return(list(BE_fdr=BE_fdr,BE_fwer=BE_fwer,EA_fdr=EA_fdr,EA_fwer=EA_fwer,
              BEA_fdr=BEA_fdr,BEA_fwer=BEA_fwer,BEEA_fdr=BEEA_fdr,BEEA_fwer=BEEA_fwer))
}
fdr_fwer_res=compute_fwer_fdr() #GTEx
# $BE_fdr
# KXD1      CERS1     SGK223      ARMC6     OR11A1      DDX49       COMP    SPAG11A    PLEKHF2 
# 0.01784353 0.01784353 0.02002136 0.02816466 0.03244586 0.03244586 0.03847389 0.04764984 0.04764984 
# PPP1R3B 
# 0.04971261 
# 
# $BE_fwer
# KXD1      CERS1 
# 0.02420133 0.03568706 
# 
# $EA_fdr
# KXD1         COMP 
# 0.0004479597 0.0357644268 
# 
# $EA_fwer
# KXD1 
# 0.0004479597 
# 
# $BEA_fdr
# NULL
# 
# $BEA_fwer
# NULL
# 
# $BEEA_fdr
# KXD1        CERS1         COMP        DDX49       SGK223        ARMC6        OR2Z1       KLHL26 
# 7.218799e-05 3.502100e-03 4.079742e-03 6.137244e-03 7.944033e-03 7.944033e-03 1.084968e-02 1.084968e-02 
# C8orf49        CRTC1      SPAG11A      DEFB136      PLEKHF2        GATA4         CTSB 
# 2.593635e-02 3.252995e-02 3.302720e-02 4.080626e-02 4.080626e-02 4.466331e-02 4.466331e-02 
# 
# $BEEA_fwer
# KXD1        CERS1         COMP        DDX49       SGK223        ARMC6 
# 7.218799e-05 7.004200e-03 1.223923e-02 2.454898e-02 4.244696e-02 4.766420e-02 
fdr_fwer_res=compute_fwer_fdr(dat=skat_min2) #outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_PC4"
# $BE_fdr
# DDX49      MEF2B 
# 0.02023329 0.02042122 
# 
# $BE_fwer
# DDX49      MEF2B 
# 0.02023329 0.04084244 
# 
# $EA_fdr
# DDX49      MEF2B 
# 0.02929032 0.02929032 
# 
# $EA_fwer
# MEF2B 
# 0.04738151 
# 
# $BEA_fdr
# NULL
# 
# $BEA_fwer
# NULL
# 
# $BEEA_fdr
# DDX49        MEF2B 
# 0.0009690309 0.0009690309 
# 
# $BEEA_fwer
# DDX49       MEF2B 
# 0.001311554 0.001938062 

updategenenames=function(genes)
{
  genemap=data.frame(old="LASS1",new="CERS1",stringsAsFactors = F)
  genemap=rbind.data.frame(genemap,data.frame(old="C19orf50",new="KXD1",stringsAsFactors = F))
  #add rows if necessary
  
  res=genes
  genes1=intersect(genes,genemap$old)
  if (length(genes1)>0)
  {
    idx=match(genes1, genes)
    idx1=match(genes1,genemap$old)
    genes[idx]=genemap$new[idx1]
  }
  return(genes)
}


library(GenomicRanges)
computeeQTL1=function(selectedsnps,snp,phenotype,covariate,gene="DDX49")
{
  res=data.frame(snp=selectedsnps,eqtl=NA,stringsAsFactors = F)
  j=which(rownames(phenotype)==gene)
  selectedsnps_alt=sapply(selectedsnps,function(x){
    tmp=unlist(strsplit(x,"_"))
    paste0(tmp[c(1,3,2)],collapse = "_")
  })
  idx3=which(selectedsnps_alt %in% rownames(snp))
  if (length(idx3)>0) selectedsnps[idx3]=selectedsnps_alt[idx3]
  idx3=which(selectedsnps %in% rownames(snp))
  if (length(j)>0 & length(idx3)>0)
  {
    Y=unlist(phenotype[j,])
    for (i in idx3)
    {
      X=t(snp[which(rownames(snp)==selectedsnps[i]),,drop=F])
      Xall=data.matrix(cbind(X,covariate))
      fit1=lm(Y~Xall)
      res$eqtl[i]=summary(fit1)$coefficients[2,4]
    }
  }
  return(res)
}

computeeQTL2=function(selectedsnps,snp,phenotype,copynumber,mutation,covariate,gene="DDX49") #TCGA
{
  res=data.frame(snp=selectedsnps,eqtl=NA,stringsAsFactors = F)
  j=which(rownames(phenotype)==gene)
  selectedsnps_alt=sapply(selectedsnps,function(x){
    tmp=unlist(strsplit(x,"_"))
    paste0(tmp[c(1,3,2)],collapse = "_")
  })
  idx3=which(selectedsnps_alt %in% rownames(snp))
  if (length(idx3)>0) selectedsnps[idx3]=selectedsnps_alt[idx3]
  idx3=which(selectedsnps %in% rownames(snp))
  if (length(j)>0 & length(idx3)>0)
  {
    Y=unlist(phenotype[j,])
    for (i in idx3)
    {
      X=t(snp[which(rownames(snp)==selectedsnps[i]),,drop=F])
      mutationV=unlist(mutation[j,])
      if (all(mutationV==0)) #no mutation
      {
        Xall=data.matrix(cbind(X,cn=unlist(copynumber[j,]),covariate))
        covariateall=data.matrix(cbind(cn=unlist(copynumber[j,]),covariate))
      }else
      {
        Xall=data.matrix(cbind(X,cn=unlist(copynumber[j,]),mutation=mutationV,covariate))
        covariateall=data.matrix(cbind(cn=unlist(copynumber[j,]),mutation=mutationV,covariate))
      }
      fit1=lm(Y~Xall)
      res$eqtl[i]=summary(fit1)$coefficients[2,4]
    }
  }
  return(res)
}

update_bcasamplenames=function(bcasamples=colnames(bcagenotype))
{
  if (grepl("_",bcasamples[1]))
  {
    bcasamples=sapply(bcasamples,function(x){
      tmp=unlist(strsplit(x,"_"))
      paste0(tmp[2:length(tmp)],collapse = "_")
    })
  }
  return(bcasamples)
}

# sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data//bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
# sampletable$site[sampletable$site=="NA"]=NA
# for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
# readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
#                         eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
#                         thesamples=colnames(bcagenotype),nskip=16,opt=1)
# {
#   tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
#   if (opt==1)
#   {
#     eigsamples=read.table(eigsampfile,stringsAsFactors = F)
#     eigsamples=eigsamples$V2
#     idx=match(thesamples,eigsamples)
#     tmp1=tmp[idx,]
#     tmp1=as.data.frame(t(tmp1))
#     colnames(tmp1)=thesamples
#   }else #don't need to change sample names
#   {
#     tmp1=as.data.frame(t(tmp))
#     colnames(tmp1)=thesamples
#   }
#   rownames(tmp1)=paste0("pc",1:nrow(tmp1))
#   return(tmp1)
# }
# eigenstratmatrix=readeigenstrat()
colnames(eigenstratmatrix)=update_bcasamplenames(colnames(eigenstratmatrix))

pvalue_assoc_selectedsnps=function(selectedsnps,bcagenotype) #outcome ~SNP
{
  comsamples=intersect(sampletable$localid,colnames(bcagenotype))
  idx=match(comsamples,colnames(bcagenotype))
  bcagenotype_=bcagenotype[,idx]
  idx=match(comsamples,sampletable$localid)
  sampletable_=sampletable[idx,]
  # idx=match(comsamples,colnames(eigenstratmatrix))
  # sampletable_=cbind.data.frame(sampletable_,t(eigenstratmatrix[1:4,idx]))
 
  res=data.frame(snp=selectedsnps,BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA)
  selectedsnps_alt=sapply(selectedsnps,function(x){
    tmp=unlist(strsplit(x,"_"))
    paste0(tmp[c(1,3,2)],collapse = "_")
  })
  idx3=which(selectedsnps_alt %in% rownames(bcagenotype_))
  if (length(idx3)>0) selectedsnps[idx3]=selectedsnps_alt[idx3]
  idx1=which(sampletable_$phenoBE_bc==2) #case
  idx2=which(sampletable_$phenoBE_bc==1)
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  for (i in 1:length(selectedsnps))
  {
    x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=sampletable_[c(idx1,idx2),],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      res$BE_p[i]=summary(fm)$coefficients[2,4]
    }
  }
  
  idx1=which(sampletable_$phenoEA_bc==2) #case
  idx2=which(sampletable_$phenoEA_bc==1)
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  for (i in 1:length(selectedsnps))
  {
    x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=sampletable_[c(idx1,idx2),],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      res$EA_p[i]=summary(fm)$coefficients[2,4]
    }
  }
  
  idx1=which(sampletable_$phenoEA_bc==2) # EAcase
  idx2=which(sampletable_$phenoBE_bc==2) #BEcase
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  for (i in 1:length(selectedsnps))
  {
    x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=sampletable_[c(idx1,idx2),],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      res$BEA_p[i]=summary(fm)$coefficients[2,4]
    }
  }
  
  idx1=which(sampletable_$phenoBE_bc==2 | sampletable_$phenoEA_bc==2) #BE or EAcase
  idx2=which(sampletable_$phenoBE_bc==1 |sampletable_$phenoEA_bc==1) #Control
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  for (i in 1:length(selectedsnps))
  {
    x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=sampletable_[c(idx1,idx2),],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      res$BEEA_p[i]=summary(fm)$coefficients[2,4]
    }
  }
  
  return(res)
}

#include GWAS BCA data
load("../result/Dong23SNPs.RData")
tmp=sapply(rownames(dat),function(x){
  tmp1=unlist(strsplit(x,"_"))
  paste0(tmp1[2:length(tmp1)],collapse = "_")
})
rownames(dat)=tmp
idx=match(allsamples,rownames(dat))
bca_gwasdat=dat[idx,1:23]

#Association snps adjust for GWAS
pvalue_assoc_selectedsnps_adjgwas=function(selectedsnps,bcagenotype,selectedsnps_gwasdist) #outcome ~SNP
{
  res=NA
  if (class(selectedsnps_gwasdist)=="data.frame")
  {
    comsamples=intersect(sampletable$localid,colnames(bcagenotype))
    idx=match(comsamples,colnames(bcagenotype))
    bcagenotype_=bcagenotype[,idx]
    idx=match(comsamples,sampletable$localid)
    sampletable_=sampletable[idx,]
    idx=match(comsamples,rownames(bca_gwasdat))
    bca_gwasdat_=bca_gwasdat[idx,]
    idxgwas=which(colnames(bca_gwasdat) %in% selectedsnps_gwasdist$gwas)
    bca_gwasdat_=bca_gwasdat_[,idxgwas,drop=F]
    res=data.frame(snp=selectedsnps,BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA)
    selectedsnps_alt=sapply(selectedsnps,function(x){
      tmp=unlist(strsplit(x,"_"))
      paste0(tmp[c(1,3,2)],collapse = "_")
    })
    idx3=which(selectedsnps_alt %in% rownames(bcagenotype_))
    if (length(idx3)>0) selectedsnps[idx3]=selectedsnps_alt[idx3]
    idx1=which(sampletable_$phenoBE_bc==2) #case
    idx2=which(sampletable_$phenoBE_bc==1)
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BE_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoEA_bc==2) #case
    idx2=which(sampletable_$phenoEA_bc==1)
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$EA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoEA_bc==2) # EAcase
    idx2=which(sampletable_$phenoBE_bc==2) #BEcase
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BEA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoBE_bc==2 | sampletable_$phenoEA_bc==2) #BE or EAcase
    idx2=which(sampletable_$phenoBE_bc==1 |sampletable_$phenoEA_bc==1) #Control
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BEEA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
  }
  
  return(res)
}

pvalue_assoc_selectedsnps_adjgwas=function(selectedsnps,bcagenotype,selectedsnps_gwasdist) #outcome ~SNP
{
  res=NA
  if (class(selectedsnps_gwasdist)=="data.frame")
  {
    comsamples=intersect(sampletable$localid,colnames(bcagenotype))
    idx=match(comsamples,colnames(bcagenotype))
    bcagenotype_=bcagenotype[,idx]
    idx=match(comsamples,sampletable$localid)
    sampletable_=sampletable[idx,]
    idx=match(comsamples,rownames(bca_gwasdat))
    bca_gwasdat_=bca_gwasdat[idx,]
    idxgwas=which(colnames(bca_gwasdat) %in% selectedsnps_gwasdist$gwas)
    bca_gwasdat_=bca_gwasdat_[,idxgwas,drop=F]
    res=data.frame(snp=selectedsnps,BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA)
    selectedsnps_alt=sapply(selectedsnps,function(x){
      tmp=unlist(strsplit(x,"_"))
      paste0(tmp[c(1,3,2)],collapse = "_")
    })
    idx3=which(selectedsnps_alt %in% rownames(bcagenotype_))
    if (length(idx3)>0) selectedsnps[idx3]=selectedsnps_alt[idx3]
    idx1=which(sampletable_$phenoBE_bc==2) #case
    idx2=which(sampletable_$phenoBE_bc==1)
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BE_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoEA_bc==2) #case
    idx2=which(sampletable_$phenoEA_bc==1)
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$EA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoEA_bc==2) # EAcase
    idx2=which(sampletable_$phenoBE_bc==2) #BEcase
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BEA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
    
    idx1=which(sampletable_$phenoBE_bc==2 | sampletable_$phenoEA_bc==2) #BE or EAcase
    idx2=which(sampletable_$phenoBE_bc==1 |sampletable_$phenoEA_bc==1) #Control
    y=c(rep(1,length(idx1)),rep(0,length(idx2)))
    for (i in 1:length(selectedsnps))
    {
      x=as.numeric(bcagenotype_[which(rownames(bcagenotype_)==selectedsnps[i]),c(idx1,idx2)])
      dat1=cbind.data.frame(y=y,x=x,sampletable_[c(idx1,idx2),c(6,7,23,24,25,26)],bca_gwasdat_[c(idx1,idx2),])
      fm=glm(y~.,data=dat1,family=binomial)
      #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
      tmp=summary(fm)$coefficients
      if ("x" %in% rownames(tmp))
      {
        res$BEEA_p[i]=summary(fm)$coefficients[2,4]
      }
    }
  }
  
  return(res)
}
#correlation matrix with GWAS
corr_selectedsnps_gwas=function(selectedsnps,bcagenotype,selectedsnps_gwasdist) #outcome ~SNP
{
  res=NA
  if (class(selectedsnps_gwasdist)=="data.frame")
  {
    comsamples=intersect(sampletable$localid,colnames(bcagenotype))
    idx=match(comsamples,colnames(bcagenotype))
    bcagenotype_=bcagenotype[,idx]
    idx=match(comsamples,rownames(bca_gwasdat))
    bca_gwasdat_=bca_gwasdat[idx,]
    idxgwas=which(colnames(bca_gwasdat) %in% selectedsnps_gwasdist$gwas)
    bca_gwasdat_=bca_gwasdat_[,idxgwas,drop=F]
    idx=match(selectedsnps,rownames(bcagenotype_))
    if (sum(is.na(idx)>0))
    {
      tmp=sapply(selectedsnps,function(x){
        tmp1=unlist(strsplit(x,"_"))
        paste0(tmp1[c(1,3,2)],collapse = "_")
      })
      idx1=tmp %in% rownames(bcagenotype_)
      selectedsnps[idx1]=tmp[idx1]
      idx=match(selectedsnps,rownames(bcagenotype_))
    }
    res=cor(data.matrix(cbind(t(bcagenotype_[idx,]),bca_gwasdat_)))
  }
  return(res)
}  

compute_p_arow_adjgwas=function(i,genemodel=genemodel1,bcagenotype=bcagenotype1,X,selectedsnps_gwasdist=selectedsnps1_gwasdist)
{
  res=NA
  if (class(selectedsnps_gwasdist)=="data.frame")
  {
    comsamples=intersect(sampletable$localid,colnames(bcagenotype))
    idx=match(comsamples,rownames(bca_gwasdat))
    bca_gwasdat_=bca_gwasdat[idx,]
    idxgwas=which(colnames(bca_gwasdat) %in% selectedsnps_gwasdist$gwas)
    bca_gwasdat_=bca_gwasdat_[,idxgwas,drop=F]
    res=compute_p_arow(i,genemodel,bcagenotype,X=cbind(X,bca_gwasdat_))
  }
  return(res)
}

#most genes came from outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4"
check_genes_in2data=function(genes=unique(c(names(fdr_fwer_res$BEEA_fdr),names(fdr_fwer_res$BE_fdr),"MEF2B")),
                             outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",
                             outfolder2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_PC4")
{
  genes=updategenenames(genes)
  load(paste0(outfolder1,"/preidiction_michigan_model.RData"))
  genemodel1=res_min
  rownames(genemodel1)=updategenenames(rownames(genemodel1))
  load(paste0(outfolder2,"/preidiction_michigan_model.RData"))
  genemodel2=res_min
  rownames(genemodel2)=updategenenames(rownames(genemodel2))
  load("../result/GTExjunctiondatafor_prediction.RData")
  phenotype1=phenotype[phenotypepos$geneid %in% proteingenes,]
  phenotypepos1=phenotypepos[phenotypepos$geneid %in% proteingenes,]
  rownames(phenotype1)=updategenenames(rownames(phenotype1))
  rownames(phenotypepos1)=updategenenames(rownames(phenotypepos1))
  snp1=snp
  snppos1=snppos
  covariate1=covariate
  load(paste0(outfolder1,"/bca_extractgenotype.RData"))
  bcagenotype1=bcagenotype
  colnames(bcagenotype1)=update_bcasamplenames(colnames(bcagenotype1))
  load("../result/TCGAdatafor_prediction_michigan.RData")
  phenotype2=phenotype
  phenotypepos2=phenotypepos
  rownames(phenotype2)=updategenenames(rownames(phenotype2))
  rownames(phenotypepos2)=updategenenames(rownames(phenotypepos2))
  snp2=snp
  snppos2=snppos
  copynumber2=copynumber
  mutation2=mutation
  covariate2=covariate
  load(paste0(outfolder2,"/bca_extractgenotype.RData"))
  bcagenotype2=bcagenotype
  colnames(bcagenotype2)=update_bcasamplenames(colnames(bcagenotype2))
  load("../result/Dong23SNPs.RData")
  gr_gwas=GRanges(seqnames = dong23snp$Chr,ranges=IRanges(start=dong23snp$Position,width = 1))
  allres=vector("list",length(genes))
  names(allres)=genes
  for (i in 1:length(genes))
  {
    #cat(i,"..")
    gene=genes[i]
    print(gene)
    selectedsnps1=selectedsnps2=selectedsnps1_gwasdist=selectedsnps2_gwasdist=selectedsnps1_eqtl1=selectedsnps1_eqtl2=selectedsnps2_eqtl1=selectedsnps2_eqtl2=assoc_selectedsnps1=assoc_selectedsnps2=assoc_selectedsnps1_adjgwas=assoc_selectedsnps2_adjgwas=corr_selectedsnps1_gwas=corr_selectedsnps2_gwas=skat_p1_adjgwas=skat_p2_adjgwas=NA
    idx1=which(rownames(genemodel1)==gene)
    if (length(idx1)>0)
    {
      if (genemodel1$glmflag[idx1]==T)
      {
        selectedsnps1=unlist(strsplit(genemodel1$selectedsnps[idx1],"|",fixed=T))
        selectedsnps1_pos=sapply(selectedsnps1,function(x){
          tmp=unlist(strsplit(x,"_"))
          as.integer(unlist(strsplit(tmp,":"))[2])
        })
        chr=unlist(strsplit(unlist(strsplit(selectedsnps1[1],"_")),":"))[1]
        gr_selectedsnps1=GRanges(seqnames = chr,ranges = IRanges(start=selectedsnps1_pos,width =1))
        tmp=distanceToNearest(gr_selectedsnps1,gr_gwas)
        selectedsnps1_gwasdist=data.frame(snp=selectedsnps1,gwas=NA,distance=NA,stringsAsFactors = F)
        idx3=match(tmp@from,1:length(selectedsnps1))
        selectedsnps1_gwasdist$gwas[idx3]=dong23snp$SNP[tmp@to]
        selectedsnps1_gwasdist$distance[idx3]=tmp@elementMetadata$distance
        skat_p1_adjgwas=compute_p_arow_adjgwas(i=which(rownames(genemodel1)[genemodel1$glmflag==1]==gene),genemodel = genemodel1,bcagenotype = bcagenotype1,X,selectedsnps_gwasdist=selectedsnps1_gwasdist)
        #print(skat_p1_adjgwas)
        #print(skat_min1[which(rownames(skat_min1)==gene),])
        selectedsnps1_eqtl1=computeeQTL1(selectedsnps = selectedsnps1,snp=snp1,phenotype=phenotype1,covariate = covariate1,gene=gene) #in GTEx data
        selectedsnps1_eqtl2=computeeQTL2(selectedsnps = selectedsnps1,snp=snp2,phenotype=phenotype2,copynumber = copynumber2,mutation=mutation2,covariate = covariate2, gene=gene) #in TCGA data
        assoc_selectedsnps1=pvalue_assoc_selectedsnps(selectedsnps = selectedsnps1,bcagenotype = bcagenotype1)
        assoc_selectedsnps1_adjgwas=pvalue_assoc_selectedsnps_adjgwas(selectedsnps = selectedsnps1,bcagenotype = bcagenotype1,selectedsnps_gwasdist=selectedsnps1_gwasdist)
        corr_selectedsnps1_gwas=corr_selectedsnps_gwas(selectedsnps = selectedsnps1,bcagenotype = bcagenotype1,selectedsnps_gwasdist=selectedsnps1_gwasdist)

      }
    }
    
    idx2=which(rownames(genemodel2)==gene)
    if (length(idx2)>0)
    {
      if (genemodel2$glmflag[idx2]==T)
      {
        selectedsnps2=unlist(strsplit(genemodel2$selectedsnps[idx2],"|",fixed=T))
        selectedsnps2_pos=sapply(selectedsnps2,function(x){
          tmp=unlist(strsplit(x,"_"))
          as.integer(unlist(strsplit(tmp,":"))[2])
        })
        chr=unlist(strsplit(unlist(strsplit(selectedsnps2[1],"_")),":"))[1]
        gr_selectedsnps2=GRanges(seqnames = chr,ranges = IRanges(start=selectedsnps2_pos,width =1))
        tmp=distanceToNearest(gr_selectedsnps2,gr_gwas)
        selectedsnps2_gwasdist=data.frame(snp=selectedsnps2,gwas=NA,distance=NA,stringsAsFactors = F)
        idx3=match(tmp@from,1:length(selectedsnps2))
        selectedsnps2_gwasdist$gwas[idx3]=dong23snp$SNP[tmp@to]
        selectedsnps2_gwasdist$distance[idx3]=tmp@elementMetadata$distance
        if (all(selectedsnps2_gwasdist$distance!=0))
        {
          skat_p2_adjgwas=compute_p_arow_adjgwas(i=which(rownames(genemodel2)[genemodel2$glmflag==1]==gene),genemodel = genemodel2,bcagenotype = bcagenotype2,X,selectedsnps_gwasdist=selectedsnps2_gwasdist)
          #print(skat_p2_adjgwas)
          #print(skat_min2[which(rownames(skat_min2)==gene),])
        }
        
        selectedsnps2_eqtl1=computeeQTL1(selectedsnps = selectedsnps2,snp=snp1,phenotype=phenotype1,covariate = covariate1,gene=gene) #in GTEx data
        selectedsnps2_eqtl2=computeeQTL2(selectedsnps = selectedsnps2,snp=snp2,phenotype=phenotype2,copynumber = copynumber2,mutation=mutation2,covariate = covariate2, gene=gene) #in TCGA data
        assoc_selectedsnps2=pvalue_assoc_selectedsnps(selectedsnps = selectedsnps2,bcagenotype = bcagenotype2)
        assoc_selectedsnps2_adjgwas=pvalue_assoc_selectedsnps_adjgwas(selectedsnps = selectedsnps2,bcagenotype = bcagenotype2,selectedsnps_gwasdist=selectedsnps2_gwasdist)
        corr_selectedsnps2_gwas=corr_selectedsnps_gwas(selectedsnps = selectedsnps2,bcagenotype = bcagenotype2,selectedsnps_gwasdist=selectedsnps2_gwasdist)

      }
    }
    
    selectsnps_com=intersect(selectedsnps1,selectedsnps2)
    if (length(selectsnps_com)==0) selectsnps_com=NA
    res=list(gene=gene,selectedsnps1=selectedsnps1,selectedsnps2=selectedsnps2,selectedsnps1_gwasdist=selectedsnps1_gwasdist,selectedsnps2_gwasdist=selectedsnps2_gwasdist,
             selectedsnps1_eqtl1=selectedsnps1_eqtl1,selectedsnps1_eqtl2=selectedsnps1_eqtl2,selectedsnps2_eqtl1=selectedsnps2_eqtl1,selectedsnps2_eqtl2=selectedsnps2_eqtl2,
             assoc_selectedsnps1=assoc_selectedsnps1,assoc_selectedsnps2=assoc_selectedsnps2,assoc_selectedsnps1_adjgwas=assoc_selectedsnps1_adjgwas,assoc_selectedsnps2_adjgwas=assoc_selectedsnps2_adjgwas,
             corr_selectedsnps1_gwas=corr_selectedsnps1_gwas,corr_selectedsnps2_gwas=corr_selectedsnps2_gwas,selectsnps_com=selectsnps_com,skat_p1_adjgwas=skat_p1_adjgwas,skat_p2_adjgwas=skat_p2_adjgwas)
   allres[[i]]=res
  }
  return(allres)
}

check_genes_in2data_res=check_genes_in2data()
save(check_genes_in2data_res,file="../result/check_genes_in2data_res.RData")



#some results
tmp=NULL
for (i in 1:nrow(genemodel2))
{
  if (genemodel2$glmflag[i]==1)
  {
    tmp1=unlist(strsplit(genemodel2$selectedsnps[i],"|",fixed=T))
    tmp=unique(c(tmp,tmp1))
  }
}
length(unique(tmp))

genemodel1_code=genemodel1[rownames(genemodel1) %in% proteingenes,]
tmp=NULL
for (i in 1:nrow(genemodel1_code))
{
  if (genemodel1_code$glmflag[i]==1)
  {
    tmp1=unlist(strsplit(genemodel1_code$selectedsnps[i],"|",fixed=T))
    tmp=unique(c(tmp,tmp1))
  }
}
length(unique(tmp))

qqplot(skat_min2$BE_p)
qqplot(skat_min2$EA_p)
qqplot(skat_min2$BEA_p)
qqplot(skat_min2$BEEA_p)

genes=rownames(genemodel2)[genemodel2$glmflag==1 & genemodel2$r2>0.05]
idx=rownames(skat_min2) %in% genes

qqplot(skat_min2$BE_p[idx])
qqplot(skat_min2$EA_p[idx])
qqplot(skat_min2$BEA_p[idx])
qqplot(skat_min2$BEEA_p[idx])

for (gene in genes)
{
  print(gene)
  idx=which(rownames(skat_min2)==gene)
  if (length(idx)>0)
    print(skat_min2[idx,])
}

for (gene in genes)
{
  print(gene)
  idx=which(rownames(skat_min1)==gene)
  if (length(idx)>0)
    print(skat_min1[idx,])
}


