#use 4pcs
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC2"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/bca_extractgenotype.RData"))
tmp=colnames(bcagenotype)
tmp=strsplit(tmp,"_")
tmp1=sapply(1:length(tmp),function(x){
  tmp1=tmp[[x]]
  paste0(tmp1[2:length(tmp1)],collapse = "_")
})
colnames(bcagenotype)=tmp1
#load(paste0(outfolder,"/bca_predict_geneexp.RData"))
library(SKAT)

if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  sampletable$site[sampletable$site=="NA"]=NA
  # geneexpsamplenames=strsplit(colnames(predict_1se)[3:ncol(predict_1se)],"_")
  # geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  #   tmp=geneexpsamplenames[[x]]
  #   paste0(tmp[2:length(tmp)],collapse = "_")
  # })
}

readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=colnames(bcagenotype),nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  if (opt==1)
  {
    eigsamples=read.table(eigsampfile,stringsAsFactors = F)
    eigsamples=eigsamples$V2
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=thesamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}

eigenstratmatrix=readeigenstrat()

removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  dat=dat[!idxconst,]
}
skat_p=function(Z,X,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  X1=data.matrix(X[c(idx1,idx2),])
  tmp=rowSums(Z1)
  if (sum(tmp!=0) !=1)
  {
    obj.s<-SKAT_Null_Model(y ~ X1,out_type="D")
    #out<-SKAT(Z1, obj.s,weights.beta=c(1,25),r.corr=0) #0.04183234
    #out<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0) #1.207545e-05
    #out=SKATBinary(Z1, obj.s, weights.beta=c(1,25),kernel = "linear.weighted",r.corr=0) #0.04183234
    out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0) #1.207545e-05
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


#genomewide:
sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
sampletable$site[sampletable$site=="NA"]=NA
eigenstratmatrix=readeigenstrat()
allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
idx=match(allsamples,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
sampletable=cbind(sampletable,t(eigenstratmatrix[1:4,]))
sampletable$sex=as.factor(sampletable$sex)
X=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3","pc4")]
rownames(X)=rownames(sampletable)
compute_p_arow=function(i,genemodel=res_min)
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
   # tmp=colnames(bcagenotype)
   # tmp=strsplit(tmp,"_")
   # tmp1=sapply(1:length(tmp),function(x){
   #   tmp1=tmp[[x]]
   #   paste0(tmp1[2:length(tmp1)],collapse = "_")
   # })
   # colnames(bcagenotype)=tmp1
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  res$BE_p=skat_p(Z,X,idx1,idx2)
  idx1=which(sampletable$phenoEA_bc==2) #case
  idx2=which(sampletable$phenoEA_bc==1)
  res$EA_p=skat_p(Z,X,idx1,idx2)
  idx1=which(sampletable$phenoBE_bc==2)
  idx2=which(sampletable$phenoEA_bc==2)
  res$BEA_p=skat_p(Z,X,idx1,idx2)
  idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  res$BEEA_p=skat_p(Z,X,idx1,idx2)
  return(res)
}

#run sequentially
compute_p_allrows=function(genemodel=res_min)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  res=data.frame(BE_p=rep(NA,length(genes)),EA_p=rep(NA,length(genes)),BEA_p=rep(NA,length(genes)),BEEA_p=rep(NA,length(genes)),stringsAsFactors = F)
  print(paste0("total genes: ",length(genes)))
  for (i in 1:length(genes))
  {
    if (i %%100 ==0) cat(i,"..")
    res[i,]=compute_p_arow(i,genemodel)
  }
  rownames(res)=genes
  return(res)
}
skat_min=compute_p_allrows()
skat_1se=compute_p_allrows(genemodel = res_1se)

library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(outfolder)
mpi.remote.exec(load(paste0(outfolder,"/preidiction_michigan_model.RData")))
mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
mpi.bcast.Robj2slave(sampletable)
mpi.bcast.Robj2slave(allsamples)
mpi.bcast.Robj2slave(X)
mpi.bcast.Robj2slave(removeconstrows)
mpi.bcast.Robj2slave(skat_p)
mpi.bcast.Robj2slave(compute_p_arow)
mpi.remote.exec(library(SKAT))

mpi_compute_p=function(genemodel=res_min)
{
  genemodel=genemodel[genemodel$glmflag==1,]
  genes=rownames(genemodel)
  rows=1:nrow(genemodel)
  res=NULL
  n=njobs
  nchunks=ceiling(length(rows)/n)
  print(paste0("number of total:",nchunks))
  for (i in 1:nchunks)
  {
    if (i %% 10==0) cat(i,"..")
    if (i<nchunks)
    {
      seq=rows[((i-1)*n+1):(i*n)]
    }else
    {
      seq=rows[((i-1)*n+1):length(rows)]
    }
    tmp=mpi.parSapply(X=seq,FUN=compute_p_arow,genemodel=genemodel,job.num=njobs)
    if (length(tmp) %% 4 !=0) print(i)
    res1=matrix(unlist(tmp),ncol=4,byrow = T)
    #rownames(res1)=rownames(genemodel)[seq]
    res=rbind(res,res1)
  }
  res=as.data.frame(res)
  rownames(res)=rownames(genemodel)
  return(res)
}
skat_1se=mpi_compute_p()
skat_min=mpi_compute_p(genemodel = res_min)
save(skat_1se,skat_min,file=paste0(outfolder,"/skat_res.RData"))
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
idx=rownames(skat_min) %in% proteingenes
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
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
#rownames(skat_min_code)=skat_min_code$gene #for GTEx
rownames(skat_min)=skat_min$gene
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
# CERS1     SGK223      CRTC1    PLEKHF2 
# 0.01203058 0.01667941 0.01667941 0.02746873 
# 
# $BE_fwer
# CERS1      CRTC1 
# 0.01203058 0.04368962 
# 
# $EA_fdr
# CERS1      OR2Z1       COMP 
# 0.02909225 0.03261448 0.03261448 
# 
# $EA_fwer
# CERS1 
# 0.02909225 
# 
# $BEA_fdr
# NULL
# 
# $BEA_fwer
# NULL
# 
# $BEEA_fdr
# CERS1        CRTC1       SGK223       KLHL26         COMP      PLEKHF2        OR2Z1 
# 0.0009422403 0.0070993153 0.0086797107 0.0086797107 0.0086797107 0.0124905720 0.0131749355 
# DDX49      C8orf49        ARMC6      SPAG11A       ISYNA1        GATA4      PPP1R3B 
# 0.0158022342 0.0192172758 0.0301337065 0.0343513080 0.0351132577 0.0392682524 0.0479425402 
# CTSB 
# 0.0479669198 
# 
# $BEEA_fwer
# CERS1        CRTC1       SGK223         COMP       KLHL26 
# 0.0009422403 0.0141986305 0.0357804014 0.0372873475 0.0433985537 
fdr_fwer_res=compute_fwer_fdr(dat=skat_min)
# $BE_fdr
# DDX49 
# 0.01039977 
# 
# $BE_fwer
# DDX49 
# 0.01039977 
# 
# $EA_fdr
# DDX49 
# 0.04825351 
# 
# $EA_fwer
# DDX49 
# 0.04825351 
# 
# $BEA_fdr
# TIPARP 
# 0.01064432 
# 
# $BEA_fwer
# TIPARP 
# 0.01064432 
# 
# $BEEA_fdr
# DDX49 
# 0.00125625 
# 
# $BEEA_fwer
# DDX49 
# 0.00125625 
