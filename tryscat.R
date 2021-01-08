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
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
library(SKAT)
# gene="DPP7"
# selectedsnps=unlist(strsplit(res_min$selectedsnps[which(rownames(res_min)==gene)],"|",fixed=T))
# idx=match(selectedsnps,rownames(bcagenotype))
# idxtocorrect=which(is.na(idx))
# if (length(idxtocorrect)>0)
# {
#   for (i in idxtocorrect)
#   {
#     tmp=unlist(strsplit(selectedsnps[i],"_"))
#     selectedsnps[i]=paste0(tmp[c(1,3,2)],collapse = "_")
#   }
#   idx=match(selectedsnps,rownames(bcagenotype))
# }
# 
# Z=t(bcagenotype[idx,])
# if (length(idxtocorrect)>0)
# {
#   Z[idxtocorrect,]=2-Z[idxtocorrect,]
# }


if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  sampletable$site[sampletable$site=="NA"]=NA
  geneexpsamplenames=strsplit(colnames(predict_1se)[3:ncol(predict_1se)],"_")
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}

readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
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
# allsamples=intersect(rownames(Z),sampletable$localid)
# idx=match(allsamples,rownames(Z))
# Z=Z[idx,]
# maf=rep(NA,ncol(Z))
# for (i in 1:ncol(Z))
# {
#   maf[i]=(sum(Z[,i]==1)+sum(Z[,i]==2)*2)/nrow(Z)/2
# }
# idx=match(allsamples,sampletable$localid)
# sampletable=sampletable[idx,]
# idx=match(allsamples,colnames(eigenstratmatrix))
# eigenstratmatrix=eigenstratmatrix[,idx]
# sampletable=cbind(sampletable,t(eigenstratmatrix[1:3,]))
# X=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3")]
removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  # idxconst=rep(F,nrow(dat))
  # for (i in 1:nrow(dat))
  # {
  #   if (i %% 10000==0) cat(i,"..")
  #   if (var(unlist(dat[i,]))==0 | is.na(var(unlist(dat[i,])))) idxconst[i]=T
  # }
  dat=dat[!idxconst,]
}
skat_p=function(Z,X,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  Z1=Z[c(idx1,idx2),,drop=F]
  X1=as.matrix(X[c(idx1,idx2),])
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

# #BE vs CO
# idx1=which(sampletable$phenoBE_bc==2) #case
# idx2=which(sampletable$phenoBE_bc==1)
# skat_p(Z,X,idx1,idx2)
# #EA vs CO
# idx1=which(sampletable$phenoEA_bc==2) #case
# idx2=which(sampletable$phenoEA_bc==1)
# skat_p(Z,X,idx1,idx2)
# #BE vs EA
# idx1=which(sampletable$phenoBE_bc==2)
# idx2=which(sampletable$phenoEA_bc==2)
# skat_p(Z,X,idx1,idx2)
# #BE EA vs CO
# idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
# idx2=which(sampletable$phenoBE_bc==1)
# skat_p(Z,X,idx1,idx2)

#genomewide:
sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
sampletable$site[sampletable$site=="NA"]=NA
eigenstratmatrix=readeigenstrat()
allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
idx=match(allsamples,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
sampletable=cbind(sampletable,t(eigenstratmatrix[1:3,]))
X=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3")]
compute_p=function(genemodel=res_1se)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  res=data.frame(gene=genes,BE_p=rep(NA,length(genes)),EA_p=rep(NA,length(genes)),BEA_p=rep(NA,length(genes)),BEEA_p=rep(NA,length(genes)),stringsAsFactors = F)
  for (i in 1:length(genes))
  {
    if (i %% 10==0) cat(i,'.')
    #cat(i,'.')
    gene=genes[i]
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
    idx=match(allsamples,rownames(Z))
    Z=Z[idx,,drop=F]
    #print(rankMatrix(Z)[[1]])
    idx1=which(sampletable$phenoBE_bc==2) #case
    idx2=which(sampletable$phenoBE_bc==1)
    res$BE_p[i]=skat_p(Z,X,idx1,idx2)
    idx1=which(sampletable$phenoEA_bc==2) #case
    idx2=which(sampletable$phenoEA_bc==1)
    res$EA_p[i]=skat_p(Z,X,idx1,idx2)
    idx1=which(sampletable$phenoBE_bc==2)
    idx2=which(sampletable$phenoEA_bc==2)
    res$BEA_p[i]=skat_p(Z,X,idx1,idx2)
    idx1=which(sampletable$phenoBE_bc==2 | sampletable$phenoEA_bc==2) #case
    idx2=which(sampletable$phenoBE_bc==1)
    res$BEEA_p[i]=skat_p(Z,X,idx1,idx2)
  }
  return(res)
}
skat_1se=compute_p()
skat_min=compute_p(genemodel = res_min)
save(skat_1se,skat_min,file=paste0(outfolder,"/skat_res.RData"))
idx=skat_min$gene %in% proteingenes
skat_min_code=skat_min[skat_min$gene %in% proteingenes,]
numsnps=res_min$numselectedsnp[match(skat_min_code$gene,rownames(res_min))]
idx1=skat_min$BEEA_p[idx]<1e-2
t.test(numsnps[idx1],numsnps[!idx1])
quantile(numsnps[!idx1],na.rm=T)
quantile(numsnps[idx1],na.rm=T)

numsnps=res_min$numselectedsnp[match(skat_min$gene,rownames(res_min))]
idx1=skat_min$BEEA_p<5e-2
t.test(numsnps[idx1],numsnps[!idx1])
quantile(numsnps[!idx1],na.rm=T)
quantile(numsnps[idx1],na.rm=T)
