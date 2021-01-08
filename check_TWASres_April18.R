#!/usr/bin/env Rscript
qqplot=function(pvalue=NULL,fwer=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (-log base 10)",
       ylab="Observed p-value (-log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
  pvalue_order=pvalue[order(pvalue)]
  if (is.null(fwer))
  {
    fwer=p.adjust(pvalue_order,method="bonferroni")
    idx=which(fwer<0.05)
    if (length(idx)>0)
    {
      points(-log((1:n)/n,base=10)[idx],-log(pvalue[order(pvalue)],base=10)[idx],pch=16,col="red")
      print(paste0("#FWER<0.05:",length(idx)))
      legend("topleft",legend = "FWER<0.05",pch=16,col="red")
    }
  }else
  {
    idx1=sum(fwer<0.05)
    if (length(idx1)>0)
    {
      idx=1:idx1
      points(-log((1:n)/n,base=10)[idx],-log(pvalue[order(pvalue)],base=10)[idx],pch=16,col="red")
      print(paste0("#FWER<0.05:",length(idx)))
      legend("topleft",legend = "FWER<0.05",pch=16,col="red")
    }
  }
}
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


load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_April18" #models stored for TCGA, came from prediction_michigan_models7_TCGA.R
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18" #models stored for GTEx, came from prediction_michigan_models6_GTEx.R
check_skatres=function(dat=skat_min1)
{
  colnames(dat)=c("BE_p","EA_p","BEA_p","BEEA_p")
  dat_code=dat[rownames(dat) %in% proteingenes,]
  tmp=compute_fwer_fdr(dat=dat_code)
  for(i in 1:length(tmp))
  {
    if (!is.null(tmp[[i]]))
    {
      print(paste0(names(tmp)[i],":"))
      print(tmp[[i]])
    }
  }
  par(mfrow=c(2,2))
  qqplot(dat_code$BE_p,main="BE")
  qqplot(dat_code$EA_p,main="EA")
  qqplot(dat_code$BEA_p,main="BE vs EA")
  qqplot(dat_code$BEEA_p,main="BEEA")
}

check_model=function(mod=res_min,r2cut=0)
{
  mod=mod[which(mod$r2>=r2cut),]
  print(nrow(mod))
  print(quantile(mod$r2))
  print(quantile(mod$numselectedsnp))
}

#GTEx------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mar=c(5,5.5,2,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
tmp=compute_fwer_fdr(dat=assoc_min_code)
par(mfrow=c(2,2))
qqplot(assoc_min_code$BE_p,main="BE")
qqplot(assoc_min_code$EA_p,main="EA")
qqplot(assoc_min_code$BEA_p,main="BEA") #EA vs BE
qqplot(assoc_min_code$BEEA_p,main="BEEA")
#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)
check_skatres(dat=skat_min3)
check_skatres(dat=skat_min4)
genes=c("HSP90AA1","KXD1","ISYNA1","UBAC1")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
genes=c("HSP90AA1","KXD1","UBAC1","FBP2")
for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",skat_min2$BE_p[idx]," BEEA-p ",skat_min2$BEEA_p[idx]))
}

genes=c("CERS1","ARMC6")
for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  #print(paste0(gene,":",skat_min2$BE_p[idx]))
  print(paste0(gene,":",skat_min2$BEEA_p[idx]))
}

#TCGA--------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_April18"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
hist(res_min$r2[res_min$glmflag==1],main="",xlab="R-squared",col="blue")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData")
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))

tmp=compute_fwer_fdr(dat=assoc_min_code)
par(mfrow=c(2,2))
qqplot(assoc_min_code$BE_p,main="BE")
qqplot(assoc_min_code$EA_p,main="EA")
qqplot(assoc_min_code$BEA_p,main="BEA") #EA vs BE
qqplot(assoc_min_code$BEEA_p,main="BEEA")
#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)
check_skatres(dat=skat_min3)
check_skatres(dat=skat_min4)
genes=c("MEF2B","DDX49","COPE","FBP2","LSM4","GDF15","PDE4C","LRRC25","SNED1")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
genes=c("FBP2")
for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",skat_min2$BE_p[idx]," BEEA-p ",skat_min2$BEEA_p[idx]))
}


#GTEx mucosa
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_April18"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
Scheck_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
tmp=compute_fwer_fdr(dat=assoc_min_code)
par(mar=c(5,5.5,3.5,2))
par(mfrow=c(2,2))
qqplot(assoc_min_code$BE_p,main="BE")
qqplot(assoc_min_code$EA_p,main="EA")
qqplot(assoc_min_code$BEA_p,main="BEA") #EA vs BE
qqplot(assoc_min_code$BEEA_p,main="BEEA")
which(rownames(phenotypepos)=="CENPQ")
phenotypepos[12130,]
#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)
check_skatres(dat=skat_min3)
check_skatres(dat=skat_min4)
which(rownames(phenotypepos)=="HEATR5B")
phenotypepos[3865,]


