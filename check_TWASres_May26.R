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


load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_May26" #models stored for GTEx, came from prediction_michigan_models6_GTExV8.R
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

scinumber=function(num)
{
  formatC(num,format="e",digits = 2)
}
#GTEx junction------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_May26"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")
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
rownames(assoc_min_code)[which.min(assoc_min_code$BEEA_p)] 

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
genes=c("HSP90AA1","ELL","ZPLD1","COMP")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
# [1] "HSP90AA1 14:102080738-102139502"
# [1] "ELL 19:18442663-18522127"
# [1] "ZPLD1 3:102099244-102479841"
# [1] "COMP 19:18782773-18791314"

for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",scinumber(skat_min2$BE_p[idx])))
  print(paste0(gene,": EA-p ",scinumber(skat_min2$EA_p[idx])))
  print(paste0(gene," BEEA-p ",scinumber(skat_min2$BEEA_p[idx])))
}
# [1] "HSP90AA1: BE-p 7.14e-07"
# [1] "HSP90AA1: EA-p 1.58e-03"
# [1] "HSP90AA1 BEEA-p 1.40e-06"
# [1] "ELL: BE-p 8.50e-06"
# [1] "ELL: EA-p 6.38e-05"
# [1] "ELL BEEA-p 4.86e-07"
# [1] "ZPLD1: BE-p 2.17e-02"
# [1] "ZPLD1: EA-p 3.66e-06"
# [1] "ZPLD1 BEEA-p 9.58e-05"
# [1] "COMP: BE-p 1.13e-04"
# [1] "COMP: EA-p 6.14e-05"
# [1] "COMP BEEA-p 5.31e-06"

#check v7 genes
genes=c("KXD1","UBAC1","ISYNA1")
for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",scinumber(skat_min2$BE_p[idx])))
  print(paste0(gene,": EA-p ",scinumber(skat_min2$EA_p[idx])))
  print(paste0(gene," BEEA-p ",scinumber(skat_min2$BEEA_p[idx])))
}

#GTEx mucosa------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_May26"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadatafor_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
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


#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)
# [1] "BE_fdr:"
# ARMC6     TMEM161A 
# 7.732076e-05 1.765087e-03 
# [1] "BE_fwer:"
# ARMC6     TMEM161A 
# 7.732076e-05 3.530175e-03 
# [1] "EA_fdr:"
# TMEM161A       ARMC6 
# 0.003450197 0.008110488 
# [1] "EA_fwer:"
# TMEM161A       ARMC6 
# 0.003450197 0.016220975 
# [1] "BEEA_fdr:"
# TMEM161A        ARMC6     SLC25A42 
# 2.952648e-07 8.806548e-06 1.521169e-02 
# [1] "BEEA_fwer:"
# TMEM161A        ARMC6     SLC25A42 
# 2.952648e-07 1.761310e-05 4.563506e-02 
check_skatres(dat=skat_min3)
check_skatres(dat=skat_min4)
genes=c("ARMC6","SLC25A42","TMEM161A")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
# [1] "ARMC6 19:19033575-19060311"
# [1] "SLC25A42 19:19063999-19112888"
# [1] "TMEM161A 19:19119169-19138513"
genes=c("ARMC6","SLC25A42","TMEM161A")
for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",formatC(skat_min2$BE_p[idx],format="e",digits = 2),
               " EA-p ",formatC(skat_min2$EA_p[idx],format="e",digits = 2),
               " BEEA-p ",formatC(skat_min2$BEEA_p[idx],format="e",digits = 2)))
}
# [1] "ARMC6: BE-p 6.50e-09 EA-p 1.36e-06 BEEA-p 1.48e-09"
# [1] "SLC25A42: BE-p 4.39e-05 EA-p 2.28e-04 BEEA-p 3.84e-06"
# [1] "TMEM161A: BE-p 2.97e-07 EA-p 2.90e-07 BEEA-p 2.48e-11"

#GTEx blood------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_May26"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8blooddatafor_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
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
gene="RAB34"
idx=which(rownames(phenotypepos)==gene)
print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))


#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)

genes=c("PGPEP1","ISYNA1")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
# [1] "PGPEP1 19:18340587-18369950"
# [1] "ISYNA1 19:18434563-18438301"

for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",formatC(skat_min2$BE_p[idx],format="e",digits = 2),
               " EA-p ",formatC(skat_min2$EA_p[idx],format="e",digits = 2),
               " BEEA-p ",formatC(skat_min2$BEEA_p[idx],format="e",digits = 2)))
}

#GTEx somatch------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_May26"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdatafor_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
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
gene="RAB34"
idx=which(rownames(phenotypepos)==gene)
print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))


#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min1)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,25)
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
colnames(skat_min3)=c("BE_p","EA_p","BEA_p","BEEA_p") #method="SKATO",weights.beta=c(1,1)
colnames(skat_min4)=c("BE_p","EA_p","BEA_p","BEEA_p") #SKAT_CommonRare
check_skatres(dat=skat_min1)
check_skatres(dat=skat_min2)

genes=c("PGPEP1","ISYNA1")
for (gene in genes)
{
  idx=which(rownames(phenotypepos)==gene)
  print(paste0(gene," ",phenotypepos$chr[idx],":",phenotypepos$s1[idx],"-",phenotypepos$s2[idx]))
}
# [1] "PGPEP1 19:18340587-18369950"
# [1] "ISYNA1 19:18434563-18438301"

for (gene in genes)
{
  idx=which(rownames(skat_min2)==gene)
  print(paste0(gene,": BE-p ",formatC(skat_min2$BE_p[idx],format="e",digits = 2),
               " EA-p ",formatC(skat_min2$EA_p[idx],format="e",digits = 2),
               " BEEA-p ",formatC(skat_min2$BEEA_p[idx],format="e",digits = 2)))
}