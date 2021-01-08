#!/usr/bin/env Rscript
library(data.table)

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

#sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
fam=read.table("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015.fam",stringsAsFactors = F)
# idx=which(sampletable[,2] %in% sampletable[which(duplicated(sampletable[,2])),2])
# View(sampletable[idx,])

idx=sampletable$localid %in% fam$V2
table(idx)

table(sampletable$phenoBE_bca[idx])
# -9    1    2 
# 2740 3204 3295 
table(sampletable$phenoEA_bca[idx])
# -9    1    2 
# 3517 3207 2515  

table(sampletable$phenoEA_bca[idx],sampletable$phenoBE_bca[idx])
table(sampletable$phenoEA_bca[idx],sampletable$phenoBE_bca[idx])

table(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1 )
# FALSE  TRUE 
# 6032  3207 

BEsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$phenoBE_bca==2)]
EAsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$phenoEA_bca==2)]
COsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & (sampletable$phenoBE_bca==1 | sampletable$phenoEA_bca==1))]

filesforGWAS=function(casesamples=BEsamples,controlsamples=COsamples,prefix="BE_CO")
{
  idx1=match(casesamples,fam$V2)
  tmp1=data.frame(FID=fam$V1[idx1],IID=fam$V2[idx1],affected=2,stringsAsFactors = F)
  idx2=match(controlsamples,fam$V2)
  tmp2=data.frame(FID=fam$V1[idx2],IID=fam$V2[idx2],affected=1,stringsAsFactors = F)
  tmp=rbind(tmp1,tmp2)
  #pcadat=read.table(pcafile,stringsAsFactors = F)
  tmp=tmp[tmp$IID %in% sampletable$localid,]
  idx=match(tmp$IID,sampletable$localid)
  tmp$pc1=sampletable$ev1_bca[idx]
  tmp$pc2=sampletable$ev2_bca[idx]
  tmp$pc3=sampletable$ev3_bca[idx]
  tmp$pc4=sampletable$ev4_bca[idx]
  idx=match(tmp$IID,sampletable$localid)
  tmp$age=sampletable$age[idx]
  tmp$sex=sampletable$sex[idx]
  #for validate_twas, ind ID changed
  tmp1=tmp
  tmp1$IID=paste0(tmp$FID,"_",tmp$IID)
  tmp1$FID=0
  write.table(tmp1[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink1.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp,file=paste0("../result/GWAS/",prefix,"_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  write.table(tmp1,file=paste0("../result/GWAS/",prefix,"_selectedsamples_pheno_plink1.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
filesforGWAS()
filesforGWAS(casesamples = EAsamples,prefix="EA_CO")
filesforGWAS(casesamples = c(BEsamples,EAsamples),prefix="BEEA_CO")

udpate_fam=function(famfile="../result/GWAS/BE_CO_19Feb2015.fam",phenofile="../result/GWAS/BE_CO_selectedsamples_pheno_plink.txt")
{
  fam=read.table(famfile,stringsAsFactors = F)
  pheno=read.table(phenofile,header=T,stringsAsFactors = F)
  idx=match(fam$V2,pheno$IID)
  fam$V5=pheno$sex[idx]
  fam$V6=pheno$affected[idx]
  write.table(fam,file=famfile,col.names=F,row.names = F,quote=F)
}
udpate_fam()
udpate_fam(famfile="../result/GWAS/EA_CO_19Feb2015.fam",phenofile="../result/GWAS/EA_CO_selectedsamples_pheno_plink.txt")
udpate_fam(famfile="../result/GWAS/BEEA_CO_19Feb2015.fam",phenofile="../result/GWAS/BEEA_CO_selectedsamples_pheno_plink.txt")
#check fam file and cov file
tmp=read.table("../result/GWAS/BEEA_CO_19Feb2015.fam")
tmp1=read.table("../result/GWAS/BEEA_CO_selectedsamples_pheno_plink.txt",header = T)
all(tmp$V2 %in% tmp1$IID) #T

BEgwas=fread("../result/GWAS/BE_CO_19Feb2015.assoc.logistic")
EAgwas=fread("../result/GWAS/EA_CO_19Feb2015.assoc.logistic")
BEEAgwas=fread("../result/GWAS/BEEA_CO_19Feb2015.assoc.logistic")
par(mfrow=c(1,2))
qqplot(BEgwas$P)
hist(BEgwas$P)
BEgwas[which.min(BEgwas$P),]
qqplot(EAgwas$P)
hist(EAgwas$P)
EAgwas[which.min(EAgwas$P),]
qqplot(BEEAgwas$P)
hist(BEEAgwas$P)
BEEAgwas[which.min(BEEAgwas$P),]
tmp1=fread("../result/GWAS/BEEA_CO_19Feb2015.ped",nrows = 100)
tmp=fread("test1.ped")
sum(tmp$V7!=0 & tmp$V8!=0) #7978
BEEAgwas$NMISS[1] #7951, 7978-7951=27 samples were removed due to missing phenotype 
#log file:Among remaining phenotypes, 5782 are cases and 2184 are controls.  (27 phenotypes are missing.)

create_sumstats=function(gwasfile="../result/GWAS/BEEA_CO_19Feb2015.assoc.logistic",
                         bimfile="../result/GWAS/BEEA_CO_19Feb2015.bim",prefix="BEEA",opt="nohm3")
{
  gwasdat=as.data.frame(fread(gwasfile))
  bimdat=as.data.frame(fread(bimfile))
  idx=which(bimdat$V1 %in% c(1:23))
  gwasdat=gwasdat[idx,]
  bimdat=bimdat[idx,]
  all(gwasdat$BP==bimdat$V4)
  if (opt=="hm3")
  {
    hm3=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/w_hm3.snplist"))
    idx=which(bimdat$V2 %in% hm3$SNP)
    prefix=paste0(prefix,"_hm3")
    gwasdat=gwasdat[idx,]
    bimdat=bimdat[idx,]
  }
  
  res=data.frame(SNP=bimdat$V2,hg19chr=bimdat$V1,bp=bimdat$V4,A1=bimdat$V5,A2=bimdat$V6,
                 or=gwasdat$OR,se=gwasdat$SE,pval=gwasdat$P,Z=gwasdat$STAT,N=gwasdat$NMISS,stringsAsFactors = F)
  write.table(res,file=paste0("../result/GWAS/",prefix,".sumstats"),col.names = T,row.names = F,sep=" ",quote=F)
  
}

create_sumstats(gwasfile="../result/GWAS/BE_CO_19Feb2015.assoc.logistic",
                         bimfile="../result/GWAS/BE_CO_19Feb2015.bim",prefix="BE")
create_sumstats(gwasfile="../result/GWAS/BE_CO_19Feb2015.assoc.logistic",
                bimfile="../result/GWAS/BE_CO_19Feb2015.bim",prefix="BE",opt="hm3")
create_sumstats(gwasfile="../result/GWAS/EA_CO_19Feb2015.assoc.logistic",
                bimfile="../result/GWAS/EA_CO_19Feb2015.bim",prefix="EA")
create_sumstats(gwasfile="../result/GWAS/EA_CO_19Feb2015.assoc.logistic",
                bimfile="../result/GWAS/EA_CO_19Feb2015.bim",prefix="EA",opt="hm3")
compute_liabilityh2=function(K=0.016,P=3287/(3287+2184),h2=0.5342)
{
  #K=pop prevalence
  #P=proportion of cases in study
  #hsq=Heritability estimate (on observed scale)
  #bigT = liability threshold
  #tau = density of gaussian
  zv <- dnorm(qnorm(K))
  h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2
  return(h2_liab)
}

#use 1000g reference
#BE
K=0.016
P=3287/(3287+2184)
h2=0.5342
zv <- dnorm(qnorm(K))
h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2 #0.345
#EA
K=0.0025
P=2495/(2495+2184)
h2=0.4334
zv <- dnorm(qnorm(K))
h2_liab <- h2 * K^2 * ( 1 - K)^2 / P / (1-P) / zv^2 #0.180

#use all snps
compute_liabilityh2(K=0.016,P=3287/(3287+2184),h2=1.4844) #0.958
compute_liabilityh2(K=0.0025,P=2495/(2495+2184),h2=1.4667) #0.608

#hm3
compute_liabilityh2(K=0.016,P=3287/(3287+2184),h2=0.6661) #0.430
compute_liabilityh2(K=0.0025,P=2495/(2495+2184),h2=0.6101) #0.253

#run gwas on selected SNPs based on impuated data, used for gwas_imputed.sh validate_twas.R
filesforGWAS1=function(casesamples=BEsamples,controlsamples=COsamples,pcafile="../result/bca_filtered_10Jan2019.pca.evec",
                       prefix="BE_CO",prefix1="TCGA",
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_",prefix1,"model/"))
{
  idx1=match(casesamples,fam$V2)
  tmp1=data.frame(FID=fam$V1[idx1],IID=fam$V2[idx1],affected=2,stringsAsFactors = F)
  idx2=match(controlsamples,fam$V2)
  tmp2=data.frame(FID=fam$V1[idx2],IID=fam$V2[idx2],affected=1,stringsAsFactors = F)
  tmp=rbind(tmp1,tmp2)
  pcadat=read.table(pcafile,stringsAsFactors = F)
  tmp=tmp[tmp$IID %in% pcadat$V1,]
  idx=match(tmp$IID,pcadat$V1)
  tmp$pc1=pcadat$V2[idx]
  tmp$pc2=pcadat$V3[idx]
  tmp$pc3=pcadat$V4[idx]
  tmp$pc4=pcadat$V5[idx]
  idx=match(tmp$IID,sampletable$localid)
  tmp$age=sampletable$age[idx]
  tmp$sex=sampletable$sex[idx]
  #for validate_twas, ind ID changed
  tmp1=tmp
  tmp1$IID=paste0(tmp$FID,"_",tmp$IID)
  tmp1$FID=0
  write.table(tmp1[,c(1,2)],file=paste0(outfolder,prefix,"_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  #write.table(tmp[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  #write.table(tmp,file=paste0("../result/GWAS/",prefix,"_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  write.table(tmp1,file=paste0(outfolder,prefix,"_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
update_bim=function(prefix="TCGA",prefix1="BE")
{
  bimfile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_",prefix1,"_CO.bim")
  bim=read.table(bimfile,stringsAsFactors = F)
  bim$V2=paste0(bim$V2,"_",bim$V5,"_",bim$V6)
  write.table(bim,file=bimfile,quote=F,row.names = F,col.names = F)
}
prefix="TCGA"
prefix="GTEx"

impsamples=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_",prefix,"model/select.fam"),stringsAsFactors = F)
impsamples$V7=NA
for (i in 1:nrow(impsamples))
{
  tmp=unlist(strsplit(impsamples$V2[i],"_"))
  impsamples$V7[i]=paste0(tmp[2:length(tmp)],collapse = "_")
}
filesforGWAS1(casesamples = intersect(BEsamples,impsamples$V7),controlsamples = intersect(COsamples,impsamples$V7),
              prefix=paste0("imp_",prefix,"_BE"),prefix1=prefix)
# tmp=read.table("../result/GWAS/imp_TCGA_BE_selectedsamples_plink1.txt",stringsAsFactors = F)
# all(tmp$V2 %in% impsamples$V2)
# tmp1=read.table("../result/GWAS/imp_TCGA_BE_CO.fam")
# tmp2=read.table("../result/GWAS/imp_TCGA_BE_selectedsamples_plink.txt",stringsAsFactors = F)
#run imputed_GWAS.R
udpate_fam(famfile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_BE_CO.fam"),
           phenofile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_BE_selectedsamples_pheno_plink.txt"))
update_bim(prefix=prefix)
#run imputed_GWAS.R

filesforGWAS1(casesamples = intersect(EAsamples,impsamples$V7),controlsamples = intersect(COsamples,impsamples$V7),
              prefix=paste0("imp_",prefix,"_EA"),prefix1=prefix)
#run imputed_GWAS.R
udpate_fam(famfile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_EA_CO.fam"),
           phenofile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_EA_selectedsamples_pheno_plink.txt"))
update_bim(prefix=prefix,prefix1="EA")
#run imputed_GWAS.R

BEEAsamples=c(BEsamples,EAsamples)
filesforGWAS1(casesamples = intersect(BEEAsamples,impsamples$V7),controlsamples = intersect(COsamples,impsamples$V7),
              prefix=paste0("imp_",prefix,"_BEEA"),prefix1=prefix)
udpate_fam(famfile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_BEEA_CO.fam"),
           phenofile=paste0("../result/GWAS/imputed_",prefix,"model/imp_",prefix,"_BEEA_selectedsamples_pheno_plink.txt"))
update_bim(prefix=prefix,prefix1="BEEA")
