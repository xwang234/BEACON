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


library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  if (length(idx)>0)
    sampletable[idx,i]=NA
}

#BCA covariate table, principle components
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
  tmp1=cbind(eigsamples,tmp)
  return(tmp1)
}

covariatetable=readeigenstrat()
rownames(covariatetable)=covariatetable$V2=gsub("SEP","",rownames(covariatetable))
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
fam=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar/select.fam")
idx=match(fam$V2,covariatetable$V2)
covariatetable=covariatetable[idx,]

BEsamples=covariatetable$V2[which(covariatetable$phenoBE_bca==2)]
EAsamples=covariatetable$V2[which(covariatetable$phenoEA_bca==2)]
COsamples=covariatetable$V2[which(covariatetable$phenoEA_bca==1 |covariatetable$phenoBE_bca==1)]
length(BEsamples) #3282
length(EAsamples) #2510
length(COsamples) #10109

filesforGWAS=function(casesamples=BEsamples,controlsamples=COsamples,prefix="BE_CO")
{
  idx1=match(casesamples,fam$V2)
  tmp1=data.frame(FID=fam$V1[idx1],IID=fam$V2[idx1],affected=2,stringsAsFactors = F)
  idx2=match(controlsamples,fam$V2)
  tmp2=data.frame(FID=fam$V1[idx2],IID=fam$V2[idx2],affected=1,stringsAsFactors = F)
  tmp=rbind(tmp1,tmp2)
  #pcadat=read.table(pcafile,stringsAsFactors = F)
  tmp=tmp[tmp$IID %in% covariatetable$V2,]
  idx=match(tmp$IID,covariatetable$V2)
  tmp$pc1=covariatetable$pc1[idx]
  tmp$pc2=covariatetable$pc2[idx]
  tmp$pc3=covariatetable$pc3[idx]
  tmp$pc4=covariatetable$pc4[idx]
  tmp$pc5=covariatetable$pc5[idx]
  tmp$pc6=covariatetable$pc6[idx]
  
  tmp$sex=as.character(covariatetable$sex[idx])
  tmp$sex[which(tmp$sex=="M")]=1 #use for plink
  tmp$sex[as.character(tmp$sex)=="F"]=0
  # #for validate_twas, ind ID changed
  # tmp1=tmp
  # tmp1$IID=paste0(tmp$FID,"_",tmp$IID)
  # tmp1$FID=0
  #write.table(tmp1[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink1.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_hrc_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp,file=paste0("../result/GWAS/",prefix,"_hrc_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  #write.table(tmp1,file=paste0("../result/GWAS/",prefix,"_hrc_selectedsamples_pheno_plink1.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
filesforGWAS()
filesforGWAS(casesamples = EAsamples,prefix="EA_CO")
filesforGWAS(casesamples = c(BEsamples,EAsamples),prefix="BEEA_CO")

udpate_fam=function(famfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar/BE_CO.fam",
                    phenofile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BE_CO_hrc_selectedsamples_pheno_plink.txt")
{
  fam1=read.table(famfile,stringsAsFactors = F)
  pheno=read.table(phenofile,header=T,stringsAsFactors = F)
  idx=match(fam1$V2,pheno$IID)
  fam1$V5=pheno$sex[idx]
  fam1$V5[fam1$V5==1]=2
  fam1$V5[fam1$V5==0]=1
  fam1$V6=pheno$affected[idx]
  write.table(fam1,file=famfile,col.names=F,row.names = F,quote=F)
}
udpate_fam(famfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/BE_CO.fam")
udpate_fam(famfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/EA_CO.fam",
           phenofile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/EA_CO_hrc_selectedsamples_pheno_plink.txt")
udpate_fam(famfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/BEEA_CO.fam",phenofile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BEEA_CO_hrc_selectedsamples_pheno_plink.txt")
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


#change snpname to rsid, chang hg38 to hg19
#load summary data
library(data.table)
summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"
#summary stat from validation
BE_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_Bonn_autosomes.txt"),header=T))
BE_Oxfordsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_oxford_autosomes.txt"),header=T))

EA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"EA_Bonn_autosomes.txt"),header=T))
BEEA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"),header=T))

tmp=intersect(BE_Bonnsummarydat$SNP,EA_Bonnsummarydat$SNP)

tmp=intersect("SNP",colnames(BE_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="pvalue")]="P"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="beta")]="BETA"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

tmp=intersect("SNP",colnames(BE_Oxfordsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="rsid")]="SNP"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="pvalue")]="P"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="beta")]="BETA"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="non-effect-allele")]="non_effect_allele"
}



tmp=intersect("SNP",colnames(EA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="pvalue")]="P"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="beta")]="BETA"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}


tmp=intersect("SNP",colnames(BEEA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="pvalue")]="P"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="beta")]="BETA"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

add_chr_tosummarydat=function(dat1=BE_Bonnsummarydat)
{
  dat2=BE_Oxfordsummarydat
  if (!"chr" %in% colnames(dat1))
  {
    dat1$chr=NA
    tmp=intersect(dat1$SNP,dat2$SNP)
    idx1=match(tmp,dat1$SNP)
    idx2=match(tmp,dat2$SNP)
    dat1$chr[idx1]=dat2$chr[idx2]
    # idx=which(is.na(dat1$chr))
    # tmp=intersect(dat1$SNP[idx],dbsnp$V2)
    # idx1=match(tmp,dat1$SNP)
    # idx2=match(tmp,dbsnp$V2)
    # dat1$chr[idx1]=dbsnp$V1[idx2]
    print(paste0(sum(is.na(dat1$chr))/nrow(dat1),"of all snps are missing chr"))
  }
  dat1$snpname=dat1$snpname1=NA
  dat1$snpname=paste0(dat1$chr,":",dat1$position,"_",dat1$effect_allele,"_",dat1$non_effect_allele)
  dat1$snpname1=paste0(dat1$chr,":",dat1$position,"_",dat1$non_effect_allele,"_",dat1$effect_allele)
  return(dat1)
}
BE_Bonnsummarydat=add_chr_tosummarydat(dat1=BE_Bonnsummarydat)
EA_Bonnsummarydat=add_chr_tosummarydat(dat1=EA_Bonnsummarydat)
BEEA_Bonnsummarydat=add_chr_tosummarydat(dat1=BEEA_Bonnsummarydat)
BE_Oxfordsummarydat=add_chr_tosummarydat(dat1=BE_Oxfordsummarydat)

library("biomaRt")
#mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
#snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp") #hg38

find_snpname1=function(gwasdat,summarydat=BE_Bonnsummarydat)
{
  tmp=unlist(strsplit(gwasdat$SNP,"_"))
  allele1=tmp[seq(2,length(tmp),3)]
  allele2=tmp[seq(3,length(tmp),3)]
  tmp1=tmp[seq(1,length(tmp),3)]
  chrs=unlist(strsplit(tmp1,":"))
  chrs=chrs[seq(1,length(chrs),2)]
  selsnps_hg19=data.frame(chr=chrs,snpname=tmp1,snphg19=gwasdat$SNP,snpname_hg38=paste0(gwasdat$CHR,":",gwasdat$BP))
  selsnps_hg19$snp=NA
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname1)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname1)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  return(selsnps_hg19)
}
getrsid=function(gwasdat)
{
  selsnpmap=find_snpname1(gwasdat = gwasdat,summarydat=BE_Oxfordsummarydat)
  for (chr in 1:22)
  {
    idx=which(selsnpmap$chr==chr & is.na(selsnpmap$snp))
    if (length(idx)>0)
    {
      cat(chr,'..')
      db153=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/db153_chr",chr,".txt"),sep="\t"))
      db153$V1=gsub("chr","",db153$V1)
      db153$snpname=paste0(db153$V1,":",db153$V2+1)
      db153=db153[which(db153$V14=="snv"),]
      #tmp=selsnpmap$snpname_hg38[which(!is.na(selsnpmap$snp) & selsnpmap$chr==chr)]
      tmp=intersect(selsnpmap$snpname_hg38[idx],db153$snpname)
      idx1=match(tmp,selsnpmap$snpname_hg38)
      idx2=match(tmp,db153$snpname)
      selsnpmap$snp[idx1]=db153$V4[idx2]
    }
  }
  
  rsid=selsnpmap$snp
  if (sum(is.na(rsid))>0) print(paste0(as.character(sum(is.na(rsid)))," snps can't find snp rsid"))
  return(selsnpmap)
}

create_sumstats=function(gwasfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/BE_CO.assoc.logistic",prefix="BCAmaf05_BE_CO")
{
  gwasdat=as.data.frame(fread(gwasfile))
  #selsnpmap=getrsid(gwasdat = gwasdat)
  if (any(gwasdat$SNP!=selsnpmap$snphg19)) print("Something is wrong with the selsnpmap")
  
  tmp=unlist(strsplit(gwasdat$SNP,"_"))
  tmp1=tmp[seq(1,length(tmp),3)]
  tmp2=unlist(strsplit(tmp1,":"))
  gwasdat$BP=as.integer(tmp2[seq(2,length(tmp2),2)])
  idx=which(!is.na(gwasdat$SNP))
  gwasdat$SNP[idx]=selsnpmap$snp[idx]
  
  write.table(gwasdat[idx,],file=paste0("../result/GWAS/",prefix,".sumstats"),col.names = T,row.names = F,sep=" ",quote=F)
  
}

gwasdat=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/BE_CO.assoc.logistic"))
#selsnpmap=getrsid(gwasdat = gwasdat)
#save(selsnpmap,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/selsnpmap.RData")
create_sumstats()

create_sumstats(gwasfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/EA_CO.assoc.logistic",prefix="BCAmaf05_EA_CO")
create_sumstats(gwasfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/BEEA_CO.assoc.logistic",prefix="BCAmaf05_BEEA_CO")

sumstatdat=as.data.frame(fread("../result/GWAS/BCAmaf05_BE_CO.sumstats"))

#get frequency file
freqdat=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/selec.frq"))
tmp=unlist(strsplit(freqdat$SNP,":"))
tmp=tmp[seq(2,length(tmp),2)]
tmp1=unlist(strsplit(tmp,"_"))
tmp2=as.integer(tmp1[seq(1,length(tmp1),3)])
all(tmp2==sumstatdat$BP)        #T        
freqdat$SNP=sumstatdat$SNP
write.table(freqdat,file="../result/GWAS/BCAmaf05.frq",col.names = T,row.names = F,sep=" ",quote=F)

sum(freqdat$A1==sumstatdat$A1) #5312259
idx=which(!freqdat$A1==sumstatdat$A1)
all(freqdat$A2[idx]==sumstatdat$A1[idx])

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
