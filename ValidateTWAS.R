
############################################################################
### This code modify Kevin's code and check eQTL prediction model and    ###
### perform TWAS analysis in BEACON data                                 ###
############################################################################

rm(list=ls())

#remove high-correlated snps
removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}


#load BCA covariate tables
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
  return(tmp)
}

covariatetable=readeigenstrat()
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

#Covariate is used for SKAT
opt="PC6"
#try different number of PCs
if (opt=="PC4")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","sex")]
if (opt=="PC6")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]





### loading summary stat data for validation ###

#beaconfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_GTExmodel/"

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

# load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
# proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

#load summary data
library(data.table)
#summary stat from validation
BE_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_Bonn_autosomes.txt"),header=T))
BE_Oxfordsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_oxford_autosomes.txt"),header=T))

EA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"EA_Bonn_autosomes.txt"),header=T))
BEEA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"),header=T))

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


#add chr to summary data (Bonn), used to match snps from BCA and summary data
#use dbsnp data
# dbsnp=NULL
# for (i in 1:22)
# {
#   tmp=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/annotation/dbsnp151_hg19/hg19_",i,".bim")))
#   dbsnp=rbind(dbsnp,tmp)
# }
# save(dbsnp,file="/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/dbsnp151.RData")
# #a large dataset to load
# load("/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/dbsnp151.RData")
# tmp1=intersect(dbsnp$V2,BE_Bonnsummarydat$SNP)
# tmp2=intersect(dbsnp$V2,BE_Oxfordsummarydat$SNP)
# tmp3=intersect(dbsnp$V2,BEEA_Bonnsummarydat$SNP)
# tmp4=intersect(dbsnp$V2,EA_Bonnsummarydat$SNP)
# tmp=unique(c(tmp1,tmp2,tmp3,tmp4))
# idx=match(tmp,dbsnp$V2)
# dbsnp=dbsnp[idx,]
#add chr based on summarydata BE_Oxfordsummarydat
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

#change cooridates of hg38 to hg19. Summarydata use hg19, BCA use hg38
hg38tohg19=function(snpnames=rownames(snp)[1:10])
{
  library(rtracklayer)
  library(GenomicRanges)
  chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
  tmp=unlist(strsplit(snpnames,":"))
  chr0=tmp[1]
  chr=paste0("chr",tmp[1])
  tmp=tmp[seq(2,length(tmp),2)]
  tmp=unlist(strsplit(tmp,"_"))
  pos=as.integer(tmp[seq(1,length(tmp),3)])
  alt=tmp[seq(2,length(tmp),3)]
  ref=tmp[seq(3,length(tmp),3)]
  gr_dat=GRanges(seqnames = chr,ranges=IRanges(start=pos,width=1))
  tmp=liftOver(gr_dat,chain)
  newsnpnames=newpos=rep(NA,length(tmp))
  for (i in 1:length(tmp))
  {
    tmp1=unlist(tmp[i])
    if (length(tmp1)==0)
    {
      warning(paste0(snpnames[i]," not transformed"))
    }else
    {
      if (length(tmp1)==1)
      {
        newpos[i]=start(tmp1)
      }else
      {
        warning(paste0(snpnames[i]," transformed to ",length(tmp1)," snps"))
        newpos[i]=start(tmp1)[1]
      }
    }
  }
  newsnpnames=paste0(chr0,":",newpos,"_",alt,"_",ref)
  res=data.frame(snphg38=snpnames,snphg19=newsnpnames,stringsAsFactors = F)
  return(res)
}

#used to find snp rsid name of BCA snps, to match with summary data
library("biomaRt")
#mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
#snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp") #hg38

#faster than using biomart to find snpnames, but can't find all the results
find_snpname=function(selsnps,summarydat=BE_Bonnsummarydat)
{
  selsnps_hg19=hg38tohg19(snpnames = selsnps)
  selsnps_hg19$selsnps=selsnps
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






#get gwas snps, check if genes are far away from these snps
dong23snp=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/Dong23snp.txt",header=T,sep="\t",stringsAsFactors = F)
extra3snp=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Extra3SNPs.txt",header = T)
dong26snp=rbind(dong23snp[,1:3],extra3snp)
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg19ToHg38.over.chain")


### Find existing GWAS hit ###
gr_dong26snp=GRanges(seqnames = paste0("chr",dong26snp$Chr),ranges=IRanges(start=dong26snp$Position,width=1))
tmp=liftOver(gr_dong26snp,chain)
dong26snp$pos38=NA
for (i in 1:length(tmp))
{
  tmp1=tmp[i]
  if (length(tmp1)>0)
    dong26snp$pos38[i]=as.numeric(start(tmp1))
}
#use hg38
gr_dong26snp=GRanges(seqnames = paste0("chr",dong26snp$Chr),ranges=IRanges(start=dong26snp$pos38,width=1))
idx=match(proteingenes,gtexv8_ge_anno$Symbol)
codinggenetable=gtexv8_ge_anno[idx,]
gr_codinggenetable=GRanges(seqnames = codinggenetable$Chromosome,ranges = IRanges(start=codinggenetable$start,end=codinggenetable$end))
tmp=distanceToNearest(gr_codinggenetable,gr_dong26snp)
codinggenetable$gwas_snp=codinggenetable$dist_to_gwas_snp=NA
codinggenetable$gwas_snp[tmp@from]=dong26snp$SNP[tmp@to]
codinggenetable$dist_to_gwas_snp[tmp@from]=tmp@elementMetadata@listData$distance

#main function starts here---
#get a result for a gene
#how to build gene model
# So I went through diagnosis: here are my choices
# 
# Use the last version of datasets with ambiguous data and : TPM+TMM+standardize
# Include covariates into model selection,
# Use standardize=T in glmnet
# Still do 100 CV to stabilize the selection
# Use correlation filtering threshold 0.9
#phenotype,phenotypepos,snp,snppos,covariate are from GTEx input data. For example, GTExV8adiposedata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData
#bcagenotype is from BCA results. For example, ../result/dist500K_GTEx_adipose_HRC_MAF005/bca_extractgenotype.RData 
library(SKAT)
TWAS_SKAT_gene <- function(genename,cvmin=1,phenotype,phenotypepos,snp,snppos,bcagenotype,covariate,verbose=0) {
  
  idx2=which(codinggenetable$Symbol==genename)
  if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  
  library(glmnet)
  library(GenomicRanges)
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  result <- data.frame(matrix(0,8,4))  
  if (verbose>0) print("get the gene model---")
  i <- which(row.names(phenotype)==genename)
  ncv=10
  distcutoff = 5e5
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)

  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  }  
  #hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  X <- removehighcorr(X1,0.9)
  #hist(ucor[colnames(X1)%in% colnames(X),1])
  #hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  #plot(cvfit)
  
  max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  cvfit$lambda.1se
  
  cverr <- matrix(NA,length(cvfit$lambda),100)
  cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  rownames(cverr)=cvfit$lambda
  for (l in 1:ncol(cverr)) {
    set.seed(l+100)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    alllambda=intersect(cvfit$lambda,fit$lambda)
    idx1=match(alllambda,cvfit$lambda)
    idx2=match(alllambda,fit$lambda)
    cverr[idx1,l] <- fit$cvm[idx2]
    cvse[idx1,l] <- fit$cvsd[idx2]
  }
  
  merr <- apply(cverr,1,mean,na.rm=T)
  mcvse <- sqrt(apply(cvse^2,1,mean))
  
  lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  lambda.min <- cvfit$lambda[which.min(merr)]
  
  #plot(log(cvfit$lambda),merr)
  #abline(v=log(lambda.1se),lty=2)
  #abline(v=log(lambda.min),lty=3)
  
  #plot(cvfit$lambda,merr)
  
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]


  mean(selsnps %in% row.names(bcagenotype))
  
  #if some imputed snps need to to flipped
  
  correctedsnps=NULL
  if (length(intersect(selsnps,rownames(bcagenotype)))<length(selsnps))
  {
    missingsnps=selsnps[!selsnps %in% rownames(bcagenotype)]
    for (j in 1:length(missingsnps))
    {
      tmp=unlist(strsplit(missingsnps[j],"_"))
      tmp1=paste0(tmp[c(1,3,2)],collapse = "_")  #change the order of allele
      idx=which(rownames(bcagenotype)==tmp1)
      if (length(idx)>0)
      {
        correctedsnps=c(correctedsnps,tmp1)
        idx1=which(selsnps==missingsnps[j])
        selsnps[idx1]=tmp1 #change the snp name to make it consistent with bca
      }
    }
  } 
  
  ### still missing snps from bca genotype data ###
  mean(selsnps %in% row.names(bcagenotype))
  
  selectedsnps <- selsnps
  
  if (verbose>0) print("compute skat result----")
  idx1=match(colnames(bcagenotype),rownames(covariatetable))
  covariatetable=covariatetable[idx1,]
  idx1=match(colnames(bcagenotype),rownames(Covariate))
  Covariate=Covariate[idx1,]
  if (any(rownames(covariatetable)!=colnames(bcagenotype))) #covariatetable keeps BCA covariates
    warning("bcagenotype doesn't match covariatetable")
  
  idx=match(selectedsnps,rownames(bcagenotype))
  ## there a few SNPs not found in bcagenotype
  if (sum(is.na(idx))>0)
    print(paste0(sum(is.na(idx))," out of ", length(idx)," selected snps not been found in genotype data"))

  selsnps <- selsnps[!is.na(idx)]
  #print(length(selsnps))
  selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  if (length(correctedsnps)>0) #flip snps
  {
    idxtocorrect=match(correctedsnps,colnames(Z))
    Z[,idxtocorrect]=2-Z[,idxtocorrect]
  }
  #colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
 
  #rownames(Z)=colnames(bcagenotype)
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[6,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[6,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[6,4] <- out$p.value
  rownames(result)[6]="BE_SKAT"
    
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- as.data.frame(cbind(Z1[,i],Covariate1))
  #   outfit <- glm(y~.,data=Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoEA_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[7,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[7,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[7,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[7,4] <- out$p.value
  rownames(result)[7]="EA_SKAT"
  
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEABE_bca==2) #case
  idx2=which(covariatetable$phenoEABE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[8,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[8,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[8,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[8,4] <- out$p.value
  rownames(result)[8]="BEEA_SKAT"
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (nrow(tmp2)>0)
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  rsid=tmp$snp
  if (sum(is.na(rsid))) print(paste0(sum(is.na(rsid)," snps can't find snp rsid")))
  library(survey)
  ### validation in Bonn EA ###
      
  tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
      
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- EA_Bonnsummarydat[EA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[1,1] <- length(tmp)
  result[1,2] <- mean(val1$P<0.05)
  result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[1]="EA_Bonn"    
      
  ### validation in Bonn EA/BE ###
      
  tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[2,1] <- length(tmp)
  result[2,2] <- mean(val1$P<0.05)
  result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[2]="BEEA_Bonn"
  ### validation in Bonn BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BE_Bonnsummarydat[BE_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[3,1] <- length(tmp)
  result[3,2] <- mean(val1$P<0.05)
  result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[3]="BE_Bonn"
  
  ### validation in Oxford BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  valsnps2 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydat$SNP %in% tmp,]
  val2 <- val2[match(rsid1,val2$SNP),]
  uu <- val2$BETA
  vv <- as.numeric(val2$se)
  ZZ <- Z1[,match(valsnps2,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  
  result[4,1] <- length(tmp)
  result[4,2] <- mean(val2$P<0.05)
  result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  #sum(lamb*(1-pchisq(Q,df=1)))
  result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[4]="BE_Oxford"
  
  #### combine the two validation datasets ####
  
  val1 <- val1[val1$SNP %in% val2$SNP,]
  valsnps1 <- valsnps1[valsnps1%in%valsnps2]
  
  val2 <- val2[val2$SNP %in% val1$SNP,]
  valsnps2 <- valsnps2[valsnps2%in%valsnps1]
  if (nrow(val1)>1)  {
    uu1 <- val1$BETA
    vv1 <- as.numeric(val1$SE)
    uu2 <- val2$BETA
    vv2 <- as.numeric(val2$se)
    
    uu <- uu1+uu2
    ZZ <- Z1[,colnames(Z1) %in% valsnps1]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV1 <- diag(vv1) %*% rr %*% diag(vv1)
    VV2 <- diag(vv2) %*% rr %*% diag(vv2)
    VV <- VV1 + VV2
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[5,1] <- nrow(val1)
    result[5,2] <- mean(val1$P<0.05 & val2$P <0.05)
    result[5,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[5,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
    rownames(result)[5]="BE_Combined"
  } 

 return(result)  
}
#to run an example using mucosa HRC data:
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData") #snp,snppos,phenotype,phenotypepose,covariate
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005/bca_extractgenotype.RData")#bca_genotype
res=TWAS_SKAT_gene(genename="SSBP4",phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate=covariate,verbose=1)

#search each organ for validation (HRC models) 
organs=c("mucosa","junction","stomach","muscularis","adipose","blood")
outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005")

#organidx is the index of organs, 1:mucosa  
validate_organs=function(outfolders,opt="HRC")
{
  allresults=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL)
  allresults1=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL) #r2-cutoff 0.05
  allresults2=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL) #r2-cutoff 0.1
  ii=jj=kk=0 #index of significant skat genes
  for (organidx in 1:6) #search organs
  {
    print(paste0("work on ",organs[organidx],"----------------------"))
    outfolder=outfolders[organidx]
    
    load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
    load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
    #load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
    load(paste0(outfolder,"/skat_res.RData")) #saved skat res, skat_min2_pc6,skat_min2
    #load GTEx gene expression data, snp,phenotype
    if (opt=="HRC")
    {
      load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
    }else #1000genome
    {
      load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_for_prediction.RData"))
    }
        
    colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
    types=c("BE","EA","BEA","BEEA")
    skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]
    skat_min2_code_fwer=matrix(NA,nrow=nrow(skat_min2_code),ncol=ncol(skat_min2_code))

    for (i in 1:4)
    {
      skat_min2_code_fwer[,i]=p.adjust(skat_min2_code[,i],method = "bonferroni")
    }
    
    #if a gene was checked in the same organ, no need to do it again
    organgenes=NULL
    for (i in 1:4)
    {
      genenames=rownames(skat_min2_code)[which(skat_min2_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults$genenames=c(allresults$genenames,genename)
          allresults$types=c(allresults$types,types[i])
          allresults$organs=c(allresults$organs,organs[organidx])
          allresults$gwassnp=c(allresults$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults$gwassnpdistance=c(allresults$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          ii=ii+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_SKAT_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          #tmp=TWAS_SKAT_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",ii,", ",genename," may be validated in ",organs[organidx]))
            allresults$validated=c(allresults$validated,T)
          }else
          {
            allresults$validated=c(allresults$validated,F)
          }
          allresults$results[[ii]]=tmp
        }
      }
    }
    
    #r2-cutoff of 0.05
    print("work on r2-cutoff of 0.05........")
    tmp=rownames(res_min)[which(res_min$r2>0.05)]
    skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes & rownames(skat_min2_pc6) %in% tmp,]
    skat_min2_code_fwer=matrix(NA,nrow=nrow(skat_min2_code),ncol=ncol(skat_min2_code))

    for (i in 1:4)
    {
      skat_min2_code_fwer[,i]=p.adjust(skat_min2_code[,i],method = "bonferroni")
    }

    for (i in 1:4)
    {
      genenames=rownames(skat_min2_code)[which(skat_min2_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults1$genenames=c(allresults1$genenames,genename)
          allresults1$types=c(allresults1$types,types[i])
          allresults1$organs=c(allresults1$organs,organs[organidx])
          allresults1$gwassnp=c(allresults1$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults1$gwassnpdistance=c(allresults1$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          
          jj=jj+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_SKAT_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one in the same organ
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",jj,", ",genename," may be validated in ",organs[organidx]))
            allresults1$validated=c(allresults1$validated,T)
          }else
          {
            allresults1$validated=c(allresults1$validated,F)
          }
          allresults1$results[[jj]]=tmp
        }
      }
    }

    
    #r2-cutoff of 0.1
    print("work on r2-cutoff of 0.1........")
    tmp=rownames(res_min)[which(res_min$r2>0.1)]
    skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes & rownames(skat_min2_pc6) %in% tmp,]
    skat_min2_code_fwer=matrix(NA,nrow=nrow(skat_min2_code),ncol=ncol(skat_min2_code))

    for (i in 1:4)
    {
      skat_min2_code_fwer[,i]=p.adjust(skat_min2_code[,i],method = "bonferroni")
    }

    for (i in 1:4)
    {
      genenames=rownames(skat_min2_code)[which(skat_min2_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults2$genenames=c(allresults2$genenames,genename)
          allresults2$types=c(allresults2$types,types[i])
          allresults2$organs=c(allresults2$organs,organs[organidx])
          allresults2$gwassnp=c(allresults2$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults2$gwassnpdistance=c(allresults2$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          
          kk=kk+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_SKAT_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one in the same organ
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          
          #tmp=TWAS_SKAT_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype)
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",kk,", ",genename," may be validated in ",organs[organidx]))
            allresults2$validated=c(allresults2$validated,T)
          }else
          {
            allresults2$validated=c(allresults2$validated,F)
          }
          allresults2$results[[kk]]=tmp
        }
      }
    }
    
  }
  return(list(allresults=allresults,allresults1=allresults1,allresults2=allresults2))
  #return(allresults)
}
  
HRCvalidres=validate_organs(outfolders=outfolders)
save(HRCvalidres,file="../result/ValidateTWASres.RData")

# [1] "work on mucosa..."
# [1] "work on BE: SSBP4. Its closest gwas snp: rs10419226, distance:0.26MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: COX7A2. Its closest gwas snp: rs76014404, distance:13.56MB"
# [1] "work on BEEA: SSBP4. Its closest gwas snp: rs10419226, distance:0.26MB"
# [1] "work on BEEA: TMEM161A. Its closest gwas snp: rs10423674, distance:0.41MB"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on BE: LRPAP1. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on BE: LRPAP1. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BE: TMEM260. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on junction----------------------"
# [1] "work on BE: FILIP1. Its closest gwas snp: rs76014404, distance:13.61MB"
# [1] "work on BE: UBA52. Its closest gwas snp: rs10419226, distance:0.11MB"
# [1] "idx 7, UBA52 may be validated in junction"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: UBA52. Its closest gwas snp: rs10419226, distance:0.11MB"
# [1] "idx 9, UBA52 may be validated in junction"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on BEEA: BORCS8. Its closest gwas snp: rs10423674, distance:0.47MB"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on BEEA: BORCS8. Its closest gwas snp: rs10423674, distance:0.47MB"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: CRLF1. Its closest gwas snp: rs10419226, distance:0.08MB"
# [1] "idx 5, CRLF1 may be validated in junction"
# [1] "work on stomach----------------------"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "idx 12, SLC25A42 may be validated in stomach"
# [1] "work on BEEA: COPE. Its closest gwas snp: rs10423674, distance:0.19MB"
# [1] "idx 13, COPE may be validated in stomach"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on BEEA: TMEM30A. Its closest gwas snp: rs76014404, distance:13.57MB"
# [1] "idx 6, TMEM30A may be validated in stomach"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on muscularis----------------------"
# [1] "work on EA: AHR. Its closest gwas snp: rs11765529, distance:35.51MB"
# [1] "work on BEEA: FOXF1. Its closest gwas snp: rs1979654, distance:0.15MB"
# [1] "work on BEEA: TMEM59L. Its closest gwas snp: rs10419226, distance:0.07MB"
# [1] "idx 16, TMEM59L may be validated in muscularis"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on EA: FOXF1. Its closest gwas snp: rs1979654, distance:0.15MB"
# [1] "work on BEEA: FOXF1. Its closest gwas snp: rs1979654, distance:0.15MB"
# [1] "work on BEEA: LRRC25. Its closest gwas snp: rs10419226, distance:0.29MB"
# [1] "work on BEEA: TMEM59L. Its closest gwas snp: rs10419226, distance:0.07MB"
# [1] "idx 10, TMEM59L may be validated in muscularis"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on EA: FOXF1. Its closest gwas snp: rs1979654, distance:0.15MB"
# [1] "work on BEEA: EXOC3. Its closest gwas snp: rs9918259, distance:0.19MB"
# [1] "idx 8, EXOC3 may be validated in muscularis"
# [1] "work on BEEA: FOXF1. Its closest gwas snp: rs1979654, distance:0.15MB"
# [1] "work on adipose----------------------"
# [1] "work on BE: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on BE: CRTC1. Its closest gwas snp: rs10419226, distance:0MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "idx 19, LDAH may be validated in adipose"
# [1] "work on BEEA: ALDH1A2. Its closest gwas snp: rs66725070, distance:0MB"
# [1] "work on BEEA: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on BEEA: CRTC1. Its closest gwas snp: rs10419226, distance:0MB"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on BE: ALDH1A2. Its closest gwas snp: rs66725070, distance:0MB"
# [1] "work on BE: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on BE: CRTC1. Its closest gwas snp: rs10419226, distance:0MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: ALDH1A2. Its closest gwas snp: rs66725070, distance:0MB"
# [1] "work on BEEA: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on BEEA: LRRC25. Its closest gwas snp: rs10419226, distance:0.29MB"
# [1] "idx 17, LRRC25 may be validated in adipose"
# [1] "work on BEEA: TMEM59L. Its closest gwas snp: rs10419226, distance:0.07MB"
# [1] "idx 18, TMEM59L may be validated in adipose"
# [1] "work on BEEA: CRTC1. Its closest gwas snp: rs10419226, distance:0MB"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on BE: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on EA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: AQP9. Its closest gwas snp: rs2464469, distance:0.07MB"
# [1] "work on blood----------------------"
# [1] "work on BE: PGPEP1. Its closest gwas snp: rs10419226, distance:0.32MB"
# [1] "idx 23, PGPEP1 may be validated in blood"
# [1] "work on BE: LSM4. Its closest gwas snp: rs10419226, distance:0.37MB"
# [1] "work on BE: KLHL26. Its closest gwas snp: rs10419226, distance:0.02MB"
# [1] "work on EA: IL2RB. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: COX7A2. Its closest gwas snp: rs76014404, distance:13.56MB"
# [1] "work on BEEA: PGPEP1. Its closest gwas snp: rs10419226, distance:0.32MB"
# [1] "idx 29, PGPEP1 may be validated in blood"
# [1] "work on BEEA: LSM4. Its closest gwas snp: rs10419226, distance:0.37MB"
# [1] "work on BEEA: GDF15. Its closest gwas snp: rs10419226, distance:0.3MB"
# [1] "idx 31, GDF15 may be validated in blood"
# [1] "work on BEEA: KLHL26. Its closest gwas snp: rs10419226, distance:0.02MB"
# [1] "work on BEEA: HOMER3. Its closest gwas snp: rs10423674, distance:0.22MB"
# [1] "idx 33, HOMER3 may be validated in blood"
# [1] "work on BEEA: COMP. Its closest gwas snp: rs10423674, distance:0.08MB"
# [1] "idx 34, COMP may be validated in blood"
# [1] "work on r2-cutoff of 0.05........"
# [1] "work on BE: HSP90AA1. Its closest gwas snp: NA, distance:NAMB"
# [1] "idx 20, HSP90AA1 may be validated in blood"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: SENP6. Its closest gwas snp: rs76014404, distance:13.92MB"
# [1] "idx 22, SENP6 may be validated in blood"
# [1] "work on BEEA: HSP90AA1. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: COMP. Its closest gwas snp: rs10423674, distance:0.08MB"
# [1] "idx 24, COMP may be validated in blood"
# [1] "work on r2-cutoff of 0.1........"
# [1] "work on BEEA: GDF7. Its closest gwas snp: rs3072, distance:0MB"


#check FILIP1, junction, BE
which(HRCvalidres$allresults$genenames=="FILIP1")
HRCvalidres$allresults$results[6]
# X1           X2           X3           X4
# EA_Bonn     22 9.090909e-02 1.216351e-01 1.157589e-01
# BEEA_Bonn   22 9.090909e-02 9.813946e-02 9.459537e-02
# BE_Bonn     22 9.090909e-02 2.202139e-01 2.088936e-01
# BE_Oxford   25 0.000000e+00 6.054553e-01 5.970953e-01
# BE_Combined 22 0.000000e+00 2.779844e-01 2.637301e-01
# BE_SKAT     26 2.484606e-06 1.032852e-05 7.460313e-05
# EA_SKAT     26 2.219994e-01 3.683174e-01 2.341559e-01
# BEEA_SKAT   26 4.871763e-06 1.823809e-05 1.844613e-04

#New genes been validated (HRC)--------------
#TMEM30A in stomach, r2-cutoff 0.05, BEEA
HRCvalidres$allresults1$results[6]
# EA_Bonn     82 7.317073e-02 7.757477e-03 0.0101408496
# BEEA_Bonn   82 1.097561e-01 5.094947e-03 0.0071391385
# BE_Bonn     83 9.638554e-02 1.688875e-01 0.1627855676
# BE_Oxford   89 1.123596e-02 9.251316e-01 0.9358895612
# BE_Combined 83 1.204819e-02 2.402131e-01 0.2299347796
# BE_SKAT     89 2.883422e-05 8.812155e-05 0.0002517309
# EA_SKAT     89 2.515922e-02 4.818494e-02 0.0880188081
# BEEA_SKAT   89 4.888336e-06 1.750600e-05 0.0001561675

#HSP90AA1 in blood, r2-cutoff 0.05,BE BEEA
# HRCvalidres$allresults1$results[20]
# EA_Bonn     32 9.375000e-02 9.329312e-02 0.0899604007
# BEEA_Bonn   32 3.125000e-02 2.160169e-01 0.2024249755
# BE_Bonn     32 1.562500e-01 1.875095e-02 0.0212780130
# BE_Oxford   33 2.121212e-01 3.009490e-02 0.0322849234
# BE_Combined 32 3.125000e-02 2.212620e-04 0.0006556874
# BE_SKAT     33 1.023311e-05 2.478727e-05 0.0126478030
# EA_SKAT     33 3.483814e-02 6.221111e-02 0.0149643446
# BEEA_SKAT   33 7.591476e-06 1.903788e-05 0.0015690568

#SENP6 in blood, r2-cutoff 0.05, BEEA
# EA_Bonn     35 8.571429e-02 7.771989e-02 0.0765877075
# BEEA_Bonn   35 5.714286e-02 2.992269e-02 0.0318517705
# BE_Bonn     35 8.571429e-02 2.885195e-01 0.2754339564
# BE_Oxford   37 0.000000e+00 9.025914e-01 0.9148757458
# BE_Combined 35 0.000000e+00 5.735766e-01 0.5634619984
# BE_SKAT     37 1.417286e-05 4.863982e-05 0.0004864302
# EA_SKAT     37 3.251049e-02 6.261970e-02 0.2789702810
# BEEA_SKAT   37 3.738383e-06 1.438743e-05 0.0005258976

#1000genome result
outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_Jan19")
oneKvalidres=validate_organs(outfolders=outfolders,opt="100genome")
save(HRCvalidres,oneKvalidres,file="../result/ValidateTWASres.RData")
# [1] "work on mucosa..."
# [1] "work on BE: SSBP4. Its closest gwas snp: rs10419226, distance:0.26MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: JUND. Its closest gwas snp: rs10419226, distance:0.41MB"
# [1] "work on BEEA: SSBP4. Its closest gwas snp: rs10419226, distance:0.26MB"
# [1] "work on BEEA: TMEM161A. Its closest gwas snp: rs10423674, distance:0.41MB"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on r2-cutoff of 0.1--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on junction..."
# [1] "work on BE: FILIP1. Its closest gwas snp: rs76014404, distance:13.61MB"
# [1] "idx 6, FILIP1 may be validated in junction"
# [1] "work on BE: UBA52. Its closest gwas snp: rs10419226, distance:0.11MB"
# [1] "idx 7, UBA52 may be validated in junction"
# [1] "work on EA: TMEM161A. Its closest gwas snp: rs10423674, distance:0.41MB"
# [1] "idx 8, TMEM161A may be validated in junction"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "1 out of 46 selected snps not been found in genotype data"
# [1] "idx 9, LDAH may be validated in junction"
# [1] "work on BEEA: COX7A2. Its closest gwas snp: rs76014404, distance:13.56MB"
# [1] "idx 10, COX7A2 may be validated in junction"
# [1] "work on BEEA: FILIP1. Its closest gwas snp: rs76014404, distance:13.61MB"
# [1] "idx 11, FILIP1 may be validated in junction"
# [1] "work on BEEA: UBA52. Its closest gwas snp: rs10419226, distance:0.11MB"
# [1] "idx 12, UBA52 may be validated in junction"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on BEEA: TMEM161A. Its closest gwas snp: rs10423674, distance:0.41MB"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on r2-cutoff of 0.1--------------"
# [1] "work on BE: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: FKBP8. Its closest gwas snp: rs10419226, distance:0.15MB"
# [1] "idx 4, FKBP8 may be validated in junction"
# [1] "work on BEEA: SLC25A42. Its closest gwas snp: rs10423674, distance:0.36MB"
# [1] "work on stomach..."
# [1] "work on BEEA: COPE. Its closest gwas snp: rs10423674, distance:0.19MB"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on r2-cutoff of 0.1--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "idx 6, LDAH may be validated in stomach"
# [1] "work on muscularis..."
# [1] "work on BEEA: COPE. Its closest gwas snp: rs10423674, distance:0.19MB"
# [1] "1 out of 48 selected snps not been found in genotype data"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on r2-cutoff of 0.1--------------"
# [1] "work on BE: ASAH2. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on adipose..."
# [1] "work on BEEA: COPE. Its closest gwas snp: rs10423674, distance:0.19MB"
# [1] "idx 17, COPE may be validated in adipose"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "idx 5, LDAH may be validated in adipose"
# [1] "work on r2-cutoff of 0.1--------------"
# [1] "work on BE: ASAH2. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on blood..."
# [1] "work on BE: PGPEP1. Its closest gwas snp: rs10419226, distance:0.32MB"
# [1] "idx 18, PGPEP1 may be validated in blood"
# [1] "work on BE: KLHL26. Its closest gwas snp: rs10419226, distance:0.02MB"
# [1] "work on EA: IL2RB. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: PGPEP1. Its closest gwas snp: rs10419226, distance:0.32MB"
# [1] "idx 22, PGPEP1 may be validated in blood"
# [1] "work on BEEA: LSM4. Its closest gwas snp: rs10419226, distance:0.37MB"
# [1] "work on BEEA: GDF15. Its closest gwas snp: rs10419226, distance:0.3MB"
# [1] "idx 24, GDF15 may be validated in blood"
# [1] "work on BEEA: KLHL26. Its closest gwas snp: rs10419226, distance:0.02MB"
# [1] "work on BEEA: COMP. Its closest gwas snp: rs10423674, distance:0.08MB"
# [1] "idx 26, COMP may be validated in blood"
# [1] "work on BEEA: HOMER3. Its closest gwas snp: rs10423674, distance:0.22MB"
# [1] "idx 27, HOMER3 may be validated in blood"
# [1] "work on r2-cutoff of 0.05--------------"
# [1] "work on BE: HSP90AA1. Its closest gwas snp: NA, distance:NAMB"
# [1] "idx 6, HSP90AA1 may be validated in blood"
# [1] "work on BEEA: LDAH. Its closest gwas snp: rs7255, distance:0MB"
# [1] "work on BEEA: HSP90AA1. Its closest gwas snp: NA, distance:NAMB"
# [1] "work on BEEA: COMP. Its closest gwas snp: rs10423674, distance:0.08MB"
# [1] "idx 9, COMP may be validated in blood"
# [1] "work on r2-cutoff of 0.1--------------"

#New genes been validated (1000genome)--------------
#COX7A2 in junction was validated, BEEA
oneKvalidres$allresults$results[10]
# X1           X2           X3           X4
# EA_Bonn     22 9.090909e-02 6.180768e-04 0.0011433202
# BEEA_Bonn   22 1.363636e-01 4.861973e-04 0.0009255241
# BE_Bonn     22 1.818182e-01 7.845778e-02 0.0753211797
# BE_Oxford   23 4.347826e-02 5.854575e-01 0.5770689779
# BE_Combined 22 0.000000e+00 8.534918e-02 0.0822651721
# BE_SKAT     23 5.504085e-05 1.592340e-04 0.0601697640
# EA_SKAT     23 5.129187e-03 1.083102e-02 0.0900800532
# BEEA_SKAT   23 2.252681e-06 8.797407e-06 0.0105436994

#FILIP1 in junction was validated, BEEA,BE
oneKvalidres$allresults$results[11]
# X1           X2           X3           X4
# EA_Bonn     19 1.052632e-01 3.165626e-02 3.261061e-02
# BEEA_Bonn   19 1.052632e-01 7.056089e-03 9.605639e-03
# BE_Bonn     19 1.052632e-01 1.083348e-01 9.844089e-02
# BE_Oxford   22 0.000000e+00 6.084927e-01 5.996947e-01
# BE_Combined 19 0.000000e+00 2.238063e-01 2.020592e-01
# BE_SKAT     23 1.014777e-06 5.269450e-06 6.949034e-05
# EA_SKAT     23 5.865351e-02 1.087750e-01 3.360826e-01
# BEEA_SKAT   23 7.873484e-07 4.197923e-06 1.145102e-04

#HSP90AA1 in blood, r2-cutoff 0.05, was validated, BE
oneKvalidres$allresults1$results[6]
# X1           X2           X3           X4
# EA_Bonn     30 3.333333e-02 5.292673e-01 5.146141e-01
# BEEA_Bonn   30 3.333333e-02 5.054913e-01 4.897769e-01
# BE_Bonn     30 2.000000e-01 3.144782e-03 4.813146e-03
# BE_Oxford   31 2.258065e-01 4.546762e-03 6.647730e-03
# BE_Combined 30 3.333333e-02 1.712827e-06 2.458087e-05
# BE_SKAT     32 6.744045e-06 1.695537e-05 1.437565e-02
# EA_SKAT     32 2.229023e-02 4.026058e-02 3.557269e-02
# BEEA_SKAT   32 4.122048e-06 1.119404e-05 2.592488e-03

#test regular-TWAS gene------------
#x is the predict_exp
pvalue_assoc=function(x,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  x=as.numeric(x[c(idx1,idx2)])
  BE_p=NA
  if (sum(is.na(x))<0.6*length(x))
  {
    #fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=covariatetable[c(idx1,idx2),],family=binomial)
    fm=glm(y~x+pc1+pc2+pc3+pc4+pc5+pc6+sex,data=covariatetable[c(idx1,idx2),],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariatetable[c(idx1,idx2),],family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      BE_p=summary(fm)$coefficients[2,4]
    }
  }
  return(BE_p)
}

TWAS_Regular_gene <- function(genename,cvmin=1,phenotype,phenotypepos,snp,snppos,bcagenotype,covariate,verbose=0) {
  
  idx2=which(codinggenetable$Symbol==genename)
  if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  
  library(glmnet)
  library(GenomicRanges)
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  result <- data.frame(matrix(0,8,4))  
  if (verbose>0) print("get the gene model---")
  i <- which(row.names(phenotype)==genename)
  ncv=10
  distcutoff = 5e5
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  
  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  }  
  #hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  X <- removehighcorr(X1,0.9)
  #hist(ucor[colnames(X1)%in% colnames(X),1])
  #hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  #plot(cvfit)
  
  max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  cvfit$lambda.1se
  
  cverr <- matrix(NA,length(cvfit$lambda),100)
  cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  rownames(cverr)=cvfit$lambda
  for (l in 1:ncol(cverr)) {
    set.seed(l+100)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    alllambda=intersect(cvfit$lambda,fit$lambda)
    idx1=match(alllambda,cvfit$lambda)
    idx2=match(alllambda,fit$lambda)
    cverr[idx1,l] <- fit$cvm[idx2]
    cvse[idx1,l] <- fit$cvsd[idx2]
  }
  
  merr <- apply(cverr,1,mean,na.rm=T)
  mcvse <- sqrt(apply(cvse^2,1,mean))
  
  lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  lambda.min <- cvfit$lambda[which.min(merr)]
  
  #plot(log(cvfit$lambda),merr)
  #abline(v=log(lambda.1se),lty=2)
  #abline(v=log(lambda.min),lty=3)
  
  #plot(cvfit$lambda,merr)
  
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  
  mean(selsnps %in% row.names(bcagenotype))
  
  #if some imputed snps need to to flipped
  
  correctedsnps=NULL
  if (length(intersect(selsnps,rownames(bcagenotype)))<length(selsnps))
  {
    missingsnps=selsnps[!selsnps %in% rownames(bcagenotype)]
    for (j in 1:length(missingsnps))
    {
      tmp=unlist(strsplit(missingsnps[j],"_"))
      tmp1=paste0(tmp[c(1,3,2)],collapse = "_")  #change the order of allele
      idx=which(rownames(bcagenotype)==tmp1)
      if (length(idx)>0)
      {
        correctedsnps=c(correctedsnps,tmp1)
        idx1=which(selsnps==missingsnps[j])
        selsnps[idx1]=tmp1 #change the snp name to make it consistent with bca
      }
    }
  } 
  
  ### still missing snps from bca genotype data ###
  mean(selsnps %in% row.names(bcagenotype))
  
  selectedsnps <- selsnps
  
  if (verbose>0) print("compute assoc result----")
  idx1=match(colnames(bcagenotype),rownames(covariatetable))
  covariatetable=covariatetable[idx1,]
  #idx1=match(colnames(bcagenotype),rownames(Covariate))
  #Covariate=Covariate[idx1,]
  if (any(rownames(covariatetable)!=colnames(bcagenotype))) #covariatetable keeps BCA covariates
    warning("bcagenotype doesn't match covariatetable")
  
  idx=match(selectedsnps,rownames(bcagenotype))
  ## there a few SNPs not found in bcagenotype
  if (sum(is.na(idx))>0)
    print(paste0(sum(is.na(idx))," out of ", length(idx)," selected snps not been found in genotype data"))
  
  selsnps <- selsnps[!is.na(idx)]
  #print(length(selsnps))
  selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  if (length(correctedsnps)>0) #flip snps
  {
    idxtocorrect=match(correctedsnps,colnames(Z))
    Z[,idxtocorrect]=2-Z[,idxtocorrect]
  }
  #colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
  predict_exp=Z %*% selcoeff
  
  all(rownames(Z)==rownames(covariatetable))
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)

  #y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]

  result[6,1] <- ncol(Z1)

  result[6,2] <- pvalue_assoc(x=predict_exp[,1],idx1=idx1,idx2=idx2)
  
  rownames(result)[6]="BE"
  
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- as.data.frame(cbind(Z1[,i],Covariate1))
  #   outfit <- glm(y~.,data=Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoEA_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  #Covariate1=as.data.frame(Covariate[c(idx1,idx2),])

  result[7,1] <- ncol(Z1)
  result[7,2]=pvalue_assoc(x=predict_exp[,1],idx1=idx1,idx2=idx2)
  rownames(result)[7]="EA"
  
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEABE_bca==2) #case
  idx2=which(covariatetable$phenoEABE_bca==1)
  
  Z1=Z[c(idx1,idx2),,drop=F]
  result[8,1] <- ncol(Z1)
  result[8,2]=pvalue_assoc(x=predict_exp[,1],idx1=idx1,idx2=idx2)
  rownames(result)[8]="BEEA"
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (nrow(tmp2)>0)
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  rsid=tmp$snp
  if (sum(is.na(rsid))) print(paste0(sum(is.na(rsid)," snps can't find snp rsid")))
  library(survey)
  ### validation in Bonn EA ###
  
  tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- EA_Bonnsummarydat[EA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[1,1] <- length(tmp)
  result[1,2] <- mean(val1$P<0.05)
  result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[1]="EA_Bonn"    
  
  ### validation in Bonn EA/BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[2,1] <- length(tmp)
  result[2,2] <- mean(val1$P<0.05)
  result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[2]="BEEA_Bonn"
  ### validation in Bonn BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BE_Bonnsummarydat[BE_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[3,1] <- length(tmp)
  result[3,2] <- mean(val1$P<0.05)
  result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[3]="BE_Bonn"
  
  ### validation in Oxford BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  valsnps2 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydat$SNP %in% tmp,]
  val2 <- val2[match(rsid1,val2$SNP),]
  uu <- val2$BETA
  vv <- as.numeric(val2$se)
  ZZ <- Z1[,match(valsnps2,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  
  result[4,1] <- length(tmp)
  result[4,2] <- mean(val2$P<0.05)
  result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  #sum(lamb*(1-pchisq(Q,df=1)))
  result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[4]="BE_Oxford"
  
  #### combine the two validation datasets ####
  
  val1 <- val1[val1$SNP %in% val2$SNP,]
  valsnps1 <- valsnps1[valsnps1%in%valsnps2]
  
  val2 <- val2[val2$SNP %in% val1$SNP,]
  valsnps2 <- valsnps2[valsnps2%in%valsnps1]
  if (nrow(val1)>1)  {
    uu1 <- val1$BETA
    vv1 <- as.numeric(val1$SE)
    uu2 <- val2$BETA
    vv2 <- as.numeric(val2$se)
    
    uu <- uu1+uu2
    ZZ <- Z1[,colnames(Z1) %in% valsnps1]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV1 <- diag(vv1) %*% rr %*% diag(vv1)
    VV2 <- diag(vv2) %*% rr %*% diag(vv2)
    VV <- VV1 + VV2
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[5,1] <- nrow(val1)
    result[5,2] <- mean(val1$P<0.05 & val2$P <0.05)
    result[5,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[5,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
    rownames(result)[5]="BE_Combined"
  } 
  
  return(result)  
}

validate_regular_organs=function(outfolders,opt="HRC")
{
  allresults=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL)
  allresults1=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL) #r2-cutoff 0.05
  allresults2=list(types=NULL,genenames=NULL,organs=NULL,results=list(),gwassnp=NULL,gwassnpdistance=NULL,validated=NULL) #r2-cutoff 0.1
  ii=jj=kk=0 #index of significant skat genes
  for (organidx in 1:6) #search organs
  {
    print(paste0("work on ",organs[organidx],"----------------------"))
    outfolder=outfolders[organidx]
    
    load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
    load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
    #load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
    load(paste0(outfolder,"/bca_assoc.RData")) #saved assoc res, assoc_min_pc6
    #load GTEx gene expression data, snp,phenotype
    if (opt=="HRC")
    {
      load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
    }else #1000genome
    {
      load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_for_prediction.RData"))
    }
    if (ncol(assoc_min_pc6)==8) assoc_min_pc6=assoc_min_pc6[,c(1,3,5,7)]
    colnames(assoc_min_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
    types=c("BE","EA","BEA","BEEA")
    assoc_min_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]
    assoc_min_code_fwer=matrix(NA,nrow=nrow(assoc_min_code),ncol=ncol(assoc_min_code))
    
    for (i in 1:4)
    {
      assoc_min_code_fwer[,i]=p.adjust(assoc_min_code[,i],method = "bonferroni")
    }
    
    #if a gene was checked in the same organ, no need to do it again
    organgenes=NULL
    for (i in 1:4)
    {
      genenames=rownames(assoc_min_code)[which(assoc_min_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults$genenames=c(allresults$genenames,genename)
          allresults$types=c(allresults$types,types[i])
          allresults$organs=c(allresults$organs,organs[organidx])
          allresults$gwassnp=c(allresults$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults$gwassnpdistance=c(allresults$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          ii=ii+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_Regular_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          #tmp=TWAS_Regular_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",ii,", ",genename," may be validated in ",organs[organidx]))
            allresults$validated=c(allresults$validated,T)
          }else
          {
            allresults$validated=c(allresults$validated,F)
          }
          allresults$results[[ii]]=tmp
        }
      }
    }
    
    #r2-cutoff of 0.05
    print("work on r2-cutoff of 0.05........")
    tmp=rownames(res_min)[which(res_min$r2>0.05)]
    assoc_min_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes & rownames(assoc_min_pc6) %in% tmp,]
    assoc_min_code_fwer=matrix(NA,nrow=nrow(assoc_min_code),ncol=ncol(assoc_min_code))
    
    for (i in 1:4)
    {
      assoc_min_code_fwer[,i]=p.adjust(assoc_min_code[,i],method = "bonferroni")
    }
    
    for (i in 1:4)
    {
      genenames=rownames(assoc_min_code)[which(assoc_min_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults1$genenames=c(allresults1$genenames,genename)
          allresults1$types=c(allresults1$types,types[i])
          allresults1$organs=c(allresults1$organs,organs[organidx])
          allresults1$gwassnp=c(allresults1$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults1$gwassnpdistance=c(allresults1$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          
          jj=jj+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_Regular_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one in the same organ
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",jj,", ",genename," may be validated in ",organs[organidx]))
            allresults1$validated=c(allresults1$validated,T)
          }else
          {
            allresults1$validated=c(allresults1$validated,F)
          }
          allresults1$results[[jj]]=tmp
        }
      }
    }
    
    
    #r2-cutoff of 0.1
    print("work on r2-cutoff of 0.1........")
    tmp=rownames(res_min)[which(res_min$r2>0.1)]
    assoc_min_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes & rownames(assoc_min_pc6) %in% tmp,]
    assoc_min_code_fwer=matrix(NA,nrow=nrow(assoc_min_code),ncol=ncol(assoc_min_code))
    
    for (i in 1:4)
    {
      assoc_min_code_fwer[,i]=p.adjust(assoc_min_code[,i],method = "bonferroni")
    }
    
    for (i in 1:4)
    {
      genenames=rownames(assoc_min_code)[which(assoc_min_code_fwer[,i]<0.05)]
      if (length(genenames)>0)
      {
        for (j in 1:length(genenames))
        {
          genename=genenames[j]
          idx2=which(codinggenetable$Symbol==genename)
          print(paste0("work on ",types[i],": ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
          allresults2$genenames=c(allresults2$genenames,genename)
          allresults2$types=c(allresults2$types,types[i])
          allresults2$organs=c(allresults2$organs,organs[organidx])
          allresults2$gwassnp=c(allresults2$gwassnp,codinggenetable$gwas_snp[idx2])
          allresults2$gwassnpdistance=c(allresults2$gwassnpdistance,codinggenetable$dist_to_gwas_snp[idx2])
          
          kk=kk+1
          if (!genename %in% organgenes) #new gene
          {
            tmp=TWAS_Regular_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate = covariate)
          }else
          {
            idx1=which(allresults$genenames==genename)[length(idx1)] #copy the last one in the same organ
            tmp=allresults$results[[idx1]]
          }
          organgenes=unique(c(organgenes,genename))
          
          #tmp=TWAS_Regular_gene(genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype)
          if (any(tmp[1:4,3:4]<0.05))
          {
            print(paste0("idx ",kk,", ",genename," may be validated in ",organs[organidx]))
            allresults2$validated=c(allresults2$validated,T)
          }else
          {
            allresults2$validated=c(allresults2$validated,F)
          }
          allresults2$results[[kk]]=tmp
        }
      }
    }
    
  }
  return(list(allresults=allresults,allresults1=allresults1,allresults2=allresults2))
  #return(allresults)
}

outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005")
HRCregvalidres=validate_regular_organs(outfolders=outfolders)
#save(HRCvalidres,file="../result/ValidateTWASres.RData")

outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_Jan19")
oneKregvalidres=validate_regular_organs(outfolders=outfolders,opt="100genome")
save(HRCregvalidres,oneKregvalidres,file="../result/Validate_REGTWASres.RData")

#check HOXB8, muscularis,1000genome,(pc4) BEEA
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19"
load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
#load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
#load(paste0(outfolder,"/bca_assoc.RData")) #saved assoc res, assoc_min_pc6
#load GTEx gene expression data, snp,phenotype
load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8","muscularis","data_ambiguous_TPM_addcontrols_for_prediction.RData"))
res=TWAS_Regular_gene(genename="HOXB8",phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate=covariate,verbose=1)

