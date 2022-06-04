
###############################################################################
##### Discover and Validate several NEW genes in summary stat data      #######
##### Aug 2021                                                          #######
###############################################################################
#for rmhighcor_allcovar results, use 1000 Genome European to calculate correlation in validation.


removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  ## first remove redundant SNPs corr=1 ##
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= 1))]
  ## next remove high correlation SNPs but not exactly same ##
  tmp <- cor(datnew)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  dlist <- (1:ncol(datnew))[apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  rlist <- (1:ncol(datnew))[!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  clist <- NULL
  for (l in dlist) {
    clist <- c(clist, which(abs(tmp[,l])>=corcutoff))
  }  
  ## this is to make sure every deleted SNP has one tagSNP retained
  l<-1
  while (l<=length(dlist)) {
    ll <- dlist[l]
    tt <- which(abs(tmp[,ll])>=corcutoff)
    if (sum(tt %in% rlist)==0) {
        tcount <- rep(0,length(tt))
        for (k in 1:length(tt)) {
          tcount[k] <- sum(clist==tt[k])
        } 
        rlist <- c(rlist,tt[which.max(tcount)])
        dlist <- dlist[-which(dlist==tt[which.max(tcount)])]
    }  
    l <- l+1
  }  
  rlist <- rlist[order(rlist)]
  datnew <- datnew[,rlist]
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

#work on validation data
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


TWAS_SKAT_gene <- function(genename,cvmin=1,phenotype,phenotypepos,snp,snppos,bcagenotype,covariate,verbose=0) {
  
  #idx2=which(codinggenetable$Symbol==genename)
  #if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  library(SKAT)
  library(glmnet)
  library(GenomicRanges)
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  result <- data.frame(matrix(NA,8,4))  
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
  dim(X1)
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
  penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
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
  cvfit$lambda.min
  
  
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
  
  #lambda.min <- cvfit$lambda.min
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1,drop=F]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  ## check the marginal eQTL association ##
  outp <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat <- data.matrix(cbind(X[,colnames(X)==selsnps[i]],covariate))
    fit <- lm(Y~Xmat)
    outp[i] <- summary(fit)$coef[2,4]
  }  
  #outp
  
  mean(selsnps %in% row.names(bcagenotype))
  ## here we should have a function to get the genotypes for SNPs being selected in selsnps ##
  #bcagenotype=extractBCAgenotype(selsnps = selsnps) #takes 2 minutes for 38 snps
  
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
  #result[6,]
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
  #result[7,]
  #for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
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
  #result[8,]
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
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
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
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
  if (sum(is.na(rsid))>0) print(paste0(as.character(sum(is.na(rsid)))," snps can't find snp rsid"))
  library(survey)
  
  ### validation in Bonn EA ###
  val1=val2=NULL
  tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  idx1=match(tmp,EA_Bonnsummarydat$SNP)
  chr=EA_Bonnsummarydat$chr[idx1][1]
  #read data from 1000 genome v3 EUR
  thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
  refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
  refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
  refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
  idx1=match(refsnps,refbim$V2)
  Z1=as.data.frame(refraw[idx1,7:ncol(refraw)])
  rownames(Z1)=refsnps
  Z1=as.data.frame(t(Z1))
  tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  tmp=intersect(tmp,refsnps) #snps should be available in 1000G ref as well
  
  if (length(tmp)>1)
  {
    valsnps1 <- selsnps[which(rsid %in% tmp)]
    rsid1 <- rsid[which(rsid %in% tmp)]
    val1 <- EA_Bonnsummarydat[EA_Bonnsummarydat$SNP %in% tmp,]
    val1 <- val1[match(rsid1,val1$SNP),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z1[,match(rsid1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[1,1] <- length(tmp)
    result[1,2] <- mean(val1$P<0.05)
    result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1) #1snp
  {
    result[1,3]=EA_Bonnsummarydat$P[which(EA_Bonnsummarydat$SNP==tmp)]
  }
  
  rownames(result)[1]="EA_Bonn"  
  result[1,]
  
  ### validation in Bonn EA/BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  tmp=intersect(tmp,refsnps)
  if (length(tmp)>1)
  {
    valsnps1 <- selsnps[which(rsid %in% tmp)]
    rsid1 <- rsid[which(rsid %in% tmp)]
    val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydat$SNP %in% tmp,]
    val1 <- val1[match(rsid1,val1$SNP),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z1[,match(rsid1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[2,1] <- length(tmp)
    result[2,2] <- mean(val1$P<0.05)
    result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[2,3]=BEEA_Bonnsummarydat$P[which(BEEA_Bonnsummarydat$SNP==tmp)]
  }
  
  rownames(result)[2]="BEEA_Bonn"
  result[2,]
  ### validation in Bonn BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  tmp=intersect(tmp,refsnps)
  if (length(tmp)>1)
  {
    valsnps1 <- selsnps[which(rsid %in% tmp)]
    rsid1 <- rsid[which(rsid %in% tmp)]
    val1 <- BE_Bonnsummarydat[BE_Bonnsummarydat$SNP %in% tmp,]
    val1 <- val1[match(rsid1,val1$SNP),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z1[,match(rsid1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[3,1] <- length(tmp)
    result[3,2] <- mean(val1$P<0.05)
    result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[3,3]=BE_Bonnsummarydat$P[which(BE_Bonnsummarydat$SNP==tmp)]
  }
  
  rownames(result)[3]="BE_Bonn"
  result[3,]
  ### validation in Oxford BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  tmp=intersect(tmp,refsnps)
  if (length(tmp)>1)
  {
    valsnps2 <- selsnps[which(rsid %in% tmp)]
    rsid1 <- rsid[which(rsid %in% tmp)]
    val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydat$SNP %in% tmp,]
    val2 <- val2[match(rsid1,val2$SNP),]
    uu <- val2$BETA
    vv <- as.numeric(val2$se)
    ZZ <- Z1[,match(rsid1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    
    result[4,1] <- length(tmp)
    result[4,2] <- mean(val2$P<0.05)
    result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    #sum(lamb*(1-pchisq(Q,df=1)))
    result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[4,3]=BE_Oxfordsummarydat$P[which(BE_Oxfordsummarydat$SNP==tmp)]
  }
  
  rownames(result)[4]="BE_Oxford"
  result[4,]
  #### combine the two validation datasets ####
  if (!is.null(val1) & !is.null(val2))
  {
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
      ZZ <- Z1[,colnames(Z1) %in% val1$SNP]
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
    }
  }
  rownames(result)[5]="BE_Combined"
  return(result)  
}

#change validation part, use position to match snps, it can include a few more snps
EA_Bonnsummarydatstr=paste0(EA_Bonnsummarydat$chr,":",EA_Bonnsummarydat$position)
BE_Bonnsummarydatstr=paste0(BE_Bonnsummarydat$chr,":",BE_Bonnsummarydat$position)
BEEA_Bonnsummarydatstr=paste0(BEEA_Bonnsummarydat$chr,":",BEEA_Bonnsummarydat$position)
BE_Oxfordsummarydatstr=paste0(BE_Oxfordsummarydat$chr,":",BE_Oxfordsummarydat$position)
#Load GTEx saved model to speed up. 
TWAS_SKAT_gene1 <- function(genename,res_min,cvmin=NULL,phenotype=NULL,phenotypepos=NULL,snp=NULL,snppos=NULL,bcagenotype=NULL,covariate=NULL,verbose=0) {
  
  #idx2=which(codinggenetable$Symbol==genename)
  #if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  library(SKAT)
  # library(glmnet)
  # library(GenomicRanges)
  # snppos$chr[snppos$chr==23]="X"
  # phenotypepos$chr[phenotypepos$chr==23]="X"
  # gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  # gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  # 
  # result <- data.frame(matrix(NA,8,4))  
  # if (verbose>0) print("get the gene model---")
  # i <- which(row.names(phenotype)==genename)
  # ncv=10
  # distcutoff = 5e5
  # Y=unlist(phenotype[i,]) #geneexp
  # r2=NA
  # glmflag=0 #if glm selected variables
  # tmp=distance(gr_snp,gr_pos[i])
  # idx=which(tmp<distcutoff)
  # tmp=rowSums(data.matrix(snp[idx,]))
  # idx=idx[tmp!=0] #remove all 0 genotypes
  # numvar=0 #number of snp selected by glmnet
  # selectedsnps=NA
  # selectedsnps_coeff=NA
  # p_gender=NA
  # p_age=NA
  # tmp=quantile(Y,probs=c(0.15,0.85))
  # if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  # 
  # X1=t(snp[idx,])
  # dim(X1)
  # ucor <- matrix(0,ncol(X1),2)
  # for (l in 1:ncol(X1)){
  #   ucov<- data.matrix(cbind(X1[,l],covariate))
  #   ufit <- lm(Y~ucov)
  #   ucor[l,1] <- summary(ufit)$coef[2,4]
  #   ucor[l,2] <- cor(Y,X1[,l])
  # }
  # #hist(ucor[,1])
  # #X<- X[,ucor[,1]<0.2]
  # #pcor<- ucor[ucor[,1]<0.2,1]
  # pcor <- ucor[,1]
  # X1 <- X1[,order(pcor,decreasing=T)]
  # X <- removehighcorr(X1,0.9)
  # #hist(ucor[colnames(X1)%in% colnames(X),1])
  # #hist(ucor[colnames(X1)%in% colnames(X),2])
  # dim(X)
  # Xall=data.matrix(cbind(X,covariate))
  # covariateall=covariate
  # 
  # penalty=rep(1,ncol(Xall))
  # penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  # set.seed(i+10000)
  # cvfit=tryCatch(
  #   {
  #     cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #   },
  #   error=function(e)
  #   {
  #     return(F)
  #   }
  # )
  # #plot(cvfit)
  # 
  # max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  # cvfit$lambda.1se
  # cvfit$lambda.min
  # 
  # 
  # cverr <- matrix(NA,length(cvfit$lambda),100)
  # cvse  <- matrix(NA,length(cvfit$lambda),100)
  # rownames(cverr)=cvfit$lambda
  # for (l in 1:ncol(cverr)) {
  #   set.seed(l+100)
  #   fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #   alllambda=intersect(cvfit$lambda,fit$lambda)
  #   idx1=match(alllambda,cvfit$lambda)
  #   idx2=match(alllambda,fit$lambda)
  #   cverr[idx1,l] <- fit$cvm[idx2]
  #   cvse[idx1,l] <- fit$cvsd[idx2]
  # }
  # 
  # merr <- apply(cverr,1,mean,na.rm=T)
  # mcvse <- sqrt(apply(cvse^2,1,mean))
  # 
  # lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  # lambda.min <- cvfit$lambda[which.min(merr)]
  # 
  # #plot(log(cvfit$lambda),merr)
  # #abline(v=log(lambda.1se),lty=2)
  # #abline(v=log(lambda.min),lty=3)
  # #plot(cvfit$lambda,merr)
  # #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  # 
  # #lambda.min <- cvfit$lambda.min
  # if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  # 
  # 
  # sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  # selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1,drop=F]
  # selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  # 
  idx1=which(rownames(res_min)==genename)
  selsnps=unlist(strsplit(res_min$selectedsnps[idx1],"|",fixed=T))
  ## check the marginal eQTL association ##
  # outp <- rep(0,length(selsnps))
  # for (i in 1:length(selsnps)) {
  #   Xmat <- data.matrix(cbind(X[,colnames(X)==selsnps[i]],covariate))
  #   fit <- lm(Y~Xmat)
  #   outp[i] <- summary(fit)$coef[2,4]
  # }  
  #outp
  
  mean(selsnps %in% row.names(bcagenotype))
  ## here we should have a function to get the genotypes for SNPs being selected in selsnps ##
  #bcagenotype=extractBCAgenotype(selsnps = selsnps) #takes 2 minutes for 38 snps
  
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
  #selcoeff <- selcoeff[!is.na(idx)]
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
  
  #if (verbose>0) print("compute skat BE result")
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[6,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[6,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[6,4] <- out$p.value
  rownames(result)[6]="BE_SKAT"
  #result[6,]
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- as.data.frame(cbind(Z1[,i],Covariate1))
  #   outfit <- glm(y~.,data=Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  #if (verbose>0) print("compute skat EA result")
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
  # out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[7,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[7,4] <- out$p.value
  rownames(result)[7]="EA_SKAT"
  #result[7,]
  #for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  #if (verbose>0) print("compute skat BEEA result")
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
  # out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[8,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[8,4] <- out$p.value
  rownames(result)[8]="BEEA_SKAT"
  #result[8,]
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  # #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
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
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )

      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
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
  # selsnpinfo=tmp
  selsnpstr=rep(length(selsnps)) #selsnps are in hg38, selsnpstr are in hg19
  for (j in 1:length(selsnps))
  {
    selsnpstr[j]=unlist(strsplit(tmp$snphg19[j],"_"))[1]
  }
  
  library(survey)
  
  ### validation in Bonn EA ###
  
  val1=val2=NULL
  tmp=intersect(selsnpstr,EA_Bonnsummarydatstr)
  tmp1=sum(EA_Bonnsummarydatstr %in% selsnpstr)
  if (tmp1>length(tmp)) warning(paste0(genename,":there are multiple snps in Bonn for a snp, check it"))
  idx1=match(tmp,EA_Bonnsummarydatstr)
  chr=EA_Bonnsummarydat$chr[idx1][1]
  #read data from 1000 genome v3 EUR
  thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
  refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
  #refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
  refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
  refbimstr=paste0(refbim$V1,":",refbim$V4)
  tmp=intersect(selsnpstr,refbimstr)
  idx1=match(tmp,refbimstr)
  Z1=as.data.frame(refraw[idx1,7:ncol(refraw)])
  rownames(Z1)=tmp
  Z1=as.data.frame(t(Z1))
  # tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps) #snps should be available in 1000G ref as well
  tmp1=intersect(selsnpstr,EA_Bonnsummarydatstr)
  tmp=intersect(tmp1,refbimstr)
  if (length(tmp)<length(tmp1)) print(paste0(genename,",EA:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    valsnps1 <- selsnps[which(selsnpstr %in% tmp)]
    selsnpstr1 <- selsnpstr[which(selsnpstr %in% tmp)]
    val1 <- EA_Bonnsummarydat[EA_Bonnsummarydatstr %in% tmp,]
    var1str=paste0(val1$chr,":",val1$position)
    val1 <- val1[match(selsnpstr1,var1str),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z1[,match(selsnpstr1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[1,1] <- length(tmp)
    result[1,2] <- mean(val1$P<0.05)
    result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1) #1snp
  {
    result[1,3]=EA_Bonnsummarydat$P[which(EA_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[1]="EA_Bonn"  
  result[1,]
  
  ### validation in Bonn EA/BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  tmp1=intersect(selsnpstr,BEEA_Bonnsummarydatstr)
  tmp=intersect(tmp1,refbimstr)
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",BEEA:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    # valsnps1 <- selsnps[which(rsid %in% tmp)]
    # rsid1 <- rsid[which(rsid %in% tmp)]
    # val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydat$SNP %in% tmp,]
    # val1 <- val1[match(rsid1,val1$SNP),]
    valsnps1 <- selsnps[which(selsnpstr %in% tmp)]
    selsnpstr1 <- selsnpstr[which(selsnpstr %in% tmp)]
    val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydatstr %in% tmp,]
    var1str=paste0(val1$chr,":",val1$position)
    val1 <- val1[match(selsnpstr1,var1str),]
    
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    #ZZ <- Z1[,match(rsid1,colnames(Z1))]
    ZZ <- Z1[,match(selsnpstr1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[2,1] <- length(tmp)
    result[2,2] <- mean(val1$P<0.05)
    result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[2,3]=BEEA_Bonnsummarydat$P[which(BEEA_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[2]="BEEA_Bonn"
  result[2,]
  ### validation in Bonn BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  tmp1=intersect(selsnpstr,BE_Bonnsummarydatstr)
  tmp=intersect(tmp1,refbimstr)
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",BE:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  
  if (length(tmp)>1)
  {
    # valsnps1 <- selsnps[which(rsid %in% tmp)]
    # rsid1 <- rsid[which(rsid %in% tmp)]
    # val1 <- BE_Bonnsummarydat[BE_Bonnsummarydat$SNP %in% tmp,]
    # val1 <- val1[match(rsid1,val1$SNP),]
    valsnps1 <- selsnps[which(selsnpstr %in% tmp)]
    selsnpstr1 <- selsnpstr[which(selsnpstr %in% tmp)]
    val1 <- BE_Bonnsummarydat[BE_Bonnsummarydatstr %in% tmp,]
    var1str=paste0(val1$chr,":",val1$position)
    val1 <- val1[match(selsnpstr1,var1str),]
    
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    #ZZ <- Z1[,match(rsid1,colnames(Z1))]
    ZZ <- Z1[,match(selsnpstr1,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[3,1] <- length(tmp)
    result[3,2] <- mean(val1$P<0.05)
    result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[3,3]=BE_Bonnsummarydat$P[which(BE_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[3]="BE_Bonn"
  result[3,]
  ### validation in Oxford BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  tmp1=intersect(selsnpstr,BE_Oxfordsummarydatstr)
  tmp=intersect(tmp1,refbimstr)
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",Oxford:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    # valsnps2 <- selsnps[which(rsid %in% tmp)]
    # rsid1 <- rsid[which(rsid %in% tmp)]
    # val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydat$SNP %in% tmp,]
    # val2 <- val2[match(rsid1,val2$SNP),]
    valsnps2 <- selsnps[which(selsnpstr %in% tmp)]
    selsnpstr2 <- selsnpstr[which(selsnpstr %in% tmp)]
    val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydatstr %in% tmp,]
    var2str=paste0(val2$chr,":",val2$position)
    val2 <- val2[match(selsnpstr2,var2str),]
    uu <- val2$BETA
    vv <- as.numeric(val2$se)
    #ZZ <- Z1[,match(rsid1,colnames(Z1))]
    ZZ <- Z1[,match(selsnpstr2,colnames(Z1))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    
    result[4,1] <- length(tmp)
    result[4,2] <- mean(val2$P<0.05)
    result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    #sum(lamb*(1-pchisq(Q,df=1)))
    result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[4,3]=BE_Oxfordsummarydat$P[which(BE_Oxfordsummarydatstr==tmp)]
  }
  
  rownames(result)[4]="BE_Oxford"
  result[4,]
  #### combine the two validation datasets ####
  if (!is.null(val1) & !is.null(val2))
  {
    val1 <- val1[val1$SNP %in% val2$SNP,]
    val1str=paste0(val1$chr,":",val1$position)
    valsnps1 <- valsnps1[valsnps1%in%valsnps2]
    
    val2 <- val2[val2$SNP %in% val1$SNP,]
    valsnps2 <- valsnps2[valsnps2%in%valsnps1]
    if (nrow(val1)>1)  {
      uu1 <- val1$BETA
      vv1 <- as.numeric(val1$SE)
      uu2 <- val2$BETA
      vv2 <- as.numeric(val2$se)
      
      uu <- uu1+uu2
      ZZ <- Z1[,colnames(Z1) %in% val1str]
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
    }
  }
  rownames(result)[5]="BE_Combined"
  return(result)  
}

EA_Bonnsummarydatstr1=paste0(EA_Bonnsummarydat$chr,":",EA_Bonnsummarydat$position,"_",EA_Bonnsummarydat$effect_allele,"_",EA_Bonnsummarydat$non_effect_allele)
BE_Bonnsummarydatstr1=paste0(BE_Bonnsummarydat$chr,":",BE_Bonnsummarydat$position,"_",BE_Bonnsummarydat$effect_allele,"_",BE_Bonnsummarydat$non_effect_allele)
BEEA_Bonnsummarydatstr1=paste0(BEEA_Bonnsummarydat$chr,":",BEEA_Bonnsummarydat$position,"_",BEEA_Bonnsummarydat$effect_allele,"_",BEEA_Bonnsummarydat$non_effect_allele)
BE_Oxfordsummarydatstr1=paste0(BE_Oxfordsummarydat$chr,":",BE_Oxfordsummarydat$position,"_",BE_Oxfordsummarydat$effect_allele,"_",BE_Oxfordsummarydat$non_effect_allele)

# EA_Bonnsummarydatstr2=paste0(EA_Bonnsummarydat$chr,":",EA_Bonnsummarydat$position,"_",EA_Bonnsummarydat$non_effect_allele,"_",EA_Bonnsummarydat$effect_allele)
# BE_Bonnsummarydatstr2=paste0(BE_Bonnsummarydat$chr,":",BE_Bonnsummarydat$position,"_",BE_Bonnsummarydat$non_effect_allele,"_",BE_Bonnsummarydat$effect_allele)
# BEEA_Bonnsummarydatstr2=paste0(BEEA_Bonnsummarydat$chr,":",BEEA_Bonnsummarydat$position,"_",BEEA_Bonnsummarydat$non_effect_allele,"_",BEEA_Bonnsummarydat$effect_allele)
# BE_Oxfordsummarydatstr2=paste0(BE_Oxfordsummarydat$chr,":",BE_Oxfordsummarydat$position,"_",BE_Oxfordsummarydat$non_effect_allele,"_",BE_Oxfordsummarydat$effect_allele)

#allelle allignment summarydat and refgenotype
TWAS_SKAT_gene2 <- function(genename,res_min,cvmin=NULL,phenotype=NULL,phenotypepos=NULL,snp=NULL,snppos=NULL,bcagenotype=NULL,covariate=NULL,verbose=0) {
  
  #idx2=which(codinggenetable$Symbol==genename)
  #if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  library(SKAT)
  # library(glmnet)
  # library(GenomicRanges)
  # snppos$chr[snppos$chr==23]="X"
  # phenotypepos$chr[phenotypepos$chr==23]="X"
  # gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  # gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  # 
   result <- data.frame(matrix(NA,8,4))  
  # if (verbose>0) print("get the gene model---")
  # i <- which(row.names(phenotype)==genename)
  # ncv=10
  # distcutoff = 5e5
  # Y=unlist(phenotype[i,]) #geneexp
  # r2=NA
  # glmflag=0 #if glm selected variables
  # tmp=distance(gr_snp,gr_pos[i])
  # idx=which(tmp<distcutoff)
  # tmp=rowSums(data.matrix(snp[idx,]))
  # idx=idx[tmp!=0] #remove all 0 genotypes
  # numvar=0 #number of snp selected by glmnet
  # selectedsnps=NA
  # selectedsnps_coeff=NA
  # p_gender=NA
  # p_age=NA
  # tmp=quantile(Y,probs=c(0.15,0.85))
  # if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  # 
  # X1=t(snp[idx,])
  # dim(X1)
  # ucor <- matrix(0,ncol(X1),2)
  # for (l in 1:ncol(X1)){
  #   ucov<- data.matrix(cbind(X1[,l],covariate))
  #   ufit <- lm(Y~ucov)
  #   ucor[l,1] <- summary(ufit)$coef[2,4]
  #   ucor[l,2] <- cor(Y,X1[,l])
  # }
  # #hist(ucor[,1])
  # #X<- X[,ucor[,1]<0.2]
  # #pcor<- ucor[ucor[,1]<0.2,1]
  # pcor <- ucor[,1]
  # X1 <- X1[,order(pcor,decreasing=T)]
  # X <- removehighcorr(X1,0.9)
  # #hist(ucor[colnames(X1)%in% colnames(X),1])
  # #hist(ucor[colnames(X1)%in% colnames(X),2])
  # dim(X)
  # Xall=data.matrix(cbind(X,covariate))
  # covariateall=covariate
  # 
  # penalty=rep(1,ncol(Xall))
  # penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  # set.seed(i+10000)
  # cvfit=tryCatch(
  #   {
  #     cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #   },
  #   error=function(e)
  #   {
  #     return(F)
  #   }
  # )
  # #plot(cvfit)
  # 
  # max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  # cvfit$lambda.1se
  # cvfit$lambda.min
  # 
  # 
  # cverr <- matrix(NA,length(cvfit$lambda),100)
  # cvse  <- matrix(NA,length(cvfit$lambda),100)
  # rownames(cverr)=cvfit$lambda
  # for (l in 1:ncol(cverr)) {
  #   set.seed(l+100)
  #   fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #   alllambda=intersect(cvfit$lambda,fit$lambda)
  #   idx1=match(alllambda,cvfit$lambda)
  #   idx2=match(alllambda,fit$lambda)
  #   cverr[idx1,l] <- fit$cvm[idx2]
  #   cvse[idx1,l] <- fit$cvsd[idx2]
  # }
  # 
  # merr <- apply(cverr,1,mean,na.rm=T)
  # mcvse <- sqrt(apply(cvse^2,1,mean))
  # 
  # lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  # lambda.min <- cvfit$lambda[which.min(merr)]
  # 
  # #plot(log(cvfit$lambda),merr)
  # #abline(v=log(lambda.1se),lty=2)
  # #abline(v=log(lambda.min),lty=3)
  # #plot(cvfit$lambda,merr)
  # #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  # 
  # #lambda.min <- cvfit$lambda.min
  # if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  # 
  # 
  # sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  # selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1,drop=F]
  # selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  # 
  idx1=which(rownames(res_min)==genename)
  selsnps=unlist(strsplit(res_min$selectedsnps[idx1],"|",fixed=T))
  ## check the marginal eQTL association ##
  # outp <- rep(0,length(selsnps))
  # for (i in 1:length(selsnps)) {
  #   Xmat <- data.matrix(cbind(X[,colnames(X)==selsnps[i]],covariate))
  #   fit <- lm(Y~Xmat)
  #   outp[i] <- summary(fit)$coef[2,4]
  # }  
  #outp
  
  mean(selsnps %in% row.names(bcagenotype))
  ## here we should have a function to get the genotypes for SNPs being selected in selsnps ##
  #bcagenotype=extractBCAgenotype(selsnps = selsnps) #takes 2 minutes for 38 snps
  
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
  #selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  # if (length(correctedsnps)>0) #flip snps
  # {
  #   idxtocorrect=match(correctedsnps,colnames(Z))
  #   Z[,idxtocorrect]=2-Z[,idxtocorrect]
  # }
  #colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
  
  #rownames(Z)=colnames(bcagenotype)
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  #if (verbose>0) print("compute skat BE result")
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[6,2] <- out$p.value
  #out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[6,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[6,4] <- out$p.value
  rownames(result)[6]="BE_SKAT"
  #result[6,]
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- as.data.frame(cbind(Z1[,i],Covariate1))
  #   outfit <- glm(y~.,data=Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  #if (verbose>0) print("compute skat EA result")
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
  # out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[7,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[7,4] <- out$p.value
  rownames(result)[7]="EA_SKAT"
  #result[7,]
  #for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  #if (verbose>0) print("compute skat BEEA result")
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
  # out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  # result[8,3] <- out$p.value
  # out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  # result[8,4] <- out$p.value
  rownames(result)[8]="BEEA_SKAT"
  #result[8,]
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  #outp1
  
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  # #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
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
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
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
  # selsnpinfo=tmp
  selsnpstr0=rep(nrow(tmp)) #selsnps are in hg38, selsnpstr are in hg19, only use position
  for (j in 1:nrow(tmp))
  {
    selsnpstr0[j]=unlist(strsplit(tmp$snphg19[j],"_"))[1]
  }
  selsnpmap=tmp
  selsnpstr=tmp$snphg19 #should use id in summary data,use position+allelells
  library(survey)
  
  ### validation in Bonn EA ###
  
  chr=unlist(strsplit(selsnps[1],":"))[1]
  #read data from 1000 genome v3 EUR
  thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
  refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
  #refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
  refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
  refbimstr1=paste0(refbim$V1,":",refbim$V4,"_",refbim$V5,"_",refbim$V6)
  refbimstr2=paste0(refbim$V1,":",refbim$V4,"_",refbim$V6,"_",refbim$V5)
  
  val1=val2=NULL
  tmp=intersect(selsnpstr0,EA_Bonnsummarydatstr)
  idx1=match(tmp,EA_Bonnsummarydatstr) #use id from summary data
  BonnEAsnpstr=EA_Bonnsummarydatstr1[idx1]
  
  Z1=Z2=NULL
  tmp1=intersect(BonnEAsnpstr,refbimstr1)
  tmp2=intersect(BonnEAsnpstr,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  # tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps) #snps should be available in 1000G ref as well
  tmp1=intersect(selsnpstr0,EA_Bonnsummarydatstr)
  tmp=colnames(Z3)
  if (length(tmp)<length(tmp1)) print(paste0(genename,",EA:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    valsnps1 <- colnames(Z3)
    selsnpstr1 <- colnames(Z3)
    val1 <- EA_Bonnsummarydat[match(tmp,EA_Bonnsummarydatstr1),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z3[,match(selsnpstr1,colnames(Z3))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[1,1] <- length(tmp)
    result[1,2] <- mean(val1$P<0.05)
    result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1) #1snp
  {
    result[1,3]=EA_Bonnsummarydat$P[which(EA_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[1]="EA_Bonn"  
  #result[1,]
  
  ### validation in Bonn EA/BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  # tmp1=intersect(selsnpstr,BEEA_Bonnsummarydatstr)
  # tmp=intersect(tmp1,refbimstr)
  
  tmp=intersect(selsnpstr0,BEEA_Bonnsummarydatstr)
  idx1=match(tmp,BEEA_Bonnsummarydatstr) #use id from summary data
  BonnBEEAsnpstr=BEEA_Bonnsummarydatstr1[idx1]
  
  Z1=Z2=NULL
  tmp1=intersect(BonnBEEAsnpstr,refbimstr1)
  tmp2=intersect(BonnBEEAsnpstr,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  tmp1=intersect(selsnpstr0,BEEA_Bonnsummarydatstr)
  tmp=colnames(Z3)
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",BEEA:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    valsnps1 <- colnames(Z3)
    selsnpstr1 <- colnames(Z3)
    val1 <- BEEA_Bonnsummarydat[match(tmp,BEEA_Bonnsummarydatstr1),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z3[,match(selsnpstr1,colnames(Z3))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[2,1] <- length(tmp)
    result[2,2] <- mean(val1$P<0.05)
    result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[2,3]=BEEA_Bonnsummarydat$P[which(BEEA_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[2]="BEEA_Bonn"
  result[2,]
  ### validation in Bonn BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  # tmp1=intersect(selsnpstr,BE_Bonnsummarydatstr)
  # tmp=intersect(tmp1,refbimstr)
  tmp=intersect(selsnpstr0,BE_Bonnsummarydatstr)

  idx1=match(tmp,BE_Bonnsummarydatstr) #use id from summary data
  BonnBEsnpstr=BE_Bonnsummarydatstr1[idx1]
  
  Z1=Z2=NULL
  tmp1=intersect(BonnBEsnpstr,refbimstr1)
  tmp2=intersect(BonnBEsnpstr,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  tmp1=intersect(selsnpstr0,BE_Bonnsummarydatstr)
  tmp=colnames(Z3)
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",BE:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  
  if (length(tmp)>1)
  {
    valsnps1 <- colnames(Z3)
    selsnpstr1 <- colnames(Z3)
    val1 <- BE_Bonnsummarydat[match(tmp,BE_Bonnsummarydatstr1),]
    uu <- val1$BETA
    vv <- as.numeric(val1$SE)
    ZZ <- Z3[,match(selsnpstr1,colnames(Z3))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    result[3,1] <- length(tmp)
    result[3,2] <- mean(val1$P<0.05)
    result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[3,3]=BE_Bonnsummarydat$P[which(BE_Bonnsummarydatstr==tmp)]
  }
  
  rownames(result)[3]="BE_Bonn"
  #result[3,]
  ### validation in Oxford BE ###
  
  # tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  # tmp=intersect(tmp,refsnps)
  # tmp1=intersect(selsnpstr,BE_Oxfordsummarydatstr)
  # tmp=intersect(tmp1,refbimstr)
  tmp=intersect(selsnpstr0,BE_Oxfordsummarydatstr)
  idx1=match(tmp,BE_Oxfordsummarydatstr) #use id from summary data
  OxfordBEsnpstr=BE_Oxfordsummarydatstr1[idx1]
  
  Z1=Z2=NULL
  tmp1=intersect(OxfordBEsnpstr,refbimstr1)
  tmp2=intersect(OxfordBEsnpstr,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  tmp1=intersect(selsnpstr0,BE_Oxfordsummarydatstr)
  tmp=colnames(Z3)
  
  if (length(tmp)<length(tmp1)) warning(paste0(genename,",Oxford:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))
  if (length(tmp)>1)
  {
    # valsnps2 <- selsnps[which(selsnpstr %in% tmp)]
    # selsnpstr2 <- selsnpstr[which(selsnpstr %in% tmp)]
    # val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydatstr %in% tmp,]
    # var2str=paste0(val2$chr,":",val2$position)
    # val2 <- val2[match(selsnpstr2,var2str),]
    # uu <- val2$BETA
    # vv <- as.numeric(val2$se)
    # #ZZ <- Z1[,match(rsid1,colnames(Z1))]
    # ZZ <- Z1[,match(selsnpstr2,colnames(Z1))]
    valsnps2 <- colnames(Z3)
    selsnpstr2 <- colnames(Z3)
    val2 <- BE_Oxfordsummarydat[match(tmp,BE_Oxfordsummarydatstr1),]
    uu <- val2$BETA
    vv <- as.numeric(val2$se)
    ZZ <- Z3[,match(selsnpstr2,colnames(Z3))]
    rr <- cor(ZZ,use="pairwise.complete.obs")
    VV <- diag(vv) %*% rr %*% diag(vv)
    lamb <- eigen(VV)$values
    Q <- drop(t(uu) %*% uu)
    
    result[4,1] <- length(tmp)
    result[4,2] <- mean(val2$P<0.05)
    result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    #sum(lamb*(1-pchisq(Q,df=1)))
    result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  }else if (length(tmp)==1)
  {
    result[4,3]=BE_Oxfordsummarydat$P[which(BE_Oxfordsummarydatstr==tmp)]
  }
  
  rownames(result)[4]="BE_Oxford"
  result[4,]
  #### combine the two validation datasets ####
  # if (!is.null(val1) & !is.null(val2))
  # {
  #   val1 <- val1[val1$SNP %in% val2$SNP,]
  #   val1str=paste0(val1$chr,":",val1$position)
  #   valsnps1 <- valsnps1[valsnps1%in%valsnps2]
  #   
  #   val2 <- val2[val2$SNP %in% val1$SNP,]
  #   valsnps2 <- valsnps2[valsnps2%in%valsnps1]
  #   if (nrow(val1)>1)  {
  #     uu1 <- val1$BETA
  #     vv1 <- as.numeric(val1$SE)
  #     uu2 <- val2$BETA
  #     vv2 <- as.numeric(val2$se)
  #     
  #     uu <- uu1+uu2
  #     ZZ <- Z1[,colnames(Z1) %in% val1str]
  #     rr <- cor(ZZ,use="pairwise.complete.obs")
  #     VV1 <- diag(vv1) %*% rr %*% diag(vv1)
  #     VV2 <- diag(vv2) %*% rr %*% diag(vv2)
  #     VV <- VV1 + VV2
  #     lamb <- eigen(VV)$values
  #     Q <- drop(t(uu) %*% uu)
  #     result[5,1] <- nrow(val1)
  #     result[5,2] <- mean(val1$P<0.05 & val2$P <0.05)
  #     result[5,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  #     result[5,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  #   }
  # }
  rownames(result)[5]="BE_Combined"
  return(result)  
}

extractBCAgenotype=function(selsnps=NULL)
{
  allgenotypedat=NULL
  dataset="merge_beacon_cambridge_hrc_maf001_snp" #for HRC maf001
  chr=unlist(strsplit(selsnps[1],":"))[1]
  genotypefile=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/",dataset,"/","chr",chr,"_filter_hg19tohg38_flip.traw")
  bimfile=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/",dataset,"/","chr",chr,"_filter_hg19tohg38_flip.bim")
  bim=as.data.frame(data.table::fread(bimfile))
  columnname=as.data.frame(data.table::fread(genotypefile,nrows=1))
  #sample name
  samplename=colnames(columnname)[7:ncol(columnname)]
  samplename=unlist(strsplit(samplename,"_"))
  samplename=samplename[seq(2,length(samplename),2)] 
  samplename=gsub("SEP","",samplename)
  for (i in 1:length(selsnps))
  {
    mypos=unlist(strsplit(selsnps[i],"_"))[1]
    mypos=unlist(strsplit(mypos,":"))[2]
    #where is the snp
    nline=which(bim$V4==mypos)
    
    if (length(nline)>0) #read data
    {
      genotypedat=as.data.frame(data.table::fread(genotypefile,nrows=length(nline),skip=nline))
      rownames(genotypedat)=paste0(genotypedat[,1],":",genotypedat[,4],"_",genotypedat[,5],"_",genotypedat[,6])
      genotypedat=genotypedat[,7:ncol(genotypedat)]
      colnames(genotypedat)=samplename
      allgenotypedat=rbind(allgenotypedat,genotypedat)
    }
  }
  
  return(allgenotypedat)
}

#check a list of genes
genelist=read.csv("../result/SKAT_TWAS_newpc6_genes_allcovar_r2_01_.csv")
genelist=read.csv("../result/SKAT_TWAS_newpc6_genes_allcovar_heritability.csv")

genelist=genelist[genelist$cutoff=="FWER",]
length(unique(genelist$gene)) 
genelist=genelist[genelist$cutoff=="FDR",]
length(unique(genelist$gene))
#merge the same gene in the same tissue
tmp=unique(paste0(genelist$gene,"_",genelist$tissue))
genelist1=NULL
for (i in 1:length(tmp))
{
  genename=unlist(strsplit(tmp[i],"_"))[1]
  organ=unlist(strsplit(tmp[i],"_"))[2]
  idx=which(genelist$gene==genename & genelist$tissue==organ)
  tmp1=genelist[idx[1],,drop=F]
  if (length(idx)>0)
  {
    tmp1$type=paste0(unique(c(genelist$type[idx])),collapse=";")
    tmp1$cutoff=paste0(unique(c(genelist$cutoff[idx])),collapse=";")
  }
  genelist1=rbind(genelist1,tmp1)
}
genelist=genelist1
nrow(genelist) #34
length(unique(genelist$gene)) #21
genelist$chr=gsub("chr","",genelist$chr)
genelist$chr=as.numeric(genelist$chr)
genelist=genelist[order(genelist$chr,genelist$position),]
#genelist$H_pvalue=genelist$H=NA
genelist$skatpBEEA=genelist$skatpEA=genelist$skatpBE=genelist$skat_numsnp=NA
genelist$Bonn_BEpsaddlepoint=genelist$Bonn_BEpsatterthwaite=genelist$Bonn_BEp005=genelist$Bonn_numsnp=NA
genelist$Bonn_EApsaddlepoint=genelist$Bonn_EApsatterthwaite=genelist$Bonn_EAp005=NA
genelist$Bonn_BEEApsaddlepoint=genelist$Bonn_BEEApsatterthwaite=genelist$Bonn_BEEAp005=NA
genelist$Oxford_BEpsaddlepoint=genelist$Oxford_BEpsatterthwaite=genelist$Oxford_BEp005=genelist$Oxford_numsnp=NA
genelist$Comb_BEpsaddlepoint=genelist$Comb_BEpsatterthwaite=genelist$Comb_BEp005=genelist$Comb_numsnp=NA
organs=unique(genelist$tissue)
set.seed(1000)
for (ii in 1:length(organs))
{
  organ=organs[ii]
  print(paste0(organ,"--------------------------------"))
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  #no need to load bcagenotype data
  load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
  load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
  load(paste0(outfolder,"/skat_res.RData")) #saved skat res, skat_min2_pc6,skat_min2
  load(paste0("../result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")) #snp,phenotype...
  outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  
  idxs=which(genelist$tissue==organ)
  for (jj in idxs)
  {
    genename=genelist$gene[jj]
    #if (is.na(genelist$skat_numsnp[jj]))
    {
      print(genename)
      myres=TWAS_SKAT_gene(genename=genename,cvmin=1,phenotype=phenotype,phenotypepos=phenotypepos,snp=snp,snppos=snppos,bcagenotype=bcagenotype,covariate=covariate,verbose=1)
      idx=which(rownames(heritability)==genename)
      if (length(idx)>0)
      {
         genelist$H[jj]=heritability$V1[idx]
         genelist$H_pvalue[jj]=heritability$V3[idx]
      }
      genelist$skat_numsnp[jj]=myres[6,1]
      genelist$skatpBE[jj]=myres[6,2]
      genelist$skatpEA[jj]=myres[7,2]
      genelist$skatpBEEA[jj]=myres[8,2]
      genelist$Bonn_numsnp[jj]=myres[1,1]
      genelist$Bonn_BEp005[jj]=myres[3,2]
      genelist$Bonn_BEpsatterthwaite[jj]=myres[3,3]
      genelist$Bonn_BEpsaddlepoint[jj]=myres[3,4]
      genelist$Bonn_EAp005[jj]=myres[1,2]
      genelist$Bonn_EApsatterthwaite[jj]=myres[1,3]
      genelist$Bonn_EApsaddlepoint[jj]=myres[1,4]
      genelist$Bonn_BEEAp005[jj]=myres[2,2]
      genelist$Bonn_BEEApsatterthwaite[jj]=myres[2,3]
      genelist$Bonn_BEEApsaddlepoint[jj]=myres[2,4]
      genelist$Oxford_numsnp[jj]=myres[4,1]
      genelist$Oxford_BEp005[jj]=myres[4,2]
      genelist$Oxford_BEpsatterthwaite[jj]=myres[4,3]
      genelist$Oxford_BEpsaddlepoint[jj]=myres[4,4]
      genelist$Comb_numsnp[jj]=myres[5,1]
      genelist$Comb_BEp005[jj]=myres[5,2]
      genelist$Comb_BEpsatterthwaite[jj]=myres[5,3]
      genelist$Comb_BEpsaddlepoint[jj]=myres[5,4]
    }
  }
}
genelist$validated=F
idx=which(genelist$Bonn_BEpsaddlepoint<0.05 | genelist$Bonn_BEpsaddlepoint<0.05 | genelist$Bonn_EApsaddlepoint<0.05 | genelist$Bonn_EApsaddlepoint<0.05 |genelist$Bonn_BEEApsaddlepoint<0.05 | genelist$Bonn_BEEApsaddlepoint<0.05
          |genelist$Oxford_BEpsatterthwaite<0.05 | genelist$Oxford_BEpsaddlepoint<0.05)
genelist$validated[idx]=T
idx=order(genelist$validated,decreasing = T)
# write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_r2_01_FWER.csv",quote=F,row.names = F)
# write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_r2_01_FDR.csv",quote=F,row.names = F)
#in validation calculate correlation using 1000 genome v3 EUR
write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_r2_01_FWER_1000G.csv",quote=F,row.names = F)
write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_r2_01_FDR_1000G.csv",quote=F,row.names = F)

# write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_heritability_FWER.csv",quote=F,row.names = F)
# write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_heritability_FDR.csv",quote=F,row.names = F)
write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_heritability_FWER_1000G.csv",quote=F,row.names = F)
write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_heritability_FDR_1000G.csv",quote=F,row.names = F)

#gene_list based on allFDR

#merge the same gene in the same tissue
genelist=read.csv("../result/SKAT_TWAS_newpc6_genes_allcovar_AllFDR05.csv")
set.seed(10000)
tmp=unique(paste0(genelist$gene,"_",genelist$tissue))
genelist1=NULL
for (i in 1:length(tmp))
{
  genename=unlist(strsplit(tmp[i],"_"))[1]
  organ=unlist(strsplit(tmp[i],"_"))[2]
  idx=which(genelist$gene==genename & genelist$tissue==organ)
  tmp1=genelist[idx[1],,drop=F]
  if (length(idx)>0)
  {
    tmp1$type=paste0(unique(c(genelist$type[idx])),collapse=";")
    tmp1$cutoff=paste0(unique(c(genelist$cutoff[idx])),collapse=";")
    tmp1$p=paste0(unique(c(genelist$p[idx])),collapse=";")
    tmp1$fdr=paste0(unique(c(genelist$fdr[idx])),collapse=";")
    tmp1$fwer=paste0(unique(c(genelist$fwer[idx])),collapse=";")
    tmp1$FDR=paste0(unique(c(genelist$FDR[idx])),collapse=";")
    tmp1$FWER=paste0(unique(c(genelist$FWER[idx])),collapse=";")
  }
  genelist1=rbind(genelist1,tmp1)
}
colnames(genelist1)[which(colnames(genelist1)=="FWER")]="AllFWER"
colnames(genelist1)[which(colnames(genelist1)=="FDR")]="AllFDR"
genelist=genelist1
nrow(genelist) #24
length(unique(genelist$gene)) #21
genelist$chr=gsub("chr","",genelist$chr)
genelist$chr=as.numeric(genelist$chr)
genelist=genelist[order(genelist$chr,genelist$position),]
#genelist$H_pvalue=genelist$H=NA
genelist$skatpBEEA=genelist$skatpEA=genelist$skatpBE=genelist$skat_numsnp=NA
genelist$Bonn_BEpsaddlepoint=genelist$Bonn_BEpsatterthwaite=genelist$Bonn_BEp005=genelist$Bonn_numsnp=NA
genelist$Bonn_EApsaddlepoint=genelist$Bonn_EApsatterthwaite=genelist$Bonn_EAp005=NA
genelist$Bonn_BEEApsaddlepoint=genelist$Bonn_BEEApsatterthwaite=genelist$Bonn_BEEAp005=NA
genelist$Oxford_BEpsaddlepoint=genelist$Oxford_BEpsatterthwaite=genelist$Oxford_BEp005=genelist$Oxford_numsnp=NA
genelist$Comb_BEpsaddlepoint=genelist$Comb_BEpsatterthwaite=genelist$Comb_BEp005=genelist$Comb_numsnp=NA
organs=unique(genelist$tissue)
for (ii in 1:length(organs))
{
  organ=organs[ii]
  print(paste0(organ,"--------------------------------"))
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  #no need to load bcagenotype data
  load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
  # load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
  load(paste0(outfolder,"/skat_res.RData")) #saved skat res, skat_min2_pc6,skat_min2
  # load(paste0("../result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")) #snp,phenotype...
  outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  
  idxs=which(genelist$tissue==organ)
  for (jj in idxs)
  {
    genename=genelist$gene[jj]
    #if (is.na(genelist$skat_numsnp[jj]))
    {
      print(genename)
      myres=TWAS_SKAT_gene2(genename=genename,res_min = res_min,bcagenotype = bcagenotype,covariate=covariate,verbose = 1)
      idx=which(rownames(heritability)==genename)
      if (length(idx)>0)
      {
        genelist$H[jj]=heritability$V1[idx]
        genelist$H_pvalue[jj]=heritability$V3[idx]
      }
      genelist$skat_numsnp[jj]=myres[6,1]
      genelist$skatpBE[jj]=myres[6,2]
      genelist$skatpEA[jj]=myres[7,2]
      genelist$skatpBEEA[jj]=myres[8,2]
      genelist$Bonn_numsnp[jj]=myres[1,1]
      genelist$Bonn_BEp005[jj]=myres[3,2]
      genelist$Bonn_BEpsatterthwaite[jj]=myres[3,3]
      genelist$Bonn_BEpsaddlepoint[jj]=myres[3,4]
      genelist$Bonn_EAp005[jj]=myres[1,2]
      genelist$Bonn_EApsatterthwaite[jj]=myres[1,3]
      genelist$Bonn_EApsaddlepoint[jj]=myres[1,4]
      genelist$Bonn_BEEAp005[jj]=myres[2,2]
      genelist$Bonn_BEEApsatterthwaite[jj]=myres[2,3]
      genelist$Bonn_BEEApsaddlepoint[jj]=myres[2,4]
      genelist$Oxford_numsnp[jj]=myres[4,1]
      genelist$Oxford_BEp005[jj]=myres[4,2]
      genelist$Oxford_BEpsatterthwaite[jj]=myres[4,3]
      genelist$Oxford_BEpsaddlepoint[jj]=myres[4,4]
      genelist$Comb_numsnp[jj]=myres[5,1]
      genelist$Comb_BEp005[jj]=myres[5,2]
      genelist$Comb_BEpsatterthwaite[jj]=myres[5,3]
      genelist$Comb_BEpsaddlepoint[jj]=myres[5,4]
    }
  }
}
#SENP6,EA:there are 1 not found in 1000G ref.
#AQP9,EA:there are 1 not found in 1000G ref.
#EXOC3,EA:there are 1 not found in 1000G ref
genelist$validated=F
idx=which(genelist$Bonn_BEpsaddlepoint<0.05 | genelist$Bonn_BEpsaddlepoint<0.05 | genelist$Bonn_EApsaddlepoint<0.05 | genelist$Bonn_EApsaddlepoint<0.05 |genelist$Bonn_BEEApsaddlepoint<0.05 | genelist$Bonn_BEEApsaddlepoint<0.05
          |genelist$Oxford_BEpsatterthwaite<0.05 | genelist$Oxford_BEpsaddlepoint<0.05)
genelist$validated[idx]=T
idx=order(genelist$validated,decreasing = T)
#write.csv(genelist[idx,],file="../result/NEWGenesValidate_allcovar_ALLFDR.csv",quote=F,row.names = F)
#update AllFDR
genelist=read.csv("../result/NEWGenesValidate_allcovar_ALLFDR.csv")
for (i in 1:nrow(genelist))
{
  if (grepl(";",genelist$type[i]))
  {
    tmp=unlist(strsplit(genelist$type[i],";"))
    tmp1=paste0(tmp,"_",genelist$tissue[i],"_",genelist$gene[i])
    tmp2=paste0(allres$type,"_",allres$tissue,"_",allres$gene)
    idx=match(tmp1,tmp2)
    genelist$FDR[i]=paste0(formatC(allres$FDR[idx],format="e",digits = 2),collapse = ";")
    genelist$FWER[i]=paste0(formatC(allres$FWER[idx],format="e",digits = 2),collapse = ";")
  }else
  {
    genelist$FDR[i]=paste0(formatC(genelist$FDR[i],format="e",digits = 2),collapse = ";")
    genelist$FWER[i]=paste0(formatC(genelist$FWER[i],format="e",digits = 2),collapse = ";")
  }
}
write.csv(genelist,file="../result/NEWGenesValidate_allcovar_ALLFDR.csv",quote=F,row.names = F)
#add pvaules/fwer and check if SKAT gene is significant after adjusting for gwas snps
update_genelist=function(genelist=genelist_newpc6,adjustopt="PC6",r2cutoff=0)
{
  #genelist$fwer_addgwas=genelist$pvalue_addgwas=genelist$BEEA_fwer=genelist$BEEA_pvalue=genelist$BEA_fwer=genelist$BEA_pvalue=genelist$EA_fwer=genelist$EA_pvalue=genelist$BE_fwer=genelist$BE_pvalue=NA
  genelist$pvalue_addgwas=NA
  for (i in 1:nrow(genelist))
  {
    if (genelist$tissue[i]=="adipose")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="blood")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="junction")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="mucosa")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="muscularis")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="stomach")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor_allcovar"
    }
    #work on pvalues
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/skat_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    # skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    # skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    # skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    # skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    # skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    # skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    # idx=which(rownames(skat_min2_code)==genelist$gene[i])
    # genelist$BEEA_pvalue[i]=skat_min2_code$BEEA_pvalue[idx]
    # genelist$BEEA_fwer[i]=skat_min2_code$BEEA_fwer[idx]
    # genelist$BEA_pvalue[i]=skat_min2_code$BEA_pvalue[idx]
    # genelist$BEA_fwer[i]=skat_min2_code$BEA_fwer[idx]
    # genelist$EA_pvalue[i]=skat_min2_code$EA_pvalue[idx]
    # genelist$EA_fwer[i]=skat_min2_code$EA_fwer[idx]
    # genelist$BE_pvalue[i]=skat_min2_code$BE_pvalue[idx]
    # genelist$BE_fwer[i]=skat_min2_code$BE_fwer[idx]

    #work on gwas snp adjusted pvalues
    load(paste0(outfolder,"/skat_gwas_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    genesr2=rownames(res_min)[which(res_min$r2>r2cutoff)]
    skat_min2_code=skat_min2_code[rownames(skat_min2_code) %in% genesr2,]
    skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    idx=which(rownames(skat_min2_code)==genelist$gene[i])
    tmp=unlist(strsplit(genelist$type[i],";"))
    tmp1=rep(NA,length(tmp))
    for (j in 1:length(tmp))
    {
      if ("EAvsBE"==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BEA_pvalue[idx],digits = 3)
      }
      if ("EAvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$EA_pvalue[idx],digits = 3)
      }
      if ("BEvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BE_pvalue[idx],digits = 3)
      }
      if ("BEEAvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BEEA_pvalue[idx],digits = 3)
      }
    }
    genelist$pvalue_addgwas[i]=paste0(tmp1,collapse = ";")
    # if ("EAvsBE" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BEA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BEA_fwer[idx]
    # }
    # if ("EAvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$EA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$EA_fwer[idx]
    # }
    # if ("BEvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BE_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BE_fwer[idx]
    # }
    # if ("BEEAvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BEEA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BEEA_fwer[idx]
    # }
  }
  return(genelist)
}
genelist=read.csv("../result/NEWGenesValidate_allcovar_r2_01_FWER.csv")
genelist1=update_genelist(genelist=genelist)
#gwas data (24 snps available)
allgwassnps=read.table("../result/dong26snp_addcontrols_genotypedat.txt")
organidx <- 2
opt <- "HRC"

organs=c("mucosa","junction","stomach","muscularis","adipose","blood")
outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005")

outfolder=outfolders[organidx]

#load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
#no need to load bcagenotype data
#load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
#load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
#load(paste0(outfolder,"/skat_res.RData")) #saved skat res, skat_min2_pc6,skat_min2
#load GTEx gene expression data, snp,phenotype
if (opt=="HRC")
{
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
}


genename <- "IL2RB"
genename <- "COX7A2"
genename <- "FILIP1"
genename <- "HSP90AA1"
genename <- "KRTAP5-8"
genename <- "FOXF1"
genename <- "LDAH"
genename <- "ISYNA1"
genename <- "UBAC1"

cvmin=1
verbose=T
eqtl <- vector("list", 6)
for (organidx in 1:6){
  
  outfolder=outfolders[organidx]
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
  
  
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
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
  #ucor <- matrix(0,ncol(X1),2)
  #for (l in 1:ncol(X1)){
  #  ucov<- data.matrix(cbind(X1[,l],covariate))
  #  ufit <- lm(Y~ucov)
  #  ucor[l,1] <- summary(ufit)$coef[2,4]
  #  ucor[l,2] <- cor(Y,X1[,l])
  #}  
  #hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  #pcor <- ucor[,1]
  #X1 <- X1[,order(pcor,decreasing=T)]
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
  plot(cvfit)
  
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  eqtl[[organidx]] <- selsnps
  print(organs[organidx])
  print(length(Y))
  print(eqtl[[organidx]])
}  
  

