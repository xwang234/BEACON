
############################################################################
### This code modify Kevin's code and check eQTL prediction model and    ###
### perform TWAS analysis in BEACON data                                 ###
############################################################################

rm(list=ls())


removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}

library(SKAT)

### loading genotype data ###

#For mucosa
#prefix="dist500K_GTEx_mucosa_May26"
#For junction
prefix="dist500K_GTEx_blood_June11"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(file=paste0(outfolder,"/bca_extractgenotype.RData"))
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


### loading phenotype data ###


library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)

sampletable <- data.frame(sampletable)
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


### loading summary stat data for validation ###

beaconfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_GTExmodel/"

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
proteingenesV8=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

# allsnps=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt"))
# allsnps=allsnps[allsnps$`#chrom` %in% paste0("chr",c(1:22,"X","Y")),]
# allsnps$`#chrom`=gsub("chr","",allsnps$`#chrom`)

library("biomaRt")
#variation = useEnsembl(biomart="snp")
#tmp1=listDatasets(variation)
#hg38
snpmart_hg38 = useEnsembl(biomart="snp", dataset="hsapiens_snp")

summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_autosomes.txt")

summaryfile2=paste0(summaryfolder,"EA_Bonn_autosomes.txt")
summaryfile3=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt")

#beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic")

library(data.table)
#summary stat from validation
summarydat=as.data.frame(fread(summaryfile,header=T))
summarydat1=as.data.frame(fread(summaryfile1,header=T))

summarydat2=as.data.frame(fread(summaryfile2,header=T))
summarydat3=as.data.frame(fread(summaryfile3,header=T))

tmp=intersect("SNP",colnames(summarydat))
if (length(tmp)==0)
{
  colnames(summarydat)[which(colnames(summarydat)=="rsid")]="SNP"
  colnames(summarydat)[which(colnames(summarydat)=="pvalue")]="P"
  colnames(summarydat)[which(colnames(summarydat)=="beta")]="BETA"
  colnames(summarydat)[which(colnames(summarydat)=="effect-allele")]="effect_allele"
  colnames(summarydat)[which(colnames(summarydat)=="non-effect-allele")]="non_effect_allele"
}

tmp=intersect("SNP",colnames(summarydat1))
if (length(tmp)==0)
{
  colnames(summarydat1)[which(colnames(summarydat1)=="rsid")]="SNP"
  colnames(summarydat1)[which(colnames(summarydat1)=="pvalue")]="P"
  colnames(summarydat1)[which(colnames(summarydat1)=="beta")]="BETA"
  colnames(summarydat1)[which(colnames(summarydat1)=="effect-allele")]="effect_allele"
  colnames(summarydat1)[which(colnames(summarydat1)=="non-effect-allele")]="non_effect_allele"
}



tmp=intersect("SNP",colnames(summarydat2))
if (length(tmp)==0)
{
  colnames(summarydat2)[which(colnames(summarydat2)=="rsid")]="SNP"
  colnames(summarydat2)[which(colnames(summarydat2)=="pvalue")]="P"
  colnames(summarydat2)[which(colnames(summarydat2)=="beta")]="BETA"
  colnames(summarydat2)[which(colnames(summarydat2)=="effect-allele")]="effect_allele"
  colnames(summarydat2)[which(colnames(summarydat2)=="non-effect-allele")]="non_effect_allele"
}


tmp=intersect("SNP",colnames(summarydat3))
if (length(tmp)==0)
{
  colnames(summarydat3)[which(colnames(summarydat3)=="rsid")]="SNP"
  colnames(summarydat3)[which(colnames(summarydat3)=="pvalue")]="P"
  colnames(summarydat3)[which(colnames(summarydat3)=="beta")]="BETA"
  colnames(summarydat3)[which(colnames(summarydat3)=="effect-allele")]="effect_allele"
  colnames(summarydat3)[which(colnames(summarydat3)=="non-effect-allele")]="non_effect_allele"
}


#summary from beacon
#beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))




#summary from beacon
#beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))

##load GTEx gene expression data
#For mucosa
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadatafor_prediction.RData")
#For junction: three version of junction datasets, V7, V8GTEx,V8selfprocessed

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8blooddata_ambiguous_TPM_for_prediction.RData")

library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp


#compute
#test=compute_cor_arow(i=31038,opt="min",ncv=10,distcutoff = 5e5)
which(row.names(phenotype)=="THAP6")
#[1] 5657
which(row.names(phenotype)=="CERS1")
#[1] 31036
which(row.names(phenotype)=="KXD1")
#i <- 31038
which(row.names(phenotype)=="HSP90AA1")
#i <- 25105
which(row.names(phenotype)=="UBAC1")
#i <- 17803
which(row.names(phenotype)=="FBP2")
#i <- 17243

which(row.names(phenotype)=="HSP90AA1")
genename <- "HSP90AA1"
genename <- "ZPLD1"

genename <- "PLA2G4C"

genename <- "KLHL26"

genename <- "PGPEP1"

genename <- "LSM4"
genename <- "GDF15"
genename <- "HOMER3"

TWAS_SKAT <- function(genename,cvmin=1) 
{
  result <- matrix(0,8,4)  
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
  hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  dim(X1)
  X <- removehighcorr(X1,0.9)
  hist(ucor[colnames(X1)%in% colnames(X),1])
  hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      #I used standardize
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  cvfit$lambda.1se
  
  cverr <- matrix(NA,length(cvfit$lambda),10)
  cvse  <- matrix(NA,length(cvfit$lambda),10)
  
  rownames(cverr)=cvfit$lambda
  for (l in 1:10) {
    set.seed(l+10)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    alllambda=intersect(cvfit$lambda,fit$lambda)
    idx1=match(alllambda,cvfit$lambda)
    idx2=match(alllambda,fit$lambda)
    cverr[idx1,l] <- fit$cvm[idx2]
    cvse[idx1,l] <- fit$cvsd[idx2]
  }
  
  merr <- apply(cverr,1,mean,na.rm=T)
  mcvse <- sqrt(apply(cvse^2,1,mean))/10
  
  lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  lambda.min <- cvfit$lambda[which.min(merr)]
  
  plot(log(cvfit$lambda),merr)
  abline(v=log(lambda.1se),lty=2)
  abline(v=log(lambda.min),lty=3)
  
  #plot(cvfit$lambda,merr)
  glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) 
  #glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.1se))
  #if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
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
  ## there a few SNPs not found in bcagenotype
  sum(is.na(idx))
  
  selsnps <- selsnps[!is.na(idx)]
  selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
  #if (length(idxtocorrect)>0)
  #{
  #  Z[idxtocorrect,]=2-Z[idxtocorrect,]
  #}
  
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bca==2) #case
  idx2=which(sampletable$phenoBE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[6,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[6,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[6,4] <- out$p.value
  
  
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoEA_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[7,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[7,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,4] <- out$p.value
  
  
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
  ##############################
  ### compare EA vs BE       ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoEA_bca==2) #case
  idx2=which(is.na(sampletable$phenoEA_bca))
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[7,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[7,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,4] <- out$p.value
  
  outp <- matrix(0,ncol(Z1),1)
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    outp[i] <- summary(outfit)$coef[2,4]
    #print(summary(outfit)$coef[2,4])
  }  
  outp <- data.frame(outp)
  rownames(outp) <- colnames(Z1)
  
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoEABE_bca==2) #case
  idx2=which(sampletable$phenoEABE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[8,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[8,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[8,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[8,4] <- out$p.value
  
  
  
  selsnps <- colnames(Z1)
  tmp=unlist(strsplit(selsnps[1],":"))
  chr=tmp[1]
  pos=rsid=allele1=allele2=p_uni=p_mul=cor_uni=rep(NA,length(selsnps))
  #selsnps="19:19259262_C_A"
  for (j in 1:length(selsnps))
  {
    cat(j,"\n")
    tmp=unlist(strsplit(selsnps[j],":"))
    tmp=unlist(strsplit(tmp[2],"_"))
    allalleles=tmp[2:3]
    allele1[j]=tmp[2] #minor allele
    allele2[j]=tmp[3]
    pos[j]=as.numeric(tmp[1])
    #find rsid based on position and allles
    #tmp=listAttributes(snpmart_hg38)
    tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart_hg38)
    if (nrow(tmp)>0)
    {
      for (k in 1:nrow(tmp))
      {
        myalleles=unlist(strsplit(tmp$allele[k],"/",fixed=T))
        if (all(allalleles %in% myalleles))
        {
          rsid[j]=tmp$refsnp_id[k]
          break
        }
      }
    }
  }
  
  library(survey)
  
  ### validation in Bonn EA ###
  
  tmp=intersect(rsid[!is.na(rsid)],summarydat2$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[rsid %in% tmp]
  val1 <- summarydat2[summarydat2$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,colnames(Z1) %in% valsnps1]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[1,1] <- length(tmp)
  result[1,2] <- mean(val1$P<0.05)
  result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  
  
  ### validation in Bonn EA/BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],summarydat3$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[rsid %in% tmp]
  val1 <- summarydat3[summarydat3$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,colnames(Z1) %in% valsnps1]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[2,1] <- length(tmp)
  result[2,2] <- mean(val1$P<0.05)
  result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  
  ### validation in Bonn BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],summarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[rsid %in% tmp]
  val1 <- summarydat[summarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,colnames(Z1) %in% valsnps1]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[3,1] <- length(tmp)
  result[3,2] <- mean(val1$P<0.05)
  result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  
  ### validation in Oxford BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],summarydat1$SNP)
  valsnps2 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[rsid %in% tmp]
  val2 <- summarydat1[summarydat1$SNP %in% tmp,]
  val2 <- val2[match(rsid1,val2$SNP),]
  uu <- val2$BETA
  vv <- as.numeric(val2$se)
  ZZ <- Z1[,colnames(Z1) %in% valsnps2]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  
  result[4,1] <- length(tmp)
  result[4,2] <- mean(val2$P<0.05)
  result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  #sum(lamb*(1-pchisq(Q,df=1)))
  result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  
  
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
  } 
  
  return(result)  
        
}
        
        
out <- TWAS_SKAT(genename="THAP6",cvmin=1) 
# [,1]       [,2]        [,3]         [,4]
# [1,]   21 0.04761905 0.395937664 3.785297e-01
# [2,]   21 0.00000000 0.315019932 2.980214e-01
# [3,]   21 0.04761905 0.262445326 2.470835e-01
# [4,]   21 0.14285714 0.005026368 6.933387e-03
# [5,]   21 0.04761905 0.806214402 8.142514e-01
# [6,]   22 0.03810708 0.020821147 3.978476e-05
# [7,]   22 0.02433690 0.115928856 5.505900e-03
# [8,]   22 0.01955131 0.018361642 4.100227e-05

        
        
out <- TWAS_SKAT(genename="HSP90AA1",cvmin=1)




        
out <- TWAS_SKAT(genename="HSP90AA1",cvmin=0)


  

out <- TWAS_SKAT(genename="KXD1",cvmin=1)





out <- TWAS_SKAT(genename="UBAC1",cvmin=1)



out <- TWAS_SKAT(genename="UBAC1",cvmin=0)

out <- TWAS_SKAT(genename="ISYNA1",cvmin=1)



out <- TWAS_SKAT(genename="ISYNA1",cvmin=1)





#functions to compare snpid in V7 and V8, used in compare SNPs selected in V7 and V8--------
#snpid hg38 --->hg19---
#input: rownames(genotypeV8) are a set of SNPIDs such as 4:1234_T_C, coordiates are in hg38
#output: data.frame(snphg38=snpnames,snphg19=newsnpnames),snphg38 is the original snpid,snphg19 is the snpid in hg19
hg38tohg19=function(snpnames=rownames(genotypeV8))
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

#check overlap of V8 snps and V7 snps, put V8 first---
overlap_snp1_snp2=function(snpnames1=rownames(genotypeV8),snpnames2=rownames(genotypeV7))
{
  #transform hg38id to hg19id
  newsnpnames1=hg38tohg19(snpnames=snpnames1)
  newsnpnames1=newsnpnames1[!is.na(newsnpnames1$snphg19),]
  #use hg19id to merge with snpnames2 (which is based on hg19)
  newsnpnames1_hg19=newsnpnames1$snphg19
  sum(newsnpnames1_hg19 %in% snpnames2)
  #some thimes the above not working, they may sit on different strands (A/G and T/C), to use positon to match
  #use position name (posname,such as 4:1234) to match 2 datasets in hg19
  posnames1=rep(NA,length(newsnpnames1_hg19))
  posnames2=rep(NA,length(snpnames2))
  get_posname=function(snpnames=newsnpnames1_hg19)
  {
    tmp=unlist(strsplit(snpnames,"_"))
    posnames=tmp[seq(1,length(tmp),3)]
    return(posnames)
  }
  posnames1=get_posname()
  posnames2=get_posname(snpnames=snpnames2)
  newsnpnames1$posnames=posnames1
  newsnpnames2=data.frame(snphg19=snpnames2,posnames=posnames2,stringsAsFactors = F)
  comsnps=intersect(posnames1,posnames2)
  print(paste0("snp1: ",length(snpnames1),", snp2: ",length(snpnames2),", common snp: ",length(comsnps)))
  newsnpnames1$overlap=F
  idx=match(comsnps,posnames1)
  newsnpnames1$overlap[idx]=T
  idx=match(comsnps,posnames2)
  newsnpnames2$overlap=F
  newsnpnames2$overlap[idx]=T
  return(list(snpnames1=newsnpnames1,snpnames2=newsnpnames2))
}






