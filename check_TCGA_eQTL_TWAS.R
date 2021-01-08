
############################################################################
### This code modify Kevin's code and check eQTL prediction model and    ###
### perform TWAS analysis in BEACON data                                 ###
############################################################################

rm(list=ls())

#!/usr/bin/env Rscript
#USE 87 TCGA EAC samples and imputed genotypes
# #get all the data
# # #load copynumber
#use the code from James 4/18/2020

load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGA_EAC_imp_michigan_genotype.RData") #filter_michigan_impuation.R in EAprogression
#geneexp
phenotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_norm.txt"
phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
covariatefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_COVA_GE_PEER_clinical_pc4_used.txt"
phenotype=fread(phenotypefile,header=T,sep="\t")
phenotype=as.data.frame(phenotype)
rownames(phenotype)=phenotype[,1]
phenotype=phenotype[,-1]
phenotypepos=fread(phenotypeposfile,header=T,sep="\t")
phenotypepos=as.data.frame(phenotypepos)
rownames(phenotypepos)=phenotypepos[,1]
#all(rownames(phenotypepos)==rownames(phenotype)) #T
covariate=as.data.frame(fread(covariatefile))
rownames(covariate)=covariate[,1]
covariate=covariate[,-1]
covariate=t(covariate)
covariate=as.data.frame(covariate)
covariate$stage[is.na(covariate$stage)]=3 #glmnet not allow NA
idx=match(rownames(covariate),tcgaclinical$bcr_patient_barcode)
covariate$gender=tcgaclinical$gender[idx]
covariate$gender[covariate$gender=="MALE"]=1
covariate$gender[covariate$gender=="FEMALE"]=0
covariate$disease=NA #"Adenomas"
tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Adenomas_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
idx=match(tmp$submitter_id,rownames(covariate))
covariate$disease[idx]=0 #"Adenomas"
tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Squmous_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
idx=match(tmp$submitter_id,rownames(covariate))
covariate$disease[idx]=1 #"Squamous"
#1 sample "TCGA-2H-A9GG" doesn't have copynumber data
cosample=intersect(colnames(snp),colnames(copynumber_new))
cosample=intersect(colnames(phenotype),cosample)
cosample=intersect(rownames(covariate)[covariate$disease==0],cosample)
idx=match(cosample,colnames(snp))
snp=snp[,idx]
idx=match(cosample,colnames(phenotype))
phenotype=phenotype[,idx]
idx=match(cosample,colnames(copynumber_new))
copynumber=copynumber_new[,idx]
idx=match(cosample,rownames(covariate))
covariate=covariate[idx,]
covariate=covariate[,-ncol(covariate)]
#sum(rownames(covariate)==colnames(phenotype))
#sum(colnames(phenotype)==colnames(snp))
#sum(colnames(copynumber)==colnames(snp))
#sum(rownames(copynumber)==rownames(phenotype))
idx=match(cosample,colnames(mutation))
mutation=mutation[,idx]
sum(rownames(mutation) %in% rownames(phenotype))
rownames(mutation)[!rownames(mutation) %in% rownames(phenotype)][1:3]
mutation1=mutation[rownames(mutation) %in% rownames(phenotype),]
mutation_new=data.frame(matrix(0,nrow=nrow(phenotype),ncol=ncol(phenotype)))
rownames(mutation_new)=rownames(phenotype)
colnames(mutation_new)=colnames(phenotype)
idx=match(rownames(mutation1),rownames(mutation_new))
mutation_new[idx,]=mutation1
mutation=mutation_new
all(colnames(mutation_new)==colnames(phenotype)) #T
#save(snp,snppos,phenotype,phenotypepos,copynumber,mutation,covariate,file="../result/TCGAdatafor_prediction_michigan.RData")

removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
  return(datnew)
}




library(SKAT)

outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_April18" #models stored for TCGA, came from prediction_michigan_models7_TCGA.R
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18" #models stored for GTEx, came from prediction_michigan_models6_GTEx.R
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_April18" #models stored for GTEx mucosa

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


beaconfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_GTExmodel/"

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

# allsnps=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt"))
# allsnps=allsnps[allsnps$`#chrom` %in% paste0("chr",c(1:22,"X","Y")),]
# allsnps$`#chrom`=gsub("chr","",allsnps$`#chrom`)

library("biomaRt")
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")

summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_autosomes.txt")
summaryfile2=paste0(summaryfolder,"EA_Bonn_autosomes.txt")
summaryfile3=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt")


beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic")

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



load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_allsnp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_allpos=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
allsnp=snp
allsnppos=snppos
allphenotype=phenotype
allphenotypepos=phenotypepos
allcopynumber=copynumber
allmutation=mutation


#compute
#test=compute_cor_arow(i=31038,opt="min",ncv=10,distcutoff = 5e5)
which(row.names(phenotype)=="DDX49")
#[1] 4777
which(row.names(phenotype)=="CERS1")
#[1] 31036
which(row.names(phenotype)=="KXD1")
#i <- 31038
which(row.names(phenotype)=="HSP90AA1")
#i <- 25105
which(row.names(phenotype)=="UBAC1")
#i <- 17803
which(row.names(phenotype)=="FBP2")
#i <- 6235




TWAS_SKAT <- function(genename,cvmin=1) {
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
p_cn=NA
p_mutation=NA
p_stage=NA
p_gender=NA
p_age=NA
tmp=quantile(Y,probs=c(0.15,0.85))
if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)

  #too many highly correlated SNPs are included for glmnet
  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  } 
  #rownames(ucor)=colnames(X1)
  #hist(ucor[,1])
  pcor <- ucor[,1]
  ## I order the snps based on their correlation with gene expression, first delete those low correlation-with-gene SNPs
  X1 <- X1[,order(pcor,decreasing=T)]
  X <- removehighcorr(X1,0.8)
  if (class(X)=="numeric") #only 1 snp left
  {
    X=matrix(X,nrow=length(X),ncol=1)
  }
  #sum(colnames(X)=="19:18817903_A_C")
  #hist(ucor[colnames(X1)%in% colnames(X),1])
  #hist(ucor[colnames(X1)%in% colnames(X),2])
  #dim(X)
  
  mutationV=unlist(mutation[i,])
  if (all(mutationV==0)) #no mutation
  {
    Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),covariate))
    covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),covariate))
  }else
  {
    Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
    covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
  }
  
  penalty=rep(1,ncol(Xall))
  penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      ### I change alpha to 0.5 for better variable selection when highly correlated features
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5)
    },
    error=function(e)
    {
      return(F)
    }
  )
  
  plot(cvfit)
  
  if (is.list(cvfit)) {
    ## do 100 times cv.glmnet and take average for cverr
    ## the number of cvfit$lambda may be less than 100 sometimes even you specified 100
    cverr <- matrix(NA,length(cvfit$lambda),100)
    rownames(cverr)=cvfit$lambda
    for (l in 1:100) {
      set.seed(l+100)
      fit=tryCatch(
        {
          ### I change alpha to 0.5 for better variable selection when highly correlated features
          cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5)
        },
        error=function(e)
        {
          return(F)
        }
      )
      if (is.list(fit)) #Error in predmat[which, seq(nlami)] = preds : replacement has length zero
      {
        #fit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
        alllambda=intersect(cvfit$lambda,fit$lambda)
        idx1=match(alllambda,cvfit$lambda)
        idx2=match(alllambda,fit$lambda)
        cverr[idx1,l] <- fit$cvm[idx2]
      }
    }
    merr <- apply(cverr,1,mean,na.rm=T)
    plot(log(cvfit$lambda),merr)
    lambda.best <- cvfit$lambda[which.min(merr)]
    
    fit=glmnet(as.matrix(Xall),Y, penalty.factor=penalty,alpha=0.5)
    glmcoeff=as.matrix(coef(fit,s=lambda.best))
    sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
    #glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
    glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
    idx1=match(glmleftsnp,rownames(glmcoeff))
    glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snp by effect size
    numvar=length(glmleftsnp)
    selcoeff <- glmcoeff[idx1,1]
    selcoeff <- selcoeff[order(abs(selcoeff),decreasing=T)]
  }  
    
  selsnps <- glmleftsnp
  
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
  
  ##############################
  ### compare BE vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  
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
  
  
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
      
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
          #tmp=listAttributes(snpmart)
          tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
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
      if (nrow(val1)>1)  {
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
      } else  {
        result[1,1] <- length(tmp)
        result[1,2] <- mean(val1$P<0.05)
        result[1,3] <- val1$P
        result[1,4] <- val1$P
      }
      
      ### validation in Bonn EA/BE ###
      
      tmp=intersect(rsid[!is.na(rsid)],summarydat3$SNP)
      valsnps1 <- selsnps[which(rsid %in% tmp)]
      rsid1 <- rsid[rsid %in% tmp]
      val1 <- summarydat3[summarydat3$SNP %in% tmp,]
      val1 <- val1[match(rsid1,val1$SNP),]
      if (nrow(val1)>1)  {
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
      }else  {
        result[2,1] <- length(tmp)
        result[2,2] <- mean(val1$P<0.05)
        result[2,3] <- val1$P
        result[2,4] <- val1$P
      }
      ### validation in Bonn BE ###
      
      tmp=intersect(rsid[!is.na(rsid)],summarydat$SNP)
      valsnps1 <- selsnps[which(rsid %in% tmp)]
      rsid1 <- rsid[rsid %in% tmp]
      val1 <- summarydat[summarydat$SNP %in% tmp,]
      val1 <- val1[match(rsid1,val1$SNP),]
      if (nrow(val1)>1)  {
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
      }else  {
        result[3,1] <- length(tmp)
        result[3,2] <- mean(val1$P<0.05)
        result[3,3] <- val1$P
        result[3,4] <- val1$P
      }
      
      
      
      ### validation in Oxford BE ###
      
      library(survey)
      
      tmp=intersect(rsid[!is.na(rsid)],summarydat1$SNP)
      valsnps2 <- selsnps[which(rsid %in% tmp)]
      rsid1 <- rsid[rsid %in% tmp]
      val2 <- summarydat1[summarydat1$SNP %in% tmp,]
      val2 <- val2[match(rsid1,val2$SNP),]
      if (nrow(val2)>1)  {
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
      }else  {
        result[4,1] <- length(tmp)
        result[4,2] <- mean(val2$P<0.05)
        result[4,3] <- val2$P
        result[4,4] <- val2$P
      }
      
      
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
      result[5,2] <- mean(val1$P<0.05)
      result[5,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
      result[5,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
      } 
      return(result)  
      
}

out <- TWAS_SKAT(genename="DDX49",cvmin=0)   

        
#[,1]         [,2]         [,3]         [,4]
#[1,]    1 0.000000e+00 5.660090e-01 5.660090e-01
#[2,]    1 0.000000e+00 5.721390e-01 5.721390e-01
#[3,]    1 0.000000e+00 6.686600e-01 6.686600e-01
#[4,]    1 1.000000e+00 2.039240e-03 2.039240e-03
#[5,]    1 3.589674e-06 3.589674e-06 3.589674e-06
#[6,]    1 3.303583e-05 3.303583e-05 3.303583e-05
#[7,]    1 2.693097e-07 2.693097e-07 2.693097e-07        
      


out <- TWAS_SKAT(genename="MEF2B",cvmin=0)   

        
#[,1]         [,2]         [,3]       [,4]
#[1,]   29 6.896552e-02 6.158139e-01 0.59871391
#[2,]   29 6.896552e-02 5.595357e-01 0.53195499
#[3,]   29 0.000000e+00 7.760890e-01 0.79464297
#[4,]   29 1.034483e-01 4.566042e-01 0.43079065
#[5,]   29 4.664422e-06 8.661454e-07 0.01636555
#[6,]   29 1.008548e-05 2.179656e-06 0.06428328
#[7,]   29 1.270885e-08 4.236283e-09 0.03113452       
        
out <- TWAS_SKAT(genename="FBP2",cvmin=0)   

#[,1]         [,2]         [,3]       [,4]
#[1,]   21 0.000000e+00 4.940105e-01 0.46779890
#[2,]   21 0.000000e+00 1.809853e-01 0.16616466
#[3,]   21 4.761905e-02 2.205167e-01 0.20133741
#[4,]   21 9.523810e-02 5.142105e-02 0.05284884
#[5,]   24 3.538334e-04 1.494097e-04 0.02434512
#[6,]   24 3.790426e-04 1.624341e-04 0.01186162
#[7,]   24 2.386209e-05 8.170888e-06 0.00613826



out <- TWAS_SKAT(genename="COPE",cvmin=0)   

#[,1]         [,2]         [,3]       [,4]
#[1,]   21 4.761905e-02 6.728062e-01 0.67484578
#[2,]   21 0.000000e+00 8.664114e-01 0.91345791
#[3,]   21 0.000000e+00 8.767791e-01 0.92575634
#[4,]   21 9.523810e-02 8.686593e-01 0.95744637
#[5,]   24 7.174821e-06 3.859688e-06 0.04364718
#[6,]   24 9.357740e-05 4.303845e-05 0.01778153
#[7,]   24 1.090448e-06 3.908016e-07 0.01541469



out <- TWAS_SKAT(genename="GDF15",cvmin=0)   

#[,1]         [,2]         [,3]         [,4]
#[1,]    8 0.000000e+00 7.416098e-01 0.7577637674
#[2,]    8 0.000000e+00 5.092974e-01 0.4885700559
#[3,]    8 0.000000e+00 3.605494e-01 0.3289163750
#[4,]    8 1.250000e-01 4.910167e-01 0.4625381537
#[5,]    9 2.128220e-04 1.059368e-04 0.0027411960
#[6,]    9 3.155392e-04 1.521150e-04 0.0003822332
#[7,]    9 1.593454e-05 6.983828e-06 0.0001977101


  





which(row.names(phenotype)=="CERS1")
#[1] 31036
which(row.names(phenotype)=="PLEKHF2")
#[1] 15931
which(row.names(phenotype)=="GATA4")

which(row.names(phenotype)=="CTSB")

which(row.names(phenotype)=="SGK223")

which(row.names(phenotype)=="KXD1")

which(row.names(phenotype)=="PPP1R3B")

which(row.names(phenotype)=="COMP")

which(row.names(phenotype)=="OR11A1")

which(row.names(phenotype)=="SPAG11A")

which(row.names(phenotype)=="MEF2B")
