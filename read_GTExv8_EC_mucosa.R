
#GTEx method
#https://gtexportal.org/home/documentationPage#staticTextAnalysisMethods
#https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
#gene expression use TPM+TMM. In the first verion, GTEx processed TMM was used

#Install peer on R3.6
# ml CMake/3.0.0-foss-2014b
# R CMD INSTALL -l /home/xwang234/R/x86_64-pc-linux-gnu-library/3.6 R_peer_source_1.3.tgz

library(data.table)

gtexv8_samplestable=fread("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt",sep="\t")
gtexv8_samplestable=as.data.frame(gtexv8_samplestable)
gtexv8_samplestable$SUBJID=NA
for (i in 1:nrow(gtexv8_samplestable))
{
  tmp=unlist(strsplit(gtexv8_samplestable$SAMPID[i],"-"))
  gtexv8_samplestable$SUBJID[i]=paste0(tmp[1:2],collapse = "-")
}
gtexv8_samples=gtexv8_samplestable$SAMPID[gtexv8_samplestable$SMTS=="Esophagus" & gtexv8_samplestable$SMTSD=="Esophagus - Mucosa"]


find_subjectid=function(sampleid=gtexv8_samples)
{
  subjectid=rep(NA,length(sampleid))
  for (i in 1:length(sampleid))
  {
    if (grepl("-",sampleid[i]))
    {
      tmp=unlist(strsplit(sampleid[i],"-"))
      subjectid[i]=paste0(tmp[1:2],collapse = "-")
    }else
    {
      subjectid[i]=sampleid[i]
    }
  }
  return(subjectid)
}
gtexv8_subjects=find_subjectid(gtexv8_samples)

subjectablev8=read.table("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt",header=T,fill=T,quote="",sep="\t",stringsAsFactors = F)


tmp=gtexv8_subjects[gtexv8_subjects %in% subjectablev8$SUBJID[subjectablev8$RACE==3]] #336
length(unique(tmp)) #511

#processed expression---
# gtexv8_ge_anno=as.data.frame(fread("../data/GTEx/V8/gencode.v26.GRCh38.genes.gtf"))
# tmp=strsplit(gtexv8_ge_anno$V9,";",fixed = T)
# gtexv8_ge_anno$Probe_Id=NA #gene_id is the same as transcript_id
# gtexv8_ge_anno$Symbol=NA
# gtexv8_ge_anno$gene_type=NA
# 
# for (i in 1:nrow(gtexv8_ge_anno))
# {
#   if (i %% 10000==0) cat(i,"..")
#   tmp1=tmp[[i]]
#   tmp2=unlist(strsplit(tmp1[1],'\"',fixed=T))
#   tmp3=unlist(strsplit(tmp1[3],'\"',fixed=T))
#   tmp4=unlist(strsplit(tmp1[4],'\"',fixed=T))
#   gtexv8_ge_anno$Probe_Id[i]=tmp2[2]
#   gtexv8_ge_anno$Symbol[i]=tmp4[2]
#   gtexv8_ge_anno$gene_type[i]=tmp3[2]
# }
# colnames(gtexv8_ge_anno)[1]="Chromosome"
# colnames(gtexv8_ge_anno)[4]="start"
# colnames(gtexv8_ge_anno)[5]="end"
# save(gtexv8_ge_anno,file="../data/GTEx/gtexv8_ge_anno.RData",version=2)
load("../data/GTEx/gtexv8_ge_anno.RData")
#including all samples
# test=fread("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
#            sep="\t",nrows = 3)
# 
# gtexv8_GE=fread("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
#                 sep="\t",select=c(colnames(test)[1:2],gtexv8_samples))
# gtexv8_GEsamples=colnames(gtexv8_GE)[3:ncol(gtexv8_GE)]
# gtexv8_GEsubjects=rep(NA,length(gtexv8_GEsamples))
# for (i in 1:length(gtexv8_GEsamples))
# {
#   tmp=unlist(strsplit(gtexv8_GEsamples[i],"-"))
#   gtexv8_GEsubjects[i]=paste0(tmp[1:2],collapse = "-")
# }

#find sampleid in specific organ, there could be duplicated samples
find_sampleid=function(subjectid=gtexv8_GEsubjects)
{
  sampleid=NULL
  subjid=NULL
  for (i in 1:length(subjectid))
  {
    idx=which(gtexv8_subjects==subjectid[i])
    sampleid=c(sampleid,gtexv8_samples[idx])
    subjid=c(subjid,rep(subjectid[i],length(idx)))
  }
  return(list(sampleid=sampleid,subjectid=subjid))
}
#gtexv8_GEsamples=find_sampleid()
# idx=match(gtexv8_GEsamples,gtexv8_samplestable$SAMPID)
# table(gtexv8_samplestable$SMTSD[idx])

gtexv8_GE_norm=fread('gunzip -cq ../data/GTEx/V8/expression_matrices/Esophagus_Mucosa.v8.EUR.normalized_expression.bed.gz')
gtexv8_GE_norm=as.data.frame(gtexv8_GE_norm)
#gtexv8_GE_norm=fread('gunzip -cq ../data/GTEx/V8/GTEx_Analysis_v8_eQTL_expression_matrices/Esophagus_Gastroesophageal_Junction.v8.normalized_expression.bed.gz')
gtexv8_ge_anno1=gtexv8_ge_anno[gtexv8_ge_anno$V3=="gene",]
idx=which(duplicated(gtexv8_ge_anno1$Symbol))
gtexv8_ge_anno1=gtexv8_ge_anno1[-idx,]
gtexv8_GE_norm=gtexv8_GE_norm[gtexv8_GE_norm$gene_id %in% gtexv8_ge_anno1$Probe_Id,]
gtexv8_GE_normsubjects=colnames(gtexv8_GE_norm)[5:ncol(gtexv8_GE_norm)]
gtexv8_GE_IDs=find_sampleid(gtexv8_GE_normsubjects)
gtexv8_GE_normsamples=gtexv8_GE_IDs$sampleid
idx=match(gtexv8_GE_norm$gene_id,gtexv8_ge_anno1$Probe_Id)
#table(gtexv8_ge_anno1$gene_type[idx])
rownames(gtexv8_GE_norm)=gtexv8_ge_anno1$Symbol[idx]
#read TPM
gtexv8_GE_TPM=as.data.frame(fread("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",select =c("Name","Description",gtexv8_GE_normsamples )))
#align TPM with norm data
idx=match(gtexv8_GE_norm$gene_id,gtexv8_GE_TPM$Name)
if (sum(is.na(idx)>0)) warning("some genes are missing!")
gtexv8_GE_TPM=gtexv8_GE_TPM[idx,]
colnames(gtexv8_GE_TPM)=find_subjectid(colnames(gtexv8_GE_TPM))
idx=match(colnames(gtexv8_GE_norm)[5:ncol(gtexv8_GE_norm)],colnames(gtexv8_GE_TPM))
if (sum(is.na(idx)>0)) warning("some samples are missing!")
gtexv8_GE_TPM=gtexv8_GE_TPM[,idx]
rownames(gtexv8_GE_TPM)=rownames(gtexv8_GE_norm)
#do TMM
library(edgeR)
dgList=DGEList(gtexv8_GE_TPM)
scale.factors <- calcNormFactors(dgList, method = "TMM")
norm.data <- as.data.frame(t(t(scale.factors$counts)/(scale.factors$samples$norm.factors*scale.factors$samples$lib.size)))
normdata=cpm(scale.factors)
cor(norm.data[,100],normdata[,100])#1
normdata=t(scale(t(normdata)))
cor(as.numeric(norm.data[1,]),as.numeric(normdata[1,]))
boxplot(normdata[,1:10])
phenotype=as.data.frame(normdata)
#phenotype=gtexv8_GE_norm[,5:ncol(gtexv8_GE_norm)] #the first version
#rownames(phenotype)=rownames(gtexv8_GE_norm)
#idx=match(rownames(gtexv8_GE_norm),gtexv8_ge_anno1$Symbol)
idx=match(rownames(phenotype),gtexv8_ge_anno1$Symbol)
phenotypepos=data.frame(geneid=rownames(gtexv8_GE_norm),chr=gtexv8_ge_anno1$Chromosome[idx],s1=gtexv8_ge_anno1$start[idx],s2=gtexv8_ge_anno1$end[idx],stringsAsFactors = F)
phenotypepos$chr=gsub("chr","",phenotypepos$chr)
#rownames(phenotypepos)=rownames(gtexv8_GE_norm)
rownames(phenotypepos)=rownames(phenotype)


gtexv8_GE_norm_covariate=as.data.frame(fread('../data/GTEx/V8/expression_covariates/Esophagus_Mucosa.v8.EUR.covariates.txt'))
gtexv8_GE_norm_covariatesubjects=colnames(gtexv8_GE_norm_covariate)[2:ncol(gtexv8_GE_norm_covariate)]
gtexv8_GE_norm_covariateIDs=find_sampleid(gtexv8_GE_norm_covariatesubjects)
#all(gtexv8_GE_normsubjects == gtexv8_GE_norm_covariatesubjects) #T
all(colnames(phenotype) == gtexv8_GE_norm_covariatesubjects) #T
peer_number=function(dat=GE,numfactors=15,maxit=100)
{
  
  totvar=0
  for (i in 1:nrow(dat)) totvar=totvar+var(unlist(dat[i,]))
  
  library(peer)
  
  model=PEER()
  PEER_setNk(model,numfactors)
  PEER_setPhenoMean(model, as.matrix(t(dat)))
  PEER_setNmax_iterations(model, maxit)
  PEER_update(model)
  resi=PEER_getResiduals(model)
  resivar=0
  for (j in 1:ncol(resi)) resivar=resivar+var(resi[,j])
  print((totvar-resivar)/totvar)
  varexplain=(totvar-resivar)/totvar
  
  print(paste0("number of factors: ",numfactors))
  print(paste0("variance explained: ",round(varexplain,digits = 3)))
  #the factors
  x=PEER_getX(model)
  x=as.data.frame(t(x))
  colnames(x)=colnames(dat)
  res=cbind.data.frame(id=paste0("factor",1:numfactors),x)
  return(res)
}
#
# phenotype_peer=peer_number(dat=gtexv8_GE_norm[,5:ncol(gtexv8_GE_norm)]) #use the GTEx normalized data 
# peer1=phenotype_peer[,2:ncol(phenotype_peer)]
# peer2=gtexv8_GE_norm_covariate[6:20,2:ncol(gtexv8_GE_norm_covariate)]
# tmp=cor(t(peer1),t(peer2)) #peer factors are consistent
# colnames(tmp)=c(1:ncol(tmp))
# heatmap(abs(tmp),Rowv=NA,Colv=NA)
# diag(tmp)
# # [1]  0.9999439  0.9998904  0.9997950  0.9997611  0.9998623  0.9996186  0.9997615 -0.9988768 -0.9981203
# # [10] -0.9986696  0.9987593 -0.9947988 -0.9691366  0.9894929  0.9118300
phenotype_peer1=peer_number(dat=phenotype) #based on TPM+TMM data
peer3=phenotype_peer1[,2:ncol(phenotype_peer1)]
# tmp=cor(t(peer3),t(peer2))
# colnames(tmp)=c(1:ncol(tmp))
# heatmap(abs(tmp),Rowv=NA,Colv=NA)
# diag(tmp)
# # [1] -0.99791280 -0.97531602  0.96793411 -0.95000428  0.98084499 -0.96572729 -0.97101021  0.93415751 -0.88788750
# # [10]  0.83424600 -0.80271821  0.88997924 -0.32199838 -0.83959153 -0.06239708

covariate=data.frame(matrix(NA,nrow=ncol(phenotype),ncol=4+15+4))
rownames(covariate)=colnames(phenotype)
all(rownames(covariate)==colnames(gtexv8_GE_norm_covariate)[2:ncol(gtexv8_GE_norm_covariate)]) #T
all(rownames(covariate)==colnames(phenotype)) #T
colnames(covariate)=c(paste0("pc",1:4),paste0("factor",1:15),"gender","age","pcr","platform")
covariate[,1:4]=t(gtexv8_GE_norm_covariate[1:4,2:ncol(gtexv8_GE_norm_covariate)])
#covariate[,5:19]=t(gtexv8_GE_norm_covariate[6:20,2:ncol(gtexv8_GE_norm_covariate)])
all(rownames(covariate)==colnames(peer3))#T
covariate[,5:19]=t(peer3)
covariate$gender=as.factor(t(gtexv8_GE_norm_covariate[which(gtexv8_GE_norm_covariate$ID=="sex"),2:ncol(gtexv8_GE_norm_covariate)]))
idx=match(rownames(covariate),subjectablev8$SUBJID)
covariate$age=subjectablev8$AGE[idx]
covariate$pcr=as.factor(t(gtexv8_GE_norm_covariate[which(gtexv8_GE_norm_covariate$ID=="pcr"),2:ncol(gtexv8_GE_norm_covariate)]))
covariate$platform=as.factor(t(gtexv8_GE_norm_covariate[which(gtexv8_GE_norm_covariate$ID=="platform"),2:ncol(gtexv8_GE_norm_covariate)]))

# Covariates
# Top 5 genotyping principal components.
# A set of covariates identified using the Probabilistic Estimation of Expression Residuals (PEER) method (Stegle et al., PLoS Comp. Biol., 2010 ), calculated for the normalized expression matrices (described below). For eQTL analyses, the number of PEER factors was determined as function of sample size (N): 15 factors for N<150, 30 factors for 150≤ N<250, 45 factors for 250≤ N<350, and 60 factors for N≥350, as a result of optimizing for the number of eGenes discovered. For sQTL analyses, 15 PEER factors were computed for each tissue.
# Sequencing platform (Illumina HiSeq 2000 or HiSeq X).
# Sequencing protocol (PCR-based or PCR-free).
# Sex.


#genotype
#wgssamples=read.table("../data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.fam")

#the code is not used, was used to use rsid to merge genotype data, but a lot genotype can't find in dbsnp
# dbsnp151_hg19=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180423.bim",drop="V3",sep="\t"))
# dbsnp151_hg38=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180418.bim",drop="V3",sep="\t"))
# tmp=intersect(dbsnp151_hg19$V2,dbsnp151_hg38$V2)
# idx1=match(tmp,dbsnp151_hg19$V2)
# idx2=match(tmp,dbsnp151_hg38$V2)
# #most rsid presverves allles from hg19 to hg38 (99.7%), use this subset
# sum((dbsnp151_hg19$V5[idx1]==dbsnp151_hg38$V5[idx2] & dbsnp151_hg19$V6[idx1]==dbsnp151_hg38$V6[idx2])) #535875322
# sum((dbsnp151_hg19$V5[idx1]==dbsnp151_hg38$V5[idx2] & dbsnp151_hg19$V6[idx1]==dbsnp151_hg38$V6[idx2]) |(dbsnp151_hg19$V5[idx1]==dbsnp151_hg38$V6[idx2] & dbsnp151_hg19$V6[idx1]==dbsnp151_hg38$V5[idx2])) #536217948
# sum(dbsnp151_hg19$V4[idx1]==dbsnp151_hg38$V4[idx2]) #4317225
# idx=which(dbsnp151_hg19$V5[idx1]==dbsnp151_hg38$V5[idx2] & dbsnp151_hg19$V6[idx1]==dbsnp151_hg38$V6[idx2])
# dbsnp151_hg19=dbsnp151_hg19[idx1[idx],]
# dbsnp151_hg38=dbsnp151_hg38[idx2[idx],]
# tmp=data.frame(snp=dbsnp151_hg19$V2,stringsAsFactors = F)
# fwrite(tmp,file="/fh/fast/stanford_j/Xiaoyu/Tools/annotation/dbsnp151_hg19_hg38_snpid.txt",col.names = F,row.names = F,quote=F)
# #run read_GTExv8_WGS.sh to split on chrs.
##

#find the common set of snps (GTEx and BCA, GTEx was liftover from hg38 to hg19)
#to remove ambigous SNPs
library(plink2R) #raw file in Tools/Rpackages
# snp=snppos=NULL
# for (chr in 1:22)
# {
#   cat(chr,'..')
#   # dat_bca=read_plink(paste0("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_noambiguous_chr",chr),impute="avg")
#   # bim1=dat_bca$bim
#   # traw1=dat_bca$bed
#   # tmp=unlist(strsplit(rownames(traw1),":"))
#   # rownames(traw1)=tmp[seq(1,length(tmp),2)]
#   # colnames(traw1)=paste0(chr,":",bim1$V4,"_",bim1$V5,"_",bim1$V6)
#   # traw1=as.data.frame(t(traw1))
#   # idx1=match(gtexv8_GE_normsubjects,colnames(traw1))
#   # if (sum(is.na(idx1))>0) warning(paste0("some RNAseq samples are missing in genotype data,",chr))
#   # traw1=traw1[,idx1]
#   
#   bim1file=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_noambiguous_chr",chr,".bim")
#   bim1=fread(bim1file)
#   bim1=as.data.frame(bim1)
#   bim2file=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr",chr,"_filter_noambiguous_hg19tohg38_flip.bim")
#   bim2=fread(bim2file)
#   bim2=as.data.frame(bim2)
#   str11=paste0(bim1$V4,"_",bim1$V5,"_",bim1$V6)
#   str12=paste0(bim1$V4,"_",bim1$V6,"_",bim1$V5)
#   str2=paste0(bim2$V4,"_",bim2$V5,"_",bim2$V6)
#   idx=str11 %in% str2 | str12 %in% str2
#   traw1file=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_noambiguous_chr",chr,".traw")
#   traw1=fread(traw1file)
#   traw1=as.data.frame(traw1)
#   tmp=paste0(traw1$CHR,":",traw1$POS,"_",traw1$COUNTED,"_",traw1$ALT)
#   rownames(traw1)=tmp
#   traw1=traw1[,7:ncol(traw1)]
#   tmp=unlist(strsplit(colnames(traw1),"_"))
#   colnames(traw1)=tmp[seq(1,length(tmp),2)]
#   idx1=match(gtexv8_GE_normsubjects,colnames(traw1))
#   if (sum(is.na(idx1))>0) warning(paste0("some RNAseq samples are missing in genotype data,",chr))
#   traw1=traw1[,idx1]
#   snp_=traw1[idx,]
#   snppos_=data.frame(chr=bim1$V1[idx],pos=bim1$V4[idx],minor=bim1$V5[idx],major=bim1$V6[idx],stringsAsFactors = F)
#   rownames(snppos_)=paste0(chr,":",snppos_$pos,"_",snppos_$minor,"_",snppos_$major)
#   snp=rbind.data.frame(snp,snp_)
#   snppos=rbind.data.frame(snppos,snppos_)
# }
# 
# # tmp=rowSums(snp)
# # snp=snp[!is.na(tmp),]
# # idx=which(is.na(tmp))
# # snppos=snppos[!is.na(tmp),]
# k <- which(is.na(snp), arr.ind=TRUE)
# snp[k] <- rowMeans(snp, na.rm=TRUE)[k[,1]]
# all(rownames(snppos)==rownames(snp))
# if (any(colnames(snp)!=colnames(phenotype))) warning("samples name not match!")
# all(colnames(snp)==colnames(phenotype))
# save(snp,snppos,phenotype,phenotypepos,covariate,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")
# 

snp=snppos=NULL
for (chr in 1:22)
{
  cat(chr,'..')
  
  bim1file=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter_flip_chr",chr,".bim")
  bim1=fread(bim1file)
  bim1=as.data.frame(bim1)
  bim2file=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr",chr,"_filter_hg19tohg38_flip.bim")
  bim2=fread(bim2file)
  bim2=as.data.frame(bim2)
  str11=paste0(bim1$V4,"_",bim1$V5,"_",bim1$V6)
  str12=paste0(bim1$V4,"_",bim1$V6,"_",bim1$V5)
  str2=paste0(bim2$V4,"_",bim2$V5,"_",bim2$V6)
  idx=str11 %in% str2 | str12 %in% str2
  traw1file=paste0("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter_flip_chr",chr,".traw")
  traw1=fread(traw1file)
  traw1=as.data.frame(traw1)
  tmp=paste0(traw1$CHR,":",traw1$POS,"_",traw1$COUNTED,"_",traw1$ALT)
  rownames(traw1)=tmp
  traw1=traw1[,7:ncol(traw1)]
  tmp=unlist(strsplit(colnames(traw1),"_"))
  colnames(traw1)=tmp[seq(1,length(tmp),2)]
  idx1=match(colnames(phenotype),colnames(traw1))
  if (sum(is.na(idx1))>0) warning(paste0("some RNAseq samples are missing in genotype data,",chr))
  traw1=traw1[,idx1]
  snp_=traw1[idx,]
  snppos_=data.frame(chr=bim1$V1[idx],pos=bim1$V4[idx],minor=bim1$V5[idx],major=bim1$V6[idx],stringsAsFactors = F)
  rownames(snppos_)=paste0(chr,":",snppos_$pos,"_",snppos_$minor,"_",snppos_$major)
  snp=rbind.data.frame(snp,snp_)
  snppos=rbind.data.frame(snppos,snppos_)
}

k <- which(is.na(snp), arr.ind=TRUE)
length(k)/nrow(snp)/ncol(snp)
snp[k] <- rowMeans(snp, na.rm=TRUE)[k[,1]]
all(rownames(snppos)==rownames(snp))
if (any(colnames(snp)!=colnames(phenotype))) warning("samples name not match!")
all(colnames(snp)==rownames(covariate))
save(snp,snppos,phenotype,phenotypepos,covariate,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")







