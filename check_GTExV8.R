library(data.table)

outfolderV8="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_May26"
outfolderV7="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_April18"


#load results----
#genemodel
load(paste0(outfolderV8,"/preidiction_michigan_model.RData"))
res_min_V8=res_min
load(paste0(outfolderV7,"/preidiction_michigan_model.RData"))
res_min_V7=res_min
#all data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadatafor_prediction.RData")
snpV8=snp
snpposV8=snppos
phenotypeV8=phenotype
phenotypeposV8=phenotypepos
covariateV8=covariate

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExmucosadatafor_prediction.RData")
snpV7=snp
snpposV7=snppos
phenotypeV7=as.data.frame(phenotype)
phenotypeposV7=phenotypepos
covariateV7=covariate

#compare gene models in V8 vs V7---
quantile(res_min_V8$numtotalsnp[res_min_V8$glmflag==1])
# 0%     25%     50%     75%    100% 
# 3.00 1274.25 1628.50 1984.00 7952.00 
quantile(res_min_V7$numtotalsnp[res_min_V7$glmflag==1])
# 0%     25%     50%     75%    100% 
# 6.0  1697.0  2117.5  2571.0 17212.0

#check sizes of common gene models---
comgenes=intersect(rownames(res_min_V7)[res_min_V7$glmflag==1],rownames(res_min_V8)[res_min_V8$glmflag==1])
idx1=match(comgenes,rownames(res_min_V7))
idx2=match(comgenes,rownames(res_min_V8))
table(res_min_V7$numselectedsnp[idx1]<res_min_V8$numselectedsnp[idx2])
#18% of models in V8 are larger
# FALSE  TRUE 
# 7630  1715
quantile(res_min_V7$numselectedsnp[idx1])
# 0%  25%  50%  75% 100% 
# 1   20   56  110  302
quantile(res_min_V8$numselectedsnp[idx2])
# 0%  25%  50%  75% 100% 
# 1    6   12   21  131 

quantile(res_min_V7$r2[idx1])
# 0%          25%          50%          75%         100% 
# 1.161863e-09 2.336548e-02 5.287927e-02 1.014025e-01 8.763441e-01
quantile(res_min_V8$r2[idx2])
# 0%          25%          50%          75%         100% 
# 3.222742e-08 1.345866e-02 3.443940e-02 7.614716e-02 7.496287e-01
quantile(res_min_V8$r2,na.rm=T)
# 0%          25%          50%          75%         100% 
# 1.328076e-08 1.034921e-02 2.720609e-02 6.030286e-02 7.903352e-01 

#genes to check, last three are for mucosa------
genes=c("HSP90AA1","UBAC1","ISYNA1","UPF1","THAP6","FOXF1")
#pick a gene to check
gene=genes[6] # an example V8 selected no SNPs
gene=genes[5] # and example V8 selected fewer SNPs

#gene model---
get_genemodel=function(gene="HSP90AA1",dat=res_min_V8)
{
  idx=which(rownames(dat)==gene)
  genemodel=dat[idx,]
  return(genemodel)
}

genemodelV8=get_genemodel(gene)
genemodelV8$numselectedsnp #22
genemodelV7=get_genemodel(gene,dat=res_min_V7)
genemodelV7$numselectedsnp #157



#check all genotype cis to a gene---
library(GenomicRanges)
get_cisSNP=function(phenotypepos=phenotypeposV8,snppos=snpposV8,dat=snpV8)
{
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2))
  gr_snp=GRanges(seqnames=snppos$chr,ranges=IRanges(start=snppos$pos,width=1))
  idx=which(rownames(phenotypepos)==gene)
  tmp=distance(gr_snp,gr_pos[idx])
  idx=which(tmp<5e5)
  return(dat[idx,])
}
cisSNPV8=get_cisSNP()
cisSNPV7=get_cisSNP(phenotypepos=phenotypeposV7,snppos=snpposV7,dat=snpV7)
readcorr=function(dat=t(cisSNPV8))
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  tmp=c(tmp)
  tmp=tmp[tmp!=0]
  return(tmp)
}
corV8=readcorr()
corV7=readcorr(dat=t(cisSNPV7))
quantile(abs(corV8))
# 0%          25%          50%          75%         100% 
# 4.041680e-21 3.347456e-02 7.219006e-02 1.453869e-01 1.000000e+00
quantile(abs(corV7))
# 0%          25%          50%          75%         100% 
# 6.437218e-22 3.490303e-02 7.336674e-02 1.506525e-01 1.000000e+00


#check selected genotype data---
get_genotype=function(dat=snpV8,genemodel=genemodelV8)
{
  selectedsenps=genemodel$selectedsnps
  selectedsenps=unlist(strsplit(selectedsenps,"|",fixed=T))
  idx=match(selectedsenps,rownames(dat))
  genotype=dat[idx,]
  return(genotype)
}

genotypeV8=get_genotype()
genotypeV7=get_genotype(dat=snpV7,genemodel=genemodelV7)

#hg38 --->hg19---
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
  newsnpnames1=hg38tohg19(snpnames=snpnames1)
  newsnpnames1=newsnpnames1[!is.na(newsnpnames1$snphg19),]
  newsnpnames1_hg19=newsnpnames1$snphg19
  sum(newsnpnames1_hg19 %in% snpnames2)
  #some thimes the above not working, they may sit on different strands (A/G and T/C), to use positon to match
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
#for selected snps
olapsnpnames=overlap_snp1_snp2()
#[1] "snp1: 22, snp2: 157, common snp: 10"
head(olapsnpnames$snpnames1,n=1)
#snphg38 is the original snpid on hg38 (V8), can be used to extract genotype data
# snphg38        snphg19   posnames overlap
# 1 4:75238084_C_T 4:76163294_C_T 4:76163294   FALSE
tmp=olapsnpnames$snpnames2
tmp$alt=tmp$ref=NA
for (i in 1:nrow(tmp))
{
  tmp1=unlist(strsplit(tmp$snphg19[i],"_"))
  tmp$alt[i]=tmp1[2]
  tmp$ref[i]=tmp1[3]
}
sum(tmp$ref %in% c("A","T") & tmp$alt %in% c("A","T"))+sum(tmp$ref %in% c("C","G") & tmp$alt %in% c("C","G"))
#[1] 19 are ambiguous SNPs can't be found in V8
head(olapsnpnames$snpnames2,n=1)
#snphg19 is the original snpid on hg19 (V7)
# snphg19   posnames overlap
# 1 4:76121093_T_G 4:76121093   FALSE

#overlap for all cisSNPs--
olapcissnpnames=overlap_snp1_snp2(snpnames1=rownames(cisSNPV8),snpnames2=rownames(cisSNPV7))
# [1] "snp1: 2158, snp2: 2857, common snp: 2126"

#extract overlap snp data, put V8 first,used to compare snp data---
extract_overlapsnpdat=function(snpnames1=rownames(genotypeV8),snpnames2=rownames(genotypeV7))
{
  olapsnpnames=overlap_snp1_snp2(snpnames1,snpnames2)
  idx=which(olapsnpnames$snpnames1$overlap==T)
  idx1=match(olapsnpnames$snpnames1$snphg38[idx],rownames(snpV8))
  snpdat1=snpV8[idx1,]
  rownames(snpdat1)=olapsnpnames$snpnames1$posnames[idx]
  idx=which(olapsnpnames$snpnames2$overlap==T)
  idx1=match(olapsnpnames$snpnames2$snphg19[idx],rownames(snpV7))
  snpdat2=snpV7[idx1,]
  rownames(snpdat2)=olapsnpnames$snpnames2$posnames[idx]
  idx=match(rownames(snpdat1),rownames(snpdat2))
  snpdat2=snpdat2[idx,]
  tmp=intersect(colnames(snpdat1),colnames(snpdat2))
  idx1=match(tmp,colnames(snpdat1))
  idx2=match(tmp,colnames(snpdat2))
  snpdat1=snpdat1[,idx1]
  snpdat2=snpdat2[,idx2]
  return(list(snpdat1=snpdat1,snpdat2=snpdat2))
}
allsnpdat=extract_overlapsnpdat()
#genotype data of V8 and V7 are in agreement---
sum(as.numeric(allsnpdat$snpdat1[1,])==as.numeric(allsnpdat$snpdat2[1,])) #T
which(as.numeric(allsnpdat$snpdat1[10,])!=as.numeric(allsnpdat$snpdat2[10,])) #1 difference (0.33 vs 1)

#check SKAT-TWAS----
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenesV7=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")
proteingenesV8=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
sum(proteingenesV7 %in% proteingenesV8)
load(paste0(outfolderV8,"/skat_res.RData"))
skat_minV8=skat_min2[rownames(skat_min2) %in% proteingenesV8,]
colnames(skat_minV7)=colnames(skat_minV8)=c("BE_p","EA_p","BEA_p","BEEA_p")
check_skattwas=function(skatres=skat_minV8)
{
  idx=which(rownames(skatres)==gene)
  print(skatres[idx,])
}
check_skattwas()
# BE_p      EA_p     BEA_p     BEEA_p
# THAP6 0.02151914 0.1159289 0.8319972 0.01879405
check_skattwas(skatres=skat_minV7)
# BE_p         EA_p     BEA_p       BEEA_p
# THAP6 0.0005539082 0.0003572668 0.5262568 6.060771e-05


#check geneexp---
get_geneexp=function(dat=phenotypeV8)
{
  idx=which(rownames(dat)==gene)
  geneexp=dat[idx,]
  return(geneexp)
}

geneexpV8=get_geneexp()
geneexpV7=get_geneexp(dat=phenotypeV7)
sum(colnames(geneexpV8) %in% colnames(geneexpV7)) #294
comsamples=intersect(colnames(geneexpV8),colnames(geneexpV7))
idx11=match(comsamples,colnames(geneexpV8))
idx22=match(comsamples,colnames(geneexpV7))
cor(as.numeric(geneexpV8[1,idx11]),as.numeric(geneexpV7[1,idx22])) #0.95 for THAP6,0.84 for FOXF1
plot(as.numeric(geneexpV8[1,idx11]),as.numeric(geneexpV7[1,idx22]))
idx=colnames(geneexpV8) %in% comsamples
boxplot(as.numeric(geneexpV8[1,])~idx)

#check cor across all the genes---
comgenes=intersect(rownames(phenotypeV7),rownames(phenotypeV8))
idx1=match(comgenes,rownames(phenotypeV8))
idx2=match(comgenes,rownames(phenotypeV7))
geneexpV8=phenotypeV8[idx1,idx11]
geneexpV7=phenotypeV7[idx2,idx22]
gneexpcor=rep(NA,length(comgenes))
for (i in 1:length(comgenes))
{
  gneexpcor[i]=cor(as.numeric(geneexpV8[i,]),as.numeric(geneexpV7[i,]))
}
quantile(gneexpcor,na.rm=T)
# 0%        25%        50%        75%       100% 
# -0.1136685  0.8614022  0.9190779  0.9487955  0.9926722
sum(abs(gneexpcor)<0.5,na.rm=T) #461
hist(gneexpcor)















