#to check numbers in the flowchart

# wc -l ../data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.bim
# 46569704 ../data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.bim
#check bca genotypes

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  if (length(idx)>0)
    sampletable[idx,i]=NA
}

#add covariate table (pc1-pc4)
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

impdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"
dataset="cambridgewtccc_qc_hrc_maf005_snp"
famfile=paste0(impdir,dataset,"/","chr7_filter_hg19tohg38_flip.fam")
fam=read.table(famfile)
fam$V2=gsub("SEP","",fam$V2)
all(fam$V2 %in% rownames(covariatetable))
idx=match(fam$V2,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
#   1    2 
# 3408  881 
table(covariatetable$phenoEA_bca[idx])
#   1    2 
# 3408 1003 
table(covariatetable$phenoEABE_bca[idx])
#   1    2 
# 3408 1884
cambridge_count=0
for (chr in 1:22)
{
  tmp=fread(paste0(impdir,dataset,"/","chr",chr,"_filter_hg19tohg38_flip.bim"))
  cambridge_count=cambridge_count+nrow(tmp)
}

dataset="beacondbgapcontrol_qc_hrc_maf005_snp"
famfile=paste0(impdir,dataset,"/","chr7_filter_hg19tohg38_flip.fam")
fam=read.table(famfile)
fam$V2=gsub("SEP","",fam$V2)
all(fam$V2 %in% rownames(covariatetable))
idx=match(fam$V2,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
#  1    2 
# 6700 2401 
table(covariatetable$phenoEA_bca[idx])
#  1    2 
# 6701 1507 
table(covariatetable$phenoEABE_bca[idx])
#   1    2 
# 6701 3908 
beacon_count=0
for (chr in 1:22)
{
  tmp=fread(paste0(impdir,dataset,"/","chr",chr,"_filter_hg19tohg38_flip.bim"))
  beacon_count=beacon_count+nrow(tmp)
}

dataset="merge_beacon_cambridge_hrc_maf005_snp"
famfile=paste0(impdir,dataset,"/","chr7_filter_hg19tohg38_flip.fam")
fam=read.table(famfile)
fam$V2=gsub("SEP","",fam$V2)
all(fam$V2 %in% rownames(covariatetable))
idx=match(fam$V2,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
#   1     2 
# 10108  3282 
table(covariatetable$phenoEA_bca[idx])
#   1     2 
# 10109  2510 
table(covariatetable$phenoEABE_bca[idx])
#   1     2 
# 10109  5792 
merge_count=0
for (chr in 1:22)
{
  tmp=fread(paste0(impdir,dataset,"/","chr",chr,"_filter_hg19tohg38_flip.bim"))
  merge_count=merge_count+nrow(tmp)
}