library("gdata")
library("data.table")
sampletable=read.csv("../data/bc_AT_JD1.csv",stringsAsFactors = F)
cambridgetable=sampletable[which(sampletable$site==30),c(2,1)]
cambridgetable$localid=as.integer(cambridgetable$localid)
write.table(cambridgetable,file="../result/cambridge_samples.txt",col.names = F,row.names = F,sep="\t") #for plink to select cambridge samples
fam=read.table("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015.fam",stringsAsFactors = F)
sum(cambridgetable$fam %in% fam$V1) #1939
sum(cambridgetable$localid %in% fam$V2) #1939
length(unique(cambridgetable$fam)) #1935
which(duplicated(cambridgetable$fam)) #[1]  321 1283 1605 1939
idx=match(cambridgetable$localid,fam$V2)
sum(fam$V3[idx]==0) #1938
which(fam$V3[idx]!=0)
fam[idx[1232],] #one father 302347 was included
#     V1     V2     V3 V4 V5 V6
#  3904 25 302286 302347  0  1 -9
sum(fam$V4[idx]==0) #1939
snps=c("rs11765529","rs7798662","rs4930068","rs2687201","rs10419226")
genotype=read.table("../result/cambridge_genotype.ped",sep="\t",stringsAsFactors = F,fill=T)
colnames(genotype)[7:11]=snps
geno2dosage=function(geno=genotype$V11)
{
  res=rep(NA,length(geno))
  res[geno=="A A"]=0
  res[geno=="A B" | geno=="B A"]=1
  res[geno=="B B"]=2
  res[geno=="0 0"]=NA
  return(res)
}
for (i in 7:11)
{
  genotype[,i]=geno2dosage(geno = genotype[,i])
}
idx=match(genotype$V2,sampletable$localid)
phenotype=sampletable[idx,]
table(phenotype$phenoBE_bc)
dat=cbind(phenotype,genotype[,7:11])
save(dat,file="../result/cambridge_dat.RData")

#task2----------------------

#update snp allel coding---
updateallele=fread("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/upall1.txt")
updateallele=as.data.frame(updateallele)
updateallele1=updateallele[!is.na(updateallele$V4) | !is.na(updateallele$V5),]
sum(bim$V2 %in% updateallele1$V1)
idx=match(bim$V2,updateallele$V1)
tmp=updateallele[idx,]
tmp=tmp[!is.na(tmp$V4),]
write.table(tmp,file="../result/updateallele.txt",col.names = F,row.names = F,sep="\t",quote=F)

#SNPs coding with A B after updating snp alleles
abbim=fread("../result/bca_filtered_30Nov2018.bim")
abbim=as.data.frame(abbim)
excludeABsnps=abbim$V2[abbim$V5=="B"]
excludeABsnps=c(excludeABsnps,abbim$V2[abbim$V6=="B"])
excludeABsnps=unique(excludeABsnps)
write.table(excludeABsnps,file="../result/excludeABsnps.txt",col.names = F,row.names = F,quote=F)

dong23snp=read.table("../data/Dong23snp.txt",header=T,sep="\t",stringsAsFactors = F)
bim=read.table("../result/bca_filtered_30Nov2018.bim",stringsAsFactors = F,sep="\t",fill=T)
sum(dong23snp$SNP %in% bim$V2) #8
dong23snp$genotype=0
for (i in 1:nrow(dong23snp))
{
  idx=which(bim$V1==dong23snp$Chr[i] & bim$V4==dong23snp$Position[i])
  if (length(idx)>0) dong23snp$genotype[i]=1
}

#check if the snps are in 1000genome
legendir="/fh/fast/stanford_j/Xiaoyu/Tools/impute_v2.3.2_x86_64_static/1000GP_Phase3/"
allchrs=unique(dong23snp$Chr)
dong23snp$g1000snp=NA
dong23snp$a0=NA
dong23snp$a1=NA
for(i in 1:length(allchrs))
{
  idx=which(dong23snp$Chr==allchrs[i])
  legend=fread(paste0("gunzip -cq ",legendir,"1000GP_Phase3_chr",allchrs[i],".legend.gz"))
  for (j in idx)
  {
    idx1=which(legend$position==dong23snp$Position[j] |legend$id==dong23snp$SNP[j] )
    if (length(idx1)>0)
    {
      dong23snp$g1000snp[j]=legend$id[idx1[1]]
      dong23snp$a0[j]=legend$a0[idx1[1]]
      dong23snp$a1[j]=legend$a1[idx1[1]]
    }
    if (length(idx1)>1) warnings(paste0("snp",j,"has multiple matches!"))
  }
}
sum(dong23snp$genotype==1 | !is.na(dong23snp$g1000snp)) #13
# impbim=read.table("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/Impute/bcaImpute_Filtered_19Feb2015.bim",stringsAsFactors = F,sep="\t",fill=T)
# sum(dong23snp$SNP %in% impbim$V2)

#select subsets of snps around those 23 for imputation
library(GenomicRanges)
subsetsnps=rep("",length(allchrs))
numsnps=rep(0,length(allchrs))
for (i in 1:length(allchrs))
{
  idx=which(dong23snp$Chr==allchrs[i])
  minidx=which(bim$V1==allchrs[i])[1]
  maxidx=which(bim$V1==allchrs[i])[sum(bim$V1==allchrs[i])]
  startidxs=endidxs=rep(0,length(idx))
  gr=GRanges()
  for (j in 1:length(idx))
  {
    idx1=which(bim$V1==allchrs[i] & bim$V4<=dong23snp$Position[idx[j]]-5e6)
    if (length(idx1)>0)
    {
      startidxs[j]=max(idx1)
    }else
    {
      startidxs[j]=minidx
    }
    idx1=which(bim$V1==allchrs[i] & bim$V4>=dong23snp$Position[idx[j]]+5e6)
    if (length(idx1)>0)
    {
      endidxs[j]=min(idx1)
    }else
    {
      endidxs[j]=maxidx
    }
    gr=append(gr,GRanges(seqnames = allchrs[i],ranges=IRanges(start=startidxs[j],end=endidxs[j])))
  }
  gr=reduce(gr)
  tmp=rep(NA,length(gr))
  for(j in 1:length(tmp))
  {
    tmp[j]=paste0(bim$V2[start(gr)[j]],"-",bim$V2[end(gr)[j]])
  }
  numsnps[i]=sum(end(gr)-start(gr))
  subsetsnps[i]=paste0(tmp,collapse = ",")
}
write.table(data.frame(snps=subsetsnps),file="../result/subsetsnp_5M.txt",row.names = F,col.names = F,quote=F)

# snpanno=read.csv("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/BC_SNP_annotation.csv")
# sum(bim$V2 %in% snpanno$rsID)
# bim$V2[!bim$V2 %in% snpanno$rsID]

#after imputation
impdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_snp23_1000g/"
allchrs=unique(dong23snp$Chr)
dong23snp$imp_g1000=NA
dong23snp$imp_g1000_r2=NA
for(i in 1:length(allchrs))
{
  idx=which(dong23snp$Chr==allchrs[i])
  imp=fread(paste0("gunzip -cq ",impdir,"chr",allchrs[i],".info.gz"))
  tmp=unlist(strsplit(imp$SNP,":"))
  imp$position=as.integer(tmp[seq(2,length(tmp),2)])
  for (j in idx)
  {
    idx1=which(imp$position==dong23snp$Position[j] )
    if (length(idx1)>0)
    {
      dong23snp$imp_g1000[j]=imp$SNP[idx1[1]]
      dong23snp$imp_g1000_r2[j]=as.numeric(imp$Rsq[idx1[1]])
    }
    if (length(idx1)>1) warnings(paste0("snp",j,"has multiple matches!"))
  }
}

impdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_snp23_phase1/"
allchrs=unique(dong23snp$Chr)
dong23snp$imp_phase1=NA
dong23snp$imp_phase1_r2=NA
for(i in 1:length(allchrs))
{
  idx=which(dong23snp$Chr==allchrs[i])
  imp=fread(paste0("gunzip -cq ",impdir,"chr",allchrs[i],".info.gz"))
  tmp=unlist(strsplit(imp$SNP,":"))
  imp$position=as.integer(tmp[seq(2,length(tmp),2)])
  for (j in idx)
  {
    idx1=which(imp$position==dong23snp$Position[j] )
    if (length(idx1)>0)
    {
      dong23snp$imp_phase1[j]=imp$SNP[idx1[1]]
      dong23snp$imp_phase1_r2[j]=as.numeric(imp$Rsq[idx1[1]])
    }
    if (length(idx1)>1) warnings(paste0("snp",j,"has multiple matches!"))
  }
}

impdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_snp23_hrc/"
allchrs=unique(dong23snp$Chr)
dong23snp$imp_hrc=NA
dong23snp$imp_hrc_r2=NA
for(i in 1:length(allchrs))
{
  idx=which(dong23snp$Chr==allchrs[i])
  imp=fread(paste0("gunzip -cq ",impdir,"chr",allchrs[i],".info.gz"))
  tmp=unlist(strsplit(imp$SNP,":"))
  imp$position=as.integer(tmp[seq(2,length(tmp),2)])
  for (j in idx)
  {
    idx1=which(imp$position==dong23snp$Position[j] )
    if (length(idx1)>0)
    {
      dong23snp$imp_hrc[j]=imp$SNP[idx1[1]]
      dong23snp$imp_hrc_r2[j]=as.numeric(imp$Rsq[idx1[1]])
    }
    if (length(idx1)>1) warnings(paste0("snp",j,"has multiple matches!"))
  }
}

idx1=which(!is.na(dong23snp$imp_g1000))
idx2=which(!is.na(dong23snp$imp_phase1))
idx3=which(!is.na(dong23snp$imp_hrc))
tmp=c(dong23snp$imp_g1000_r2[idx1],dong23snp$imp_phase1_r2[idx2],dong23snp$imp_hrc_r2[idx3])
boxplot(outline=F,tmp~c(rep("phase3",21),rep("phase1",23),rep("HRC",18)),ylab="quality score",cex.axis=1.3,cex.lab=1.3,ylim=c(0.65,1))
stripchart(tmp~c(rep("phase3",21),rep("phase1",23),rep("HRC",18)),vertical = TRUE, pch=1,method = "jitter",add = TRUE)
write.table(dong23snp,file="../result/Dong23SNPs_impscore.txt",row.names = F,col.names = T,sep="\t",quote=F)
write.table(data.frame(snp=paste0(dong23snp$Chr,":",dong23snp$Position,"-",dong23snp$Position),chr=dong23snp$Chr),file="../result/Dong23SNPs_postion.txt",row.names = F,col.names = F,quote=F)
dong23snppos=read.table("../result/Dong23SNPs_postion.txt",stringsAsFactors = F)
addsnps=c("rs11765529", "rs4930068", "rs2687201")
for (i in 1:length(addsnps))
{
  cmd=paste0("grep -w ",addsnps[i]," /fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt")
  system(cmd)
}
# chr7	52922212	52922213	rs11765529	+	single	C,T,
# chr11	2297614	2297615	rs4930068	+	single	A,C,
# chr3	70928929	70928930	rs2687201	+	single	A,C,
tmp_chr=c(7,11,3)
tmp_pos=c(52922213,2297615,70928930)
tmp=data.frame(V1=paste0(tmp_chr,":",tmp_pos,"-",tmp_pos),V2=tmp_chr)
# dong26snppos=rbind(dong23snppos,tmp)
# write.table(dong26snppos,file="../result/Dong26SNPs_postion.txt",row.names = F,col.names = F,sep="\t",quote=F)
write.table(tmp,file="../result/Extra3SNPs_postion.txt",row.names = F,col.names = F,sep="\t",quote=F)
tmp=data.frame(SNP=addsnps,Chr=tmp_chr,Position=tmp_pos,stringsAsFactors = F)
write.table(tmp,file="../result/Extra3SNPs.txt",row.names = F,col.names = T,sep="\t",quote=F)

#extract genotype data
#base on phase1,generated by extract_genotype.sh
impfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Dong23SNPs_genotype.txt"
tmp=as.data.frame(fread(impfile,sep="\t"))
tmp1=data.frame(matrix(NA,nrow=9228+4,ncol=23))
colnames(tmp1)=dong23snp$SNP
tmp1[1,]=tmp[,1]
tmp1[2,]=tmp[,2]
tmp1[3,]=tmp[,4]
tmp1[4,]=tmp[,5]
for (i in 1:23)
{
  idx=which(grepl("0/0",tmp[i,10:ncol(tmp)]))
  tmp1[idx+4,i]=0
  idx=which(grepl("0/1",tmp[i,10:ncol(tmp)]))
  tmp1[idx+4,i]=1
  idx=which(grepl("1/1",tmp[i,10:ncol(tmp)]))
  tmp1[idx+4,i]=2
}
rownames(tmp1)[1:4]=c("chr","position","ref","alt")

tmp=fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_snp23_phase1/chr9.dose.vcf_head.txt",skip=13,nrow=1)
rownames(tmp1)[5:nrow(tmp1)]=colnames(tmp)[10:ncol(tmp)]
effect_alt=rep(NA,23) #if effect allele is the alt allele (T)
sum(tmp1[3,] %in% c("A","T") & tmp1[4,] %in% c("A","T")) #0
sum(tmp1[3,] %in% c("C","G") & tmp1[4,] %in% c("C","G")) #2, need to take extra steps for these 2
for (i in 1:23)
{
  if (dong23snp$Effect.allele[i] %in% c("T","A") & tmp1[3,i] %in% c("T","A"))
  {
    effect_alt[i]=F
  }
  if (dong23snp$Effect.allele[i] %in% c("C","G") & tmp1[3,i] %in% c("C","G"))
  {
    effect_alt[i]=F
  }
  if (dong23snp$Effect.allele[i] %in% c("T","A") & tmp1[4,i] %in% c("T","A"))
  {
    effect_alt[i]=T
  }
  if (dong23snp$Effect.allele[i] %in% c("C","G") & tmp1[4,i] %in% c("C","G"))
  {
    effect_alt[i]=T
  }
}
idx=which(tmp1[3,] %in% c("C","G") & tmp1[4,] %in% c("C","G"))
for (i in idx)
{
  if (dong23snp$Effect.allele[i] ==tmp1[3,i])
  {
    effect_alt[i]=F
  }else
  {
    effect_alt[i]=T
  }
}
table(effect_alt)
# FALSE  TRUE 
# 7    16 
tmp1$PRS=NA
for (i in 5:nrow(tmp1))
{
  if (i %%1000==0) cat(i,"..")
  #compute numer of effect alleles
  tmp2=rep(0,23)
  idx=tmp1[i,1:23]==1
  tmp2[idx]=1
  idx=effect_alt==T
  idx=tmp1[i,1:23]==2 & effect_alt==T
  tmp2[idx]=2
  idx=tmp1[i,1:23]==0 & effect_alt==F
  tmp2[idx]=2
  tmp1$PRS[i]=sum(log(dong23snp$Published.OR)*tmp2)/23
}
#merge with phenotype data
sampleid=paste0(sampletable$fam,"_",sampletable$localid)
idx=match(rownames(tmp1)[5:nrow(tmp1)],sampleid)
missingsamples=sampleid[!sampleid %in% rownames(tmp1)[5:nrow(tmp1)]]
idx=match(missingsamples,sampleid)
sampletable[idx,1:8]
# localid       fam phenoBE_bc phenoEA_bc phenoEABE_bc sex age site
# 2211  151021 113321423          1          1            1   1  69   15
# 3016  171295 111273007         -9         -9           -9   1  61   17
# 7000  302106 115275511          2         -9            2   2  73   30
availablesamples=intersect(sampleid,rownames(tmp1)[5:nrow(tmp1)])
allgenotypedat=tmp1
idx=match(availablesamples,rownames(allgenotypedat))
dat=allgenotypedat[idx,]
idx=match(availablesamples,sampleid)
dat=cbind(dat,sampletable[idx,])
#PCA
eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca"
eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind"
PC=read.table(eigfile,skip=nskip,stringsAsFactors = F)
colnames(PC)=paste0("pc",1:ncol(PC))
eigsamples=read.table(eigsampfile,stringsAsFactors = F)
eigsamples=paste0(eigsamples$V1,"_",eigsamples$V2)
rownames(PC)=eigsamples

#add two PRS
all(dong23snp$SNP==colnames(allgenotypedat)[1:23]) #T
allgenotypedat$PRS_BE=allgenotypedat$PRS_EA=NA
rmsnps=c("rs2687202", "rs17749155", "rs10108511", "rs7255")
idxkeep1=!dong23snp$SNP %in% rmsnps
rmsnps1=c("rs17749155", "rs10108511")
idxkeep2=!dong23snp$SNP %in% rmsnps1
for (i in 5:nrow(allgenotypedat))
{
  if (i %%1000==0) cat(i,"..")
  #compute numer of effect alleles
  tmp2=rep(0,23)
  idx=allgenotypedat[i,1:23]==1
  tmp2[idx]=1
  idx=effect_alt==T
  idx=allgenotypedat[i,1:23]==2 & effect_alt==T
  tmp2[idx]=2
  idx=allgenotypedat[i,1:23]==0 & effect_alt==F
  tmp2[idx]=2
  allgenotypedat$PRS_BE[i]=sum(log(dong23snp$Published.OR[idxkeep1])*tmp2[idxkeep1])/19
  allgenotypedat$PRS_EA[i]=sum(log(dong23snp$Published.OR[idxkeep2])*tmp2[idxkeep2])/21
}

#add rs11765529, rs4930068, rs2687201
addsnps=c("rs11765529", "rs4930068", "rs2687201")
snps=c("rs11765529","rs7798662","rs4930068","rs2687201","rs10419226")
genotype=read.table("../result/cambridgeall_genotype.ped",sep="\t",stringsAsFactors = F,fill=T)
colnames(genotype)[7:11]=snps

for (i in 7:11)
{
  genotype[,i]=geno2dosage(geno = genotype[,i])
}
idx=match(addsnps,colnames(genotype))
idx1=match(rownames(allgenotypedat),paste0(genotype$V1,"_",genotype$V2))
tmp=genotype[idx1,idx]
allgenotypedat=cbind(allgenotypedat,tmp)
idx=match(rownames(dat),rownames(allgenotypedat))
dat=cbind(dat,allgenotypedat[idx,25:29])
save(missingsamples,PC,allgenotypedat,dat,dong23snp,file="../result/Dong23SNPs.RData")


#task3----------------------
#check samples included in each chrom, some samples failed --mind 0.05 --geno 0.05 --hwe 0.000001 --maf 0.05
allsamples=comsamples=NULL
for (i in 1:23)
{
  tmp=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_19Dec2018_chr",i,".fam"),stringsAsFactors = F)
  cat(i)
  print(nrow(tmp))
  if (is.null(comsamples))
  {
    comsamples=allsamples=tmp$V2
  }else
  {
    allsamples=unique(c(allsamples,tmp$V2))
    comsamples=intersect(comsamples,tmp$V2)
  }
}
removedsamples=allsamples[!allsamples %in% comsamples]
removedsamples1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_michigan_imp_statistics_lowcallratesamples.txt",stringsAsFactors = F)
sum(removedsamples %in% removedsamples1$V2)
tmp=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_19Dec2018_chr",22,".fam"),stringsAsFactors = F)
idx=match(removedsamples,tmp$V2)
removedsamples2=data.frame(V1=tmp$V1[idx],V2=tmp$V2[idx])
removedsamples=rbind(removedsamples1,removedsamples2)
write.table(removedsamples,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_michigan_imp_removedsamples.txt",row.names = F,col.names = F,sep="\t",quote=F)

