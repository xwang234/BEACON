
#!/usr/bin/env Rscript
rm(list=ls())
library(data.table)
library(GenomicRanges)
library(PMA)

removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}

removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  dat=dat[!idxconst,]
}
processsemicolon=function(txt="5'UTR;1stExon;1stExon;1stExon;5'UTR;5'UTR;5'UTR;1stExon;5'UTR")
{
  txt=as.character(txt)
  res=txt
  for (i in 1:length(txt))
  {
    if (grepl(";",txt[i]))
    {
      tmp=unlist(strsplit(txt[i],";"))
      res[i]=paste0(unique(tmp),collapse="|")
    }
  }
  
  return(res)
}

#prepare data matrice------
load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp_methylation.RData")
#gene transcript position, not used in mwas
phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
phenotypepos=fread(phenotypeposfile,header=T,sep="\t")
phenotypepos=as.data.frame(phenotypepos)
rownames(phenotypepos)=phenotypepos[,1]

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr) 
anno$UCSC_RefGene_Name=processsemicolon(anno$UCSC_RefGene_Name)
anno$UCSC_RefGene_Group=processsemicolon(anno$UCSC_RefGene_Group)
#order anno
anno$chr=factor(anno$chr,levels = c(1:22,"X","Y"))
idx=order(anno$chr,anno$pos)
anno=as.data.frame(anno)
anno=anno[idx,]
idx=match(rownames(methydata),anno$Name)
methydata=methydata[!is.na(idx),]
idx=match(rownames(methydata),anno$Name)
gr_meth=GRanges(seqnames = anno$chr[idx],ranges=IRanges(start=anno$pos[idx],width=1))
#gr_geneexp=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp

#try genotyped SNP
snp6_anno <- fread(file="/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",skip=18,header=T,sep=",",stringsAsFactors=F)
snp6_anno=as.data.frame(snp6_anno)
snp6_anno$`Physical Position`=as.numeric(snp6_anno$`Physical Position`)

#get gene symbol of SNP6
# snp6_anno$gene=NA
# library("biomaRt")
# mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
# tmp=listAttributes(mart=mart)
# for (i in 1:nrow(snp6_anno))
# {
#   tmp=unlist(strsplit(snp6_anno$`Associated Gene`[i]," "))[1]
#   if (i %% 1000==0) cat(i,'..')
#   if (tmp != "---")
#   {
#     my_gene <- getBM(attributes='hgnc_symbol',
#                        filters = 'ensembl_transcript_id',
#                        values =tmp,
#                        mart = mart)
#     if (nrow(my_gene)>0)
#     {
#       snp6_anno$gene[i]=my_gene[1,1]
#     }
#   }
# }
# save(snp6_anno,file="../result/snp6_anno.RData")

all(rownames(genotypedata) %in% snp6_anno$`Probe Set ID`) #T
idx=match(rownames(genotypedata),snp6_anno$`Probe Set ID`)
genotypedata=genotypedata[snp6_anno$Chromosome[idx] %in% c(1:22,"X","Y"),]
idx=match(rownames(genotypedata),snp6_anno$`Probe Set ID`)
gr_snp=GRanges(seqnames = snp6_anno$Chromosome[idx], ranges = IRanges(start=snp6_anno$`Physical Position`[idx],width=1))

idx=match(colnames(genotypedata),colnames(methydata))
methydata=methydata[,idx]
all(colnames(genotypedata)==colnames(methydata))
#only include EAC
tmp=tcgaclinical$bcr_patient_barcode[tcgaclinical$histological_type=="Esophagus Adenocarcinoma, NOS"]
genotypedata=genotypedata[,colnames(genotypedata) %in% tmp]
methydata=methydata[,colnames(methydata) %in% tmp]
idx=match(colnames(methydata),colnames(geneexpdata))
geneexpdata=geneexpdata[,idx]

##CCA analysis----
genes=c("MEF2B","DDX49","KXD1","CERS1","COMP","ARMC6","SGK223","TMEM161A","CRTC1","HOMER3")
#pick a gene
gene=genes[2]
#idx=which(phenotypepos$geneid==gene)
#chr=phenotypepos$chr[idx]
# startloc=phenotypepos$s1[idx]-5e5
# startloc=max(1,startloc)
# endloc=phenotypepos$s2[idx]+5e5

#to select CpGs associated with the gene in annoation
idx <- grep(gene,anno$UCSC_RefGene_Name)


#idx=unique(c(which(anno$UCSC_RefGene_Name==gene),
#             grep(paste0("^",gene,";"),anno$UCSC_RefGene_Name),
#             grep(paste0(";",gene,";"),anno$UCSC_RefGene_Name),
#             grep(paste0(";",gene,"$"),anno$UCSC_RefGene_Name)))
chr=anno$chr[idx[1]]
#to select GpGs in nearby intergenic region
tmp=anno$UCSC_RefGene_Name[(min(idx)-1):(min(idx)-500)]
len=0
for (i in 1:500)
{
  if (tmp[i]=="")
  {
    len=len+1
  }else
  {
    break
  }
}
startcpgidx=min(idx)-len
tmp=anno$UCSC_RefGene_Name[(max(idx)+1):(max(idx)+500)]
len=0
for (i in 1:500)
{
  if (tmp[i]=="")
  {
    len=len+1
  }else
  {
    break
  }
}
endcpgidx=max(idx)+len
startcpgloc=anno$pos[startcpgidx]
endcpgloc=anno$pos[endcpgidx]
#add 500KB to both ends
startsnploc=max(1,startcpgloc-5e5)
endsnploc=endcpgloc+5e5
gr_snpregion=GRanges(seqnames = chr,ranges = IRanges(start=startsnploc,end=endsnploc))
tmp=distance(gr_snp,gr_snpregion)
sum(tmp==0,na.rm=T) #161 snps
idx=which(tmp==0)
CCAX=genotypedata[idx,]
CCAX=removeconstrows(CCAX)

gr_methregion=GRanges(seqnames = chr,ranges=IRanges(start=startcpgloc,end=endcpgloc))
tmp=distance(gr_meth,gr_methregion)
sum(tmp==0,na.rm=T) #46 CpGs
idx=which(tmp==0)
CCAZ=methydata[idx,,drop=F]
#check NA, Z should have no NA
# tmp=sapply(1:ncol(CCAZ),function(x){
#   sum(is.na(CCAZ[,x]))
# })
#remove CpGs having NA
tmp=sapply(1:nrow(CCAZ),function(x){
  sum(is.na(CCAZ[x,]))
})
CCAZ=CCAZ[tmp==0,]

#order CCAX and CCAZ
idx=match(rownames(CCAX),snp6_anno$`Probe Set ID`)
pos=snp6_anno$`Physical Position`[idx]
tmp=sort(pos,index.return=T)
CCAX=CCAX[tmp$ix,]
idx=match(rownames(CCAZ),anno$Name)
pos=anno$pos[idx]
tmp=sort(pos,index.return=T)
CCAZ=CCAZ[tmp$ix,]

#run CCA
#Sys.time()
#no permutation needed for ordered data
# remove the redundant SNPs 
CCAX <- t(removehighcorr(t(CCAX)))
perm.out <- CCA.permute(t(CCAX),t(CCAZ),typex="standard",typez="standard",nperms=1000,trace=F,standardize = F)
#Sys.time()
CCAout <- CCA(t(CCAX),t(CCAZ),typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init,standardize=F)
#CCAout <- CCA(t(CCAX),t(CCAZ),typex="ordered",typez="ordered",K=1)


pcpg <- drop(t(CCAZ)%*% CCAout$v)
gid <- grep(gene,row.names(geneexpdata))[1]
plot(pcpg,geneexpdata[gid,])
summary(glm(as.numeric(geneexpdata[gid,])~pcpg))

psnp <- drop(t(CCAX)%*%CCAout$u)
plot(psnp,log(geneexpdata[gid,]))
summary(glm(as.numeric(geneexpdata[gid,])~psnp))

cor(psnp,pcpg)
cor(pcpg,as.numeric(geneexpdata[gid,]))
#check result
table(CCAout$u!=0)
# FALSE  TRUE 
#150  7
table(CCAout$v!=0)
# FALSE  TRUE 
# 28    12

az <- cor(t(CCAZ[CCAout$v!=0,]),t(CCAX[CCAout$u!=0,]))

tmp=rownames(CCAX)[CCAout$u!=0]
idx=match(tmp,snp6_anno$`Probe Set ID`)
selsnp=snp6_anno$`dbSNP RS ID`[idx]
#"rs4641872"  "rs1363119"  "rs6512257"  "rs9636202"  "rs7251764"  "rs10854166" "rs8103648"
#table(snp6_anno$gene[idx])
selsnppos=snp6_anno$`Physical Position`[idx]
#[1] 18440942 18444809 18447551 18449238 18648860 18652844 18655434

gwassnp=c("rs10419226","rs10423674")
idx=match(gwassnp,snp6_anno$`dbSNP RS ID`) #NA NA
gwaspos=c(18803172,18817903)

selCpG=rownames(CCAZ)[CCAout$v!=0]
idx=match(selCpG,anno$Name)
selCpGpos=anno$pos[idx]
table(anno$Regulatory_Feature_Group[idx])
table(anno$UCSC_RefGene_Group[idx])
table(anno$UCSC_RefGene_Name[idx])


