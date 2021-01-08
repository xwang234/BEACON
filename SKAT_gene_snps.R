#!/usr/bin/env Rscript

#The code is used to label SNPs within 50KB of a gene
library(data.table)
ucsc_refseq=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/ucsc_refseqgenes.txt",header = T,stringsAsFactors = F,comment.char = "")
ucsc_refseq=ucsc_refseq[ucsc_refseq$cdsStartStat=="cmpl",]
knowngenes=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",header=F,stringsAsFactors = F,sep="\t")
knongenes_utr3=read.table("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/UCSC_gene_FiveUTR.bed")
#The 5′ UTR begins at the transcription start site
#the three prime untranslated region (3'-UTR) is the section of messenger RNA (mRNA) that immediately follows the translation termination codon
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
gtex_ge_anno=gtex_ge_anno[gtex_ge_anno$V3=="gene" & gtex_ge_anno$gene_type=="protein_coding",]
gtex_ge_anno=gtex_ge_anno[gtex_ge_anno$Chromosome %in% c(1:22),]
gtex_ge_utr5=read.table("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/GENCODE27_gene_fiveUTR.bed",stringsAsFactors = F)
gtex_ge_utr3=read.table("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/GENCODE27_gene_threeUTR.bed",stringsAsFactors = F)
idx=which(duplicated(gtex_ge_anno$Symbol))
idx=which(gtex_ge_anno$Symbol %in% gtex_ge_anno$Symbol[idx])
View(gtex_ge_anno[idx,])
idx=duplicated(gtex_ge_anno$Symbol)
gtex_ge_anno=gtex_ge_anno[!idx,]
idx=order(gtex_ge_anno$Chromosome,gtex_ge_anno$start)
gtex_ge_anno=gtex_ge_anno[idx,]
library(GenomicRanges)
gr_gene=GRanges(seqnames = gtex_ge_anno$Chromosome,ranges=IRanges(start=gtex_ge_anno$start-5e4,end=gtex_ge_anno$end+5e4))
bimdat=as.data.frame(fread("../../BeagessCambridgeAmos/bca_filteredMAF_20Feb2015.bim"))
gr_bim=GRanges(seqnames = bimdat$V1,ranges = IRanges(start=bimdat$V4,width=1))
genesnptable=data.frame(gene=gtex_ge_anno$Symbol,snp=NA,numsnp=0,stringsAsFactors = F)
snpgenetable=data.frame(snp=bimdat$V2,gene="_",numgene=0,stringsAsFactors = F)
for (i in 1:nrow(genesnptable))
{
  tmp=distance(gr_gene[i],gr_bim)
  idx=which(tmp==0)
  if (length(idx)>0)
  {
    if (i %% 500==0) cat(i,'..')
    genesnptable$snp[i]=paste0(bimdat$V2[idx],collapse=",")
    snpgenetable$gene[idx]=paste0(snpgenetable$gene[idx],",",gtex_ge_anno$Symbol[i])
    genesnptable$numsnp[i]=length(idx)
    snpgenetable$numgene[idx]=snpgenetable$numgene[idx]+1
  }
}
genesnptable=genesnptable[!is.na(genesnptable$snp),]
snpgenetable$gene=gsub("_,","",snpgenetable$gene)
snpgenetable$gene=gsub("_","",snpgenetable$gene)
all(bimdat$V2==snpgenetable$snp) #T
idx=which(snpgenetable$gene=="")
snpgenetable$gene[idx]=NA
snpgenetable$chr=bimdat$V1
snpgenetable$position=bimdat$V4
save(genesnptable,snpgenetable,file="../result/skat_gene_snp_bca_filteredMAF_19Feb2015.RData")
quantile(genesnptable$numsnp)

quantile(snpgenetable$numgene[snpgenetable$numgene>0])
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

for (l in 1:nrow(genesnptable))
{
  if (l %% 1000==0) cat(l,'..')
  #if (!is.na(genesnptable[l,2]))
  #{
    selectedsnps <- unlist(strsplit(genesnptable[l,2],",",fixed=T))   
    idx=match(selectedsnps,rownames(allgenotype))
    #selectedsnps[idx]
    if (sum(is.na(idx))>0)
    {
      print(l)
    }
  #}
}

