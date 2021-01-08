#!/usr/bin/env Rscript
library(data.table)
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log10)",
       ylab="Observed p-value (log10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

getrsid=function(oldid=mqtl_23highrisk_cis$SNP)
{
  res=sapply(oldid,function(x){
    unlist(strsplit(x,":"))[1]
  })
  return(res)
}

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4" #junction
outfolder2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTExmucosa" #mucosa

#compare skat
compare_2col=function(dat1=skat_min1,dat2=skat_min2,n=1,cutoff=0.05)
{
  res=NA
  tmp1=p.adjust(dat1[,n],method="fdr")
  tmp2=p.adjust(dat2[,n],method="fdr")
  idx1=which(tmp1<cutoff)
  idx2=which(tmp2<cutoff)
  genes=unique(c(rownames(dat1)[idx1],rownames(dat2)[idx2]))
  if (length(genes)>0)
  {
    res=data.frame(p1=rep(NA,length(genes)),p2=rep(NA,length(genes)),fdr1=rep(NA,length(genes)),fdr2=rep(NA,length(genes)),stringsAsFactors = F)
    rownames(res)=genes
    for (i in 1:length(genes))
    {
      idx1=which(rownames(dat1)==genes[i])
      if (length(idx1)>0)
      {
        res$p1[i]=dat1[idx1,n]
        res$fdr1[i]=tmp1[idx1]
      }
      idx2=which(rownames(dat2)==genes[i])
      if (length(idx2)>0)
      {
        res$p2[i]=dat2[idx2,n]
        res$fdr2[i]=tmp2[idx2]
      }
    }
  }
  return(res)
}
compare_skat=function(outfolder1,outfolder2,cutoff=0.05,
                      eqtlfile1="/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtex_23highrisk_cn_mutation_cis",
                      eqtlfile2="/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtex_mucosa_23highrisk_cn_mutation_cis")
{
  load(paste0(outfolder1,"/skat_res.RData"))
  skat_min1=skat_min[rownames(skat_min) %in% proteingenes,]
  load(paste0(outfolder1,"/preidiction_michigan_model.RData"))
  model1=res_min[rownames(res_min) %in% proteingenes,]
  load(paste0(outfolder2,"/skat_res.RData"))
  skat_min2=skat_min[rownames(skat_min) %in% proteingenes,]
  sum(rownames(skat_min1) %in% rownames(skat_min2)) #6316
  load(paste0(outfolder2,"/preidiction_michigan_model.RData"))
  model2=res_min[rownames(res_min) %in% proteingenes,]
  #compare BE
  BE=compare_2col(dat1=skat_min1,dat2=skat_min2,n=1,cutoff = cutoff)
  EA=compare_2col(dat1=skat_min1,dat2=skat_min2,n=2,cutoff = cutoff)
  BEA=compare_2col(dat1=skat_min1,dat2=skat_min2,n=3,cutoff = cutoff)
  BEEA=compare_2col(dat1=skat_min1,dat2=skat_min2,n=4,cutoff = cutoff)
  eqtl1=fread(eqtlfile1,header = T)
  eqtl1=eqtl1[eqtl1$gene %in% proteingenes,]
  eqtl1$SNP=getrsid(eqtl1$SNP)
  eqtl2=fread(eqtlfile2,header = T)
  eqtl2=eqtl2[eqtl2$gene %in% proteingenes,]
  eqtl2$SNP=getrsid(eqtl2$SNP)
  pair1=paste0(eqtl1$SNP,":",eqtl1$gene)
  pair2=paste0(eqtl2$SNP,":",eqtl2$gene)
  idx1=which(eqtl1$qvalue<cutoff)
  idx2=which(eqtl2$qvalue<cutoff)
  pair=unique(c(pair1[idx1],pair2[idx2]))
  EQTL=NA
  if (length(pair)>0)
  {
    tmp=unlist(strsplit(pair,":"))
    EQTL=data.frame(snp=tmp[seq(1,length(tmp),2)],gene=tmp[seq(2,length(tmp),2)],p1=rep(NA,length(pair)),p2=rep(NA,length(pair)),fdr1=rep(NA,length(pair)),fdr2=rep(NA,length(pair)))
    for (i in 1:length(pair))
    {
      idx1=which(pair1==pair[i])
      if (length(idx1)>0)
      {
        EQTL$p1[i]=eqtl1$`p-value`[idx1]
        EQTL$fdr1[i]=eqtl1$qvalue[idx1]
      }
      idx2=which(pair2==pair[i])
      if (length(idx2)>0)
      {
        EQTL$p2[i]=eqtl2$`p-value`[idx2]
        EQTL$fdr2[i]=eqtl2$qvalue[idx2]
      }
    }
  }
  
  return(list(BE=BE,EA=EA,BEA=BEA,BEEA=BEEA,EQTL=EQTL))
}

compare_skat_res=compare_skat(outfolder1 = outfolder1,outfolder2 = outfolder2)

