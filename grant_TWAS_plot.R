#!/usr/bin/env Rscript

qqplotnew=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL,genes=rownames(skat_min_code),legend=NULL)
{
  library(RColorBrewer)
  
  par("mar"=c(6,5.5,5,3.5))
  idx=order(pvalue)
  pvalue=pvalue[idx]
  genes=genes[idx]
  fwer=p.adjust(pvalue,method="bonferroni")
  idx=which(fwer<0.05)
  n=length(pvalue)
  if (is.null(xlim)) xlim=6.7
  if (is.null(ylim)) ylim=-log10(min(pvalue,na.rm = T))+0.5
  plot(-log((1:n)/n,base=10),-log(pvalue,base=10),xlab="Expected p-value (log base 10)",
       ylab="Observed p-value (log base 10)",main=main,xlim=c(0,xlim),ylim=c(0,ylim),cex.lab=2.2,cex.axis=2,cex.main=2.5)
  abline(0,1,lty=2)
  
  if (length(idx)>0)
  {
    # abline(h=-log10(0.05/length(pvalue)),lty=2,col="green")
    # text(x=0.2*(par("usr")[1]+par("usr")[2]),y=-log10(0.05/length(pvalue))+0.2,labels = "FWER=0.05",col="green",cex=1.2)
    idx=rev(idx)
    color2=colorRampPalette(c("blue", "red"))(length(idx))
    if (length(color2)==1) color2="red"
    points(-log10(idx/n),-log10(pvalue[idx]),col=color2[1:length(idx)],pch=19,cex=1.5)
    tmp=-log10(idx/n)
    len=nchar(genes[idx])
    # tmp[seq(1,length(idx),2)]=tmp[seq(1,length(idx),2)]-0.08*len[seq(1,length(idx),2)]
    # tmp[seq(2,length(idx),2)]=tmp[seq(2,length(idx),2)]+0.08*len[seq(2,length(idx),2)]
    rt=(par("usr")[2]-par("usr")[1])/5
    leftmove=rep(0.95*rt,length(len))
    leftmove[len==3]=0.53*rt
    leftmove[len==4]=0.7*rt
    leftmove[len==5]=0.8*rt
    leftmove[len==6]=1*rt
    rightmove=rep(0.9*rt,length(len))
    rightmove[len==3]=0.53*rt
    rightmove[len==4]=0.7*rt
    rightmove[len==5]=0.8*rt
    rightmove[len==6]=1*rt
    tmp[seq(1,length(idx),2)]=tmp[seq(1,length(idx),2)]-leftmove[seq(1,length(idx),2)]
    if (length(idx)>1)
    {
      tmp[seq(2,length(idx),2)]=tmp[seq(2,length(idx),2)]+rightmove[seq(2,length(idx),2)]
    }
    rty=(par("usr")[4]-par("usr")[3])/10
    tmpy=rep(0,length(idx))
    if (length(idx)>2)
    {
      k=length(idx)-3+1
      for (j in 1:length(idx))
      {
        # if (j<k)
        # {
        #   tmpy[j]=tmpy[j]+rty*0.1
        # }
        if (j==k)
        {
          tmpy[j]=tmpy[j]-rty*0.2 #move #3 down
        }
        if (j>k)
        {
          tmpy[j]=tmpy[j]+rty*0.3
        }
      }
    }
    
    text(x=tmp,y=-log10(pvalue[idx])+tmpy,labels = genes[idx],col=color2[1:length(idx)],cex=2)
  }
  if (!is.null(legend)) legend("bottomright",legend=legend,cex=2.5)
}

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
col1=c(5,3,7) #regular twas 
col2=c(3,2,4)
names=c("GTEx junction","GTEx mucosa","GTEx whole blood")
prefixs=c("dist500K_GTEx_June11","dist500K_GTEx_mucosa_June11","dist500K_GTEx_blood_June11")
mains=c("EAC vs BE","EAC vs control","EAC/BE vs control")
r2cutoffs=c(0.1,0,0)
xlims=c(5.5,6.6,6.5)
ylims=c(5.5,6.6,7)
pdf("../result/TWAS_plot.pdf",height = 9,width = 15)
par(mfrow=c(2,3))
#Regular TWAS
for (i in 1:length(names))
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefixs[i])
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  load(paste0(outfolder,"/bca_assoc.RData"))
  genesr2=rownames(res_min_code)[res_min_code$r2>r2cutoffs[i]]
  assoc_min_code=assoc_min[rownames(assoc_min) %in% intersect(proteingenes,genesr2),]
  qqplotnew(pvalue=assoc_min_code[,col1[i]],genes=rownames(assoc_min_code),
            legend = names[i],main=mains[i],xlim=xlims[i],ylim=ylims[i])
  
}
#SKAT
for (i in 1:length(names))
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefixs[i])
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p")
  genesr2=rownames(res_min_code)[res_min_code$r2>r2cutoffs[i]]
  skat_min2_code=skat_min2[rownames(skat_min2) %in% intersect(proteingenes,genesr2),]
  qqplotnew(pvalue=skat_min2_code[,col2[i]],genes=rownames(skat_min2_code),
            legend = names[i],main=mains[i],xlim=xlims[i],ylim=ylims[i])
}
dev.off()

#CCA on chr19

library(GenomicRanges)
library(PMA)
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

genes=c("MEF2B","DDX49","KXD1","CERS1","COMP","ARMC6","TMEM161A","CRTC1","HOMER3","LSM4","KLHL26","PGPEP1","GDF15")
gene="COMP"
validatedgenes=c("THEM161A","KLHL26","PGPEP1","HOMER3")
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8blooddata_ambiguous_TPM_for_prediction.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")

#load(paste0(outfolder,"/preidiction_michigan_model.RData"))

idx=which(rownames(phenotype) %in% proteingenes,)
phenotype=phenotype[idx,]
phenotypepos=phenotypepos[idx,]
all(genes %in% rownames(phenotype))
genes=intersect(genes,rownames(phenotype))

idx=which(rownames(phenotype)==gene)
chr=phenotypepos$chr[idx]
startloc=phenotypepos$s1[idx]-5e5
startloc=max(1,startloc)
endloc=phenotypepos$s2[idx]+5e5
gr_region=GRanges(seqnames = chr,ranges = IRanges(start=startloc,end=endloc))
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
tmp=distance(gr_pos,gr_region)
sum(tmp==0,na.rm=T) #29 genes
all(genes %in% phenotypepos$geneid[tmp==0])

idx=which(tmp==0)
CCAZ=phenotype[idx,,drop=F]

idx=match(rownames(CCAZ),phenotypepos$geneid)
idx1=order(phenotypepos$s1[idx])
CCAZ1=CCAZ[idx1,]
library(gplots)
# cormat <- cor(t(CCAZ1),use="pairwise.complete.obs")
# col<- colorRampPalette(c("blue", "white", "red"))(10)
# heatmap(cormat, Rowv=NA, Colv=NA,col=col, symm=TRUE)
heatmap(cormat,col=col, symm=TRUE)

heatmap.2(cor(t(CCAZ1)), col=colorRampPalette(c("blue", "white", "red"))(10), keysize = 1.2,key.title = "",cexRow = 0.8,cexCol = 0.8,
          density.info="none", trace="none", dendrogram="none", 
          symm=T,symkey=T,symbreaks=T, scale="none")

library(ComplexHeatmap)
library(GetoptLong)
library(circlize)

ht_global_opt(
  # heatmap_legend_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  # heatmap_legend_labels_gp = gpar(fontsize = 10,fontface = "bold"), 
  heatmap_column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_row_names_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_column_title_gp = gpar(fontsize = 10,fontface = "bold"),
  heatmap_row_title_gp = gpar(fontsize = 10,fontface = "bold")
)


mat2=cor(t(CCAZ1))
pdf("../result/CHR19_corr_blood1.pdf",width = 5,height = 8)
ht_list = Heatmap(mat2, name = "Correlation",
                  #col = greenred(5),# colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  #col=colorRampPalette(c("blue", "white", "red"))(30),
                  col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                  #col=colorpanel(8,"green","orange","red"),
                  #column_dend_height = unit(4, "cm"),
                  cluster_rows = F,
                  cluster_columns=F,
                  column_dend_reorder=T,
                  show_row_dend=F,
                  show_column_dend=F,
                  row_dend_width=unit(2.5,"cm"),
                  show_heatmap_legend=F,
                  heatmap_height=unit(12,"cm"),
                  #top_annotation = ha,
                  show_column_names = T, show_row_names=T)
ht_list = draw(ht_list, heatmap_legend_side = "right")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
exprlgd=Legend(col_fun = col_fun, title = "Correlation", at = c(-1,0,1),direction = "horizontal",
               labels_gp=gpar(fontsize = 14,fontface = "bold"),
               legend_gp=gpar(fontsize = 14,fontface = "bold"),
               title_gp=gpar(fontsize = 14,fontface = "bold"))
draw(exprlgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
dev.off()

idx=match(rownames(CCAZ1),rownames(phenotypepos))
View(phenotypepos[idx,])
