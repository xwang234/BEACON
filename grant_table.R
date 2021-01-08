
#First is a table below, including TCGA/GTEx junction/GTEx mucosal/ cis-eQTL mapping result, using 500KB and FDR control.
library(GenomicRanges)
updatesnpname=function(oldnames=eqtlres$SNP)
{
  newnames=sapply(oldnames,function(x){
    unlist(strsplit(x,":"))[1]
  })
  return(newnames)
}

snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_26highrisk_SNP_POS.txt"
snppos=read.table(snpposfile,header = T,stringsAsFactors = F)
snppos$snp=updatesnpname(snppos$snp)
gr_snppos=GRanges(seqnames = snppos$chr,ranges = IRanges(start=snppos$position,width = 1))
cytobandtable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt",header=T,stringsAsFactors = F)
cytobandtable$chrom=gsub("chr","",cytobandtable$chrom)
gr_cytobandtable=GRanges(seqnames = cytobandtable$chrom,ranges=IRanges(start=cytobandtable$start,end=cytobandtable$end))
snppos$cytoband=NA
for (i in 1:nrow(snppos))
{
  tmp=distanceToNearest(gr_snppos[i],gr_cytobandtable)
  if (tmp@elementMetadata$distance==0)
  {
    snppos$cytoband[i]=paste0(cytobandtable$chrom[tmp@to],unlist(strsplit(cytobandtable$name[tmp@to],".",fixed=T))[1])
  }
}
snppos$gene=NA
genetable=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",header = F,stringsAsFactors = F,sep="\t")
genetable=genetable[genetable$V2 %in% paste0("chr",c(1:22,"X","Y")),]
genetable$V2=gsub("chr","",genetable$V2)
genetable=genetable[genetable$V11!="",]
gr_tanscript=GRanges(seqnames = genetable$V2,ranges = IRanges(start=genetable$V4,end=genetable$V5))
for (i in 1:nrow(snppos))
{
  tmp=distance(gr_tanscript,gr_snppos[i])
  idx=which(tmp==0)
  if (length(idx)==0) #intergenic
  {
    idx=which(genetable$V2==snppos$chr[i] & genetable$V4<snppos$position[i])
    tmp=idx[length(idx)]
    snppos$gene[i]=paste0("Intergenic;",genetable$V14[tmp],"/",genetable$V14[tmp+1])
  }else
  {
    idx=idx[length(idx)]
    tmpstarts=as.integer(unlist(strsplit(genetable$V9[idx],",")))
    tmpends=as.integer(unlist(strsplit(genetable$V10[idx],",")))
    gr_tmp=GRanges(seqnames = snppos$chr[i],ranges = IRanges(start=tmpstarts,end=tmpends))
    tmp=distance(gr_tmp,gr_snppos[i])
    if (sum(tmp==0,na.rm=T)>0)
    {
      snppos$gene[i]=genetable$V14[idx]
    }else
    {
      snppos$gene[i]=paste0("Intronic;",genetable$V14[idx])
    }
  }
}
#mqtlres=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_26highrisk_cn_mutation_cis",header=T,stringsAsFactors = F)
mqtlres=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_26highrisk_cn_mutation_nopeers_cis",header=T,stringsAsFactors = F) #no 15 PCs
mqtlres$SNP=updatesnpname(oldnames = mqtlres$SNP)
load("/fh/fast/dai_j/CancerGenomics/prostate_methylation/annotation.RData")
anno$UCSC_RefGene_Name=as.character(anno$UCSC_RefGene_Name)
idx=match(mqtlres$gene,anno$IlmnID)
mqtlres$genename=processsemicolon(anno$UCSC_RefGene_Name[idx])
tcgaeqtl=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_EAC_26highrisk_cn_mutation_nopeers_cis",header=T,stringsAsFactors = F)
tcgaeqtl$SNP=updatesnpname(tcgaeqtl$SNP)
gtexjunctioneqtl=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtexcoding_26highrisk_cn_mutation_cis",header = T,stringsAsFactors = F)
gtexjunctioneqtl$SNP=updatesnpname(gtexjunctioneqtl$SNP)
gtexmucosaeqtl=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtexcoding_mucosa_26highrisk_cn_mutation_cis",header = T,stringsAsFactors = F)
gtexmucosaeqtl$SNP=updatesnpname(gtexmucosaeqtl$SNP)
table1=function(cutoff=0.2)
{
  res=data.frame(snp=snppos$snp,region=snppos$gene,chr=snppos$cytoband,TCGA_EAC="-",TCGA_eGenes="-",TCGA_mGenes="-",GTEx_junction="-",junction_eGenes="-",GTEx_mucosa="-",mucosa_eGenes="-",stringsAsFactors = F)
  idx1=which(tcgaeqtl$qvalue<cutoff)
  idx2=which(gtexjunctioneqtl$qvalue<cutoff)
  idx3=which(gtexmucosaeqtl$qvalue<cutoff)
  idx4=which(mqtlres$qvalue<cutoff)
  for (i in 1:nrow(res))
  {
    if (snppos$snp[i] %in% tcgaeqtl$SNP[idx1] & !snppos$snp[i] %in% mqtlres$SNP[idx4])
    {
      res$TCGA_EAC[i]="cis-eQTL"
    }
    if (snppos$snp[i] %in% tcgaeqtl$SNP[idx1] & snppos$snp[i] %in% mqtlres$SNP[idx4])
    {
      res$TCGA_EAC[i]="cis-eQTL,cis-meQTL"
    }
    if (!snppos$snp[i] %in% tcgaeqtl$SNP[idx1] & snppos$snp[i] %in% mqtlres$SNP[idx4])
    {
      res$TCGA_EAC[i]="cis-meQTL"
    }
    
    if (snppos$snp[i] %in% gtexjunctioneqtl$SNP[idx2])
    {
      res$GTEx_junction[i]="cis-eQTL"
    }
    if (snppos$snp[i] %in% gtexmucosaeqtl$SNP[idx3])
    {
      res$GTEx_mucosa[i]="cis-eQTL"
    }
    #egenes
    if (snppos$snp[i] %in% tcgaeqtl$SNP[idx1])
    {
      res$TCGA_eGenes[i]=paste0(tcgaeqtl$gene[idx1][tcgaeqtl$SNP[idx1]==snppos$snp[i]],collapse = ";")
    }
    if (snppos$snp[i] %in% gtexjunctioneqtl$SNP[idx2])
    {
      res$junction_eGenes[i]=paste0(gtexjunctioneqtl$gene[idx2][gtexjunctioneqtl$SNP[idx2]==snppos$snp[i]],collapse = ";")
    }
    if (snppos$snp[i] %in% gtexmucosaeqtl$SNP[idx3])
    {
      res$mucosa_eGenes[i]=paste0(gtexmucosaeqtl$gene[idx3][gtexmucosaeqtl$SNP[idx3]==snppos$snp[i]],collapse = ";")
    }
    #mGenes
    if (snppos$snp[i] %in% mqtlres$SNP[idx4])
    {
      idx11=which(mqtlres$SNP==snppos$snp[i] & mqtlres$qvalue<cutoff)
      mprobes=mqtlres$gene[idx11]
      mgenes=rep("",length(idx11))
      for (j in 1:length(idx11))
      {
        tmp=paste0(mqtlres$genename[idx11[j]][mqtlres$genename[idx11[j]]==snppos$snp[i]],collapse = ";")
        if (tmp == "") #Intergenic/Intronic
        {
          idx22=which(anno$IlmnID==mprobes[j])
          tmp1=processsemicolon(anno$UCSC_RefGene_Name[(idx22-1000):idx22])
          tmp1=tmp1[tmp1!=""]
          tmp2=processsemicolon(anno$UCSC_RefGene_Name[idx22:(idx22+1000)])
          tmp2=tmp2[tmp2!=""]
          if (tmp1[length(tmp1)]!=tmp2[1])
          {
            mgenes[j]=paste0("Intergeneic;",tmp1[length(tmp1)],"/",tmp2[1])
          }else
          {
            mgenes[j]=paste0("Intronic;",tmp1[length(tmp1)])
          }
        }else
        {
          mgenes[j]=tmp
        }
      }
      tmp=NULL
      for (j in 1:length(mprobes))
      {
        tmp=c(tmp,mprobes[j],"(")
        tmp=c(tmp,mgenes[j],")")
        if (length(mprobes)>1 & j<length(mprobes))
        {
          tmp=c(tmp,"|")
        }
      }
      res$TCGA_mGenes[i]=paste0(tmp,collapse = "")
    }
  }
  return(res)
}
res=table1()
write.csv(res,file="../result/Ciseqtllist_genes_qvalue005_mqtlnopc.csv",row.names = F)
# gtexjunctioneqtl=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtexcoding_26highrisk_cn_mutation_nopeers_cis",header = T,stringsAsFactors = F)
# gtexjunctioneqtl$SNP=updatesnpname(gtexjunctioneqtl$SNP)
# gtexmucosaeqtl=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtexcoding_mucosa_26highrisk_cn_mutation_nopeers_cis",header = T,stringsAsFactors = F)
# gtexmucosaeqtl$SNP=updatesnpname(gtexmucosaeqtl$SNP)
# res1=table1()

#draw qqplot
#skat-based results
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_PC4"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTExmucosa"
load(paste0(outfolder,"/skat_res.RData"))
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
qqplotnew=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL,genes=rownames(skat_min_code))
{
  library(RColorBrewer)
  
  par("mar"=c(5.1,4.1,2.1,2.1))
  idx=order(pvalue)
  pvalue=pvalue[idx]
  genes=genes[idx]
  fwer=p.adjust(pvalue,method="bonferroni")
  idx=which(fwer<0.05)
  n=length(pvalue)
  if (is.null(xlim)) xlim=5
  if (is.null(ylim)) ylim=-log10(min(pvalue,na.rm = T))+0.5
  plot(-log((1:n)/n,base=10),-log(pvalue,base=10),xlab="Expected p-value (log base 10)",
       ylab="Observed p-value (log base 10)",main=main,xlim=c(0,xlim),ylim=c(0,ylim),cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
  
  if (length(idx)>0)
  {
    abline(h=-log10(0.05/length(pvalue)),lty=2,col="green")
    text(x=0.3*(par("usr")[1]+par("usr")[2]),y=-log10(0.05/length(pvalue))+0.2,labels = "FWER=0.05",col="green")
    idx=rev(idx)
    color2=colorRampPalette(c("blue", "red"))(length(idx))
    if (length(color2)==1) color2="red"
    points(-log10(idx/n),-log10(pvalue[idx]),col=color2[1:length(idx)],pch=19)
    tmp=-log10(idx/n)
    len=nchar(genes[idx])
    # tmp[seq(1,length(idx),2)]=tmp[seq(1,length(idx),2)]-0.08*len[seq(1,length(idx),2)]
    # tmp[seq(2,length(idx),2)]=tmp[seq(2,length(idx),2)]+0.08*len[seq(2,length(idx),2)]
    rt=(par("usr")[2]-par("usr")[1])/5
    leftmove=rep(0.65*rt,length(len))
    leftmove[len==3]=0.32*rt
    leftmove[len==4]=0.37*rt
    leftmove[len==5]=0.42*rt
    leftmove[len==6]=0.48*rt
    rightmove=rep(0.65*rt,length(len))
    rightmove[len==3]=0.32*rt
    rightmove[len==4]=0.37*rt
    rightmove[len==5]=0.42*rt
    rightmove[len==6]=0.44*rt
    tmp[seq(1,length(idx),2)]=tmp[seq(1,length(idx),2)]-leftmove[seq(1,length(idx),2)]
    if (length(idx)>1)
    {
      tmp[seq(2,length(idx),2)]=tmp[seq(2,length(idx),2)]+rightmove[seq(2,length(idx),2)]
    }
    text(x=tmp,y=-log10(pvalue[idx]),labels = genes[idx],col=color2[1:length(idx)],cex=1)
  }
}
#TCGA
qqplotnew(pvalue=skat_min$BEEA_p,genes=rownames(skat_min))
#GTEx
qqplotnew(pvalue=skat_min_code$BEEA_p)

#for regular (not using skat) results


#TCGA
prefix="dist500K_EAC2"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_assoc.RData"))
qqplotnew(pvalue=assoc_min$BEEA_p,genes=rownames(assoc_min))
#GTEx junction
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx"
load(paste0(outfolder,"/bca_assoc.RData"))
qqplotnew(pvalue=assoc_min$BEEA_p,genes=rownames(assoc_min))
#GTEx mucosa
prefix="dist500K_GTExmucosa"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_assoc.RData"))
qqplotnew(pvalue=assoc_min$BEEA_p,genes=rownames(assoc_min))

#check the modified skat (coeff)
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_PC4"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTExmucosa"
rm(skat_min)
load(paste0(outfolder,"/skat_coef_abs_res.RData"))
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,] #GTEx
#TCGA
qqplotnew(pvalue=skat_min$BEEA_p,genes=rownames(skat_min))
#GTEx
qqplotnew(pvalue=skat_min_code$BEEA_p)


#for Fig7, rs17749155, meQTL and CpG~BMI
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr)  
anno <- anno[anno$chr!="X"&anno$chr!="Y",]
anno$chr <- as.numeric(anno$chr)
anno$pos <- as.numeric(anno$pos)
anno=as.data.frame(anno)
idx=order(anno$chr,anno$pos)
anno=anno[idx,]
#GSE data
load("../data/GSE89181.RData")
library(limma)
idx=complete.cases(GSEclinical[,c("bmi","type","gender","age")])
quantile(GSEclinical$bmi[idx])
# 0%      25%      50%      75%     100% 
# 20.63956 24.95216 29.98314 30.98900 42.76745 
designmat=model.matrix(~I(bmi>30)+type+gender+age,data=GSEclinical[idx,])
#fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
fit <- lmFit(as.matrix(GSEmethall[,idx]), designmat)
result=eBayes(fit)
colnames(result$p.value)
#qqplot(result$p.value[,2])
tmp=p.adjust(result$p.value[,2],method="bonferroni")
sum(tmp<0.2) #1
min(tmp) #[1] 0.008789511
rownames(result$p.value)[which.min(tmp)] #"cg13601977"
anno$UCSC_RefGene_Name[which(anno$Name=="cg13601977")] #""
which(anno$Name=="cg13601977") #242668
unique(anno$UCSC_RefGene_Name[(242668-20):242668]) #"PHF2"
unique(anno$UCSC_RefGene_Name[242668:(242668+30)]) #"BARX1"
# start1=which(anno$UCSC_RefGene_Name=="PHF2")
# end1=which(anno$UCSC_RefGene_Name=="BARX1")
# idx=(start1[length(start1)]+1):(end1[1]-1)
idx=which(grepl("MSRA",anno$UCSC_RefGene_Name))
min(idx) #220948
max(idx) #221034
unique(anno$UCSC_RefGene_Name[(220948-50):220948]) #MIR124-1" "MSRA;MSRA"
unique(anno$UCSC_RefGene_Name[221034:(221034+30)]) #T-SP1
start1=which(grepl("MIR124-1",anno$UCSC_RefGene_Name))
end1=which(anno$UCSC_RefGene_Name=="T-SP1")
idxMSRA=(start1[length(start1)]+1):(end1[1]-1)
MSRAprobes=anno$Name[idxMSRA]
idx=match(rownames(result$p.value),anno$Name)
result$p.value=as.data.frame(result$p.value)
result$p.value$gene=anno$UCSC_RefGene_Name[idx]
#idx=which(grepl("MSRA",result$p.value$gene))
#MSRAresult=result$p.value[idx,]
MSRAresult=result$p.value[rownames(result$p.value) %in% MSRAprobes,]
MSRAresult$gene=processsemicolon(MSRAresult$gene)
#qqplot(result$p.value$`I(bmi > 30)TRUE`[idx])
quantile(MSRAresult$`I(bmi > 30)TRUE`,na.rm=T)
#     0%         25%         50%         75%        100% 
# 0.007928466 0.135953940 0.455220264 0.697733963 0.969974741 
tmp=p.adjust(MSRAresult$`I(bmi > 30)TRUE`,method="fdr")
quantile(tmp)
which.min(tmp)
idx=match(rownames(MSRAresult),anno$Name)
MSRAresult$enhancer=anno$Enhancer[idx]
MSRAresult$promoter=anno$Regulatory_Feature_Group[idx]
MSRAresult$group=anno$UCSC_RefGene_Group[idx]
MSRAresult$group=processsemicolon(MSRAresult$group)
MSRAresult$position=anno$pos[idx]
tmp=mqtlres$gene[mqtlres$p.value<0.05 & mqtlres$SNP=="rs17749155"]
tmp1=rownames(MSRAresult)[MSRAresult$`I(bmi > 30)TRUE`<0.05]
tmp1=rownames(result$p.value)[result$p.value$`I(bmi > 30)TRUE`<0.05]
tmp2=intersect(tmp,tmp1) #0 overlap
mqtlres_MSRA=mqtlres[mqtlres$gene %in% rownames(MSRAresult),]
idx=match(mqtlres_MSRA$gene,rownames(MSRAresult))
mqtlres_MSRA$bmipvalue=MSRAresult$`I(bmi > 30)TRUE`[idx]
mqtlres_MSRA$genename=MSRAresult$gene[idx]
mqtlres_MSRA$enhancer=MSRAresult$enhancer[idx]
mqtlres_MSRA$promoter=MSRAresult$promoter[idx]
mqtlres_MSRA$group=MSRAresult$group[idx]
which(mqtlres_MSRA$gene=="cg17575721") #53
mqtlres_MSRA[53,]
# SNP       gene    p.value    t.stat        beta     snp_r2 betacopy pcopy copy_r2 betamutation pmutation
# 3343 rs17749155 cg17575721 0.02618801 -2.264415 -0.03701827 0.05885135       NA    NA      NA           NA        NA
# mutation_r2       FDR    qvalue  bmipvalue genename enhancer promoter group
# 3343          NA 0.8128994 0.7275781 0.05369424     MSRA     TRUE           Body
medat=fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_EAC_ME_rmprobes.txt",header = T)
medat=as.data.frame(medat)
rownames(medat)=medat$id
medat=medat[,-1]
snpdat=fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_EAC_26highrisk_SNP_ME.txt")
snpdat=as.data.frame(snpdat)
rownames(snpdat)=snpdat$id
snpdat=snpdat[,-1]
all(colnames(medat)==colnames(snpdat)) #T
rownames(snpdat)[which(grepl("rs17749155",rownames(snpdat)))]
plots_snp_meth_bmi=function(snp="rs17749155",probe1="cg17575721",probe2="cg01511457",
                           #snpdat=snpdat,
                           colors=c("red","green","blue"),
                           pchs=c(1,2,3),
                           #medat=medat,
                           ylim=c(0.6,1),
                           ylim1=c(0.15,0.85),
                           tx=2,ty=0.99,
                           Aale="G",Bale="A",plotflag=T)
{
  gt1=paste0(Aale,'/',Aale)
  gt2=paste0(Aale,'/',Bale)
  gt3=paste0(Bale,'/',Bale)
  # n=length(snpdat)
  # nsamples=rep(0,n)
  # for (i in 1:n)
  # {
  #   nsamples[i]= ncol(snpdat[[i]])
  # }
  # 
  
  idx=which(grepl(snp,rownames(snpdat)))
  dosage=unlist(snpdat[idx,])
  x=convert_dosage2_gt(gt1,gt2,gt3,dosage)
  
    
  idx=which(rownames(medat)==probe1)
  y=matrix(unlist(medat[idx,]),ncol=1)
  
  x=factor(x,levels = c(gt1,gt2,gt3))
  x=as.character(x)
  x[x=="G/A"]="G/A or A/A"
  x[x=="A/A"]="G/A or A/A"
  x=factor(x,levels=c("G/G","G/A or A/A"))
  ymin=min(y,na.rm=T)
  ymax=max(y,na.rm=T)
  if (is.null(ylim)) ylim=c(ymin,ymax)
  #legend("top",inset=0,legend=legends,col=colors,pch=16,horiz=T,cex=1.3,bty="n")
  par(mfrow=c(1,2))
  par(mar=c(3,4,2,1))
  plot(0,xlim=c(0,4),ylim=ylim,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',axes=F,type='n')
  boxplot(y~x,outline=F,xlim=c(0,4),ylab='DNA methylation beta value',cex.lab=1.2,cex.axis=1.2,at=c(1,3),add=T)
  stripchart(y~x,vertical = TRUE, method = "jitter",pch = pchs[1],col=colors[1],at=c(1,3),add = TRUE)
  fm=glm(y~dosage)
  print(summary(fm)$coefficients)
  text(tx,ty,labels=paste0('p-value=',round(summary(fm)$coefficients[2,4],3)),cex=1.2)
  idx=which(rownames(GSEmeth)==probe2)
  idx1=GSEclinical$bmi>30
  y1=unlist(GSEmeth[idx,!is.na(idx1)])
  ymin1=min(y1,na.rm=T)
  ymax1=max(y1,na.rm=T)
  if (is.null(ylim1)) ylim1=c(ymin1,ymax1)
  boxplot(y1~!idx1[!is.na(idx1)],ylim=ylim1,names=c("BMI<=30","BMI>30"),cex.axis=1.2,cex.lab=1.2,outline=F,ylab='DNA methylation beta value',yaxt='n')
  stripchart(y1~!idx1[!is.na(idx1)],vertical = TRUE, method = "jitter",pch = pchs[c(2,2)],col=colors[c(3,3)],add = TRUE)
  #axis(1,at=c(1,2),cex.lab=1.2,labels = c("BMI<=30","BMI>30"))
  axis(2,at=seq(0.1,0.9,0.1),cex.axis=1.2)
  fm=glm(y1~idx1[!is.na(idx1)])
  text(1.5,0.83,labels=paste0('p-value=',round(summary(fm)$coefficients[2,4],3)),cex=1.2)
}
pdf("../result/fig7.pdf",width=6,height=4)
plots_snp_meth_bmi()
dev.off()
anno[which(anno$Name=="cg17575721"),]

plots_snp_meth_bmi(probe1="cg00958676",probe2="cg01511457")
plots_snp_meth_bmi1=function(snp="rs17749155",probe1="cg17575721",probe2="cg17575721",
                            #snpdat=snpdat,
                            colors=c("red","green","blue"),
                            pchs=c(1,2,3),
                            #medat=medat,
                            ylim=c(0.6,1),
                            ylim1=c(0.1,1),
                            tx1=2,tx2=5.5,ty=0.99,
                            Aale="G",Bale="A",plotflag=T)
{
  gt1=paste0(Aale,'/',Aale)
  gt2=paste0(Aale,'/',Bale)
  gt3=paste0(Bale,'/',Bale)
  # n=length(snpdat)
  # nsamples=rep(0,n)
  # for (i in 1:n)
  # {
  #   nsamples[i]= ncol(snpdat[[i]])
  # }
  # 
  
  idx=which(grepl(snp,rownames(snpdat)))
  dosage=unlist(snpdat[idx,])
  x=convert_dosage2_gt(gt1,gt2,gt3,dosage)
  
  
  idx=which(rownames(medat)==probe1)
  y=matrix(unlist(medat[idx,]),ncol=1)
  
  x=factor(x,levels = c(gt1,gt2,gt3))
  x=as.character(x)
  x[x=="G/A"]="G/A or A/A"
  x[x=="A/A"]="G/A or A/A"
  x=factor(x,levels=c("G/G","G/A or A/A"))
  
  ymin=min(y,na.rm=T)
  ymax=max(y,na.rm=T)
  if (is.null(ylim)) ylim=c(ymin,ymax)
  #legend("top",inset=0,legend=legends,col=colors,pch=16,horiz=T,cex=1.3,bty="n")
  par(mfrow=c(1,1))
  par(mar=c(3,4,2,1))
  plot(0,xlim=c(0,7),ylim=ylim,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',axes=F,type='n')
  boxplot(y~x,outline=F,xlim=c(0,4),ylab='DNA methylation beta value',cex.lab=1.2,cex.axis=1.2,at=c(1,3),add=T)
  stripchart(y~x,ylim=ylim,vertical = TRUE, method = "jitter",pch = pchs[i],col=colors[i],at=c(1,3),add = TRUE)
  fm=glm(y~dosage)
  print(summary(fm)$coefficients)
  text(tx,ty,labels=paste0('p-value=',round(summary(fm)$coefficients[2,4],3)),cex=1.2)
  idx=which(rownames(GSEmeth)==probe2)
  idx1=GSEclinical$bmi>30
  y1=unlist(GSEmeth[idx,!is.na(idx1)])
  ymin1=min(y1,na.rm=T)
  ymax1=max(y1,na.rm=T)
  if (is.null(ylim1)) ylim1=c(ymin1,ymax1)
  boxplot(y1~!idx1[!is.na(idx1)],ylim=ylim1,names=c("BMI<=30","BMI>30"),cex.axis=1.2,cex.lab=1.2,outline=F,ylab='DNA methylation beta value',at=c(5,6),add=T)
  stripchart(y1~!idx1[!is.na(idx1)],vertical = TRUE, method = "jitter",pch = pchs[c(2,2)],col=colors[c(3,3)],add = TRUE,at=c(5,6))
  #axis(1,at=c(1,2),cex.lab=1.2,labels = c("BMI<=30","BMI>30"))
 
  fm=glm(y1~idx1[!is.na(idx1)])
  text(tx2,ty,labels=paste0('p-value=',round(summary(fm)$coefficients[2,4],3)),cex=1.2)
}  
plots_snp_meth_bmi1()

which(grepl("MSRA",rownames(geneexpdata))) #11255
idx=match(colnames(medat),colnames(geneexpdata))
which(rownames(medat)=="cg17575721") #208805
cor(unlist(medat[208805,]),unlist(geneexpdata[11255,idx]))
tmp1=unlist(geneexpdata[11255,idx])
tmp2=unlist(medat[208805,])
idx=tmp2>0.6 & tmp1 <quantile(tmp1,0.95)
cor(tmp1[idx],tmp2[idx])
cor.test(tmp1[idx],tmp2[idx])

plot(tmp2[idx],tmp1[idx])
fm=glm(tmp1[idx]~tmp2[idx])
abline(fm)

#for fig5
genes=c("MEF2B","DDX49","CERS1","COMP","ARMC6","TMEM161A","CRTC1","HOMER3")
prefix="dist500K_GTEx_PC4"
prefix="dist500K_TCGA_PC4"
prefix="dist500K_GTExmucosa"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data//bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  sampletable$site[sampletable$site=="NA"]=NA
  geneexpsamplenames=strsplit(colnames(predict_1se)[3:ncol(predict_1se)],"_")
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}

readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  if (opt==1)
  {
    eigsamples=read.table(eigsampfile,stringsAsFactors = F)
    eigsamples=eigsamples$V2
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=thesamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}

eigenstratmatrix=readeigenstrat()

effectsize_arow=function(i,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  x=as.numeric(predict_bcageneexp[i,c(idx1,idx2)])
  x=scale(x)
  effectsize=NA
  names(effectsize)=rownames(predict_bcageneexp)[i]
  if (sum(is.na(x))<0.6*length(x))
  {
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=covariates,family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      effectsize=summary(fm)$coefficients[2,1]
    }
  }
  return(effectsize)
}
predict_bcageneexp=predict_min
predict_bcageneexp=predict_bcageneexp[,3:ncol(predict_bcageneexp)]
geneexpsamplenames=colnames(predict_bcageneexp)
geneexpsamplenames=strsplit(geneexpsamplenames,"_")
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_bcageneexp)=geneexpsamplenames
geneexpsamplenames=intersect(sampletable$localid,geneexpsamplenames)
idx=match(geneexpsamplenames,colnames(predict_bcageneexp))
predict_bcageneexp=predict_bcageneexp[,idx]
idx=match(geneexpsamplenames,sampletable$localid)
sampletable1=sampletable[idx,]
for (i in 1:ncol(sampletable1)) sampletable1[,i][sampletable1[,i]==-9]=NA
idx=match(geneexpsamplenames,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:4,]))
idx1=which(sampletable1$phenoBE_bc==2 | sampletable1$phenoEA_bc==2) #BE or EAcase
length(idx1) #5802
idx2=which(sampletable1$phenoBE_bc==1 |sampletable1$phenoEA_bc==1) #Control
length(idx2) #2176
length(intersect(idx1,idx2))
covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                      bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                      pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],pc4=sampletable1$pc4[c(idx1,idx2)],
                      stringsAsFactors = F)
if ("LASS1" %in% rownames(predict_bcageneexp))
{
  rownames(predict_bcageneexp)[which(rownames(predict_bcageneexp)=="LASS1")]="CERS1"
}
if ("LASS1" %in% rownames(res_min))
{
  rownames(res_min)[which(rownames(res_min)=="LASS1")]="CERS1"
}
sum(genes %in% rownames(predict_bcageneexp))
sum(genes %in% rownames(res_min))
genes[!genes %in% rownames(res_min)]
genes[genes %in% rownames(predict_bcageneexp)]

res1=data.frame(gene=genes,effectsize=NA,stringsAsFactors = F)
for (i in 1:nrow(res1))
{
  if (genes[i] %in% rownames(predict_bcageneexp))
  {
    idx=which(rownames(predict_bcageneexp)==genes[i])
    res1$effectsize[i]=effectsize_arow(idx,idx1,idx2)
  }
}
res_junction=res1
res_tcga=res1
res_mucosa=res1
res=cbind(res_tcga,junction=res_junction$effectsize,mucosa=res_mucosa$effectsize)
colnames(res)[2]="tcga"
write.csv(res,file="../result/effectsize.csv",row.names = F)

#compute the beta difference in DMR
load("../data/GSE89181.RData")
load("../result/methylation_riskfactors_DMRres.RData")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr)  
anno <- anno[anno$chr!="X"&anno$chr!="Y",]
anno$chr <- as.numeric(anno$chr)
anno$pos <- as.numeric(anno$pos)
anno$UCSC_RefGene_Name=processsemicolon(anno$UCSC_RefGene_Name)
idx=match(rownames(GSEmethall),anno$Name)
anno2=anno[idx,]
idx=order(anno2$chr,anno2$pos)
idx1=GSEclinical$bmi>30
betadiff=function(dat=DMRs_bmi$table,i=1)
{
  idx2=dat$indexStart[i]:dat$indexEnd[i]
  probes=anno2$Name[idx[idx2]]
  idx3=match(probes,rownames(GSEmethall))
  res=data.frame(probe=probes,betadiff=NA,std=NA,stringsAsFactors = F)
  for (i in 1:length(idx3))
  {
    res$beta1[i]=round(mean(as.numeric(GSEmethall[idx3[i],which(idx1==T)])),3)
    res$beta2[i]=round(mean(as.numeric(GSEmethall[idx3[i],which(idx1==F)])),3)
    res$betadiff[i]=round(mean(as.numeric(GSEmethall[idx3[i],which(idx1==T)]))-mean(as.numeric(GSEmethall[idx3[i],which(idx1==F)])),3)
    res$std[i]=round(sd(as.numeric(GSEmethall[idx3[i],which(!is.na(idx1))])),3)
    res$mean[i]=round(mean(as.numeric(GSEmethall[idx3[i],which(!is.na(idx1))])),3)
    
  }
  idx3=match(res$probe,anno$Name)
  res$gene=anno$UCSC_RefGene_Name[idx3]
  return(res)
}
betadiff()
#BMI:
#   probe       betadiff  std beta1  beta2 mean gene
# 1 cg04131969    0.372 0.353 0.657 0.285 0.462 MYADML
#Smoke:
idx1=rep(NA,nrow(GSEclinical))
idx1[GSEclinical$smoke=="Yes"]=T
idx1[GSEclinical$smoke=="No"]=F
betadiff(dat=DMRs_smoke$table,i=1)
#     probe     betadiff  std  beta1 beta2  mean  gene
# 1  cg20265733    0.110 0.240 0.604 0.494 0.576 GATA5
# 2  cg02484469    0.093 0.227 0.608 0.515 0.584 GATA5
# 3  cg14980983    0.099 0.229 0.598 0.499 0.572 GATA5
# 4  cg24320612    0.131 0.186 0.514 0.383 0.480 GATA5
# 5  cg16714055    0.117 0.223 0.589 0.472 0.559 GATA5
# 6  cg11982072    0.133 0.222 0.589 0.456 0.555 GATA5
# 7  cg24500900    0.124 0.230 0.475 0.351 0.443 GATA5
# 8  cg08568720    0.134 0.244 0.574 0.441 0.540 GATA5
# 9  cg18278265    0.119 0.199 0.458 0.339 0.427 GATA5
# 10 cg02146001    0.115 0.188 0.547 0.432 0.517 GATA5
# 11 cg23770904    0.140 0.219 0.601 0.461 0.565 GATA5
# 12 cg09339194    0.158 0.189 0.571 0.413 0.530 GATA5
# 13 cg25988317    0.127 0.175 0.395 0.268 0.362 GATA5
# 14 cg00677866    0.144 0.194 0.489 0.345 0.452 GATA5
# 15 cg18872321    0.150 0.190 0.377 0.227 0.338 GATA5
# 16 cg25667841    0.141 0.176 0.333 0.192 0.296 GATA5
# 17 cg14388488    0.180 0.229 0.519 0.339 0.472 GATA5
# 18 cg07914822    0.137 0.173 0.350 0.213 0.314 GATA5
# 19 cg03777459    0.152 0.231 0.611 0.459 0.571 GATA5
# 20 cg04742780    0.171 0.214 0.598 0.427 0.554 GATA5
# 21 cg20525917    0.184 0.197 0.716 0.533 0.669 GATA5
betadiff(dat=DMRs_smoke$table,i=2)
#      probe     betadiff  std beta1 beta2  mean    gene
# 1  cg04337096    0.186 0.236 0.405 0.219 0.357 PPP2R2B
# 2  cg20560075    0.229 0.243 0.363 0.134 0.304 PPP2R2B
# 3  cg03225002    0.144 0.202 0.273 0.129 0.236 PPP2R2B
# 4  cg21938261    0.160 0.183 0.318 0.158 0.276 PPP2R2B
# 5  cg03693911    0.107 0.193 0.317 0.210 0.289 PPP2R2B
# 6  cg14693112    0.180 0.235 0.346 0.167 0.300 PPP2R2B
# 7  cg08774231    0.134 0.212 0.277 0.143 0.242 PPP2R2B
# 8  cg04583285    0.149 0.243 0.291 0.141 0.252 PPP2R2B
# 9  cg09528265    0.136 0.221 0.285 0.149 0.250 PPP2R2B
# 10 cg11826826    0.153 0.202 0.266 0.113 0.226 PPP2R2B
# 11 cg03851984    0.175 0.222 0.284 0.109 0.238 PPP2R2B
# 12 cg25149751    0.146 0.200 0.339 0.193 0.301 PPP2R2B
# 13 cg03659005    0.174 0.196 0.418 0.245 0.373 PPP2R2B
# 14 cg25363885    0.188 0.210 0.399 0.211 0.351 PPP2R2B
# 15 cg13971892    0.152 0.193 0.439 0.287 0.400 PPP2R2B
# 16 cg15927927    0.166 0.194 0.421 0.255 0.378 PPP2R2B
# 17 cg25021259    0.236 0.231 0.569 0.333 0.508 PPP2R2B
# 18 cg13983063    0.170 0.208 0.438 0.269 0.394 PPP2R2B
betadiff(dat=DMRs_smoke$table,i=3)
#      probe     betadiff std beta1 beta2  mean   gene
# 1  cg14378848    0.233 0.273 0.451 0.218 0.391 ELOVL5
# 2  cg25514947    0.239 0.264 0.486 0.246 0.424 ELOVL5
# 3  cg27599958    0.276 0.243 0.495 0.219 0.423 ELOVL5
# 4  cg13529912    0.191 0.215 0.389 0.198 0.339 ELOVL5
# 5  cg18564099    0.160 0.168 0.323 0.164 0.282 ELOVL5
# 6  cg10369451    0.186 0.168 0.351 0.165 0.303 ELOVL5
# 7  cg03998871    0.227 0.212 0.456 0.229 0.398 ELOVL5
# 8  cg21195414    0.161 0.204 0.517 0.356 0.475 ELOVL5
# 9  cg13636408    0.203 0.235 0.545 0.342 0.492 ELOVL5
# 10 cg15133313    0.185 0.215 0.519 0.335 0.471 ELOVL5
# 11 cg00024396    0.120 0.208 0.531 0.411 0.500 ELOVL5
# 12 cg24524396    0.113 0.179 0.535 0.423 0.506 ELOVL5
# 13 cg20419072    0.115 0.176 0.515 0.400 0.485 ELOVL5
# 14 cg01865606    0.161 0.189 0.378 0.217 0.336 ELOVL5
# 15 cg20397543    0.155 0.191 0.449 0.293 0.409 ELOVL5
# 16 cg10410213    0.193 0.214 0.385 0.192 0.335 ELOVL5
betadiff(dat=DMRs_smoke$table,i=4)
#        probe betadiff   std beta1 beta2  mean  gene
# 1 cg18391209   -0.244 0.303 0.369 0.612 0.432 CAPN8
betadiff(dat=DMRs_smoke$table,i=5)
#        probe betadiff   std beta1 beta2  mean   gene
# 1 cg22623223    0.252 0.239 0.307 0.056 0.242 PTPRN2
# 2 cg07602744    0.236 0.245 0.344 0.108 0.283 PTPRN2
betadiff(dat=DMRs_smoke$table,i=6)
#        probe betadiff   std beta1 beta2  mean   gene
# 1 cg20312821    0.212 0.271 0.478 0.266 0.423 IGDCC3
# 2 cg05054468    0.181 0.223 0.466 0.285 0.419 IGDCC3
# 3 cg03028236    0.235 0.248 0.475 0.240 0.414 IGDCC3
# 4 cg23168483    0.224 0.282 0.510 0.286 0.452 IGDCC3
# 5 cg01107006    0.300 0.258 0.497 0.198 0.420 IGDCC3