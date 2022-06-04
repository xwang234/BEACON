#!/usr/bin/env Rscript

# Figure 1:   a) violin plot for heritability estimates among the ones with p<0.05, for six organs; b) violin plot for r^2 estimates among the ones with p<0.05, for six organs;
# 
# Figure 2:  a) b) c) Manhattam plot for BE vs control, EA vs control, and BE/EA vs control, with a red line FDR=0.05, this is for SKAT
# 
# Supplementary figure 1:  a) the q-q plot for Figure 2
# 
# Figure 3: a) b) c) d) The locus zoom plot for the four genes validated
# 
# Figure 4. Can you plot the scatterplot of OR of individual SNPs (BCA data) vs coefficients of genetic association with gene expression (GTEx)? For the four genes above. Note these are individual SNPs effect, using a model for one snp at a time. Keep the code to pull these four genes, and plot them. Heng will use these codes to do some causal mediation tests
# 
# Table 1: Guo2020 Gastro paper: check their Table 1. We need a table like that, first for four new genes we have (all FDR<0.05). List gene name, locus, SKAT p-value, heritability, r^2, validation p-values (Bonn only), closest previous GWAS snps (distance in MB)
# 
# Table 2: like Table 1, but now list the genes who are not validated or not new.

#Figure1-------------
library("vioplot")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
organs=c("adipose","blood","junction","mucosa","muscularis","stomach")
samplesizes=c(393,558,275,411,385,260)
allheritability=allgenes=data.frame(matrix(nrow=20000,ncol=length(organs)))
colnames(allgenes)=colnames(allheritability)=organs
median_H=rep(NA,6)
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  heritability=heritability[rownames(heritability) %in% proteingenes,]
  heritability=heritability[order(rownames(heritability)),]
  idx=which(heritability$V3 <0.05)
  allheritability[1:length(idx),i]=heritability$V1[idx]
  allgenes[1:length(idx),i]=rownames(heritability)[idx]
}

# idx=complete.cases(allheritability)
# table(idx)
# # FALSE  TRUE 
# # 1697 10051 
tmp=unlist(allheritability)
tmp1=rep(organs,each=nrow(allheritability))
tmp2=sapply(1:ncol(allheritability),function (x) sum(!is.na(allheritability[,x])))
allcolors=c("black","red","blue","darkorchid1","limegreen","goldenrod1","brown4","darkseagreen","darkseagreen1")
#pdf(file="../result/TWAS_figure1a.pdf",width = 12)
png(filename="../result/TWAS_figure1a.png",width=960,type = "cairo")
par(mar=c(3,5,2,1))
vioplot(tmp~tmp1,names=organs,xlab="",ylab="",col=allcolors[1:6],main="",frame.plot=FALSE,cex=1,ylim=c(0,1.3),cex.axis=1.4,cex.main=1.3,yaxt="n")
axis(side = 1, at = 1:6, labels = organs,cex.axis=1.4,cex.lab=1.4)
axis(side = 2, at = seq(0,1,0.2), labels = seq(0,1,0.2),cex.axis=1.4,cex.lab=1.4)
mtext("Heritability",side=2,cex=1.4,line=3)
text(1:6,1.2,labels = paste0(samplesizes," samples"),cex=1.2)
text(1:6,1.1,labels = paste0(tmp2," genes"),cex=1.2)
dev.off()


allr2=allr2genes=data.frame(matrix(nrow=20000,ncol=length(organs)))
colnames(allr2)=colnames(allr2genes)=organs
r2_median=rep(NA,6)
names(r2_median)=organs
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  res_min_code=res_min_code[order(rownames(res_min_code)),]
  idx=which(res_min_code$r2>0.01)
  allr2genes[1:length(idx),i]=rownames(res_min_code)[idx]
  allr2[1:length(idx),i]=res_min_code$r2[idx]
  r2_median[i]=median(res_min_code$r2[idx])
}
r2_median=round(r2_median,3)
sum(duplicated(allr2genes$junction[!is.na(allr2genes$junction)]))
tmp=intersect(allr2genes$junction,allr2genes$mucosa)
length(tmp[!is.na(tmp)])
tmp=intersect(allr2genes$junction,allr2genes$muscularis)
length(tmp[!is.na(tmp)]) #4328
tmp=intersect(allr2genes$mucosa,allr2genes$muscularis)
length(tmp[!is.na(tmp)]) #4317
overlapgene=data.frame(matrix(NA,nrow=6,ncol=6))
rownames(overlapgene)=colnames(overlapgene)=organs
for (i in 1:6)
{
  for (j in 1:6)
  {
    tmp=intersect(allr2genes[,i],allr2genes[,j])
    overlapgene[i,j]=length(tmp[!is.na(tmp)])
  }
}

tmp=unlist(allr2)
tmp1=rep(organs,each=nrow(allr2))
tmp2=sapply(1:ncol(allr2),function (x) sum(!is.na(allr2[,x])))

tmp=unlist(allr2[idx,])
tmp1=rep(organs,each=length(idx))
#pdf(file="../result/TWAS_figure1b.pdf",width=12)
png(file="../result/TWAS_figure1b.png",width=960,type = "cairo")
par(mar=c(3,5,2,1))
vioplot(tmp~tmp1,names=organs,xlab="",ylab="",col=allcolors[1:6],main="",frame.plot=FALSE,cex=1,ylim=c(0,1.3),cex.axis=1.4,cex.main=1.3,yaxt="n")
axis(side = 1, at = 1:6, labels = organs,cex.axis=1.4,cex.lab=1.4)
axis(side = 2, at = seq(0,1,0.2), labels = seq(0,1,0.2),cex.axis=1.4,cex.lab=1.4)
mtext("R-squared",side=2,cex=1.4,line=3)
text(1:6,1.2,labels = paste0(samplesizes," samples"),cex=1.2)
text(1:6,1.1,labels = paste0(tmp2," genes"),cex=1.2)
dev.off()

#draw Venn diagram for overlap of genes across junction/mucosa/muscularis

library(VennDiagram)
library(Cairo)

area1=overlapgene[3,3] #junction
area2=overlapgene[4,4] #mucosa
area3=overlapgene[5,5] #muscularis
n12=overlapgene[3,4]
n23=overlapgene[4,5]
n13=overlapgene[3,5]
tmp=allr2genes[,3]
for (i in 4:5)
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n123=length(tmp)
CairoPNG("../result/venn_3overlapgenes_.png",width = 480,height=480,res=300,pointsize = 2.5)
grid.newpage()
draw.triple.venn(area1=area1,area2=area2,area3=area3,n12=n12,n23=n23,n13=n13,n123=n123,
                 category = c("Junction", "Mucosa","Muscularis"), lty = rep("blank", 3), fill = allcolors[3:5], scaled = T,euler.d=T,
                   cat.cex=rep(3,3),cex=rep(3,7),cat.fontface = rep("bold",3),margin=0.05,resolution = 300,
                 compression = "lzw",cat.dist =c(0.075, 0.075, 0.05),alpha = rep(0.6, 3),fontface = rep("bold", 7))
dev.off()

area1=overlapgene[3,3] #junction
area2=overlapgene[4,4] #mucosa
area3=overlapgene[5,5] #muscularis
area4=overlapgene[6,6] #stomach
n12=overlapgene[3,4]
n13=overlapgene[3,5]
n14=overlapgene[3,6]
n23=overlapgene[4,5]
n24=overlapgene[4,6]
n34=overlapgene[5,6]
tmp=allr2genes[,3]
for (i in 4:5)
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n123=length(tmp)
tmp=allr2genes[,3]
for (i in c(4,6))
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n124=length(tmp)
tmp=allr2genes[,3]
for (i in c(5,6))
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n134=length(tmp)
tmp=allr2genes[,4]
for (i in c(5,6))
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n234=length(tmp)
tmp=allr2genes[,3]
for (i in c(4,5,6))
{
  tmp=intersect(tmp,allr2genes[,i])
}
tmp=tmp[!is.na(tmp)]
n1234=length(tmp)
CairoPNG("../result/venn_4overlapgenes.png",width = 480,height=480,res=300,pointsize = 2.5)
grid.newpage()
draw.quad.venn(area1=area1,area2=area2,area3=area3,area4=area4,n12=n12,n13=n13,n14=n14,n23=n23,n24=n24,n34=n34,n123=n123,n124=n124,n134=n134,n234=n234,n1234=n1234,
                 category = c("Junction", "Mucosa","Muscularis","Stomach"), lty = rep("blank", 4), fill = allcolors[3:6], scaled = TRUE,
                 cat.cex=rep(3,4),cex=rep(3,15),cat.fontface = rep("bold",4),margin=0.05,resolution = 300,
                 compression = "lzw",cat.dist = c(0.25, 0.25, 0.13, 0.13),alpha = rep(0.6, 4),fontface = rep("bold", 15))
dev.off()

#Figure 2----------------------
#proteingenestable
load("../result/proteingenestable.RData")
proteingenestable=proteingenestable[!is.na(proteingenestable$starthg19),]
proteingenestable$Chromosome==proteingenestable$Chromosomehg19
organs=c("adipose","blood","junction","mucosa","muscularis","stomach")
allres=NULL
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  genesr2=rownames(res_min_code)[!is.na(res_min_code$r2) & res_min_code$r2>0.01]
  
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
  skat_min2_pc6=skat_min2_pc6[rownames(skat_min2_pc6) %in% genesr2,]
  res1=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="BEvsCO",p=skat_min2_pc6$BE_p,fdr=p.adjust(skat_min2_pc6$BE_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$BE_p,method="bonferroni"))
  res2=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="EAvsCO",p=skat_min2_pc6$EA_p,fdr=p.adjust(skat_min2_pc6$EA_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$EA_p,method="bonferroni"))
  res3=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="BEEAvsCO",p=skat_min2_pc6$BEEA_p,fdr=p.adjust(skat_min2_pc6$BEEA_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$BEEA_p,method="bonferroni"))
  res=rbind(res1,res2,res3)
  allres=rbind(allres,res)
}

#allres$FWER=p.adjust(allres$p,method="bonferroni")
allres$FDR=p.adjust(allres$p,method="fdr")


#qqplot based on pvaule and all FDR
qqplot=function(res=allres[allres$type=="BEEAvsCO",],main="",xlim=NULL,ylim=NULL)
{
  tmp=data.frame(pvalue=res$p,fdr=res$FDR)
  tmp=tmp[order(tmp$pvalue),]
  n=nrow(tmp)
  par(mar=c(6,6,3,1))
  plot(-log((1:n)/n,base=10),-log(tmp$pvalue[order(tmp$pvalue)],base=10),xlab="Expected p-value (-log base 10)",
       ylab="Observed p-value (-log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4)
  abline(0,1,lty=2)
  idx1=sum(tmp$fdr<0.05)
  if (length(idx1)>0)
  {
    idx=1:idx1
    points(-log((1:n)/n,base=10)[idx],-log(tmp$pvalue[order(tmp$pvalue)],base=10)[idx],pch=16,col="red")
    #print(paste0("#FDR<0.05:",length(idx)))
    legend("topleft",legend = "FDR<0.05",pch=16,col="red",cex=1.3)
  }
}
png(filename = "../result/qqplot_TWAS_BEEA.png",type="cairo")
qqplot()
dev.off()
png(filename = "../result/qqplot_TWAS_EA.png",type="cairo")
qqplot(res=allres[allres$type=="EAvsCO",])
dev.off()
png(filename = "../result/qqplot_TWAS_BE.png",type="cairo")
qqplot(res=allres[allres$type=="BEvsCO",])
dev.off()

#draw q-q plot for different organ/type
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  tmp=data.frame(pvalue=pvalue)
  tmp=tmp[order(tmp$pvalue),,drop=F]
  n=nrow(tmp)
  par(mar=c(6,6,3,1))
  plot(-log((1:n)/n,base=10),-log(tmp$pvalue[order(tmp$pvalue)],base=10),xlab="Expected p-value (-log base 10)",
       ylab="Observed p-value (-log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.4,cex.axis=1.4)
  abline(0,1,lty=2)
}

for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
  skat_min2_pc6=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]
  pdf(file=paste0("../result/",organ,"_QQplot.pdf"),width=15)
  par(mfrow=c(1,3))
  qqplot(pvalue=skat_min2_pc6$BE_p,main="BE")
  qqplot(pvalue=skat_min2_pc6$EA_p,main="EA")
  qqplot(pvalue=skat_min2_pc6$BEEA_p,main="BE/EA")
  dev.off()
}



twasplot_all=function(skatres=allres,ylim=NULL,opt1="nogenes")
{
  allshapes=c(15:18,24,25)
  idx=match(skatres$gene,proteingenestable$Symbol)
  skatres$type=gsub("vsCO","",skatres$type)
  twas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=skatres$gene,organ=skatres$tissue,type=skatres$type,P=NA,BP=proteingenestable$starthg19[idx],fdr=skatres$FDR,stringsAsFactors = F)
  twas.dat$CHR=gsub("chr","",twas.dat$CHR)
  twas.dat$P=skatres$p
  twas.dat=twas.dat[twas.dat$CHR %in% 1:22,]
  twas.dat=twas.dat[!is.na(twas.dat$P),]
  twas.dat$CHR=as.integer(twas.dat$CHR)
  idx=order(twas.dat$CHR,twas.dat$BP)
  twas.dat=twas.dat[idx,]
  
  nCHR <- length(unique(twas.dat$CHR))
  #BPcum is the x coordinate
  twas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(twas.dat$CHR)){
    nbp[i] <- max(twas.dat[twas.dat$CHR == i,]$BP)
    twas.dat[twas.dat$CHR == i,"BPcum"] <- twas.dat[twas.dat$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- twas.dat %>% 
    group_by(CHR) %>% 
    dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(twas.dat$P)))) + 0.5
  }
  sig <- 0.05/nrow(twas.dat)
  print(paste0("sig-p=",sig))
  #twas.dat$fdr=p.adjust(twas.dat$P,method="fdr")
  #twas.dat$fewr=p.adjust(twas.dat$P,method="bonferroni")
  tmp=twas.dat
  tmp=tmp[order(tmp$P,tmp$fdr),]
  fdrcutoff=NULL
  if (sum(twas.dat$fdr<0.05)>0)
  {
    idx=which(tmp$fdr>0.05)
    fdrcutoff=tmp$P[idx[1]]
    fdrcutoff1=round(tmp$fdr[idx[1]],3)
  }
  
  manhplot=NULL
  fweridx=which(twas.dat$P<=sig)
  print(twas.dat$Gene[fweridx])
  fweridx=which(twas.dat$fdr<0.05) #use FDR cutoff
  print(paste0(twas.dat$Gene[fweridx],"_",twas.dat$CHR[fweridx],"_",twas.dat$organ[fweridx],"_",twas.dat$type[fweridx],"_",formatC(twas.dat$P,format="e",digits = 2)[fweridx]))
  fweridx=fweridx[order(twas.dat$P[fweridx],decreasing = T)]
  
  fwergenes=NULL
  if (length(fweridx)>0)
  {
    nlenx=rep(NA,length(fweridx))
    fwergenes=twas.dat$Gene[fweridx]
    fwergenes=paste0(twas.dat$Gene[fweridx],"_",twas.dat$CHR[fweridx],"_",twas.dat$organ[fweridx],"_",twas.dat$type[fweridx])
    #print(fwergenes)
    nlen=sapply(fwergenes,nchar) #length of gene names
    nlen[nlen>7]=7
    fwerx=twas.dat$BPcum[fweridx]
    fwery=rep(0.2,length(fweridx))
    idx1=seq(1,length(fwerx),2)
    idx2=NULL
    if (length(fwerx)>1) idx2=seq(2,length(fwerx),2)
    nlenx[idx1]=20000000*nlen[idx1]
    fwerx[idx1]=fwerx[idx1]-nlenx[idx1]
    fwery[idx1]=fwery[idx1]-0.05
    
    if (length(idx2)>0)
    {
      nlenx[idx2]=15000000*nlen[idx2]
      fwerx[idx2]=fwerx[idx2]+nlenx[idx2]
      fwery[idx2]=fwery[idx2]+0.05
    }
  }
  
  twas.dat$BPcum=twas.dat$BPcum-min(twas.dat$BPcum)
  twas.dat$organ1=factor(twas.dat$organ,levels=organs)
  twas.dat$organ1=allshapes[twas.dat$organ1]
  twas.dat$organ1=as.factor(twas.dat$organ1)
  xmin=min(twas.dat$BPcum)
  xmax=max(twas.dat$BPcum)
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="genes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR),shape=organ1)) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "red", linetype = "dashed",size=1) + 
      #annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      annotate("text", x=120000000+xmin, y=-log10(fdrcutoff)+0.2, label= paste0("FDR=0.05"),size=6) +
      scale_x_continuous(expand = c(0,20000000),label = axis.set$CHR, breaks = axis.set$center,limits=c(xmin,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "SKAT -log10(p)") + 
      annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=2) +
      scale_shape_manual(name = "Organ",
                         labels = organs,
                         values = allshapes) + guides(colour=FALSE)+
      theme_minimal() +
      theme( 
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 5, 10, "mm"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(size = 1, color="black") , 
        axis.ticks.length = unit(.2, "cm"),
        text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust = 0.5)
      )
  }
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="nogenes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR),shape=organ1)) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "red", linetype = "dashed",size=1) + 
      #annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      annotate("text", x=120000000+xmin, y=-log10(fdrcutoff)+0.2, label= paste0("FDR=0.05"),size=6) +
      scale_x_continuous(expand = c(0,20000000),label = axis.set$CHR, breaks = axis.set$center,limits=c(xmin,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "SKAT -log10(p)") + 
      #annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=6) +
      scale_shape_manual(name = "Organ",
                         labels = organs,
                         values = allshapes) + guides(colour=FALSE)+
      theme_minimal() +
      theme( 
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 10, 5, "mm"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(size = 1, color="black") , 
        axis.ticks.length = unit(.2, "cm"),
        text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust = 0.5)
      )
  }
  return(manhplot)
}

pdf(paste0("../result/Manhattan_BE_AllFDR.pdf"),width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="BEvsCO",],opt1="nogenes")
print(manhplot)
dev.off()
pdf("../result/Manhattan_AllFDR_BE_genes.pdf",width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="BEvsCO",],opt1="genes")
print(manhplot)
dev.off()

pdf(paste0("../result/Manhattan_EA_AllFDR.pdf"),width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="EAvsCO",],opt1="nogenes")
print(manhplot)
dev.off()
pdf("../result/Manhattan_AllFDR_EA_genes.pdf",width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="EAvsCO",],opt1="genes")
print(manhplot)
dev.off()

pdf(paste0("../result/Manhattan_BEEA_AllFDR.pdf"),width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="BEEAvsCO",],opt1="nogenes")
print(manhplot)
dev.off()
pdf("../result/Manhattan_AllFDR_BEEA_genes.pdf",width = 16,onefile = F)
manhplot=twasplot_all(skatres=allres[allres$type=="BEEAvsCO",],opt1="genes")
print(manhplot)
dev.off()


#add regular TWAS into Figure 2:
twasplot=function(skatres=allres[allres$type=="BEvsCO",],ylim=NULL,opt1="nogenes")
{
  allshapes=c(15:18,24,25)
  idx=match(skatres$gene,proteingenestable$Symbol)
  skatres$type=gsub("vsCO","",skatres$type)
  twas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=skatres$gene,organ=skatres$tissue,type=skatres$type,P=NA,BP=proteingenestable$starthg19[idx],fdr=skatres$FDR,stringsAsFactors = F)
  twas.dat$CHR=gsub("chr","",twas.dat$CHR)
  twas.dat$P=skatres$p
  twas.dat=twas.dat[twas.dat$CHR %in% 1:22,]
  twas.dat=twas.dat[!is.na(twas.dat$P),]
  twas.dat$CHR=as.integer(twas.dat$CHR)
  idx=order(twas.dat$CHR,twas.dat$BP)
  twas.dat=twas.dat[idx,]
  
  nCHR <- length(unique(twas.dat$CHR))
  #BPcum is the x coordinate
  twas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(twas.dat$CHR)){
    nbp[i] <- max(twas.dat[twas.dat$CHR == i,]$BP)
    twas.dat[twas.dat$CHR == i,"BPcum"] <- twas.dat[twas.dat$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- twas.dat %>% 
    group_by(CHR) %>% 
    dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(twas.dat$P)))) + 0.5
  }
  sig <- 0.05/nrow(twas.dat)
  print(paste0("sig-p=",sig))
  #twas.dat$fdr=p.adjust(twas.dat$P,method="fdr")
  #twas.dat$fewr=p.adjust(twas.dat$P,method="bonferroni")
  tmp=twas.dat
  tmp=tmp[order(tmp$P,tmp$fdr),]
  fdrcutoff=NULL
  if (sum(twas.dat$fdr<0.05)>0)
  {
    idx=which(tmp$fdr>0.05)
    fdrcutoff=tmp$P[idx[1]]
    fdrcutoff1=round(tmp$fdr[idx[1]],3)
  }
  
  manhplot=NULL
  fweridx=which(twas.dat$P<=sig)
  print(twas.dat$Gene[fweridx])
  fweridx=which(twas.dat$fdr<0.05) #use FDR cutoff
  print(paste0(twas.dat$Gene[fweridx],"_",twas.dat$CHR[fweridx],"_",twas.dat$organ[fweridx],"_",twas.dat$type[fweridx],"_",formatC(twas.dat$P,format="e",digits = 2)[fweridx]))
  fweridx=fweridx[order(twas.dat$P[fweridx],decreasing = T)]
  
  fwergenes=NULL
  if (length(fweridx)>0)
  {
    nlenx=rep(NA,length(fweridx))
    fwergenes=twas.dat$Gene[fweridx]
    fwergenes=paste0(twas.dat$Gene[fweridx],"_",twas.dat$CHR[fweridx],"_",twas.dat$organ[fweridx],"_",twas.dat$type[fweridx])
    #print(fwergenes)
    nlen=sapply(fwergenes,nchar) #length of gene names
    nlen[nlen>7]=7
    fwerx=twas.dat$BPcum[fweridx]
    fwery=rep(0.2,length(fweridx))
    idx1=seq(1,length(fwerx),2)
    idx2=NULL
    if (length(fwerx)>1) idx2=seq(2,length(fwerx),2)
    nlenx[idx1]=20000000*nlen[idx1]
    fwerx[idx1]=fwerx[idx1]-nlenx[idx1]
    fwery[idx1]=fwery[idx1]-0.05
    
    if (length(idx2)>0)
    {
      nlenx[idx2]=15000000*nlen[idx2]
      fwerx[idx2]=fwerx[idx2]+nlenx[idx2]
      fwery[idx2]=fwery[idx2]+0.05
    }
  }
  
  twas.dat$BPcum=twas.dat$BPcum-min(twas.dat$BPcum)
  twas.dat$organ1=factor(twas.dat$organ,levels=organs)
  twas.dat$organ1=allshapes[twas.dat$organ1]
  twas.dat$organ1=as.factor(twas.dat$organ1)
  xmin=min(twas.dat$BPcum)
  xmax=max(twas.dat$BPcum)
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="genes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR),shape=organ1)) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "red", linetype = "dashed",size=1) + 
      #annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      annotate("text", x=120000000+xmin, y=-log10(fdrcutoff)+0.2, label= paste0("FDR=0.05"),size=6) +
      scale_x_continuous(expand = c(0,20000000),label = axis.set$CHR, breaks = axis.set$center,limits=c(xmin,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "SKAT -log10(p)") + 
      annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=2) +
      scale_shape_manual(name = "Organ",
                         labels = organs,
                         values = allshapes) + guides(colour=FALSE)+
      theme_minimal() +
      theme( 
        legend.position = "right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 5, 10, "mm"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(size = 1, color="black") , 
        axis.ticks.length = unit(.2, "cm"),
        text = element_text(size=20),
        axis.text.x = element_text(angle = 0, vjust = 0.5)
      )
  }
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="nogenes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR),shape=organ1)) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "red", linetype = "dashed",size=1) + 
      #annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      annotate("text", x=120000000+xmin, y=-log10(fdrcutoff)+0.4, label= paste0("FDR=0.05"),size=6) +
      scale_x_continuous(expand = c(0,20000000),label = axis.set$CHR, breaks = axis.set$center,limits=c(xmin,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "SKAT -log10(p)") + 
      #annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=6) +
      scale_shape_manual(name = "Organ",
                         labels = organs,
                         values = allshapes) + guides(colour=FALSE)+
      theme_minimal() +
      theme( 
        legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #plot.margin = margin(0, 0, 10, 5, "mm"),
        plot.margin = margin(4, 5, 1, 5, "mm"),
        axis.line.y = element_line(colour = 'black', size = 1),
        axis.ticks.y = element_line(size = 1, color="black") , 
        axis.ticks.length = unit(.2, "cm"),
        text = element_text(size=20)
        # axis.line = element_line(colour = 'black', size = 1),
        # axis.ticks = element_line(size = 1, color="black") , 
        # axis.ticks.length = unit(.2, "cm"),
        # text = element_text(size=16),
        # axis.text.x = element_text(angle = 0, vjust = 0.5)
      )
  }
  return(list(manhplot=manhplot,fdrcutoff=fdrcutoff))
}

organs=c("adipose","blood","junction","mucosa","muscularis","stomach")
allres=NULL
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  genesr2=rownames(res_min_code)[!is.na(res_min_code$r2) & res_min_code$r2>0.01]
  
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
  skat_min2_pc6=skat_min2_pc6[rownames(skat_min2_pc6) %in% genesr2,]
  res1=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="BEvsCO",p=skat_min2_pc6$BE_p,fdr=p.adjust(skat_min2_pc6$BE_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$BE_p,method="bonferroni"))
  res2=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="EAvsCO",p=skat_min2_pc6$EA_p,fdr=p.adjust(skat_min2_pc6$EA_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$EA_p,method="bonferroni"))
  res3=data.frame(tissue=rep(organ,nrow(skat_min2_pc6)),gene=rownames(skat_min2_pc6),type="BEEAvsCO",p=skat_min2_pc6$BEEA_p,fdr=p.adjust(skat_min2_pc6$BEEA_p,method="fdr"),fwer=p.adjust(skat_min2_pc6$BEEA_p,method="bonferroni"))
  res=rbind(res1,res2,res3)
  allres=rbind(allres,res)
}

#allres$FWER=p.adjust(allres$p,method="bonferroni")
allres$FDR=p.adjust(allres$p,method="BH")

#regular TWAS
allres1=NULL
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  genesr2=rownames(res_min_code)[!is.na(res_min_code$r2) & res_min_code$r2>0.01]
  
  load(paste0(outfolder,"/bca_assoc.RData"))
  assoc_min_pc6=assoc_min_pc6[rownames(assoc_min_pc6) %in% genesr2,]
  res1=data.frame(tissue=rep(organ,nrow(assoc_min_pc6)),gene=rownames(assoc_min_pc6),type="BEvsCO",p=assoc_min_pc6$BE_p,fdr=p.adjust(assoc_min_pc6$BE_p,method="fdr"),fwer=p.adjust(assoc_min_pc6$BE_p,method="bonferroni"))
  res2=data.frame(tissue=rep(organ,nrow(assoc_min_pc6)),gene=rownames(assoc_min_pc6),type="EAvsCO",p=assoc_min_pc6$EA_p,fdr=p.adjust(assoc_min_pc6$EA_p,method="fdr"),fwer=p.adjust(assoc_min_pc6$EA_p,method="bonferroni"))
  res3=data.frame(tissue=rep(organ,nrow(assoc_min_pc6)),gene=rownames(assoc_min_pc6),type="BEEAvsCO",p=assoc_min_pc6$BEEA_p,fdr=p.adjust(assoc_min_pc6$BEEA_p,method="fdr"),fwer=p.adjust(assoc_min_pc6$BEEA_p,method="bonferroni"))
  res=rbind(res1,res2,res3)
  allres1=rbind(allres1,res)
}

#allres1$FWER=p.adjust(allres1$p,method="bonferroni")
allres1$FDR=p.adjust(allres1$p,method="BH")
idx=match(allres1$gene,proteingenestable$Symbol)
allres1$chr=proteingenestable$Chromosomehg19[idx]
plot(-log10(allres1$p[allres1$type=="BEvsCO"]),-log10(allres$p[allres$type=="BEvsCO"]),main="BE",xlab="TWAS",ylab="SKAT")
abline(0,1,col="red")
plot(-log10(allres1$p[allres1$type=="EAvsCO"]),-log10(allres$p[allres$type=="EAvsCO"]),main="EA",xlab="TWAS",ylab="SKAT")
abline(0,1,col="red")
plot(-log10(allres1$p[allres1$type=="BEEAvsCO"]),-log10(allres$p[allres$type=="BEEAvsCO"]),main="BEEA",xlab="TWAS",ylab="SKAT")
abline(0,1,col="red")

table(allres$p[allres$type=="BEvsCO"]<allres1$p[allres1$type=="BEvsCO"])
table(allres$p[allres$type=="EAvsCO"]<allres1$p[allres1$type=="EAvsCO"])
table(allres$p[allres$type=="BEEAvsCO"]<allres1$p[allres1$type=="BEEAvsCO"])
sum(allres$p[allres$type=="BEvsCO"]<0.05) #3033
sum(allres1$p[allres1$type=="BEvsCO"]<0.05) #2333
sum(allres$p[allres$type=="EAvsCO"]<0.05) #2721
sum(allres1$p[allres1$type=="EAvsCO"]<0.05) #2285
sum(allres$p[allres$type=="BEEAvsCO"]<0.05) #3307
sum(allres1$p[allres1$type=="BEEAvsCO"]<0.05) #2509


idx=which(allres$p<0.05 & allres1$p<0.05)
plot(-log10(allres1[idx,]$p[allres1$type=="BEvsCO"]),-log10(allres[idx,]$p[allres$type=="BEvsCO"]),main="BE",xlab="TWAS",ylab="SKAT")
abline(0,1,col="red")

library(dplyr)
gwasplot=function(gwasres=allres1[allres1$type=="BEvsCO",],ylim=NULL,sig=NULL)
{
  allshapes=c(15:18,24,25)
  idx=match(gwasres$gene,proteingenestable$Symbol)
  gwasres$type=gsub("vsCO","",gwasres$type)
  gwas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=gwasres$gene,organ=gwasres$tissue,type=gwasres$type,P=NA,BP=proteingenestable$starthg19[idx],fdr=gwasres$FDR,stringsAsFactors = F)
  gwas.dat$CHR=gsub("chr","",gwas.dat$CHR)
  gwas.dat$P=gwasres$p
  gwas.dat=gwas.dat[gwas.dat$CHR %in% 1:22,]
  gwas.dat=gwas.dat[!is.na(gwas.dat$P),]
  gwas.dat$CHR=as.integer(gwas.dat$CHR)
  idx=order(gwas.dat$CHR,gwas.dat$BP)
  gwas.dat=gwas.dat[idx,]
  
  nCHR <- length(unique(gwas.dat$CHR))
  #BPcum is the x coordinate
  gwas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(gwas.dat$CHR)){
    nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
    gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- gwas.dat %>% 
    group_by(CHR) %>% 
    dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(gwas.dat$P)))) + 0.5
  }
  if (is.null(sig)) sig <- 5e-8
  
  #x lables not show
  idx=axis.set$CHR %in% c(19,21)
  uselabels=axis.set$CHR
  #uselabels[idx]=""
  gwas.dat$logP=-log10(gwas.dat$P)
  gwas.dat$Color=as.factor(gwas.dat$CHR)
  gwas.dat$organ1=factor(gwas.dat$organ,levels=organs)
  gwas.dat$organ1=allshapes[gwas.dat$organ1]
  gwas.dat$organ1=as.factor(gwas.dat$organ1)
  xmin=min(gwas.dat$BPcum)
  xmax=max(gwas.dat$BPcum)
  manhplotr <- ggplot(gwas.dat, aes(x = BPcum, y = logP, 
                                    color = Color, shape=organ1)) +
    geom_point() +
    #geom_rangeframe(data = data.frame(BPcum = c(axis.set$center[1], axis.set$center[22]), logP = c(ylim, 0)))+
    #geom_rangeframe(color="black")+
    #coord_cartesian(xlim=c(axis.set$center[1],max(axis.set$center)),expand = F)+
    #geom_segment(aes_all(c('x', 'y', 'xend', 'yend')),
    #            data = data.frame(x = c(0,axis.set$center[1]), xend = c(0, max(axis.set$center)), y = c(0, 0), yend = c(ylim, 0))) +
    
    #geom_hline(yintercept = -log10(sig), color = "gray", linetype = "dashed",size=1) + 
    #scale_x_continuous(label = uselabels, breaks = axis.set$center) +
    scale_x_continuous(expand = c(0,20000000),label = axis.set$CHR, breaks = axis.set$center,limits=c(xmin,xmax)) +
    scale_y_continuous(trans = "reverse",limits = c(ylim,0 )) +
    #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
    scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosome", 
         y = "TWAS -log10(p)") + 
    scale_shape_manual(name = "Organ",
                       labels = organs,
                       values = allshapes) + guides(colour=FALSE)+
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # panel.grid.major.y = element_blank(),
      # panel.grid.minor.y = element_blank(),
      plot.margin = margin(0, 5, 10, 5, "mm"),
      axis.line = element_line(colour = 'black', size = 1),
      axis.ticks = element_line(size = 1, color="black") , 
      axis.ticks.length = unit(.2, "cm"),
      text = element_text(size=20),
      axis.text.x = element_text(angle = 0, vjust = 0.5)
    )
  return(manhplotr)
}

plot2manhattan=function(skatres=allres[allres$type=="BEvsCO",],opt1="nogenes",
                        gwasres=allres1[allres1$type=="BEvsCO",],
                        prefix="BE",ylim=6.5)
{
  manhplot=twasplot(skatres=skatres,opt1=opt1,ylim=ylim)
  manhplotr=gwasplot(gwasres=gwasres,ylim=ylim,sig = manhplot$fdrcutoff)
  #png("../result/test.png",width=900,res=100,pointsize = 2)
  pdf(paste0("../result/Manhattan2_",prefix,".pdf"),width = 16,onefile = F)
  # library(grDevices)
  #cairo_ps(filename = "../result/test.eps",onefile = F,fallback_resolution = 600)
  egg::ggarrange(manhplot$manhplot, manhplotr, heights = c(0.5, 0.5))
  dev.off()
}
plot2manhattan()
plot2manhattan(skatres=allres[allres$type=="EAvsCO",],opt1="nogenes",
                        gwasres=allres1[allres1$type=="EAvsCO",],
                        prefix="EA",ylim=6.5)

plot2manhattan(skatres=allres[allres$type=="BEEAvsCO",],opt1="nogenes",
               gwasres=allres1[allres1$type=="BEEAvsCO",],
               prefix="BEEA",ylim=9.5)

#Figure 3----------------------
#code saved in locuszoom2.R


#Figure 4----------------------
#load BCA covariate tables
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  if (length(idx)>0)
    sampletable[idx,i]=NA
}

#BCA covariate table, principle components
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
#Covariate is used for SKAT

BCAcovariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]


mediationplot=function(genename="HSP90AA1",organ="blood",type="BE")
{
  #genetic effect from GTEx
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  idx=which(rownames(res_min)==genename)
  #all the selected snps for that gene
  selectedsnps=unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T))
  #GTEXdata
  load(paste0("../result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")) #snp,phenotype...
  GTEx_coeff=rep(NA,length(selectedsnps))
  names(GTEx_coeff)=selectedsnps
  all(colnames(snp)==rownames(covariate))
  all(colnames(snp)==colnames(phenotype))
  
  correctedsnps=NULL #snps to be recoding
  selectedsnps1=selectedsnps
  for (i in 1:length(selectedsnps))
  {
    idx=which(rownames(snp)==selectedsnps[i])
    tmp1=data.frame(snp=snp[idx,])
    colnames(tmp1)=colnames(snp)
    tmp2=t(covariate)
    dat=rbind(tmp1,tmp2)
    idx=which(rownames(phenotype)==genename)
    y=unlist(phenotype[idx,])
    dat1=as.data.frame(t(dat))
    for (j in c(1:20,22))
    {
      dat1[,j]=as.numeric(dat1[,j])
    }
    colnames(dat1)[1]="snp"
    fit=glm(y~.,data=dat1)
    coeffs=coefficients(fit)
    if ("snp" %in% names(coeffs))
    {
      if (coeffs[2]>0) #make sure effect size is positive
      {
        GTEx_coeff[i]=as.numeric(coeffs[2])
      }else #effect size is negative, need to recode snp
      {
        GTEx_coeff[i]=-as.numeric(coeffs[2])
        tmp1=unlist(strsplit(selectedsnps1[i],"_"))
        selectedsnps1[i]=paste(tmp1[c(1,3,2)],collapse="_")
        correctedsnps=c(correctedsnps,i)
      }
    }
  }
  
  #BCAdata
  #OR from BCA
  BCA_OR=rep(NA,length(selectedsnps)) #log odds ratio
  names(BCA_OR)=selectedsnps
  load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved BCA genotype data, bcabim,bcagenotype
  idx=match(colnames(bcagenotype),rownames(BCAcovariate))
  if (type=="BE")
  {
    y=covariatetable$phenoBE_bca[idx]
  }
  if (type=="EA")
  {
    y=covariatetable$phenoEA_bca[idx]
  }
  if (type=="BEEA")
  {
    y=covariatetable$phenoEABE_bca[idx]
  }
  y=y-1
  
  tmp1=paste0(bcabim$V1,":",bcabim$V4,"_",bcabim$V5,"_",bcabim$V6)
  tmp2=paste0(bcabim$V1,":",bcabim$V4,"_",bcabim$V6,"_",bcabim$V5) #snp coding in BCA is different, need to change back
  idx1=match(selectedsnps,tmp1)
  idx2=match(selectedsnps,tmp2)
  idx3=idx1
  idx3[is.na(idx3)]=idx2[is.na(idx3)]
  if(sum(is.na(idx3))>0) stop("Some snps are missing in BCA")
  selectedbcagenotype=bcagenotype[idx3,]
  idx4=which(is.na(idx1))
  if (length(idx4)>0) selectedbcagenotype[idx4,]=2-selectedbcagenotype[idx4,] #make sure snp coding in BCA is the same as in GTEx
  
  
  for (i in 1:length(selectedsnps))
  {
    tmp1=data.frame(snp=selectedbcagenotype[i,])
    dat2=cbind(t(tmp1),BCAcovariate[idx,])
    colnames(dat2)[1]="snp"
    if (i %in% correctedsnps)
    {
      dat2$snp=2-dat2$snp
    }
    fit=glm(y~.,data=dat2,family = "binomial")
    coeffs1=coefficients(fit)
    if ("snp" %in% names(coeffs1)) BCA_OR[i]=as.numeric(coeffs1[2])
  }
  par(mar=c(6,6,3,1))
  plot(GTEx_coeff,BCA_OR,xlab="GTEx genetic effect on expression",ylab="BCA genetic effect (Log OR) on trait",cex.lab=1.3,cex.axis=1.3,
       ylim=c(min(BCA_OR),max(BCA_OR)*1.1),main=paste0(genename,", ",type),cex.main=1.3)
  fit=lm(BCA_OR~GTEx_coeff)
  abline(fit,col="red")
  tmp=summary(fit)$coefficients
  tmp1=par("usr")
  txt=paste0("y=",round(tmp[1,1],2),"+",round(tmp[2,1],2),"x, p value for coeff x is ",round(tmp[2,4],4))
  text(tmp1[1]+0.4*(tmp1[2]-tmp1[1]),tmp1[3]+0.85*(tmp1[4]-tmp1[3]),txt,cex=1.2)
  return(list(GTEx_coeff=GTEx_coeff,BCA_OR=BCA_OR))
}
ex1=mediationplot(genename="HSP90AA1",organ="blood",type="BE")
ex1_=mediationplot(genename="HSP90AA1",organ="blood",type="BEEA")
ex2=mediationplot(genename="SENP6",organ="blood",type="BEEA")
ex2_=mediationplot(genename="SENP6",organ="blood",type="BE")
ex3=mediationplot(genename="ZNF641",organ="junction",type="BEEA")
ex3_=mediationplot(genename="ZNF641",organ="junction",type="BE")
ex4=mediationplot(genename="CFDP1",organ="stomach",type="BEEA")
ex4_=mediationplot(genename="CFDP1",organ="stomach",type="BE")
ex5=mediationplot(genename="EXOC3",organ="adipose",type="BEEA")
ex6=mediationplot(genename="JUND",organ="blood",type="BEEA")
save(ex1,ex2,ex3,ex4,file="../result/TWAS_figure4.RData")

#Table 1--------------------------
#GWAS snps
dong23snp=read.table("../data/Dong23snp.txt",header=T,sep="\t",stringsAsFactors = F)
extra3snp=read.table("../result/Extra3SNPs.txt",header = T)
dong26snp=rbind(dong23snp[,1:3],extra3snp)
#write.table(dong26snp,file="../data/Dong26snp.txt",row.names = F,col.names = T,quote=F,sep="\t")
tmp=read.table("../result/Dong26SNPs_postion.txt")
tmp$pos=NA
for (i in 1:nrow(tmp))
{
  tmp$pos[i]=unlist(strsplit(tmp$V1[i],"-"))[2]
}
all(tmp$pos==dong26snp$Position) #T
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg19ToHg38.over.chain")

gr_dong26snp=GRanges(seqnames = paste0("chr",dong26snp$Chr),ranges=IRanges(start=dong26snp$Position,width=1))
tmp=liftOver(gr_dong26snp,chain)
dong26snp$pos38=NA
for (i in 1:length(tmp))
{
  tmp1=tmp[i]
  if (length(tmp1)>0)
    dong26snp$pos38[i]=as.numeric(start(tmp1))
}

#saved in NewGenesValidate3.R
allgenelist=read.csv("../result/NEWGenesValidate_allcovar_ALLFDR.csv")
#add missing SNPs in Bonn using closest SNPs with dist<=50BP 
allgenelist=read.csv("../result/NEWGenesValidate_allcovar_ALLFDR_1213.csv")
#use cor>0.6
allgenelist=read.csv("../result/NEWGenesValidate_allcovar_ALLFDR_01_06.csv")
update_genelist=function(genelist=allgenelist,adjustopt="PC6",r2cutoff=0)
{
  #genelist$fwer_addgwas=genelist$pvalue_addgwas=genelist$BEEA_fwer=genelist$BEEA_pvalue=genelist$BEA_fwer=genelist$BEA_pvalue=genelist$EA_fwer=genelist$EA_pvalue=genelist$BE_fwer=genelist$BE_pvalue=NA
  genelist$pvalue_addgwas=NA
  for (i in 1:nrow(genelist))
  {
    if (genelist$tissue[i]=="adipose")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="blood")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="junction")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="mucosa")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="muscularis")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor_allcovar"
    }
    if (genelist$tissue[i]=="stomach")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor_allcovar"
    }
    #work on pvalues
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/skat_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    # skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    # skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    # skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    # skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    # skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    # skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    # idx=which(rownames(skat_min2_code)==genelist$gene[i])
    # genelist$BEEA_pvalue[i]=skat_min2_code$BEEA_pvalue[idx]
    # genelist$BEEA_fwer[i]=skat_min2_code$BEEA_fwer[idx]
    # genelist$BEA_pvalue[i]=skat_min2_code$BEA_pvalue[idx]
    # genelist$BEA_fwer[i]=skat_min2_code$BEA_fwer[idx]
    # genelist$EA_pvalue[i]=skat_min2_code$EA_pvalue[idx]
    # genelist$EA_fwer[i]=skat_min2_code$EA_fwer[idx]
    # genelist$BE_pvalue[i]=skat_min2_code$BE_pvalue[idx]
    # genelist$BE_fwer[i]=skat_min2_code$BE_fwer[idx]
    
    #work on gwas snp adjusted pvalues
    rm(skat_min2_pc6,skat_min2)
    load(paste0(outfolder,"/skat_gwas_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    genesr2=rownames(res_min)[which(res_min$r2>r2cutoff)]
    skat_min2_code=skat_min2_code[rownames(skat_min2_code) %in% genesr2,]
    skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    idx=which(rownames(skat_min2_code)==genelist$gene[i])
    tmp=unlist(strsplit(genelist$type[i],";"))
    tmp1=rep(NA,length(tmp))
    for (j in 1:length(tmp))
    {
      if ("EAvsBE"==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BEA_pvalue[idx],digits = 3)
      }
      if ("EAvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$EA_pvalue[idx],digits = 3)
      }
      if ("BEvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BE_pvalue[idx],digits = 3)
      }
      if ("BEEAvsCO" ==tmp[j])
      {
        tmp1[j]=signif(skat_min2_code$BEEA_pvalue[idx],digits = 3)
      }
    }
    genelist$pvalue_addgwas[i]=paste0(tmp1,collapse = ";")
    # if ("EAvsBE" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BEA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BEA_fwer[idx]
    # }
    # if ("EAvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$EA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$EA_fwer[idx]
    # }
    # if ("BEvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BE_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BE_fwer[idx]
    # }
    # if ("BEEAvsCO" %in% unlist(strsplit(genelist$type[i],";")))
    # {
    #   genelist$pvalue_addgwas[i]=skat_min2_code$BEEA_pvalue[idx]
    #   genelist$fwer_addgwas[i]=skat_min2_code$BEEA_fwer[idx]
    # }
  }
  return(genelist)
}

allgenelist1=update_genelist()

#collapse genes
update_genelist1=function(genelist=allgenelist1)
{
  allBCAtypes=c("BE","EA","BEEA")
  allres=NULL
  for (i in 1:nrow(genelist))
  {
    tmp1=unlist(strsplit(genelist$type[i],";"))
    res=NULL
    for (j in 1:length(tmp1))
    {
      res=rbind(res,genelist[i,])
    }
    
    tmp2=unlist(strsplit(genelist$AllFWER[i],";"))
    tmp3=unlist(strsplit(genelist$AllFDR[i],";"))
    tmp4=unlist(strsplit(genelist$pvalue_addgwas[i],";"))
    for (j in 1:length(tmp1))
    {
      res$type[j]=tmp1[j]
      res$AllFWER[j]=tmp2[j]
      res$AllFDR[j]=tmp3[j]
      res$pvalue_addgwas[j]=tmp4[j]
    }
    allres=rbind(allres,res)
  }
  allres$type=gsub("vsCO","",allres$type)
  #keep 1 result for 1 gene
  allres1=NULL
  allgenes=unique(allres$gene)
  for (i in 1:length(allgenes))
  {
    idx=which(allres$gene==allgenes[i])
    
    # tmp1=allBCAtypes[which.min(allres[idx[1],c("skatpBE","skatpEA","skatpBEEA")])]
    # idx1=match(tmp1,allres$type[idx]) #pick the minimum p-value for a gene
    idx1=which.min(allres$AllFWER[idx]) 
    allres1=rbind(allres1,allres[idx[idx1],])
  }
  allres1$pvalue_addgwas=as.numeric(allres1$pvalue_addgwas)
  allres1$AllFDR=as.numeric(allres1$AllFDR)
  return(allres1)
}
allgenelist2=update_genelist1()
#add regular twas p-value
add_regP=function(dat=allgenelist2)
{
  dat$TWASp=NA
  for (i in 1:nrow(dat))
  {
    organ=dat$tissue[i]
    outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
    load(paste0(outfolder,"/bca_assoc.RData"))
    idx=which(rownames(assoc_min_pc6)==dat$gene[i])
    dat$TWASp[i]=assoc_min_pc6[idx,paste0(dat$type[i],"_p")]
  }
  return(dat)
}
allgenelist3=add_regP()

#TWASp value for BCA and Bonn
library(survcomp)
add_metaTWASp=function(dat=allgenelist3)
{
  dat$BCABonnTWASp=NA
  for (i in 1:nrow(dat))
  {
    if (dat$type[i]=="BE")
     p=c(dat$skatpBE[i],dat$Bonn_BEpsatterthwaite[i])
    if (dat$type[i]=="EA")
      p=c(dat$skatpEA[i],dat$Bonn_EApsatterthwaite[i])
    if (dat$type[i]=="BEEA")
      p=c(dat$skatpBEEA[i],dat$Bonn_BEEApsatterthwaite[i])
    
    dat$BCABonnTWASp[i]=combine.test(p,weight=c(15925,6183),method="fisher")
  }
  return(dat)
}

allgenelist3=add_metaTWASp()
tmp=data.frame(allgenelist3[,c("gene","BCABonnTWASp")])
tmp
#missinggwassnp had been fixed
missinggwassnp=dong26snp$SNP[!dong26snp$SNP %in% colnames(allgwassnps)] #[1] "rs9918259"  "rs75783973"
#hg19
#5:663092-663092 5:668309-668309
#hg38
#5:662977 5:668194

#dong26snp_addcontrols_genotypedat
allgwassnps=read.table("../result/dong26snp_addcontrols_genotypedat.txt")
# cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand19.txt',header=T)
cytoband=read.table('/fh/fast/dai_j/CancerGenomics/Tools/database/other/cytoBand38.txt',fill=T)
cytoband=cytoband[cytoband$V1 %in% paste0("chr",c(1:22,"X","Y")),]
gr_cytoband=GRanges(seqnames = cytoband$V1,ranges=IRanges(start=cytoband$V2,end=cytoband$V3),cytoband=cytoband$V4)

gettable1=function(genelist=allgenelist2[allgenelist2$AllFDR<0.05 & allgenelist2$pvalue_addgwas<0.05,])
{
  
  gr_genelist=GRanges(seqnames = paste0("chr",genelist$chr),IRanges(start=genelist$position,width=1))
  allres=NULL
  for (i in 1:nrow(genelist))
  {
    tmp1=genelist$type[i]
    tmp1=gsub("vsCO","",tmp1) #type
    tmp2=genelist[i,paste0("skatp",tmp1)] #skat p-value
    tmp3=min(genelist[i,c("Bonn_BEpsatterthwaite","Bonn_BEpsatterthwaite","Bonn_EApsatterthwaite","Bonn_BEEApsatterthwaite")]) #Bonn-validated
    tmp4=genelist[i,"Oxford_BEpsatterthwaite"] #Oxford-validated
    tmp5=genelist$pvalue_addgwas[i] #p-value add gwas hits
    res=data.frame(Locus=NA,gene=genelist$gene[i],Numsnp=genelist$numsnp[i],R2=genelist$R2[i],Heritability=genelist$H[i],Leedsnp=genelist$gwassnp[i],Distance=genelist$dist_gwassnp[i],Pvalue=tmp2,FDR=genelist$AllFDR[i],Pvalue_afteradjustSNP=tmp5,BonnPvalue=tmp3,OxfordPvalue=tmp4)
    tmp=distance(gr_cytoband,gr_genelist[i])
    idx=which(tmp==0)
    res$Locus=paste0(genelist$chr[i],cytoband$V4[idx])
    allres=rbind(allres,res)
  }
  allres$Validated=F
  idx=which(allres$BonnPvalue<0.05 | allres$OxfordPvalue<0.05)
  allres$Validated[idx]=T
  idx=which(is.na(allres$Leedsnp))
  #remove adj_gwas if it is new
  allres$Pvalue_afteradjustSNP[idx]=NA
  idx=order(allres$Validated,allres$Pvalue,decreasing = c(T,F))
  
  return(allres[idx,])
}

table1=gettable1()
#write.csv(table1,file="../result/TWAS_table1.csv",row.names = F)

table2=gettable1(genelist=allgenelist2[allgenelist2$AllFDR<0.05 & allgenelist2$pvalue_addgwas>0.05,])
#write.csv(table2,file="../result/TWAS_table2.csv",row.names = F)

#Table 3, validation results
genes=c("EXOC3",
        "SENP6",
        "KRTAP5-8",
        "ZNF641",
        "HSP90AA1",
        "CFDP1",
        "CHST5",
        "BCAR1"
)
BonnP=NULL
mytypes=c("BE","EA","BEEA")
for (i in 1:length(genes))
{
  idx=which(allgenelist3$gene==genes[i])
  # idx1=which.min(allgenelist3[idx,c("Bonn_BEpsatterthwaite","Bonn_EApsatterthwaite","Bonn_BEEApsatterthwaite")])
  # mytype=mytypes[idx1]
  mytype=allgenelist3$type[idx]
  tmp1=data.frame(gene=genes[i],trait=mytype,satterthwaite=allgenelist3[idx,paste0("Bonn_",mytype,"psatterthwaite")],
                 saddlepoint=allgenelist3[idx,paste0("Bonn_",mytype,"psaddlepoint")])
  BonnP=rbind(BonnP,tmp1)
}
BonnP
#get number of p-value<0.05 and minimum p-value
BonnP1=data.frame(gene=genes,numsnp=NA,minp=NA,nsnp05=NA,minsnp=NA)
#get BE_Bonnsummarydat ... in NewGenesValidate3.R
for (i in 1:length(genes))
{
  idx0=which(allgenelist3$gene==genes[i])
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",allgenelist3$tissue[idx0],"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  idx1=which(rownames(res_min)==genes[i])
  selsnps=unlist(strsplit(res_min$selectedsnps[idx1],"|",fixed=T))
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  rsid=tmp$snp
  if (sum(is.na(rsid))>0) print(paste0(as.character(sum(is.na(rsid)))," snps can't find snp rsid"))
  library(survey)
  if (allgenelist3$type[idx0]=="BE")
  {
    dat=BE_Bonnsummarydat
  }
  if (allgenelist3$type[idx0]=="EA")
  {
    dat=EA_Bonnsummarydat
  }
  if (allgenelist3$type[idx0]=="BEEA")
  {
    dat=BEEA_Bonnsummarydat
  }
  tmp=intersect(dat$SNP,rsid[!is.na(rsid)])
  idx2=match(tmp,dat$SNP)
  BonnP1$numsnp[i]=length(tmp)
  BonnP1$minp[i]=min(dat$P[idx2],na.rm=T)
  BonnP1$minsnp[i]=dat$SNP[idx2][which.min(dat$P[idx2])]
  BonnP1$nsnp05[i]=sum(dat$P[idx2]<0.05,na.rm=T)
}
idx=match(genes,allgenelist3$gene)
tmp=allgenelist3[idx,]
table3=data.frame(gene=tmp$gene,Trait=BonnP$trait,GTEXnsnp=tmp$numsnp,Bonnsnp=tmp$Bonn_numsnp,minp=BonnP1$minp,minsnp=BonnP1$minsnp,nsnp05=BonnP1$nsnp05,Bonn_P=BonnP$satterthwaite,Hochberg_p=p.adjust(BonnP$satterthwaite,method="hochberg"))
#write.csv(table3,file="../result/TWAS_table3.csv",row.names = F)

table3$organ=c("adipose","blood","adipose","junction","blood","stomach","junction","blood")
#add pvalue in BC
table3$BCminp=NA
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
#Covariate is used for SKAT
opt="PC6"
#try different number of PCs
if (opt=="PC4")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","sex")]
if (opt=="PC6")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]


findrsid=function(selsnps=NULL,verbose=0)
{
  print(paste0(length(selsnps)," snprsid to be found"))
  chr=unlist(strsplit(selsnps[1],":"))[1]
  pos=allele1=allele2=rep(NA,length(selsnps))
  res=data.frame(selsnps=selsnps,rsid=NA,ref=NA,alt=NA,gene=NA)
  idx=which(is.na(res$rsid))
  n=0
  while(n<5 & length(idx)>0)
  {
    n=n+1
    if (n %%2==0) print(n)
    for (j in idx)
    {
      tmp1=unlist(strsplit(selsnps[j],":"))
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      res$ref[j]=allele2[j]
      res$alt[j]=allele1[j]
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand','associated_gene'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = F)
        },
        error=function(e)
        {
          return(F)
        }
      )
      Sys.sleep(1)
      if (verbose==1 & class(tmp2)[1]=="logical") print(paste0(selsnps[j]," getBM can't work!"))
      #if (class(tmp2)[1]=="logical") print(paste0(selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            res$rsid[j]=tmp2$refsnp_id[k]
            res$gene[j]=tmp2$associated_gene[k]
            break
          }
        }
      }
    }
    idx=which(is.na(res$rsid))
  }
  idx=which(is.na(res$rsid))
  if (length(idx)>0) print(paste0(length(idx)," snps not found rsid"))
  return(res)
}

all(table3$minsnp %in% BE_Oxfordsummarydat$SNP)

for (i in 1:nrow(table3))
{
  organ=table3$organ[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
  idx1=match(colnames(bcagenotype),rownames(covariatetable))
  covariatetable=covariatetable[idx1,]
  idx1=match(colnames(bcagenotype),rownames(Covariate))
  Covariate=Covariate[idx1,]
  gene=table3$gene[i]
  idx=which(rownames(res_min) ==gene)
  selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
  
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  selsnspinfo=tmp
  rsid=tmp$snp
  if (sum(is.na(rsid))>0) print(paste0(as.character(sum(is.na(rsid)))," snps can't find snp rsid"))
  library(survey)
  type=table3$Trait[i]
  if (type=="BE")
  {
    dat=BE_Bonnsummarydat
  }
  if (type=="EA")
  {
    dat=EA_Bonnsummarydat
  }
  if (type=="BEEA")
  {
    dat=BEEA_Bonnsummarydat
  }
  METALdat=data.frame(MarkerName=selsnspinfo$snphg38,snp=selsnspinfo$snp,chr=NA,position=NA,Pvalue=NA,PvalueBonn=NA,stringsAsFactors = F)
  for (j in 1:nrow(METALdat))
  {
    tmp=unlist(strsplit(selsnspinfo$snphg19[j],"_"))
    tmp1=unlist(strsplit(tmp[1],":"))
    METALdat$chr[j]=tmp1[1]
    METALdat$position[j]=as.numeric(tmp1[2])
  }
  tmp=intersect(METALdat$snp,dat$SNP)
  idx2=match(tmp,dat$SNP)
  idx1=match(tmp,METALdat$snp)
  METALdat$PvalueBonn[idx1]=dat$P[idx2]
  
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
  
  idx=match(selsnps,rownames(bcagenotype))
  ## there a few SNPs not found in bcagenotype
  if (sum(is.na(idx))>0)
    print(paste0(sum(is.na(idx))," out of ", length(idx)," selected snps not been found in genotype data"))
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  if (length(correctedsnps)>0) #flip snps
  {
    idxtocorrect=match(correctedsnps,colnames(Z))
    Z[,idxtocorrect]=2-Z[,idxtocorrect]
  }
  type=table3$Trait[i]
  if (type=="BE")
  {
    idx1=which(covariatetable$phenoBE_bca==2) #case
    idx2=which(covariatetable$phenoBE_bca==1)
  }
  if (type=="BEEA")
  {
    idx1=which(covariatetable$phenoEABE_bca==2) #case
    idx2=which(covariatetable$phenoEABE_bca==1)
  }
  if (type=="EA")
  {
    idx1=which(covariatetable$phenoEA_bca==2) #case
    idx2=which(covariatetable$phenoEA_bca==1)
  }
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  for (j in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[j]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    METALdat$Pvalue[j] <- summary(fit)$coef[2,4]
  }
  idx=which(BE_Oxfordsummarydat$SNP==table3$minsnp[i])
  tmp=paste0(BE_Oxfordsummarydat$chr[idx],":",BE_Oxfordsummarydat$position[idx])
  idx1=which(grepl(tmp,bcabim$V2)) #hg19
  tmp1=paste0(bcabim$V1[idx1],":",bcabim$V4[idx1])
  idx=which(grepl(tmp1,METALdat$MarkerName))
  table3$BCminp[i]=METALdat$Pvalue[idx]
  
  if (table3$Bonn_P[i]<0.05)
  {
    idx=order(METALdat$position)
    write.table(METALdat[idx,],file=paste0("../result/",gene,"_pvalues.txt"),row.names = F,quote=F,sep="\t")
  }
  
}

write.csv(table3,file="../result/TWAS_table3.csv",row.names = F)

#Figure 3, locuszoom (generated using PPT)
# genes=c("EXOC3","SENP6","ZNF641","HSP90AA1","CFDP1","JUND")
# types=c("BC","Bonn","GTEx")
# png("../result/Fig3_1.png",width = 600, height = 800,type="cairo")
# par(mar=c(1,1,1,1))
# for (i in 1:3)
# {
#   gene=genes[i]
#   img1=readPNG(paste0("../result/",gene,"_",types[1],".png"),native=T,info=T)
#   img2=readPNG(paste0("../result/",gene,"_",types[2],".png"),native=T,info=T)
#   plot(NA, xlim = c(0, 800), ylim = c(0, 960), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",bty='n')
#   rasterImage(img1, 0, 0, 800, 480)
#   rasterImage(img2, 0, 480, 800, 960)
#   text(0, 960, 'a', cex=1.5)
#   text(0, 480, 'b', cex=1.5)
#   dev.off()
# }

#check data
#Beacon 2,413 BE cases, 1,512 EA cases and 2,185 controls +4541 controls
#Cambridge  873 BE cases, 995 EA cases (1,868 BE/EA combined), and 3,408 unscreened controls

sampletable=readxl::read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
# 11=England-Sheffield
# 12=Kaiser
# 13=Sweden - Karolinska
# 14=Mayo
# 15=EGA-WA
# 16=Ireland - FINBAR
# 17=Australia - Queensland
# 18=Toronto
# 19=UNC
# 20=WA reflux
# 21=Canada - Nova Scotia
# 22=EGA-NJ
# 23=USC Keck
# 25=Australia-wide
# 27=WA Reid
# 30=Cambridge
# 55=AMOS
sampletable_beacon=sampletable$localid[!sampletable$site %in% c(30,55)]
sampletable_cambridge=sampletable$localid[sampletable$site %in% c(30)]
sampletable_amos=sampletable$localid[sampletable$site %in% c(55)]
idx=match(sampletable_beacon,sampletable$localid)
#sampletable has all cases included in Puya paper
table(sampletable$phenoBE_bca[idx]) #2413 BE
table(sampletable$phenoEA_bca[idx]) #1512 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2185 controls

idx=match(sampletable_cambridge,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #882 BE
table(sampletable$phenoEA_bca[idx]) #1003 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #0 controls

#genotyped samples
BEACON_CRIC1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.fam")
BEACONsamples=gsub("SEP","",BEACON_CRIC1$V2)
#imputed samples
imputedBEACON=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc_maf005_snp/chr1_filter_hg19tohg38_flip.fam")
imputedBEACONsamples=gsub("SEP","",imputedBEACON$V2)

Cambridge_WTCCC1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Cambridge_WTCCC_plinksamples.txt")
Cambridgesamples=gsub("SEP","",Cambridge_WTCCC1$V2)
imputedCambridge=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_hrc_maf005_snp/chr17_filter_hg19tohg38_flip.fam")
imputedCambridgesamples=gsub("SEP","",imputedCambridge$V2)

nrow(BEACON_CRIC1)+nrow(Cambridge_WTCCC1) #15925
nrow(imputedBEACON)+nrow(imputedCambridge) #15091
#merged imputed samples
mergesample=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/","merge_beacon_cambridge_hrc_maf005_snp","/chr",1,"_filter_hg19tohg38_flip.fam"))
mergesample=gsub("SEP","",mergesample$V2)
length(mergesample) #15901

all(BEACONsamples %in% rownames(covariatetable))
#[1] TRUE
all(Cambridgesamples %in% rownames(covariatetable))
#[1] TRUE
#check cases in genotyped data
idx=match(BEACONsamples,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
# 1    2 
# 6717 2406
table(covariatetable$phenoEA_bca[idx])
# 1    2 
# 6718 1508 
table(covariatetable$phenoEABE_bca[idx])
# 1    2 
# 6718 3914 
idx=match(Cambridgesamples,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
# 1    2 
# 3408  882 
table(covariatetable$phenoEA_bca[idx])
# 1    2 
# 3408 1003 
table(covariatetable$phenoEABE_bca[idx])
# 1    2 
# 3408 1885
#check cases in imputed data
idx=match(imputedBEACONsamples,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
# 1    2 
# 6700 2401 
table(covariatetable$phenoEA_bca[idx])
# 1    2 
# 6701 1507 
table(covariatetable$phenoEABE_bca[idx])
# 1    2 
# 6701 3908 
idx=match(imputedCambridgesamples,rownames(covariatetable))
table(covariatetable$phenoBE_bca[idx])
# 1    2 
# 3408  881 
table(covariatetable$phenoEA_bca[idx])
# 1    2 
# 3408 1003 
table(covariatetable$phenoEABE_bca[idx])
# 1    2 
# 3408 1884 

count_imputedsnps=function(impfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc_maf005_snp/")
{
  n=0
  for (i in 1:22)
  {
    tmp=fread(paste0(impfolder,"chr",i,"_filter_hg19tohg38_flip.bim"))
    n=n+nrow(tmp)
  }
  print(n)
  return(n)
}
count_imputedsnps() #5406696
count_imputedsnps(impfolder = "/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_hrc_maf005_snp/") #5384012
count_imputedsnps(impfolder = "/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/") #5312829
