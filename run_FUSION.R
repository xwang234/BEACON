#!/usr/bin/env Rscript

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

load("../data/GTEx/gtexv8_ge_anno.RData")
#load("../data/GTEx/gtex_ge_anno.RData")
genetable=gtexv8_ge_anno[gtexv8_ge_anno$V3=="gene",]
idx=duplicated(genetable$Symbol)
genetable=genetable[!idx,]
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
proteingenestable=gtexv8_ge_anno[gtexv8_ge_anno$Symbol %in% proteingenes & gtexv8_ge_anno$gene_type=="protein_coding",]
proteingenestable=proteingenestable[proteingenestable$V3=="gene",]

#update snpnames
library(data.table)
snptable=as.data.frame(fread(input="../data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"))
snptable$variant_id=gsub("chr","",snptable$variant_id)

bim=bim1=as.data.frame(fread("../data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid.bim"))
table(bim$V2 %in% snptable$variant_id)
idx=match(bim$V2,snptable$variant_id)
snptable=snptable[idx,]
idx=which(grepl("^rs",snptable$rs_id_dbSNP151_GRCh38p7))
bim$V2[idx]=snptable$rs_id_dbSNP151_GRCh38p7[idx]
idx=which(duplicated(bim$V2))
bim$V2[idx]=bim1$V2[idx]
sum(duplicated(bim$V2))
write.table(bim,file="../data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid.bim",
            row.names = F,col.names = F,sep=" ",quote=F)

hm3=data.frame(fread("../data/HM3/hapmap3_r3_b36_fwd.consensus.qc.poly.map"))
all(grepl("^rs",hm3$V2))
tmp=data.frame(snp=hm3$V2)
write.table(tmp,file="../data/HM3/snplist.txt",row.names = F,col.names = F,sep=" ",quote=F)

tmp=read.table("../data/GTEx/tmp/Esophagus_Gastroesophageal_Junction.10/Esophagus_Gastroesophageal_Junction.ENSG00000162618.13.pheno",comment.char = "")
tmp=fread("../data/GTEx/expression_matrices/Esophagus_Gastroesophageal_Junction.v8.EUR.normalized_expression.bed.gz")
tmp=fread("../data/GTEx/expression_covariates/Esophagus_Gastroesophageal_Junction.v8.EUR.covariates1.txt")
tmp=as.data.frame(fread("../result/fusion/WEIGHTS/GTEx.Whole_Blood.pos"))
idx=which(sampletable$site<30 & (sampletable$phenoBE_bca==2 | sampletable$phenoEA_bca==2 )) #3925
idx=which(sampletable$site<30 & (sampletable$phenoBE_bca==1 | sampletable$phenoEA_bca==1 )) #2185 controls
idx=which(sampletable$site<30 & sampletable$phenoBE_bca==2) #2413 BE
table(sampletable$phenoEA_bca )
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt"))
#idx=tmp$non_effect_allele %in% c("A","T","G","C") & tmp$effect_allele %in% c("A","T","G","C")
hm3=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/w_hm3.snplist"))
idx=tmp$SNP %in% hm3$SNP
tmp=tmp[idx,]
tmp$Z=NA
tmp$Z=tmp$beta/tmp$se
tmp$A1=tmp$effect_allele
tmp$A2=tmp$non_effect_allele
write.table(tmp,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N_hm3.txt",
            col.names = T,row.names = F,sep=" ",quote=F)
tmp=read.table("../result/fusion/result/Blood.22.dat",header = T)
tmp0=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt"))
gwasBEEABEACON=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N_hm3.txt"))
gwasBEEABCA=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BEEA_1.tbl"))
tmp=as.data.frame(fread(input="../result/test.sumstats.gz"))

tmp1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BEACON_BEEA_Junction.15.dat",header=T)
tmp2=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/PGC2.SCZ.22.dat",header=T)

Junction_BEACON_BEEA=NULL
for (chr in 1:22)
{
  tmp=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BEACON_BEEA_Junction.",chr,".dat"),header = T,stringsAsFactors = F)
  Junction_BEACON_BEEA=rbind(Junction_BEACON_BEEA,tmp)
}

Junction_BCA_BEEA=NULL
for (chr in 1:22)
{
  tmp=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction.",chr,".dat"),header = T,stringsAsFactors = F)
  Junction_BCA_BEEA=rbind(Junction_BCA_BEEA,tmp)
}
Junction_BCA_BEEA$P=Junction_BCA_BEEA$TWAS.P
Junction_BCA_BEEA=Junction_BCA_BEEA[!is.na(Junction_BCA_BEEA$P),]
#plot results
library(ggplot2) # 3.2.1
library(dplyr)
library(lubridate)
library(cowplot) # 1.0.0
library(egg) # 0.4.5
library(GenomicRanges)

#highrisk SNPs
highrisk=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_26highrisk_SNP_POS.txt",header=T,stringsAsFactors = F)
for (i in 1:nrow(highrisk))
{
  tmp=unlist(strsplit(highrisk$snp[i],":"))
  highrisk$snp[i]=tmp[1]
}
gr_highrisk=GRanges(seqnames = highrisk$chr,ranges = IRanges(start=highrisk$position,width=1))

#use columns CHR, position
gwasplot=function(gwas.dat=gwasBEEABEACON,ylim=NULL,ylabel="GWAS-log10(p)",optsig="gwas")
{
  gwas.dat=gwas.dat[gwas.dat$CHR %in% 1:22,]
  gwas.dat=gwas.dat[!is.na(gwas.dat$P),]
  tmp=paste0(gwas.dat$CHR,"_",gwas.dat$position)
  idx=duplicated(tmp)
  
  if (nrow(gwas.dat)<5000000)
  {
    ndownsize=5
  }else
  {
    ndownsize=50
  }
  if (nrow(gwas.dat)<50000) ndownsize=1
  gwas.dat=gwas.dat[!idx,]
  sig.dat <- gwas.dat %>% 
    subset(P < 0.05)
  #downsampling none-sig SNPs
  notsig.dat <- gwas.dat %>% 
    subset(P >= 0.05) %>%
    dplyr::slice(sample(nrow(.), nrow(.) / ndownsize))
  gwas.dat <- rbind(sig.dat,notsig.dat)
  
  gwas.dat$CHR=as.integer(gwas.dat$CHR)
  idx=order(gwas.dat$CHR,gwas.dat$position)
  gwas.dat=gwas.dat[idx,]
  
  nCHR <- length(unique(gwas.dat$CHR))
  #BPcum is the x coordinate
  gwas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(gwas.dat$CHR)){
    nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$position)
    gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"position"] + s
    s <- s + nbp[i]
  }
  
  highrisk$BPcum=NA
  for (i in 1:nrow(highrisk))
  {
    chr=highrisk$chr[i]
    highrisk$BPcum[i]=highrisk$position[i]+sum(nbp[1:(chr-1)])
  }
  
  axis.set <- gwas.dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(gwas.dat$P)))) + 0.5
  }
  if (optsig=="gwas")
  {
    sig <- 5e-8
  }else
  {
    sig <- 0.05/nrow(gwas.dat) #TWAS
  }
  
  #x lables not show
  idx=axis.set$CHR %in% c(19,21)
  uselabels=axis.set$CHR
  #uselabels[idx]=""
  gwas.dat$logP=-log10(gwas.dat$P)
  gwas.dat$Color=as.factor(gwas.dat$CHR)
  manhplotr <- ggplot(gwas.dat, aes(x = BPcum, y = logP, 
                                    color = Color, size = logP)) +
    geom_point() +
    #geom_rangeframe(data = data.frame(BPcum = c(axis.set$center[1], axis.set$center[22]), logP = c(ylim, 0)))+
    #geom_rangeframe(color="black")+
    #coord_cartesian(xlim=c(axis.set$center[1],max(axis.set$center)),expand = F)+
    #geom_segment(aes_all(c('x', 'y', 'xend', 'yend')),
    #            data = data.frame(x = c(0,axis.set$center[1]), xend = c(0, max(axis.set$center)), y = c(0, 0), yend = c(ylim, 0))) +
    
    geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
    scale_x_continuous(label = uselabels, breaks = axis.set$center) +
    scale_y_continuous(trans = "reverse",limits = c(ylim,0 )) +
    #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
    scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = "Chromosome", 
         y = ylabel) + 
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
      text = element_text(size=16),
      axis.text.x = element_text(angle = 0, vjust = 0.5)
    )
  manhplotr=manhplotr+geom_vline(xintercept = highrisk$BPcum, linetype="dotted", 
                               color = "orange", size=0.5)
  return(manhplotr)
}
twasplot=function(twas.dat=Junction_BEACON_BEEA,ylim=NULL,ylabel="TWAS-log10(p)")
{
  if (sum(colnames(twas.dat)=="P")==0 & sum(colnames(twas.dat)=="TWAS.P")==1) twas.dat$P=twas.dat$TWAS.P
 
  twas.dat=twas.dat[twas.dat$CHR %in% 1:22,]
  twas.dat=twas.dat[!is.na(twas.dat$P),]
  twas.dat$CHR=as.integer(twas.dat$CHR)
  
  if (sum(colnames(twas.dat)=="position")==0 & sum(colnames(twas.dat)=="P0")==1) twas.dat$position=twas.dat$P0
  
  idx=order(twas.dat$CHR,twas.dat$position)
  twas.dat=twas.dat[idx,]
  
  nCHR <- length(unique(twas.dat$CHR))
  #BPcum is the x coordinate
  twas.dat$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(twas.dat$CHR)){
    nbp[i] <- max(twas.dat[twas.dat$CHR == i,]$position)
    twas.dat[twas.dat$CHR == i,"BPcum"] <- twas.dat[twas.dat$CHR == i,"position"] + s
    s <- s + nbp[i]
  }
  highrisk$BPcum=NA
  for (i in 1:nrow(highrisk))
  {
    chr=highrisk$chr[i]
    highrisk$BPcum[i]=highrisk$position[i]+sum(nbp[1:(chr-1)])
  }
  axis.set <- twas.dat %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(twas.dat$P)))) + 0.5
  }
  sig <- 0.05/nrow(twas.dat)
  print(paste0("sig-p=",sig))
  manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                   color = as.factor(CHR), size = -log10(P))) +
    geom_point() +
    geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
    scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
    scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = ylabel) + 
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(4, 5, 1, 5, "mm"),
      axis.line.y = element_line(colour = 'black', size = 1),
      axis.ticks.y = element_line(size = 1, color="black") , 
      axis.ticks.length = unit(.2, "cm"),
      #axis.text.y=element_text(face = "bold"),
      text = element_text(size=16)
    )
  #add gwas highrisk snps

  manhplot=manhplot+geom_vline(xintercept = highrisk$BPcum, linetype="dotted", 
                                 color = "orange", size=0.5)

  return(manhplot)
}

plot2manhattan=function(gwas.dat=gwasBEEABEACON,
                        twas.dat=Junction_BEACON_BEEA,
                        ylabel1="TWAS-log10(p)",
                        ylabel2="GWAS-log10(p)",
                        optsig="gwas",
                        prefix="FUSION_Junction_BEACON_BEEA",ylim=8.5)
{
  
  manhplot=twasplot(twas.dat=twas.dat,ylabel=ylabel1,ylim=ylim)
  manhplotr=gwasplot(gwas.dat=gwas.dat,ylabel=ylabel2,optsig=optsig,ylim=ylim)
  pdf(paste0("../result/Manhattan_",prefix,".pdf"),width = 16,onefile = F)
  egg::ggarrange(manhplot, manhplotr, heights = c(0.5, 0.5))
  dev.off()
}

plot2manhattan(gwas.dat=gwasBEEABCA,
                        twas.dat=Junction_BCA_BEEA[Junction_BCA_BEEA$ID %in% proteingenes,],
                        prefix="FUSION_Junction_BCA_BEEA",ylim=8)


plot2manhattan_files=function(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                        twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction",
                        ylabel1="TWAS-log10(p)",
                        ylabel2="GWAS-log10(p)",
                        optsig="gwas",
                        prefix="FUSION_Junction_BCA_BEEA")
{
  gwas.dat=as.data.frame(fread(gwas.file))
  twas.dat=NULL
  for (chr in 1:22)
  {
    tmp=read.table(paste0(twas.prefix,".",chr,".dat"),header = T,stringsAsFactors = F)
    twas.dat=rbind(twas.dat,tmp)
  }
  if (sum(colnames(twas.dat)=="P")==0 & sum(colnames(twas.dat)=="TWAS.P")==1) twas.dat$P=twas.dat$TWAS.P
  tmp=min(c(gwas.dat$P,twas.dat$P),na.rm=T)
  ylim=-log10(tmp)
  ylim=(as.integer(ylim/0.5)+1)*0.5
  ylim=max(ylim,7.5) # include gwas cutoff
  manhplot=twasplot(twas.dat=twas.dat,ylabel=ylabel1,ylim=ylim)
  manhplotr=gwasplot(gwas.dat=gwas.dat,ylabel=ylabel2,optsig=optsig,ylim=ylim)
  pdf(paste0("../result/Manhattan_",prefix,".pdf"),width = 16,onefile = F)
  egg::ggarrange(manhplot, manhplotr, heights = c(0.5, 0.5))
  dev.off()
}

plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction",
                     prefix="FUSION_Junction_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Junction",
                     prefix="FUSION_Junction_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Junction",
                     prefix="FUSION_Junction_BCA_BE")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Blood",
                     prefix="FUSION_Blood_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Blood",
                     prefix="FUSION_Blood_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Blood",
                     prefix="FUSION_Blood_BCA_BE")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Stomach",
                     prefix="FUSION_Stomach_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Stomach",
                     prefix="FUSION_Stomach_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Stomach",
                     prefix="FUSION_Stomach_BCA_BE")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Mucosa",
                     prefix="FUSION_Mucosa_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Mucosa",
                     prefix="FUSION_Mucosa_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Mucosa",
                     prefix="FUSION_Mucosa_BCA_BE")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Muscularis",
                     prefix="FUSION_Muscularis_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Muscularis",
                     prefix="FUSION_Muscularis_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Muscularis",
                     prefix="FUSION_Muscularis_BCA_BE")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Adiplose",
                     prefix="FUSION_Adipose_BCA_BEEA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Adiplose",
                     prefix="FUSION_Adipose_BCA_EA")
plot2manhattan_files(gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats",
                     twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Adiplose",
                     prefix="FUSION_Adipose_BCA_BE")
#check meta gwas pvalues
EACambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt"))
EABeacon=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_BEACON_autosomes_N.txt"))
#snps in both datasets
EAmeta_both=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_EA_1.tbl"))
#snps in either datasets
EAmeta=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_EAall_1.tbl"))
#check snps having small pvalues
idx1=which(EACambridge$P<1e-3)
idx2=which(EABeacon$P<1e-3)
length(intersect(EACambridge$SNP[idx1],EABeacon$SNP[idx2])) #190
snps=unique(c(EACambridge$SNP[idx1],EABeacon$SNP[idx2]))
comsnps=intersect(EACambridge$SNP[idx1],EABeacon$SNP[idx2])
idx1=match(comsnps,EAmeta_both$MarkerName)
idx2=match(comsnps,EACambridge$SNP)
idx3=match(comsnps,EABeacon$SNP)
plot(-log10(EAmeta_both$`P-value`[idx1]))
points(1:length(idx2),-log10(EACambridge$P[idx2]),col="red")
points(1:length(idx2),-log10(EABeacon$P[idx3]),col="blue")
#they have different effect direction
checksnps=comsnps[which(EAmeta_both$`P-value`[idx1]>0.01)]
idx1=match(checksnps,EACambridge$SNP)
idx2=match(checksnps,EABeacon$SNP)
# View(EACambridge[idx1,])
# View(EABeacon[idx2,])
idx=match(EAmeta_both$MarkerName,EAmeta$MarkerName)
all(EAmeta$`P-value`[idx]==EAmeta_both$`P-value`) #T
uniqsnps=EAmeta$MarkerName[!EAmeta$MarkerName %in% EAmeta_both$MarkerName]
idx=match(uniqsnps,EAmeta$MarkerName)
quantile(EAmeta$`P-value`[idx])
# 0%       25%       50%       75%      100% 
# 4.621e-11 2.320e-01 4.824e-01 7.386e-01 1.000e+00 
quantile(EAmeta_both$`P-value`)
# 0%       25%       50%       75%      100% 
# 5.147e-08 2.398e-01 4.911e-01 7.445e-01 1.000e+00
quantile(EACambridge$P)
# 0%         25%         50%         75%        100% 
# 4.62144e-11 2.42791e-01 4.93685e-01 7.46960e-01 1.00000e+00 
quantile(EABeacon$P[EABeacon$P>0])
# 0%          25%          50%          75%         100% 
# 8.258310e-08 2.416665e-01 4.914950e-01 7.452910e-01 1.000000e+00

#check TWAS muscularis BEEA
gettwas=function(twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Muscularis")
{
  twas.dat=NULL
  for (chr in 1:22)
  {
    tmp=read.table(paste0(twas.prefix,".",chr,".dat"),header = T,stringsAsFactors = F)
    twas.dat=rbind(twas.dat,tmp)
  }
  twas.dat$adjP=p.adjust(twas.dat$TWAS.P,method="bonferroni")
  return(twas.dat)
}
twas_BEEA_muscularis=gettwas()
min(twas_BEEA_muscularis$adjP,na.rm = T) #0.039
twas_BEEA_muscularis[which.min(twas_BEEA_muscularis$adjP),]

#Regular TWAS result
plot2twas=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11",
                   twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction",
                   type="BEEA",
                   prefix="TWAS_Junction_BCA_BEEA",plotopt=T)
{
  twas.dat=NULL
  for (chr in 1:22)
  {
    tmp=read.table(paste0(twas.prefix,".",chr,".dat"),header = T,stringsAsFactors = F)
    twas.dat=rbind(twas.dat,tmp)
  }
  if (sum(colnames(twas.dat)=="P")==0 & sum(colnames(twas.dat)=="TWAS.P")==1) twas.dat$P=twas.dat$TWAS.P
  twas.dat$adjP=p.adjust(twas.dat$TWAS.P,method="bonferroni")
  
  load(paste0(outfolder,"/bca_assoc.RData"))
  regTWAS=assoc_min
  if (type=="BEEA")
  {
    regTWAS$P=regTWAS$BEEA_p
  }
  if (type=="BE")
  {
    regTWAS$P=regTWAS$BE_p
  }
  if (type=="EA")
  {
    regTWAS$P=regTWAS$EA_p
  }
  idx=match(rownames(regTWAS),genetable$Symbol)
  regTWAS$CHR=genetable$Chromosome[idx]
  regTWAS$CHR=gsub("chr","",regTWAS$CHR)
  regTWAS$position=genetable$start[idx]
  
  if (plotopt==T)
  {
    tmp1=min(min(regTWAS$P,na.rm=T),0.05/nrow(regTWAS))
    tmp2=min(min(twas.dat$TWAS.P,na.rm=T),0.05/nrow(twas.dat))
    tmp=min(tmp1,tmp2)
    ylim=(as.integer(-log10(tmp)/0.5)+1)*0.5
    plot2manhattan(gwas.dat=regTWAS,
                   twas.dat=twas.dat,
                   ylabel1 = "FUSION_TWAS-log10(p)",
                   ylabel2 = "Our_TWAS-log10(p)",
                   optsig="twas",
                   prefix=prefix,ylim=ylim)
  }else
  {
    genes=twas.dat$ID[which(twas.dat$adjP<0.2)]
    if (length(genes)>0)
    {
      print(paste0("number of fusion genes:",nrow(twas.dat)))
      for (i in 1:length(genes))
      {
        gene=genes[i]
        idx=which(twas.dat$ID==gene)
        print(twas.dat[idx,])
      }
    }
    return(list(twas.dat=twas.dat,regTWAS=regTWAS))
  }
  
}

tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction",
              type="BEEA",
              prefix="TWAS_Junction_BCA_BEEA",plotopt = F)
# [1] "number of fusion genes:4786"
# PANEL
# 1206 Esophagus_Gastroesophageal_Junction
# FILE
# 1206 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS//Esophagus_Gastroesophageal_Junction/Esophagus_Gastroesophageal_Junction.ENSG00000237125.4.wgt.RDat
# ID CHR       P0       P1    HSQ BEST.GWAS.ID BEST.GWAS.Z    EQTL.ID EQTL.R2 EQTL.Z EQTL.GWAS.Z NSNP NWGT MODEL MODELCV.R2
# 1206 HAND2-AS1   4 1.74e+08 1.75e+08 0.2025    rs4585277        3.48 rs17059549   0.109   5.47      -2.592  345    5 lasso     0.1785
# MODELCV.PV   TWAS.Z   TWAS.P        P      adjP
# 1206   8.28e-11 -4.25878 2.06e-05 2.06e-05 0.0981796
# PANEL
# 3585 Esophagus_Gastroesophageal_Junction
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Junction",
              type="BE",
              prefix="TWAS_Junction_BCA_BE",plotopt = F)

tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Stomach",
              type="BEEA",
              prefix="TWAS_Stomach_BCA_BEEA",plotopt = F)
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Mucosa",
              type="BEEA",
              prefix="TWAS_Mucosa_BCA_BEEA",plotopt = F)
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Mucosa",
              type="BE",
              prefix="TWAS_Mucosa_BCA_BE",plotopt = F)
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11",
          twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Muscularis",
          type="BEEA",
          prefix="TWAS_Muscularis_BCA_BEEA",plotopt = F)
#"number of fusion genes:7647"
# FILE
# 910 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS//Esophagus_Muscularis/Esophagus_Muscularis.ENSG00000118961.10.wgt.RDat
# ID CHR       P0      P1    HSQ BEST.GWAS.ID BEST.GWAS.Z    EQTL.ID EQTL.R2 EQTL.Z EQTL.GWAS.Z NSNP
# 910 C2orf43   2 20900000 2.1e+07 0.5015       rs7255       -5.41 rs13385191   0.305 -10.23      -3.731  512
# NWGT MODEL MODELCV.R2 MODELCV.PV  TWAS.Z   TWAS.P        P       adjP
# 910   32  enet   0.336881    1.1e-31 4.55818 5.16e-06 3.36e-06 0.02559312
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Adiplose",
              type="BEEA",
              prefix="TWAS_Adipose_BCA_BEEA",plotopt = F)
tmp=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Adiplose",
              type="BE",
              prefix="TWAS_Adipose_BCA_BE",plotopt = F)
tmp1=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11",
              twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Mucosa",
              type="BE",
              prefix="TWAS_Mucosa_BCA_BE",plotopt = F)
tmp2=plot2twas(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11",
               twas.prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Adiplose",
               type="BEEA",
               prefix="TWAS_Adipose_BCA_BEEA",plotopt = F)

#check gwas
gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats"
gwas.dat=as.data.frame(fread(gwas.file))
which(gwas.dat$P<5e-8)
gwas.dat[which(gwas.dat$P<5e-8),]
# SNP A1 A2 Weight     Z         P Direction CHR position
# 6883406 rs8031215  a  g  15908 5.671 1.417e-08        ++  15 58362207
gwas.beacon=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt"))
gwas.beacon[which(gwas.beacon$SNP=="rs181436442"),] #P=0.0532759
gwas.cambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt"))
gwas.cambridge[which(gwas.cambridge$SNP=="rs181436442"),] #P=2.5925e-12
gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats"
gwas.dat=as.data.frame(fread(gwas.file))
which(gwas.dat$P<5e-8)
gwas.dat[which(gwas.dat$P<5e-8),]
gwas.file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats"
gwas.dat=as.data.frame(fread(gwas.file))
sum(gwas.dat$P<5e-8)
gwas.dat[which(gwas.dat$P<5e-8),]
