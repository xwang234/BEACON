#!/usr/bin/env Rscript

library(data.table)
#sample table
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
proteingenestable=gtexv8_ge_anno[gtexv8_ge_anno$Symbol %in% proteingenes & gtexv8_ge_anno$gene_type=="protein_coding",]
proteingenestable=proteingenestable[proteingenestable$V3=="gene",]

#hg38->hg19
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
gr_proteingenetable=GRanges(seqnames = proteingenestable$Chromosome,ranges = IRanges(start=proteingenestable$start,end=proteingenestable$end))
proteingenestable$starthg19=proteingenestable$endhg19=proteingenestable$Chromosomehg19=NA
for (i in 1:nrow(proteingenestable))
{
  if (i %% 1000==0) cat(i,'..')
  tmp=liftOver(gr_proteingenetable[i],chain)
  if (length(unlist(start(tmp)))>0)
  {
    proteingenestable$starthg19[i]=unlist(start(tmp))[1]
    proteingenestable$endhg19[i]=unlist(end(tmp))[length(unlist(end(tmp)))]
    proteingenestable$Chromosomehg19[i]=as.character(seqnames(tmp)@unlistData@values)[1]
  }
}
save(proteingenestable,file="../result/proteingenestable.RData")
table(is.na(proteingenestable$starthg19))
proteingenestable=proteingenestable[!is.na(proteingenestable$starthg19),]
idx=proteingenestable$Chromosome==proteingenestable$Chromosomehg19
quantile(proteingenestable$start[idx]- proteingenestable$starthg19[idx])
quantile(proteingenestable$end[idx]- proteingenestable$endhg19[idx])
table(proteingenestable$Chromosome==proteingenestable$Chromosomehg19)

#GWAS result:
#gwasEA=as.data.frame(fread("../result/GWAS/Beacon_genotyped_EA_CO_gwas.assoc.logistic"))
gwasBEBA=as.data.frame(fread("../result/GWAS/BEEA_CO_19Feb2015.assoc.logistic"))
gwasEA=as.data.frame(fread("../result/GWAS/EA_CO_19Feb2015.assoc.logistic"))
gwasBE=as.data.frame(fread("../result/GWAS/BE_CO_19Feb2015.assoc.logistic"))
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"

#TWAS result:
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]



library(ggplot2) # 3.2.1
library(dplyr)
library(lubridate)
library(cowplot) # 1.0.0
library(egg) # 0.4.5

#' #' Create some data to play with. Two time series with the same timestamp.
#' df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), 
#'                  series1 = rnorm(8761), 
#'                  series2 = rnorm(8761, 100))
#' 
#' #' Create the two plots.
#' plot1 <- df %>%
#'   dplyr::select(DateTime, series1) %>%
#'   na.omit() %>%
#'   ggplot() +
#'   geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
#'   ylab("Red dots / m") +
#'   theme_minimal() +
#'   theme(axis.title.x = element_blank(),
#'         axis.text.x = element_blank())
#' 
#' plot2 <- df %>%
#'   dplyr::select(DateTime, series2) %>%
#'   na.omit() %>%
#'   ggplot() +
#'   geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
#'   ylab("Blue drops / L") +scale_y_continuous(trans = "reverse")+
#'   theme_minimal() +
#'   theme(axis.title.x = element_blank())
#' 
#' # Draw the two plot aligned vertically, with the top plot 1/3 of the height
#' # of the bottom plot
#' cowplot::plot_grid(plot1, plot2, align = "v", ncol = 1, rel_heights = c(0.25, 0.75))
#' egg::ggarrange(plot1, plot2, heights = c(0.5, 0.5))


# library(normentR)
# set.seed(2404)
# gwas.dat <- simulateGWAS(nSNPs = 1e5, nSigCols = 3)
# gwas.dat=gwas.dat[gwas.dat$CHR %in% 1:22,]
# 
# sig.dat <- gwas.dat %>% 
#   subset(P < 0.05)
# notsig.dat <- gwas.dat %>% 
#   subset(P >= 0.05) %>%
#   slice(sample(nrow(.), nrow(.) / 5))
# gwas.dat <- rbind(sig.dat,notsig.dat)
# 
# 
# tmp=paste0(gwas.dat$CHR,"_",gwas.dat$BP)
# idx=duplicated(tmp)
# gwas.dat=gwas.dat[!idx,]
# gwas.dat$CHR=as.integer(gwas.dat$CHR)
# idx=order(gwas.dat$CHR,gwas.dat$BP)
# gwas.dat=gwas.dat[idx,]
# 
# nCHR <- length(unique(gwas.dat$CHR))
# gwas.dat$BPcum <- NA
# s <- 0
# nbp <- c()
# for (i in unique(gwas.dat$CHR)){
#   nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
#   gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
#   s <- s + nbp[i]
# }
# 
# axis.set <- gwas.dat %>% 
#   group_by(CHR) %>% 
#   summarize(center = (max(BPcum) + min(BPcum)) / 2)
# ylim <- abs(floor(log10(min(gwas.dat$P)))) + 2 
# sig <- 5e-8
# 
# manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(P), 
#                                  color = as.factor(CHR), size = -log10(P))) +
#   geom_point(alpha = 0.75) +
#   geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
#   scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
#   scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "-log10(p)") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     plot.margin = margin(4, 5, 1, 5, "mm"),
#     axis.line.y = element_line(colour = 'black', size = 1),
#     axis.ticks.y = element_line(size = 1, color="black") , 
#     axis.ticks.length = unit(.2, "cm"),
#   )
# print(manhplot)
# 
# manhplotr <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(P), 
#                                  color = as.factor(CHR), size = -log10(P))) +
#   geom_point(alpha = 0.75) +
#   geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
#   scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
#   scale_y_continuous(trans = "reverse",limits = c(ylim,0 )) +
#   scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
#   scale_size_continuous(range = c(0.5,3)) +
#   labs(x = NULL, 
#        y = "-log10(p)") + 
#   theme_minimal() +
#   theme( 
#     legend.position = "none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     plot.margin = margin(0, 5, 1, 5, "mm"),
#     axis.line = element_line(colour = 'black', size = 1),
#     axis.ticks = element_line(size = 1, color="black") , 
#     axis.ticks.length = unit(.2, "cm"),
#     axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
#   )
#   
#   #annotate(x=100, xend=100, y=ylim, yend=0, colour="black", lwd=0.75, geom="segment") #+
#   #annotate(x=0, xend=max(gwas.dat$BPcum), y=0, yend=0, colour="black", lwd=0.75, geom="segment")
# print(manhplotr)
# egg::ggarrange(manhplot, manhplotr, heights = c(0.5, 0.5))

gwasplot=function(gwas.dat=gwasBEBA,ylim=NULL)
{
  gwas.dat=gwas.dat[gwas.dat$CHR %in% 1:22,]
  gwas.dat=gwas.dat[!is.na(gwas.dat$P),]
  tmp=paste0(gwas.dat$CHR,"_",gwas.dat$BP)
  idx=duplicated(tmp)
  
  gwas.dat=gwas.dat[!idx,]
  sig.dat <- gwas.dat %>% 
    subset(P < 0.05)
  #downsampling none-sig SNPs
  notsig.dat <- gwas.dat %>% 
    subset(P >= 0.05) %>%
    dplyr::slice(sample(nrow(.), nrow(.) / 5))
  gwas.dat <- rbind(sig.dat,notsig.dat)
  
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
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(gwas.dat$P)))) + 0.5
  }
  sig <- 5e-8
  
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
         y = "GWAS-log10(p)") + 
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
  return(manhplotr)
}
twasplot=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11",
                  opt="BEEA",ylim=NULL)
{
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
  skatres=skat_min2[rownames(skat_min2) %in% proteingenestable$Symbol,]
  idx=match(rownames(skatres),proteingenestable$Symbol)
  twas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=rownames(skatres),P=NA,BP=proteingenestable$starthg19[idx],stringsAsFactors = F)
  twas.dat$CHR=gsub("chr","",twas.dat$CHR)
  if (opt=="BEEA")
  {
    twas.dat$P=skatres$BEEA_p
  }
  if (opt=="EA")
  {
    twas.dat$P=skatres$EA_p
  }
  if (opt=="BE")
  {
    twas.dat$P=skatres$BE_p
  }
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
         y = "TWAS-log10(p)") + 
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
  return(manhplot)
}

plot2manhattan=function(gwas.dat=gwasBEBA,
                        outfolder="dist500K_GTEx_June11",
                        opt="BEEA",prefix="Junction",ylim=10.5)
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",outfolder)
  manhplot=twasplot(outfolder=outfolder,opt=opt,ylim=ylim)
  manhplotr=gwasplot(gwas.dat=gwas.dat,ylim=ylim)
  #png("../result/test.png",width=900,res=100,pointsize = 2)
  pdf(paste0("../result/Manhattan_",prefix,"_",opt,".pdf"),width = 16,onefile = F)
  # library(grDevices)
  #cairo_ps(filename = "../result/test.eps",onefile = F,fallback_resolution = 600)
  egg::ggarrange(manhplot, manhplotr, heights = c(0.5, 0.5))
  dev.off()
}

#BEEA
plot2manhattan()
plot2manhattan(outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach")
plot2manhattan(outfolder="dist500K_GTEx_blood_June11",prefix = "Blood")
plot2manhattan(outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa")
plot2manhattan(outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis")
plot2manhattan(outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose")

#EA
plot2manhattan(gwas.dat=gwasEA,opt="EA",ylim=7.5)
plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach",ylim=7.5)
plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_blood_June11",prefix = "Blood",ylim=7.5)
plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa",ylim=7.5)
plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis",ylim=7.5)
plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose",ylim=7.5)

#BE
plot2manhattan(gwas.dat=gwasBE,opt="BE",ylim=9)
plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach",ylim=9)
plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_blood_June11",prefix = "Blood",ylim=9)
plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa",ylim=9)
plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis",ylim=9)
plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose",ylim=9)
