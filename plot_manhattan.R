#!/usr/bin/env Rscript
#used to draw manhattan plot of TWAS result

library(data.table)
#sample table
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
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
#save(proteingenestable,file="../result/proteingenestable.RData")
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
colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_pc6_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]



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

# gwasplot=function(gwas.dat=gwasBEBA,ylim=NULL)
# {
#   gwas.dat=gwas.dat[gwas.dat$CHR %in% 1:22,]
#   gwas.dat=gwas.dat[!is.na(gwas.dat$P),]
#   tmp=paste0(gwas.dat$CHR,"_",gwas.dat$BP)
#   idx=duplicated(tmp)
#   
#   gwas.dat=gwas.dat[!idx,]
#   sig.dat <- gwas.dat %>% 
#     subset(P < 0.05)
#   #downsampling none-sig SNPs
#   notsig.dat <- gwas.dat %>% 
#     subset(P >= 0.05) %>%
#     dplyr::slice(sample(nrow(.), nrow(.) / 5))
#   gwas.dat <- rbind(sig.dat,notsig.dat)
#   
#   gwas.dat$CHR=as.integer(gwas.dat$CHR)
#   idx=order(gwas.dat$CHR,gwas.dat$BP)
#   gwas.dat=gwas.dat[idx,]
#   
#   nCHR <- length(unique(gwas.dat$CHR))
#   #BPcum is the x coordinate
#   gwas.dat$BPcum <- NA
#   s <- 0
#   nbp <- c()
#   for (i in unique(gwas.dat$CHR)){
#     nbp[i] <- max(gwas.dat[gwas.dat$CHR == i,]$BP)
#     gwas.dat[gwas.dat$CHR == i,"BPcum"] <- gwas.dat[gwas.dat$CHR == i,"BP"] + s
#     s <- s + nbp[i]
#   }
#   
#   axis.set <- gwas.dat %>% 
#     group_by(CHR) %>% 
#     dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
#   if (is.null(ylim))
#   {
#     ylim <- abs(floor(log10(min(gwas.dat$P)))) + 0.5
#   }
#   sig <- 5e-8
#   
#   #x lables not show
#   idx=axis.set$CHR %in% c(19,21)
#   uselabels=axis.set$CHR
#   #uselabels[idx]=""
#   gwas.dat$logP=-log10(gwas.dat$P)
#   gwas.dat$Color=as.factor(gwas.dat$CHR)
#   manhplotr <- ggplot(gwas.dat, aes(x = BPcum, y = logP, 
#                                     color = Color, size = logP)) +
#     geom_point() +
#     #geom_rangeframe(data = data.frame(BPcum = c(axis.set$center[1], axis.set$center[22]), logP = c(ylim, 0)))+
#     #geom_rangeframe(color="black")+
#     #coord_cartesian(xlim=c(axis.set$center[1],max(axis.set$center)),expand = F)+
#     #geom_segment(aes_all(c('x', 'y', 'xend', 'yend')),
#      #            data = data.frame(x = c(0,axis.set$center[1]), xend = c(0, max(axis.set$center)), y = c(0, 0), yend = c(ylim, 0))) +
#     
#     geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
#     scale_x_continuous(label = uselabels, breaks = axis.set$center) +
#     scale_y_continuous(trans = "reverse",limits = c(ylim,0 )) +
#     #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
#     scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
#     scale_size_continuous(range = c(0.5,3)) +
#     labs(x = "Chromosome", 
#          y = "GWAS-log10(p)") + 
#     theme_minimal() +
#     theme( 
#       legend.position = "none",
#       panel.border = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       # panel.grid.major.y = element_blank(),
#       # panel.grid.minor.y = element_blank(),
#       plot.margin = margin(0, 5, 10, 5, "mm"),
#       axis.line = element_line(colour = 'black', size = 1),
#       axis.ticks = element_line(size = 1, color="black") , 
#       axis.ticks.length = unit(.2, "cm"),
#       text = element_text(size=16),
#       axis.text.x = element_text(angle = 0, vjust = 0.5)
#     )
#   return(manhplotr)
# }

twasplot=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar",
                  opt="BEEA",ylim=NULL,opt1="nogenes")
{
  load(paste0(outfolder,"/skat_res.RData"))
  colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") 
  skatres=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenestable$Symbol,]
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
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
    dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
  if (is.null(ylim))
  {
    ylim <- abs(floor(log10(min(twas.dat$P)))) + 0.5
  }
  sig <- 0.05/nrow(twas.dat)
  print(paste0("sig-p=",sig))
  twas.dat$fdr=p.adjust(twas.dat$P,method="fdr")
  twas.dat$fewr=p.adjust(twas.dat$P,method="bonferroni")
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
  fweridx=which(twas.dat$fdr<0.05)
  print(twas.dat$Gene[fweridx])
  fweridx=fweridx[order(twas.dat$P[fweridx],decreasing = T)]
  
  fwergenes=NULL
  if (length(fweridx)>0)
  {
    nlenx=rep(NA,length(fweridx))
    fwergenes=twas.dat$Gene[fweridx]
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
  
  xmax=max(twas.dat$BPcum)*1.01
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="genes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR), size = -log10(P))) +
      geom_point() +
      geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "black", linetype = "dashed",size=1) + 
      annotate("text", x=100000, y=-log10(fdrcutoff)+0.2, label= paste0("FDR=",fdrcutoff1),size=6) +
      scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,limits=c(0,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "-log10(p)") + 
      annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=6) +
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
  }
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="nogenes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR), size = -log10(P))) +
      geom_point() +
      geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "black", linetype = "dashed",size=1) + 
      annotate("text", x=500000000, y=-log10(fdrcutoff)+0.15, label= paste0("FDR=",fdrcutoff1),size=6) +
      scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,limits=c(0,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "-log10(p)") + 
      #annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=6) +
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
  }
  return(manhplot)
}

#another way to annotate genes, use comma for close genes, doesn't work
# twasplot=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar",
#                   opt="BEEA",ylim=NULL)
# {
#   load(paste0(outfolder,"/skat_res.RData"))
#   colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") 
#   skatres=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenestable$Symbol,]
#   load(paste0(outfolder,"/preidiction_michigan_model.RData"))
#   skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
#   idx=match(rownames(skatres),proteingenestable$Symbol)
#   twas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=rownames(skatres),P=NA,BP=proteingenestable$starthg19[idx],stringsAsFactors = F)
#   twas.dat$CHR=gsub("chr","",twas.dat$CHR)
#   if (opt=="BEEA")
#   {
#     twas.dat$P=skatres$BEEA_p
#   }
#   if (opt=="EA")
#   {
#     twas.dat$P=skatres$EA_p
#   }
#   if (opt=="BE")
#   {
#     twas.dat$P=skatres$BE_p
#   }
#   twas.dat=twas.dat[twas.dat$CHR %in% 1:22,]
#   twas.dat=twas.dat[!is.na(twas.dat$P),]
#   twas.dat$CHR=as.integer(twas.dat$CHR)
#   idx=order(twas.dat$CHR,twas.dat$BP)
#   twas.dat=twas.dat[idx,]
#   
#   nCHR <- length(unique(twas.dat$CHR))
#   #BPcum is the x coordinate
#   twas.dat$BPcum <- NA
#   s <- 0
#   nbp <- c()
#   for (i in unique(twas.dat$CHR)){
#     nbp[i] <- max(twas.dat[twas.dat$CHR == i,]$BP)
#     twas.dat[twas.dat$CHR == i,"BPcum"] <- twas.dat[twas.dat$CHR == i,"BP"] + s
#     s <- s + nbp[i]
#   }
#   
#   axis.set <- twas.dat %>% 
#     group_by(CHR) %>% 
#     dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
#   if (is.null(ylim))
#   {
#     ylim <- abs(floor(log10(min(twas.dat$P)))) + 0.5
#   }
#   sig <- 0.05/nrow(twas.dat)
#   print(paste0("sig-p=",sig))
#   twas.dat$fdr=p.adjust(twas.dat$P,method="fdr")
#   twas.dat$fewr=p.adjust(twas.dat$P,method="bonferroni")
#   fdrcutoff=NULL
#   tmp=twas.dat$fdr-0.05
#   
#   idx1=which(tmp>0)
#   fdrcutoff2=min(twas.dat$P[idx1])
#   idx1=which(tmp<0)
#   if (length(idx1)>0)
#   {
#     fdrcutoff = (min(twas.dat$P[idx1])+fdrcutoff2)/2
#   }
#   manhplot=NULL
#   #fweridx=which(twas.dat$P<=sig)
#   #annote FDR genes
#   fweridx=which(twas.dat$fdr<0.05)
#   fweridx=fweridx[order(twas.dat$CHR[fweridx],twas.dat$BPcum[fweridx])]
#   fwergenes=NULL
#   annotatedat=twas.dat[fweridx,]
#   
#   diff1=c(diff(annotatedat$BPcum))
#   diff2=diff1<1e6
#   mysegment=data.frame(start=1,end=1)
#   for (i in 1:length(diff2))
#   {
#     if (diff2[i])
#     {
#       mysegment$end[nrow(mysegment)]=mysegment$end[nrow(mysegment)]+1
#     
#     }else
#     {
#       tmp=data.frame(start=mysegment$end[nrow(mysegment)]+1,end=mysegment$end[nrow(mysegment)]+1)
#       mysegment=rbind(mysegment,tmp)
#     }
#   }
#   annotatedat1=NULL
#   for (i in 1:nrow(mysegment))
#   {
#     tmp=data.frame(gene=paste0(annotatedat$Gene[mysegment$start[i]:mysegment$end[i]],collapse = ","),
#                    P=min(annotatedat$P[mysegment$start[i]:mysegment$end[i]]),
#                    BPcum=mean(annotatedat$BPcum[mysegment$start[i]:mysegment$end[i]]))
#     annotatedat1=rbind(annotatedat1,tmp)
#   }
#   
#   # if (length(fweridx)>0)
#   # {
#   #   nlenx=rep(NA,length(fweridx))
#   #   fwergenes=twas.dat$Gene[fweridx]
#   #   print(fwergenes)
#   #   nlen=sapply(fwergenes,nchar) #length of gene names
#   #   nlen[nlen>7]=7
#   #   fwerx=twas.dat$BPcum[fweridx]
#   #   fwery=rep(0.2,length(fweridx))
#   #   idx1=seq(1,length(fwerx),2)
#   #   idx2=NULL
#   #   if (length(fwerx)>1) idx2=seq(2,length(fwerx),2)
#   #   nlenx[idx1]=20000000*nlen[idx1]
#   #   fwerx[idx1]=fwerx[idx1]-nlenx[idx1]
#   #   fwery[idx1]=fwery[idx1]-0.05
#   #   
#   #   if (length(idx2)>0)
#   #   {
#   #     nlenx[idx2]=15000000*nlen[idx2]
#   #     fwerx[idx2]=fwerx[idx2]+nlenx[idx2]
#   #     fwery[idx2]=fwery[idx2]+0.05
#   #   }
#   # }
#   
#   xmax=max(twas.dat$BPcum)*1.01
#   if (!is.null(fdrcutoff) & length(fweridx)>0)
#   {
#     manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
#                                      color = as.factor(CHR), size = -log10(P))) +
#       geom_point() +
#       geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
#       geom_hline(yintercept = -log10(fdrcutoff), color = "black", linetype = "dashed",size=1) + 
#       scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,limits=c(0,xmax)) +
#       scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
#       #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
#       scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
#       scale_size_continuous(range = c(0.5,3)) +
#       labs(x = "Chromosome", 
#            y = "-log10(p)") + 
#       annotate("text", x=annotatedat1$BPcum, y=-log10(annotatedat1$P)+0.3, label= annotatedat1$gene,size=6) +
#       theme_minimal() +
#       theme( 
#         legend.position = "none",
#         panel.border = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.grid.major.y = element_blank(),
#         # panel.grid.minor.y = element_blank(),
#         plot.margin = margin(0, 5, 10, 5, "mm"),
#         axis.line = element_line(colour = 'black', size = 1),
#         axis.ticks = element_line(size = 1, color="black") , 
#         axis.ticks.length = unit(.2, "cm"),
#         text = element_text(size=16),
#         axis.text.x = element_text(angle = 0, vjust = 0.5)
#       )
#   }
#   
#   
#   return(manhplot)
# }

plotmanhattan=function(organ="junction",opt="BEEA")
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  pdf(paste0("../result/Manhattan2_",organ,"_",opt,".pdf"),width = 16,onefile = F)
  manhplot=twasplot(outfolder=outfolder,opt=opt,opt1="nogenes")
  print(manhplot)
  dev.off()
  pdf(paste0("../result/Manhattan2_",organ,"_",opt,"_genes.pdf"),width = 16,onefile = F)
  manhplot=twasplot(outfolder=outfolder,opt=opt,opt1="genes")
  print(manhplot)
  dev.off()
}


# plotmanhattan(organ="mucosa")
# plotmanhattan(organ="junction",opt="BE")
# plotmanhattan(organ="junction",opt="BEEA")
# plotmanhattan(organ="stomach",opt="BEEA")
# plotmanhattan(organ="muscularis",opt="EA")
# plotmanhattan(organ="muscularis",opt="BEEA")
# plotmanhattan(organ="adipose",opt="BE")
# plotmanhattan(organ="adipose",opt="BEEA")
# plotmanhattan(organ="blood",opt="BE")
# plotmanhattan(organ="blood",opt="EA")
# plotmanhattan(organ="blood",opt="BEEA")
plotmanhattan(organ="mucosa",opt="BE")
plotmanhattan(organ="mucosa",opt="EA")
plotmanhattan(organ="mucosa",opt="BEEA")
plotmanhattan(organ="junction",opt="BE")
plotmanhattan(organ="junction",opt="BEEA")
plotmanhattan(organ="stomach",opt="EA")
plotmanhattan(organ="stomach",opt="BEEA")

plotmanhattan(organ="adipose",opt="BE")
plotmanhattan(organ="adipose",opt="EA")
plotmanhattan(organ="adipose",opt="BEEA")
plotmanhattan(organ="blood",opt="BE")
plotmanhattan(organ="blood",opt="EA")
plotmanhattan(organ="blood",opt="BEEA")

#plot manhattan plot on combined p-values----
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

allres$FWER=NA
allres$FDR=NA
allres$FDR=p.adjust(allres$p,method="fdr")
allres$FWER=p.adjust(allres$p,method="bonferroni")
twasplot_all=function(skatres=allres,ylim=NULL,opt1="nogenes")
{
  idx=match(skatres$gene,proteingenestable$Symbol)
  skatres$type=gsub("vsCO","",skatres$type)
  twas.dat=data.frame(CHR=proteingenestable$Chromosomehg19[idx],Gene=skatres$gene,organ=skatres$tissue,type=skatres$type,P=NA,BP=proteingenestable$starthg19[idx],stringsAsFactors = F)
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
  twas.dat$fdr=p.adjust(twas.dat$P,method="fdr")
  twas.dat$fewr=p.adjust(twas.dat$P,method="bonferroni")
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
  
  xmax=max(twas.dat$BPcum)*1.01
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="genes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR), size = -log10(P))) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "black", linetype = "dashed",size=1) + 
      annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,limits=c(0,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "-log10(p)") + 
      annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=2) +
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
  }
  if (!is.null(fdrcutoff) & length(fweridx)>0 & opt1=="nogenes")
  {
    manhplot <- ggplot(twas.dat, aes(x = BPcum, y = -log10(P), 
                                     color = as.factor(CHR), size = -log10(P))) +
      geom_point() +
      #geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed",size=1) + 
      geom_hline(yintercept = -log10(fdrcutoff), color = "black", linetype = "dashed",size=1) + 
      annotate("text", x=50000, y=-log10(fdrcutoff)+0.3, label= paste0("FDR=",fdrcutoff1),size=6) +
      scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,limits=c(0,xmax)) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_color_manual(values = rep(c("navyblue", "springgreen3"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = "Chromosome", 
           y = "-log10(p)") + 
      #annotate("text", x=fwerx, y=-log10(twas.dat$P[fweridx])+fwery, label= fwergenes,size=6) +
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
  }
  return(manhplot)
}

pdf(paste0("../result/Manhattan_AllFDR.pdf"),width = 16,onefile = F)
manhplot=twasplot_all(opt1="nogenes")
print(manhplot)
dev.off()
pdf("../result/Manhattan_AllFDR_genes.pdf",width = 16,onefile = F)
manhplot=twasplot_all(opt1="genes")
print(manhplot)
dev.off()

# plot2manhattan=function(gwas.dat=gwasBEBA,
#                         outfolder="dist500K_GTEx_June11",
#                         opt="BEEA",prefix="Junction",ylim=10.5)
# {
#   outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",outfolder)
#   manhplot=twasplot(outfolder=outfolder,opt=opt,ylim=ylim)
#   manhplotr=gwasplot(gwas.dat=gwas.dat,ylim=ylim)
#   #png("../result/test.png",width=900,res=100,pointsize = 2)
#   pdf(paste0("../result/Manhattan_",prefix,"_",opt,".pdf"),width = 16,onefile = F)
#   # library(grDevices)
#   #cairo_ps(filename = "../result/test.eps",onefile = F,fallback_resolution = 600)
#   egg::ggarrange(manhplot, manhplotr, heights = c(0.5, 0.5))
#   dev.off()
# }
# 
# #BEEA
# plot2manhattan()
# plot2manhattan(outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach")
# plot2manhattan(outfolder="dist500K_GTEx_blood_June11",prefix = "Blood")
# plot2manhattan(outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa")
# plot2manhattan(outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis")
# plot2manhattan(outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose")
# 
# #EA
# plot2manhattan(gwas.dat=gwasEA,opt="EA",ylim=7.5)
# plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach",ylim=7.5)
# plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_blood_June11",prefix = "Blood",ylim=7.5)
# plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa",ylim=7.5)
# plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis",ylim=7.5)
# plot2manhattan(gwas.dat=gwasEA,opt="EA",outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose",ylim=7.5)
# 
# #BE
# plot2manhattan(gwas.dat=gwasBE,opt="BE",ylim=9)
# plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_stomach_June11",prefix = "Stomach",ylim=9)
# plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_blood_June11",prefix = "Blood",ylim=9)
# plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_mucosa_June11",prefix = "Mucosa",ylim=9)
# plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_muscularis_June11",prefix = "Muscularis",ylim=9)
# plot2manhattan(gwas.dat=gwasBE,opt="BE",outfolder="dist500K_GTEx_adipose_June11",prefix = "Adipose",ylim=9)
