hg38tohg19=function(snpnames=rownames(snp)[1:10])
{
  library(rtracklayer)
  library(GenomicRanges)
  chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
  tmp=unlist(strsplit(snpnames,":"))
  chr0=tmp[1]
  chr=paste0("chr",tmp[1])
  tmp=tmp[seq(2,length(tmp),2)]
  tmp=unlist(strsplit(tmp,"_"))
  pos=as.integer(tmp[seq(1,length(tmp),3)])
  alt=tmp[seq(2,length(tmp),3)]
  ref=tmp[seq(3,length(tmp),3)]
  gr_dat=GRanges(seqnames = chr,ranges=IRanges(start=pos,width=1))
  tmp=liftOver(gr_dat,chain)
  newsnpnames=newpos=rep(NA,length(tmp))
  for (i in 1:length(tmp))
  {
    tmp1=unlist(tmp[i])
    if (length(tmp1)==0)
    {
      warning(paste0(snpnames[i]," not transformed"))
    }else
    {
      if (length(tmp1)==1)
      {
        newpos[i]=start(tmp1)
      }else
      {
        warning(paste0(snpnames[i]," transformed to ",length(tmp1)," snps"))
        newpos[i]=start(tmp1)[1]
      }
    }
  }
  newsnpnames=paste0(chr0,":",newpos,"_",alt,"_",ref)
  res=data.frame(snphg38=snpnames,snphg19=newsnpnames,stringsAsFactors = F)
  return(res)
}

refluxtable=read.table(file="../data/Reflux_ST2.txt",fill=T,sep=" ")
idx=which(grepl("rs",refluxtable[,1]))
snps=refluxtable[idx,1]
chrs=refluxtable[idx,2]
poss=refluxtable[idx+1,1]
flxgenes=rep(NA,length(idx))
allflxgenes=NULL
for (i in idx)
{
  tmp=refluxtable[i+2,1]
  tmp1=unlist(strsplit(tmp,"|",fixed=T))
  if (length(tmp1)>1)
  {
    flxgenes[i]=tmp1[2]
    allflxgenes=c(allflxgenes,unlist(strsplit(tmp1[2],",")))
    allflxgenes=unique(allflxgenes)
  }
}

library("readxl")
table1=read_excel("../result/TWAS_tables.xlsx",1)
table1=table1[-1,]
refluxstr=paste0(chrs,":",poss)


check_reflux=function(thetable=table1)
{
  thetable$numreflxsnp=NA
  thetable$flxsnp=NA
  thetable$numflxgene=NA
  thetable$flxgene=F
  for (i in 1:nrow(thetable))
  {
    organ=tolower(thetable$Tissue[i])
    outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    gene=thetable$Gene[i]
    idx1=which(rownames(res_min)==gene)
    selsnps=unlist(strsplit(res_min$selectedsnps[idx1],"|",fixed=T))
    selsnps_hg19=hg38tohg19(snpnames = selsnps)
    selsnps_hg19str=unlist(strsplit(selsnps_hg19$snphg19,"_"))
    selsnps_hg19str=selsnps_hg19str[seq(1,length(selsnps_hg19str),3)]
    
    thetable$numreflxsnp[i]=length(intersect(selsnps_hg19str,refluxstr))
    if (thetable$numreflxsnp[i]>0)
      thetable$flxsnp[i]=intersect(selsnps_hg19str,refluxstr)
    if (gene %in% allflxgenes)
      thetable$flxgene[i]=T
  }
  return(thetable)
}
table11=check_reflux()
table11$numrefluxsnp
table11$flxgene

table2=read_excel("../result/TWAS_tables.xlsx",2)
table2=table2[-1,]
table22=check_reflux(thetable=table2)
table22$numreflxsnp
table22$flxgene
table22$Gene[which(table22$flxgene==T)] #CRCT1
table22[which(table22$numreflxsnp>0),c("Gene","flxsnp")]
#ISYNA1
which(refluxstr=="19:18449238")
snps[29] #"rs9636202"
