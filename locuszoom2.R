#for TWAS results based on include all covariates in GTEx models
#arguments of locuszoom are in /fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom.R

setwd("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/locuszoom/allcovar")
library(data.table)
#wget -t 1220 -c --waitretry=28 --wait=33 --random-wait --referer="" --user-agent="" --limit-rate=10000k -e robots=off http://csg.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
#ml Python/2.7.15-foss-2018b #this is needed
###export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin
#export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9 #this is needed

#check if you can run locuszoom:
system("/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom --help")
#if you can't run it in rstudio, use R

# Chromosome 14, DYNC1H1 and HSP90AA1
# Chromosome 6, FILIP1/COX7A2/THEM30A/SENP6/MYO6
library(data.table)
#summary stat from validation
summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"
BE_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_Bonn_autosomes.txt"),header=T))
BE_Oxfordsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_oxford_autosomes.txt"),header=T))

EA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"EA_Bonn_autosomes.txt"),header=T))
BEEA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"),header=T))

tmp=intersect(BE_Bonnsummarydat$SNP,EA_Bonnsummarydat$SNP)

tmp=intersect("SNP",colnames(BE_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="pvalue")]="P"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="beta")]="BETA"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

tmp=intersect("SNP",colnames(BE_Oxfordsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="rsid")]="SNP"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="pvalue")]="P"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="beta")]="BETA"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="non-effect-allele")]="non_effect_allele"
}



tmp=intersect("SNP",colnames(EA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="pvalue")]="P"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="beta")]="BETA"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}


tmp=intersect("SNP",colnames(BEEA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="pvalue")]="P"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="beta")]="BETA"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

#add chr based on summarydata BE_Oxfordsummarydat
add_chr_tosummarydat=function(dat1=BE_Bonnsummarydat)
{
  dat2=BE_Oxfordsummarydat
  if (!"chr" %in% colnames(dat1))
  {
    dat1$chr=NA
    tmp=intersect(dat1$SNP,dat2$SNP)
    idx1=match(tmp,dat1$SNP)
    idx2=match(tmp,dat2$SNP)
    dat1$chr[idx1]=dat2$chr[idx2]
    # idx=which(is.na(dat1$chr))
    # tmp=intersect(dat1$SNP[idx],dbsnp$V2)
    # idx1=match(tmp,dat1$SNP)
    # idx2=match(tmp,dbsnp$V2)
    # dat1$chr[idx1]=dbsnp$V1[idx2]
    print(paste0(sum(is.na(dat1$chr))/nrow(dat1),"of all snps are missing chr"))
  }
  dat1$snpname=dat1$snpname1=NA
  dat1$snpname=paste0(dat1$chr,":",dat1$position,"_",dat1$effect_allele,"_",dat1$non_effect_allele)
  dat1$snpname1=paste0(dat1$chr,":",dat1$position,"_",dat1$non_effect_allele,"_",dat1$effect_allele)
  return(dat1)
}
BE_Bonnsummarydat=add_chr_tosummarydat(dat1=BE_Bonnsummarydat)
EA_Bonnsummarydat=add_chr_tosummarydat(dat1=EA_Bonnsummarydat)
BEEA_Bonnsummarydat=add_chr_tosummarydat(dat1=BEEA_Bonnsummarydat)
BE_Oxfordsummarydat=add_chr_tosummarydat(dat1=BE_Oxfordsummarydat)


locuszoom="/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom"

#updated functions 6/6/2021
library("biomaRt")
#mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
#snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
#snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp") #hg38
snpmart= useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
findrsid=function(selsnps=NULL,verbose=1)
{
  print(paste0(length(selsnps)," snprsid to be found"))
  chr=unlist(strsplit(selsnps[1],":"))[1]
  pos=allele1=allele2=rep(NA,length(selsnps))
  res=data.frame(selsnps=selsnps,rsid=NA,ref=NA,alt=NA,gene=NA)
  idx=which(is.na(res$rsid))
  n=0
  while(n<1 & length(idx)>0)
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
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand','associated_gene'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = T)
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

zl_validaton=function(gene="DYNC1H1",prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",
                      res_min=NULL,ymax=NULL,nsnp=NULL)
{
  # idx=which(rownames(phenotypepos)==gene)
  # chr=phenotypepos$chr[idx]
  # startpos=phenotypepos$s1[idx]-5e5
  # endpos=phenotypepos$s2[idx]+5e5
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  if (!file.exists(METALfile))
  {
    summarydat=as.data.frame(fread(summaryfile,header=T))
    #for Oxford data
    colnames(summarydat)[colnames(summarydat)=="rsid"]="SNP"
    colnames(summarydat)[colnames(summarydat)=="pvalue"]="P"
    colnames(summarydat)[colnames(summarydat)=="non-effect-allele"]="non_effect_allele"
    colnames(summarydat)[colnames(summarydat)=="effect-allele"]="effect_allele"
    idx=which(rownames(res_min) ==gene)
    selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
    rsid=findrsid(selsnps = selsnps)
    if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
    
    summarydat1=summarydat[summarydat$SNP %in% rsid$rsid,]
    
    METALdat=data.frame(MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                        alt=summarydat1$effect_allele,color="red",stringsAsFactors = F)
    
    METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
    write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  }else
  {
    METALdat=read.table(METALfile,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  if (is.null(nsnp))
  {
    title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  }else
  {
    title=paste0(prefix,", ",gene,", ",nsnp," SNPs")
  }
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

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

find_snpname=function(selsnps,summarydat=BE_Bonnsummarydat)
{
  selsnps_hg19=hg38tohg19(snpnames = selsnps)
  selsnps_hg19$selsnps=selsnps
  selsnps_hg19$snp=NA
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname1)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname1)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  return(selsnps_hg19)
}


locusplot_validation=function(organ="blood",genes=c("SENP6"),prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=NULL,nsnp=NULL)
{
  prefix=paste0(organ,"_",prefix)
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  # modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
  # load(modeldat) #phenotypepos
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
  # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
  
  # load(modelres)
  # load(skatres)
  for (gene in genes)
  {
    print(gene)
    zl_validaton(gene=gene,summaryfile = summaryfile,prefix=prefix,res_min = res_min,ymax=ymax,nsnp=nsnp)
  }
}

#locusplot_validation()
locusplot_validation(organ="junction",genes="ZNF641",prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=5)
#locusplot_validation(organ="blood",genes="HSP90AA1",prefix="BE_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt")
#locusplot_validation(organ="blood",genes="HSP90AA1",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=5)
locusplot_validation(organ="blood",genes="HSP90AA1",prefix="BE_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt",ymax=5)

#locusplot_validation(organ="blood",genes="HSP90AA1",prefix="EA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Bonn_autosomes.txt")

#locusplot_validation(organ="stomach",genes="CFDP1",prefix="EA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Bonn_autosomes.txt")
#locusplot_validation(organ="blood",genes="COX7A2",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt")
#locusplot_validation(organ="stomach",genes="CFDP1",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=5)
locusplot_validation(organ="stomach",genes="CFDP1",prefix="EA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Bonn_autosomes.txt",ymax=5.5)

locusplot_validation(organ="adipose",genes="EXOC3",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=5,nsnp=26)
#locusplot_validation(organ="blood",genes="JUND",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=7)
locusplot_validation(organ="blood",genes="JUND",prefix="EA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Bonn_autosomes.txt",ymax=7,nsnp=301)

locusplot_validation(organ="blood",genes="SENP6",prefix="BEEA_Bonn",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",ymax=6)

#add missing snps within 50BP
zl_validaton1=function(gene="HSP90AA1",prefix="BE_Bonn",summarydat=BE_Bonnsummarydat,
                       res_min=NULL,ymax=NULL,nsnp=NULL)
{
  # idx=which(rownames(phenotypepos)==gene)
  # chr=phenotypepos$chr[idx]
  # startpos=phenotypepos$s1[idx]-5e5
  # endpos=phenotypepos$s2[idx]+5e5
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  # if (!file.exists(METALfile))
  # {
  #summarydat=as.data.frame(fread(summaryfile,header=T))
  #for Oxford data
  colnames(summarydat)[colnames(summarydat)=="rsid"]="SNP"
  colnames(summarydat)[colnames(summarydat)=="pvalue"]="P"
  colnames(summarydat)[colnames(summarydat)=="non-effect-allele"]="non_effect_allele"
  colnames(summarydat)[colnames(summarydat)=="effect-allele"]="effect_allele"
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
  selsnpmap=tmp
  rsid=selsnpmap$snp
  
  if (sum(is.na(rsid)>0)) print("some rsid not been found!")
  idx0=which(!rsid %in% summarydat$SNP)
  missingrsid=rsid[!rsid %in% summarydat$SNP]
  
  if (length(missingrsid)>0)
  {
    refreplacesnps1=data.frame(missingid=missingrsid,hg19pos=selsnpmap$snphg19[idx0],replaceid=NA,dist=NA)
    
    for (i in 1:length(missingrsid))
    {
      
      tmp1=unlist(strsplit(refreplacesnps1$hg19pos[i],":"))
      tmp2=as.integer(unlist(strsplit(tmp1[2],"_"))[1])
      idx=which(summarydat$chr==tmp1[1])
      idx1=which.min(abs(summarydat$position[idx]-tmp2))
      refreplacesnps1$replaceid[i]=summarydat$SNP[idx[idx1]]
      refreplacesnps1$dist[i]=min(abs(summarydat$position[idx]-tmp2))
    }
    
  }
  
  allrsid=c(rsid,refreplacesnps1$replaceid[which(refreplacesnps1$dist<=50)])
  summarydat1=summarydat[summarydat$SNP %in% allrsid,]
  
  METALdat=data.frame(MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                      alt=summarydat1$effect_allele,color="red",stringsAsFactors = F)
  
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  #}#else
  # {
  #   METALdat=read.table(METALfile,header=T)
  # }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  if (is.null(nsnp))
  {
    title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  }else
  {
    title=paste0(prefix,", ",gene,", ",nsnp," SNPs")
  }
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}
locusplot_validation1=function(organ="blood",genes=c("SENP6"),prefix="BEEA_Bonn",summarydat=BEEA_Bonnsummarydat,ymax=NULL,nsnp=NULL)
{
  prefix=paste0(organ,"_",prefix)
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  # modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
  # load(modeldat) #phenotypepos
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
  # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
  
  # load(modelres)
  # load(skatres)
  for (gene in genes)
  {
    print(gene)
    zl_validaton1(gene=gene,summarydat = summarydat,prefix=prefix,res_min = res_min,ymax=ymax,nsnp=nsnp)
  }
}

locusplot_validation1(organ="adipose",genes="EXOC3",prefix="BEEA_Bonn",summarydat=BEEA_Bonnsummarydat,ymax=5,nsnp=27)
locusplot_validation1(organ="blood",genes="HSP90AA1",prefix="BE_Bonn",summarydat=BE_Bonnsummarydat,ymax=5,nsnp=90)


#work on GTEx data
zl_discovery=function(gene="LDAH",prefix="Mucosa_GTEx",res_min=NULL,snp=NULL,covariate=NULL,phenotype=NULL,nsnp=NULL)
{
  
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  if (!file.exists(METALfile))
  {
    idx=which(rownames(res_min) ==gene)
    selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
    rsid=findrsid(selsnps = selsnps)
    if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
    METALdat=data.frame(MarkerName=rsid$rsid,Pvalue=NA,ref=rsid$ref,
                        alt=rsid$alt,color="red",stringsAsFactors = F)
    idx=which(rownames(phenotype)==gene)
    Y=as.numeric(phenotype[idx,])
    for (i in 1:length(selsnps)) {
      Xmat <- data.matrix(cbind(t(snp[rownames(snp)==selsnps[i],]),covariate))
      fit <- lm(Y~Xmat)
      METALdat$Pvalue[i] <- summary(fit)$coef[2,4]
    }  
    
    METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
    write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  }else
  {
    METALdat=read.table(METALfile,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax)
  if (is.null(nsnp))
  {
    title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  }else
  {
    title=paste0(prefix,", ",gene,", ",nsnp," SNPs")
  }
  
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

locusplot_discovery=function(organ="blood",genes=c("SENP6"),prefix="GTEx")
{
  prefix=paste0(organ,"_",prefix)
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
  load(modeldat) #phenotypepos
  # load(paste0(outfolder,"/skat_res.RData"))
  # colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") 
  # skatres=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenestable$Symbol,]
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
  # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
  
  # load(modelres)
  # load(skatres)
  for (gene in genes)
  {
    print(gene)
    zl_discovery(gene=gene,prefix=prefix,res_min = res_min,snp = snp,phenotype = phenotype,covariate = covariate)
    # idx=which(rownames(skat_min2)==gene)
    # print(skat_min2[idx,])
  }
}

# #locusplot_discovery()
# locusplot_discovery(organ="junction",genes="ZNF641")
# locusplot_discovery(organ="blood",genes="HSP90AA1")
# locusplot_discovery(organ="stomach",genes="CFDP1")
# #locusplot_discovery(organ="blood",genes="COX7A2")
# locusplot_discovery(organ="blood",genes="SENP6")
# locusplot_discovery(organ="adipose",genes="EXOC3")
# locusplot_discovery(organ="blood",genes="JUND")

#work on BEACON and Cambridge
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
opt="PC6"
#try different number of PCs
if (opt=="PC4")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","sex")]
if (opt=="PC6")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]



zl_BC=function(gene="LDAH",prefix="Mucosa_BC",res_min=NULL,Covariate=NULL,covariatetable=NULL,bcagenotype=NULL,type="BE",ymax=NULL,nsnp=NULL)
{
  
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  if (!file.exists(METALfile))
  {
    idx1=match(colnames(bcagenotype),rownames(covariatetable))
    covariatetable=covariatetable[idx1,]
    idx1=match(colnames(bcagenotype),rownames(Covariate))
    Covariate=Covariate[idx1,]
    
    idx=which(rownames(res_min) ==gene)
    selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
    rsid=findrsid(selsnps = selsnps)
    if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
    METALdat=data.frame(MarkerName=rsid$rsid,Pvalue=NA,ref=rsid$ref,
                        alt=rsid$alt,color="red",stringsAsFactors = F)
    
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
    for (i in 1:length(selsnps)) {
      Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
      fit <- glm(y~Xmat,family=binomial) 
      METALdat$Pvalue[i] <- summary(fit)$coef[2,4]
    }  
    
    METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
    write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  }else
  {
    METALdat=read.table(METALfile,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax))  ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax)
  
  if (is.null(nsnp))
  {
    title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  }else
  {
    title=paste0(prefix,", ",gene,", ",nsnp," SNPs")
  }
  
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

locusplot_BC=function(organ="blood",genes=c("SENP6"),prefix="BC",type="BEEA",ymax=NULL,nsnp=NULL)
{
  prefix=paste0(organ,"_",type,"_",prefix)
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  #no need to load bcagenotype data
  load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
  
  for (gene in genes)
  {
    print(gene)
    zl_BC(gene=gene,prefix=prefix,res_min = res_min,Covariate=Covariate,covariatetable=covariatetable,bcagenotype=bcagenotype,type=type,ymax=ymax,nsnp=nsnp)
    # idx=which(rownames(skat_min2)==gene)
    # print(skat_min2[idx,])
  }
}

#locusplot_BC()
# locusplot_BC(organ="junction",genes="ZNF641",type="BEEA")
# locusplot_BC(organ="blood",genes="HSP90AA1",type="BE")
# #locusplot_BC(organ="blood",genes="HSP90AA1",type="BEEA")
# locusplot_BC(organ="stomach",genes="CFDP1",type="BEEA")
# #locusplot_BC(organ="blood",genes="COX7A2",type="BEEA")
# locusplot_BC(organ="blood",genes="JUND",type="BEEA")
# locusplot_BC(organ="adipose",genes="EXOC3",type="BEEA",ymax=7.5)
# locusplot_BC(organ="blood",genes="SENP6",type="BEEA")

#plot 2 results together
plot_BC_Bonn=function(gene="EXOC3",organ="adipose",type1="BEEA",type2="BEEA")
{
  METALfile1=paste0("./",organ,"_",type1,"_BC_",gene,"_metalfile.txt")
  METALdat1=read.table(METALfile1,header=T)
  METALfile2=paste0("./",organ,"_",type2,"_Bonn_",gene,"_metalfile.txt")
  METALdat2=read.table(METALfile2,header=T)
  METALdat2$color="green"
  METALdat=rbind(METALdat1,METALdat2)
  METALfile=paste0("./",organ,"_BC_Bonn_",gene,"_metalfile.txt")
  write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  title=gene
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",gene," colorCol=color refsnpTextColor=transparent "," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9 condRefsnpPch=NULL condPch='1,1,1,1,1,1,1,1,1,1,1' showRefsnpAnnot=FALSE")
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",gene," colorCol=color refsnpTextColor=transparent "," signifLine=",-log10(0.05)," signifLineColor=blue title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9 condRefsnpPch=NULL")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}
# #plot_BC_Bonn()
# plot_BC_Bonn(gene="SENP6",organ="blood")
# plot_BC_Bonn(gene="ZNF641",organ="junction")
# plot_BC_Bonn(gene="HSP90AA1",organ="blood",type1="BE",type2="BE")
# plot_BC_Bonn(gene="CFDP1",organ="stomach",type2="EA")

# #work on two neighborhood genes----
# #put the most important genes first
# zl_BC2=function(genes=c("DYNC1H1","HSP90AA1"),prefix="Mucosa_BC",res_min=NULL,Covariate=NULL,covariatetable=NULL,bcagenotype=NULL,type="BE")
# {
#   genes1=paste0(genes[c(1,length(genes))],collapse = "_")
#   
#   idx1=match(colnames(bcagenotype),rownames(covariatetable))
#   covariatetable=covariatetable[idx1,]
#   idx1=match(colnames(bcagenotype),rownames(Covariate))
#   Covariate=Covariate[idx1,]
#   
#   idx=which(rownames(res_min) %in% genes)
#   selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
#   rsid=findrsid(selsnps = selsnps)
#   if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
#   METALdat=data.frame(MarkerName=rsid$rsid,Pvalue=NA,ref=rsid$ref,
#                       alt=rsid$alt,color="red",stringsAsFactors = F)
#   
#   correctedsnps=NULL
#   if (length(intersect(selsnps,rownames(bcagenotype)))<length(selsnps))
#   {
#     missingsnps=selsnps[!selsnps %in% rownames(bcagenotype)]
#     for (j in 1:length(missingsnps))
#     {
#       tmp=unlist(strsplit(missingsnps[j],"_"))
#       tmp1=paste0(tmp[c(1,3,2)],collapse = "_")  #change the order of allele
#       idx=which(rownames(bcagenotype)==tmp1)
#       if (length(idx)>0)
#       {
#         correctedsnps=c(correctedsnps,tmp1)
#         idx1=which(selsnps==missingsnps[j])
#         selsnps[idx1]=tmp1 #change the snp name to make it consistent with bca
#       }
#     }
#   } 
#   
#   idx=match(selsnps,rownames(bcagenotype))
#   ## there a few SNPs not found in bcagenotype
#   if (sum(is.na(idx))>0)
#     print(paste0(sum(is.na(idx))," out of ", length(idx)," selected snps not been found in genotype data"))
#   idx <- idx[!is.na(idx)]
#   Z=t(bcagenotype[idx,,drop=F])
#   if (length(correctedsnps)>0) #flip snps
#   {
#     idxtocorrect=match(correctedsnps,colnames(Z))
#     Z[,idxtocorrect]=2-Z[,idxtocorrect]
#   }
#   
#   if (type=="BE")
#   {
#     idx1=which(covariatetable$phenoBE_bca==2) #case
#     idx2=which(covariatetable$phenoBE_bca==1)
#   }
#   if (type=="BEEA")
#   {
#     idx1=which(covariatetable$phenoEABE_bca==2) #case
#     idx2=which(covariatetable$phenoEABE_bca==1)
#   }
#   if (type=="EA")
#   {
#     idx1=which(covariatetable$phenoEA_bca==2) #case
#     idx2=which(covariatetable$phenoEA_bca==1)
#   }
#   y=c(rep(1,length(idx1)),rep(0,length(idx2)))
#   Z1=Z[c(idx1,idx2),,drop=F]
#   Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
#   for (i in 1:length(selsnps)) {
#     Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
#     fit <- glm(y~Xmat,family=binomial) 
#     METALdat$Pvalue[i] <- summary(fit)$coef[2,4]
#   }  
#   
#   METALfile=paste0("./",prefix,"_metalfile.txt")
#   write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
#   #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
#   #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
#   #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
#   #print(cmd)
#   ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
#   #cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax)
#   title=paste0(prefix,", ",nrow(METALdat)," SNPs")
#   cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",genes[1]," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
#   
#   print(cmd)
#   system(cmd,wait=T)
#   return(0)
# }
# 
# locusplot_BC2=function(organ="muscularis",genes=c("DYNC1H1","HSP90AA1"),prefix="BC",type="BE")
# {
#   prefix=paste0(organ,"_",type,"_",prefix)
#   prefix=paste0(prefix,"_",paste0(genes[c(1,length(genes))],collapse = "_"))
#   outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
#   load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
#   #no need to load bcagenotype data
#   load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
#   
#   zl_BC2(genes=genes,prefix=prefix,res_min = res_min,Covariate=Covariate,covariatetable=covariatetable,bcagenotype=bcagenotype,type=type)
#  
# }
# 
# locusplot_BC2()
# locusplot_BC2(type="BEEA")
# locusplot_BC2(organ="blood")
# locusplot_BC2(organ="blood",type="BEEA")
# locusplot_BC2(organ="junction",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),type="BEEA")
# locusplot_BC2(organ="junction",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),type="BE")
# locusplot_BC2(organ="stomach",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),type="BEEA")
# locusplot_BC2(organ="stomach",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),type="BE")
# 
# zl_validaton2=function(genes=c("DYNC1H1","HSP90AA1"),prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt",
#                       res_min=NULL)
# {
#   # idx=which(rownames(phenotypepos)==gene)
#   # chr=phenotypepos$chr[idx]
#   # startpos=phenotypepos$s1[idx]-5e5
#   # endpos=phenotypepos$s2[idx]+5e5
#   genes1=paste0(genes[c(1,length(genes))],collapse = "_")
#   summarydat=as.data.frame(fread(summaryfile,header=T))
#   #for Oxford data
#   colnames(summarydat)[colnames(summarydat)=="rsid"]="SNP"
#   colnames(summarydat)[colnames(summarydat)=="pvalue"]="P"
#   colnames(summarydat)[colnames(summarydat)=="non-effect-allele"]="non_effect_allele"
#   colnames(summarydat)[colnames(summarydat)=="effect-allele"]="effect_allele"
#   idx=which(rownames(res_min) %in% genes)
#   selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
#   rsid=findrsid(selsnps = selsnps)
#   if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
#   summarydat1=summarydat[summarydat$SNP %in% rsid$rsid,]
#   
#   METALdat=data.frame(MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
#                       alt=summarydat1$effect_allele,color="red",stringsAsFactors = F)
#   
#   METALfile=paste0("./",prefix,"_metalfile.txt")
#   write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
#   #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
#   #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
#   #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
#   #print(cmd)
#   ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
#   title=paste0(prefix,", ",nrow(METALdat)," SNPs")
#   cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",genes[1]," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
#   #cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent")
#   print(cmd)
#   system(cmd,wait=T)
#   return(0)
# }
# 
# locusplot_validation2=function(organ="muscularis",genes=c("DYNC1H1","HSP90AA1"),prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt")
# {
#   prefix=paste0(organ,"_",prefix)
#   prefix=paste0(prefix,"_",paste0(genes[c(1,length(genes))],collapse = "_"))
#   outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
#   # modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
#   # load(modeldat) #phenotypepos
#   load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
#   # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
#   
#   # load(modelres)
#   # load(skatres)
#   
#   zl_validaton2(genes=genes,summaryfile = summaryfile,prefix=prefix,res_min = res_min)
#   
# }
# locusplot_validation2()
# locusplot_validation2(prefix="BE_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt")
# locusplot_validation2(summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_oxford_autosomes.txt",prefix="BE_Oxford")
# 
# locusplot_validation2(organ="blood")
# locusplot_validation2(organ="blood",prefix="BE_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt")
# locusplot_validation2(organ="blood",summaryfile = "/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_oxford_autosomes.txt",prefix="BE_Oxford")
# 
# locusplot_validation2(organ="junction",prefix="BEEA_Bonn",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt")
# locusplot_validation2(organ="junction",prefix="BE_Bonn",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt")
# locusplot_validation2(organ="junction",prefix="BE_Oxford",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_oxford_autosomes.txt")
# 
# locusplot_validation2(organ="stomach",prefix="BEEA_Bonn",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt")
# locusplot_validation2(organ="stomach",prefix="BE_Bonn",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt")
# locusplot_validation2(organ="stomach",prefix="BE_Oxford",genes=c("FILIP1","COX7A2","TMEM30A","SENP6","MYO6"),summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_oxford_autosomes.txt")
# 
# 
# library(GenomicRanges)
# zl_discovery2=function(genes=c("DYNC1H1","HSP90AA1"),prefix="Mucosa_GTEx",snp=NULL,res_min=NULL,snp=NULL,covariate=NULL,phenotype=NULL,phenotypepos=NULL)
# {
#   
#   genes1=paste0(genes[c(1,length(genes))],collapse = "_")
#   idx=which(rownames(res_min) %in% genes)
#   selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
#   rsid=findrsid(selsnps = selsnps)
#   if (sum(is.na(rsid$rsid))>0) print("some rsid not been found!")
#   METALdat=data.frame(MarkerName=rsid$rsid,Pvalue=NA,ref=rsid$ref,
#                       alt=rsid$alt,color="red",stringsAsFactors = F)
#   
#   idx=which(rownames(phenotypepos) %in% genes)
#   gr_genes=GRanges(seqnames = phenotypepos$chr[idx],ranges = IRanges(start=phenotypepos$s1[idx],end=phenotypepos$s2[idx]))
#   for (i in 1:length(selsnps)) {
#     Xmat <- data.matrix(cbind(t(snp[rownames(snp)==selsnps[i],]),covariate))
#     tmp=unlist(strsplit(selsnps[i],":"))
#     chr=tmp[1]
#     position=as.integer(unlist(strsplit(tmp[2],"_"))[1])
#     gr_snp=GRanges(seqnames = chr,ranges=IRanges(start=position,width=1))
#     tmp=distance(gr_genes,gr_snp) #pick the close one as the phenotype?
#     idx1=which.min(tmp)
#     Y=as.numeric(phenotype[idx[idx1],])
#     fit <- lm(Y~Xmat)
#     METALdat$Pvalue[i] <- summary(fit)$coef[2,4]
#   }  
#   
#   METALfile=paste0("./",prefix,"_","_metalfile.txt")
#   write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
#   #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
#   #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
#   #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
#   #print(cmd)
#   ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
#   #cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax)
#   title=paste0(prefix,", ",nrow(METALdat)," SNPs")
#   cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",genes[1]," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
#   
#   print(cmd)
#   system(cmd,wait=T)
#   return(0)
# }
# 
# locusplot_discovery2=function(organ="muscularis",genes=c("DYNC1H1","HSP90AA1"),prefix="GTEx")
# {
#   prefix=paste0(organ,"_",prefix)
#   prefix=paste0(prefix,"_",paste0(genes[c(1,length(genes))],collapse = "_"))
#   outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
#   modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
#   load(modeldat) #phenotypepos
#   # load(paste0(outfolder,"/skat_res.RData"))
#   # colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") 
#   # skatres=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenestable$Symbol,]
#   load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
#   # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
#   
#   # load(modelres)
#   # load(skatres)
#  
#   zl_discovery2(genes=genes,prefix=prefix,res_min = res_min,snp = snp,phenotype = phenotype,phenotypepos=phenotypepos,covariate = covariate)
#   
# }
# 
# locusplot_discovery2()

