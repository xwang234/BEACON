
setwd("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/locuszoom")
library(data.table)
#wget -t 1220 -c --waitretry=28 --wait=33 --random-wait --referer="" --user-agent="" --limit-rate=10000k -e robots=off http://csg.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
#ml Python/2.7.15-foss-2016b
#module load anaconda2/2.4
###export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin
#export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9

locuszoom="/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom"
allsnp=data.frame(fread("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/annotation/All_SNPs.txt"))
# modeldat="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData"
# modelres="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18/preidiction_michigan_model.RData"
# skatres="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18/skat_res.RData"
# gene="KXD1"
# prefix="BEEA_BEACON"
# tmp=as.character(Sys.Date())
# tmp=unlist(strsplit(tmp,"-"))
# outfilename=paste0(prefix,"_20",tmp[2],tmp[3],"_",gene,".pdf")
summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt"
#locuszoom --metal ./Kathiresan_2009_HDL.txt --refgene FADS1 --build hg19 --pop EUR --source 1000G_Nov2010 --db /app/locuszoom/1.1/data/database/locuszoom_hg19.db --flank 500k --pvalcol "P-value"

generateMETALfile=function(gene=NULL)
{
  idx=which(rownames(phenotypepos)==gene)
  chr=phenotypepos$chr[idx]
  startpos=phenotypepos$s1[idx]-5e5
  endpos=phenotypepos$s2[idx]+5e5
  summarydat=as.data.frame(fread(summaryfile,header=T))
  if (sum(colnames(summarydat)=="CHR")==0)
  {
    summarydat=summarydat[summarydat$SNP %in% allsnp$name,]
    idx=match(summarydat$SNP,allsnp$name)
    summarydat$CHR=gsub("chr","",allsnp$X.chrom[idx])
  }
  idx=which(summarydat$CHR==chr & summarydat$position>=startpos & summarydat$position<=endpos)
  summarydat1=summarydat[idx,]
  idx=which(summarydat1$position %in% snppos$pos[snppos$chr==chr]) #snps in discovery
  summarydat1=summarydat1[idx,]
  METALdat=data.frame(SNP=summarydat1$SNP,MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                      alt=summarydat1$effect_allele,modelselect=F,color="gray90",stringsAsFactors = F)
  METALdat$MarkerName1=paste0("chr",summarydat1$CHR,":",summarydat1$position)
  idx=which(rownames(res_min) ==gene)
  selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
  modeldat=data.frame(selsnp=selsnps,MarkerName=NA,MakerName1=NA,ref=NA,alt=NA,stringsAsFactors = F)
  for (i in 1:nrow(modeldat))
  {
    tmp=modeldat$selsnp[i]
    tmp=unlist(strsplit(tmp,"_"))
    modeldat$MarkerName1[i]=paste0("chr",tmp[1])
    modeldat$ref[i]=tmp[3]
    modeldat$alt[i]=tmp[2]
    idx=which(METALdat$MarkerName1==modeldat$MarkerName1[i])
    if (length(idx)>0)
    {
      if ((METALdat$ref[idx]==modeldat$ref[i] & METALdat$alt[idx]==modeldat$alt[i])|(METALdat$ref[idx]==modeldat$alt[i] & METALdat$alt[idx]==modeldat$ref[i]))
      {
        METALdat$modelselect[idx]=T
        METALdat$color[idx]="red"
      }
    }
  }
  idx=METALdat$modelselect==T
  METALdat=rbind(METALdat[!idx,],METALdat[idx,]) # put selected in the end
  #idx=which(METALdat$modelselect==T)
  #markerdat=data.frame(snp=METALdat$MarkerName[idx],string=" ",color="red",stringsAsFactors = F)
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  write.table(METALdat[,2:7],file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent")
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

locusplot=function(modeldat="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData",
                   modelres="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18/preidiction_michigan_model.RData",
                   skatres="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18/skat_res.RData",
                   prefix="BEEA_BEACON",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt",
                   genes=c("KXD1","HSP90AA1")
)
{
  load(modeldat)
  load(modelres)
  load(skatres)
  for (gene in genes)
  {
    #generateMETALfile(gene=gene)
    idx=which(rownames(skat_min2)==gene)
    print(skat_min2[idx,])
  }
}

locusplot(prefix="BEEA_Cambridge",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt")
locusplot(prefix="BEEA_Bonn",summaryfile="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt")

#updated functions 6/6/2021
generateMETALfile1=function(gene=NULL,prefix="",summaryfile=NULL)
{
  idx=which(rownames(phenotypepos)==gene)
  chr=phenotypepos$chr[idx]
  startpos=phenotypepos$s1[idx]-5e5
  endpos=phenotypepos$s2[idx]+5e5
  summarydat=as.data.frame(fread(summaryfile,header=T))
  if (sum(colnames(summarydat)=="CHR")==0)
  {
    summarydat=summarydat[summarydat$SNP %in% allsnp$name,]
    idx=match(summarydat$SNP,allsnp$name)
    summarydat$CHR=gsub("chr","",allsnp$X.chrom[idx])
  }
  idx=which(summarydat$CHR==chr & summarydat$position>=startpos & summarydat$position<=endpos)
  summarydat1=summarydat[idx,]
  idx=which(summarydat1$position %in% snppos$pos[snppos$chr==chr]) #snps in discovery
  summarydat1=summarydat1[idx,]
  METALdat=data.frame(SNP=summarydat1$SNP,MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                      alt=summarydat1$effect_allele,modelselect=F,color="gray90",stringsAsFactors = F)
  METALdat$MarkerName1=paste0("chr",summarydat1$CHR,":",summarydat1$position)
  idx=which(rownames(res_min) ==gene)
  selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
  modeldat=data.frame(selsnp=selsnps,MarkerName=NA,MakerName1=NA,ref=NA,alt=NA,stringsAsFactors = F)
  for (i in 1:nrow(modeldat))
  {
    tmp=modeldat$selsnp[i]
    tmp=unlist(strsplit(tmp,"_"))
    modeldat$MarkerName1[i]=paste0("chr",tmp[1])
    modeldat$ref[i]=tmp[3]
    modeldat$alt[i]=tmp[2]
    idx=which(METALdat$MarkerName1==modeldat$MarkerName1[i])
    if (length(idx)>0)
    {
      if ((METALdat$ref[idx]==modeldat$ref[i] & METALdat$alt[idx]==modeldat$alt[i])|(METALdat$ref[idx]==modeldat$alt[i] & METALdat$alt[idx]==modeldat$ref[i]))
      {
        METALdat$modelselect[idx]=T
        METALdat$color[idx]="red"
      }
    }
  }
  idx=METALdat$modelselect==T
  METALdat=rbind(METALdat[!idx,],METALdat[idx,]) # put selected in the end
  #idx=which(METALdat$modelselect==T)
  #markerdat=data.frame(snp=METALdat$MarkerName[idx],string=" ",color="red",stringsAsFactors = F)
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  write.table(METALdat[,2:7],file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent")
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

locusplot1=function(organ="mucosa",genes=c("KXD1","HSP90AA1"),prefix="BEEA_BC",summaryfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BEEA_1.txt")
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
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
    generateMETALfile1(gene=gene,prefix=prefix,summaryfile = summaryfile)
    idx=which(rownames(skat_min2)==gene)
    print(skat_min2[idx,])
  }
}



generateMETALfile1=function(gene=NULL,prefix=NULL,summaryfile=NULL,res_min=NULL)
{
  idx=which(rownames(phenotypepos)==gene)
  chr=phenotypepos$chr[idx]
  startpos=phenotypepos$s1[idx]-5e5
  endpos=phenotypepos$s2[idx]+5e5
  summarydat=as.data.frame(fread(summaryfile,header=T,sep="\t",fill=T))
  if (sum(colnames(summarydat)=="CHR")==0)
  {
    summarydat=summarydat[summarydat$SNP %in% allsnp$name,]
    idx=match(summarydat$SNP,allsnp$name)
    summarydat$CHR=gsub("chr","",allsnp$X.chrom[idx])
  }
  idx=which(summarydat$CHR==chr & summarydat$position>=startpos & summarydat$position<=endpos)
  summarydat1=summarydat[idx,]
  idx=which(summarydat1$position %in% snppos$pos[snppos$chr==chr]) #snps in discovery
  summarydat1=summarydat1[idx,]
  METALdat=data.frame(SNP=summarydat1$SNP,MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                      alt=summarydat1$effect_allele,modelselect=F,color="gray90",stringsAsFactors = F)
  METALdat$MarkerName1=paste0("chr",summarydat1$CHR,":",summarydat1$position)
  idx=which(rownames(res_min) ==gene)
  selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
  modeldat=data.frame(selsnp=selsnps,MarkerName=NA,MakerName1=NA,ref=NA,alt=NA,stringsAsFactors = F)
  for (i in 1:nrow(modeldat))
  {
    tmp=modeldat$selsnp[i]
    tmp=unlist(strsplit(tmp,"_"))
    modeldat$MarkerName1[i]=paste0("chr",tmp[1])
    modeldat$ref[i]=tmp[3]
    modeldat$alt[i]=tmp[2]
    idx=which(METALdat$MarkerName1==modeldat$MarkerName1[i])
    if (length(idx)>0)
    {
      if ((METALdat$ref[idx]==modeldat$ref[i] & METALdat$alt[idx]==modeldat$alt[i])|(METALdat$ref[idx]==modeldat$alt[i] & METALdat$alt[idx]==modeldat$ref[i]))
      {
        METALdat$modelselect[idx]=T
        METALdat$color[idx]="red"
      }
    }
  }
  idx=METALdat$modelselect==T
  METALdat=rbind(METALdat[!idx,],METALdat[idx,]) # put selected in the end
  #idx=which(METALdat$modelselect==T)
  #markerdat=data.frame(snp=METALdat$MarkerName[idx],string=" ",color="red",stringsAsFactors = F)
  METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
  write.table(METALdat[,2:7],file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 500kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color refsnpTextColor=transparent")
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

# updatemetafile=function(metafile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BEEA_1.tbl",
#                         metafile1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt",
#                         metafile2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt")
# {
#   metadat=as.data.frame(fread(metafile))
#   metadat1=as.data.frame(fread(metafile1))
#   metadat2=as.data.frame(fread(metafile2))
#   colnames(metadat)[which(colnames(metadat)=="MarkerName")]="SNP"
#   colnames(metadat)[which(colnames(metadat)=="P-value")]="P"
#   metadat$Allele1=toupper(metadat$Allele1)
#   metadat$Allele2=toupper(metadat$Allele2)
#   comsnps1=intersect(metadat$SNP,metadat1$SNP)
#   comsnps2=intersect(metadat$SNP,metadat2$SNP)
#   comsnps=intersect(comsnps1,comsnps2)
#   metadat=metadat[match(comsnps,metadat$SNP),]
#   metadat1=metadat1[match(comsnps,metadat1$SNP),]
#   metadat2=metadat2[match(comsnps,metadat2$SNP),]
#   metadat$non_effect_allele=metadat1$non_effect_allele
#   metadat$effect_allele=metadat1$effect_allele
#   metadat$CHR=metadat1$CHR
#   metadat$position=metadat1$position
#   outfile=gsub("tbl","txt",metafile)
#   fwrite(metadat,file=outfile,sep="\t",quote=F)
# }
# 
# updatemetafile(metafile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BE_1.tbl",
#                         metafile1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_BEACON_autosomes_N.txt",
#                         metafile2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt")
