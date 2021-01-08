#!/usr/bin/env Rscript
library(data.table)

ME=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GSE72874/GSE72874-GPL13534_series_matrix.txt",skip=63,sep="\t",header = T))
rownames(ME)=ME$ID_REF
ME=ME[,-1]
for (i in 1:ncol(ME)) ME[,i]=as.numeric(ME[,i])
sum(is.na(ME))/nrow(ME)/ncol(ME) #[1] 0.0005562037
idx=grepl("^cg",rownames(ME))
table(idx) #non-CpG probes
# FALSE   TRUE 
# 1185 372376 
ME=ME[idx,]
MEall=ME

#use subset of probes
EACprobes=read.csv("../data/supp_bgw018_Copy_of_Supplementary_Table3.csv",skip = 1,header = T,stringsAsFactors = F)
BEprobes=read.csv("../data/supp_bgw018_Copy_of_Supplementary_Table4.csv",skip = 1,header = T,stringsAsFactors = F)
#selectedprobes=intersect(EACprobes$TargetID,BEprobes$TargetID)
#length(selectedprobes) #39191
selectedprobes=unique(c(EACprobes$TargetID,BEprobes$TargetID))
length(selectedprobes) #63500
#use all probes
#selectedprobes=rownames(MEall)

ME=MEall[rownames(MEall) %in% selectedprobes,]


MEIDtable=read.table("../data/GSE72874/methylationID.txt",stringsAsFactors = F)
library(xlsx)
clinicaltable=read.xlsx("../data/Copy of Copy_of_Supplementary_Table2.xlsx",sheetIndex = 1,startRow = 2,header = T)
idx=match(clinicaltable$Methylation.ArrayID,MEIDtable$V3)
rownames(clinicaltable)=MEIDtable$V1[idx]
idx=match(colnames(ME),rownames(clinicaltable))
clinicaltable=clinicaltable[idx,]
clinicaltable$Array.batch=factor(clinicaltable$Array.batch,levels = c(1,2,3,4))
clinicaltable$NSE=F
clinicaltable$NSE[clinicaltable$Sample.type %in% c("Control","GERD","Normal")]=T
for (i in 1:ncol(clinicaltable))
{
  idx=which(clinicaltable[,i]==".")
  if (length(idx)>0)
  {
    clinicaltable[idx,i]=NA
  }
}
save(ME,MEall,clinicaltable,file="../data/KrauseData.RData")


#work on Grady's data:---------------------------------------------------------------
#to get the clinical info of an GEO:
#https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE104707
readGSE=function(GSEfile="../data/GSE89181/GSE89181-GPL13534_series_matrix.txt")
{
  cmd=paste0("grep -n ","Sample_title"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  studyno=unlist(strsplit(tmp1,"\t\""))
  studyno=gsub("\"","",studyno)
  studyno=studyno[2:length(studyno)]
  cmd=paste0("grep -n ","Sample_geo_accession"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  sampleid=unlist(strsplit(tmp1,"\t\""))
  sampleid=gsub("\"","",sampleid)
  sampleid=sampleid[2:length(sampleid)]
  cmd=paste0("grep -n ","ID_REF"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=as.integer(unlist(strsplit(tmp1,":"))[1])
  GSE=as.data.frame(fread(GSEfile,skip=tmp1-1,header = T,sep="\t",fill=T))
  idx=grepl("^cg",GSE[,1])
  GSE=GSE[idx,]
  rownames(GSE)=GSE[,1]
  GSE=GSE[,-1]
  if(all(colnames(GSE)==sampleid)) colnames(GSE)=studyno
  # histology=NA
  # cmd=paste0("grep -n ","histology:"," ",GSEfile)
  # tmp1=system(cmd,intern = TRUE)
  # tmp1=tmp1 [length(tmp1)]
  # if (!is.na(tmp1))
  # {
  #   histology=unlist(strsplit(tmp1,"\t\""))
  #   histology=gsub("\"","",histology)
  #   histology=histology[2:length(histology)]
  #   histology=gsub("histology:","",histology)
  #   histology=gsub(" ","",histology)
  # }
  # age=NA
  # cmd=paste0("grep -n ","age:"," ",GSEfile)
  # tmp1=system(cmd,intern = TRUE)
  # tmp1=tmp1 [length(tmp1)]
  # if (!is.na(tmp1))
  # {
  #   age=unlist(strsplit(tmp1,"\t\""))
  #   age=gsub("\"","",age)
  #   age=age[2:length(age)]
  #   age=gsub("age:","",age)
  #   age=gsub(" ","",age)
  # }
 return(GSE)
}
GSE89181_GPL13534=readGSE()
GSE89181_GPL18809=readGSE(GSEfile = "../data/GSE89181/GSE89181-GPL18809_series_matrix.txt")
all(rownames(GSE89181_GPL13534)==rownames(GSE89181_GPL18809)) #T
sum(colnames(GSE89181_GPL13534) %in% colnames(GSE89181_GPL18809)) #0
table_GPL13534=read.csv("../data/GSE89181/GSE89181_GPL13534.csv")
for (i in 1:ncol(table_GPL13534))
{
  idx=which(table_GPL13534[,i]=="")
  if (length(idx)>0) table_GPL13534[idx,i]=NA
}
table_GPL18809=read.csv("../data/GSE89181/GSE89181_GPL18809.csv")
for (i in 1:ncol(table_GPL18809))
{
  idx=which(table_GPL18809[,i]=="")
  if (length(idx)>0) table_GPL18809[idx,i]=NA
}
all(colnames(GSE89181_GPL13534)==table_GPL13534$Title) #T
all(colnames(GSE89181_GPL18809)==table_GPL18809$Title) #T
idx=match(colnames(table_GPL13534),colnames(table_GPL18809))
table_GPL18809=table_GPL18809[,idx]
all(colnames(table_GPL13534)==colnames(table_GPL18809))
clinicaltable=rbind(table_GPL13534,table_GPL18809)
rownames(clinicaltable)=clinicaltable$Title
clinicaltable$Histology=as.character(clinicaltable$Histology)
clinicaltable$Histology[which(clinicaltable$Histology=="HGD/Ca")]="HGD"
clinicaltable$Bmi[which(clinicaltable$Bmi==0)]=NA
clinicaltable$Alcoholicrinksperweek[which(clinicaltable$Alcoholicrinksperweek=="Unknown")]=NA
clinicaltable$GPL18809=T
clinicaltable$GPL18809[1:nrow(table_GPL13534)]=F
MEall=cbind(GSE89181_GPL13534,GSE89181_GPL18809)
all(colnames(MEall)==clinicaltable$Title) #T


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
table(rownames(GSEmeth) %in% anno$Name)
# TRUE 
# 453444 
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr)  
anno <- anno[anno$chr!="X"&anno$chr!="Y",]
anno$chr <- as.numeric(anno$chr)
anno$pos <- as.numeric(anno$pos)
anno=as.data.frame(anno)
anno=anno[order(anno$chr,anno$pos),]
table(rownames(MEall) %in% anno$Name)
# FALSE   TRUE 
# 225 453219 
MEall=MEall[rownames(MEall) %in% anno$Name,]

#get CpGs close to 26 SNPs
snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_26highrisk_SNP_POS.txt"
snppos=read.table(snpposfile,header = T,stringsAsFactors = F)
library(GenomicRanges)
idx=match(rownames(MEall),anno$Name)
gr_MEall=GRanges(seqnames = anno$chr[idx],ranges = IRanges(start=anno$pos[idx],width=1))
gr_snppos=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$position,width=1))
selectedCpGs=NULL
distcutoff=5e5
for (i in 1:length(gr_snppos))
{
  tmp=distance(gr_MEall,gr_snppos[i])
  selectedCpGs=unique(c(selectedCpGs,rownames(MEall)[which(tmp<distcutoff)]))
}
idx=match(selectedCpGs,rownames(MEall))
ME=MEall[idx,]
dim(ME)
# [1] 7017   81

save(MEall,ME,clinicaltable,file="../data/KazData.RData")

readGSE1=function(GSEfile="../data/GSE89181/GSE89181-GPL13534_series_matrix.txt")
{
  cmd=paste0("grep -n ","Sample_title"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  studyno=unlist(strsplit(tmp1,"\t\""))
  studyno=gsub("\"","",studyno)
  studyno=studyno[2:length(studyno)]
  cmd=paste0("grep -n ","Sample_geo_accession"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  sampleid=unlist(strsplit(tmp1,"\t\""))
  sampleid=gsub("\"","",sampleid)
  sampleid=sampleid[2:length(sampleid)]
  cmd=paste0("grep -n ","ID_REF"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=as.integer(unlist(strsplit(tmp1,":"))[1])
  GSE=as.data.frame(fread(GSEfile,skip=tmp1-1,header = T,sep="\t",fill=T))
  idx=grepl("^cg",GSE[,1])
  GSE=GSE[idx,]
  rownames(GSE)=GSE[,1]
  GSE=GSE[,-1]
  all(colnames(GSE)==sampleid)
  #if(all(colnames(GSE)==sampleid)) colnames(GSE)=studyno
  type=rep(NA,length(studyno))
  for (i in 1:length(type))
  {
    tmp=unlist(strsplit(studyno[i]," "))
    type[i]=tmp[length(tmp)]
  }
  gender=NA
  cmd=paste0("grep -n ","gender:"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    gender=unlist(strsplit(tmp1,"\t\""))
    gender=gsub("\"","",gender)
    gender=gender[2:length(gender)]
    gender=gsub("gender:","",gender)
    gender=gsub(" ","",gender)
  }
  gender[which(!gender %in% c("female","male"))]=NA
  age=NA
  cmd=paste0("grep -n ","age:"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    age=unlist(strsplit(tmp1,"\t\""))
    age=gsub("\"","",age)
    age=age[2:length(age)]
    age=gsub("age:","",age)
    age=gsub(" ","",age)
  }
  age=as.numeric(age)
  clinicaltable=data.frame(type=type,age=age,gender=gender,stringsAsFactors = F)
  rownames(clinicaltable)=sampleid
  return(res=list(GSE=GSE,clinicaltable=clinicaltable))
}
GSE81334=readGSE1(GSEfile = "../data/GSE81334/GSE81334_series_matrix.txt")
MEall=GSE81334$GSE
clinicaltable=GSE81334$clinicaltable
all(colnames(MEall)==rownames(clinicaltable))
save(MEall,clinicaltable,file="../data/GSE81334.RData")
