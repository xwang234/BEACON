library(data.table)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684

library(data.table)


readGSE=function(GSEfile="../../../CancerGenomics/Colorectal_Cancer/Data/GSE48684_series_matrix.txt")
{
  cmd=paste0("grep -n ","ID_REF"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=as.integer(unlist(strsplit(tmp1,":"))[1])
  cmd=paste0("head -",tmp1," ",GSEfile," >",GSEfile,".header.txt")
  system(cmd)
  cmd=paste0("grep -n ","Sample_title"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  sample=unlist(strsplit(tmp1,"\t\""))
  sample=gsub("\"","",sample)
  sample=sample[2:length(sample)]
  sample=gsub("Genomic DNA from ","",sample)
  sample=gsub("Genomic DNA of ","",sample)
  sample=gsub("\"","",sample)
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
  #if(all(colnames(GSE)==sampleid)) colnames(GSE)=sample
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
  stage=NA
  cmd=paste0("grep -n ","Stage:"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    stage=unlist(strsplit(tmp1,"\t\""))
    stage=gsub("\"","",stage)
    stage=stage[2:length(stage)]
    stage=gsub("Stage:","",stage)
    stage=gsub("\"","",stage)
    #stage=as.integer(gsub(" ","",stage))
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
    gender[which(gender=="female")]="Female"
    gender[which(gender=="male")]="Male"
  }
  status=NA
  cmd=paste0("grep -n ","\"","disease status:","\""," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    status=unlist(strsplit(tmp1,"\t\""))
    status=gsub("\"","",status)
    status=status[2:length(status)]
    status=gsub("disease status:","",status)
    status=gsub(" ","",status)
  }
  region=NA
  cmd=paste0("grep -n ","\"","colon region:","\""," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    region=unlist(strsplit(tmp1,"\t\""))
    region=gsub("\"","",region)
    region=region[2:length(region)]
    region=gsub("colon region:","",region)
    region=gsub(" ","",region)
  }
  clinicaltable=data.frame(sample=sample,stage=stage,gender=gender,status=status,region=region,stringsAsFactors = F)
  rownames(clinicaltable)=sampleid
  return(list(GSE=GSE,clinicaltable=clinicaltable))
}

GSE48684=readGSE()

MEall=GSE48684$GSE
dim(MEall)
# 482421    147

clinicaltable=GSE48684$clinicaltable
colnames(clinicaltable)[1]="sample"
save(MEall,clinicaltable,file="../../../CancerGenomics/Colorectal_Cancer/Data/GSE48684.RData")

