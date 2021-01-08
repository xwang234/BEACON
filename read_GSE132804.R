library(data.table)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132804

library(data.table)

readGSE=function(GSEfile="../data/GSE132804/GSE132804-GPL21145_series_matrix.txt")
{
  cmd=paste0("grep -n ","ID_REF"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=as.integer(unlist(strsplit(tmp1,":"))[1])
  cmd=paste0("head -",tmp1," ",GSEfile," >",GSEfile,".header.txt")
  system(cmd)
  cmd=paste0("grep -n ","Sample_title"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  studyno=unlist(strsplit(tmp1,"\t\""))
  studyno=gsub("\"","",studyno)
  studyno=studyno[2:length(studyno)]
  studyno=gsub("genomic DNA from ","",studyno)
  studyno=gsub("normal colon sample ","",studyno)
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
  #if(all(colnames(GSE)==sampleid)) colnames(GSE)=studyno
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
    age=as.integer(gsub(" ","",age))
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
  source=NA
  cmd=paste0("grep -n ","source:"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    source=unlist(strsplit(tmp1,"\t\""))
    source=gsub("\"","",source)
    source=source[2:length(source)]
    source=gsub("source:","",source)
    source=gsub(" ","",source)
  }
  risk=NA
  cmd=paste0("grep -n ","risk:"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=tmp1 [length(tmp1)]
  if (!is.na(tmp1))
  {
    risk=unlist(strsplit(tmp1,"\t\""))
    risk=gsub("\"","",risk)
    risk=risk[2:length(risk)]
    risk=gsub("crc risk:","",risk)
    risk=gsub(" ","",risk)
  }
  clinicaltable=data.frame(studyno=studyno,risk=risk,age=age,gender=gender,source=source,stringsAsFactors = F)
  rownames(clinicaltable)=sampleid
  return(list(GSE=GSE,clinicaltable=clinicaltable))
}

GSE132804_GPL21145=readGSE()
GSE132804_GPL13534=readGSE(GSEfile="../data/GSE132804/GSE132804-GPL13534_series_matrix.txt")
all(rownames(GSE132804_GPL21145$GSE)==rownames(GSE132804_GPL13534$GSE))#T
sum(colnames(GSE132804_GPL21145$GSE) %in% colnames(GSE132804_GPL13534$GSE)) #0
MEall=cbind.data.frame(GSE132804_GPL21145$GSE,GSE132804_GPL13534$GSE)
dim(MEall)
# [1] 397969    334
GSE132804_GPL21145$clinicaltable$platform="EPIC"
GSE132804_GPL13534$clinicaltable$platform="HM450"
clinicaltable=rbind(GSE132804_GPL21145$clinicaltable,GSE132804_GPL13534$clinicaltable)
save(MEall,clinicaltable,file="../data/GSE132804.RData")

