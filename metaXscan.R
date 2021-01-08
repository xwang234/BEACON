#!/usr/bin/env Rscript

#check example inputs:
#read the model
library("RSQLite")
sqlfile="/fh/fast/dai_j/CancerGenomics/Tools/MetaXcan/software/data/DGN-WB_0.5.db"
con <- dbConnect(drv=RSQLite::SQLite(), dbname=sqlfile)

## list all tables
tables <- dbListTables(con)

## exclude sqlite_sequence (contains table information)
tables <- tables[tables != "sqlite_sequence"]

lDataFrames <- vector("list", length=length(tables))

## create a data.frame for each table
for (i in seq(along=tables)) {
  lDataFrames[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
}
gwasmodel=lDataFrames[[4]]
#read GWAS result:
gwasfolder="/fh/fast/dai_j/CancerGenomics/Tools/MetaXcan/software/data/GWAS/"
gwasresult=NULL
for (i in 1:22)
{
  cat(i,'..')
  tmp=read.table(gzfile(paste0(gwasfolder,"chr",i,".assoc.dosage.gz")),header = T,stringsAsFactors = F)
  tmp$chr=i
  gwasresult=rbind(gwasresult,tmp)
}

table(gwasmodel$rsid %in% gwasresult$SNP)
# FALSE   TRUE 
# 8 331417 
gwasmodel[!gwasmodel$rsid %in% gwasresult$SNP,]
table(gwasresult$SNP %in% gwasmodel$rsid)
# FALSE    TRUE 
# 3212790  249691 

#check reference covariance matrix
library(data.table)
covfile="/fh/fast/dai_j/CancerGenomics/Tools/MetaXcan/software/data/covariance.DGN-WB_0.5.txt.gz"
gwascov=fread(paste0("gunzip ",covfile))


