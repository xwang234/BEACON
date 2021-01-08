#/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}
library(data.table)
inputinfo=args[1]
oldinfo=as.data.frame(fread(inputinfo,stringsAsFactors = F))
str11=paste0(oldinfo$SNP,"_",oldinfo$`REF(0)`,"_",oldinfo$`ALT(1)`)
oldinfo$SNP=paste0(oldinfo$SNP,"_",oldinfo$`ALT(1)`,"_",oldinfo$`REF(0)`)
inputbim=args[2]
oldbim=read.table(inputbim,stringsAsFactors = F)
oldbim$V2=paste0(oldbim$V1,":",oldbim$V4,"_",oldbim$V5,"_",oldbim$V6)
tmp=intersect(oldbim$V2, str11)
if (length(tmp)>0)
{
  idx=match(tmp,str11)
  oldinfo$SNP[idx]=tmp
}
if (any(!oldbim$V2 %in% oldinfo$SNP))
{
  warnings(paste0("some variants not found in ",inputinfo))
}
write.table(oldinfo,file=paste0(inputinfo,".newID"),sep=" ",row.names = F,col.names = T,quote=F)
