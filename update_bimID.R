#/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

inputbim=args[1]
oldbim=read.table(inputbim,stringsAsFactors = F)
oldbim$V2=paste0(oldbim$V1,":",oldbim$V4,"_",oldbim$V5,"_",oldbim$V6)

write.table(oldbim,file=inputbim,sep=" ",row.names = F,col.names = F,quote=F)
