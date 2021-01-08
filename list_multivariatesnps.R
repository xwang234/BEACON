#/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

inputbim=args[1]
outfile=args[2]
if (!grepl(".gz",inputbim))
{
  oldbim=read.table(inputbim,stringsAsFactors = F)
  idx=duplicated(oldbim$V2)
  excludelist=data.frame(snp=unique(oldbim$V2[idx]))
}else
{
  oldinfo=read.table(gzfile(inputbim),header = T,stringsAsFactors = F)
  idx=duplicated(oldinfo[,1])
  excludelist=data.frame(snp=unique(oldinfo[idx,1]))
}


write.table(excludelist,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
