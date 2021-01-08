#/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

inputfam=args[1]
oldfam=read.table(inputfam,stringsAsFactors = F)
tmp1=tmp2=rep(NA,nrow(oldfam))
for (i in 1:nrow(oldfam))
{
  tmp=unlist(strsplit(oldfam$V2[i],"_"))
  tmp1[i]=tmp[1]
  tmp2[i]=tmp[2]
}
oldfam$V1=tmp1
oldfam$V2=tmp2
write.table(oldfam,file=inputfam,sep=" ",row.names = F,col.names = F,quote=F)
