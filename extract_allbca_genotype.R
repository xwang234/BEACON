#!/usr/bin/env Rscript
#to add allels into variant ID, to solve the problem that after lifting many alt allels are missing.
library(data.table)
outdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/"
for (chr in 1:22)
{
  bimfile=paste0(outdir,"chr",chr,"_all_noambiguous.bim")
  bimdat=as.data.frame(fread(bimfile))
  bimdat$V2=paste0(bimdat$V2,"_",bimdat$V5,"_",bimdat$V6)
  fwrite(bimdat,file=bimfile,col.names = F,row.names = F,quote=F,sep="\t")
}