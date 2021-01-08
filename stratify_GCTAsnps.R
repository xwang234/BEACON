#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
prefix=args[1]
library(data.table)
#stratify the SNPs by LD scores of individual SNPs in R
#https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata
lds_seg = as.data.frame(fread(paste0(prefix,".score.ld"),header=T))
if (sum(colnames(lds_seg)=="ldscore_region")==0) lds_seg$ldscore_region=lds_seg$ldscore
quartiles=summary(lds_seg$ldscore_region)

lb1 = which(lds_seg$ldscore_region <= quartiles[2])
lb2 = which(lds_seg$ldscore_region > quartiles[2] & lds_seg$ldscore_region <= quartiles[3])
lb3 = which(lds_seg$ldscore_region > quartiles[3] & lds_seg$ldscore_region <= quartiles[5])
lb4 = which(lds_seg$ldscore_region > quartiles[5])

lb1_snp = data.frame(snp=lds_seg$SNP[lb1])
lb2_snp = data.frame(snp=lds_seg$SNP[lb2])
lb3_snp = data.frame(snp=lds_seg$SNP[lb3])
lb4_snp = data.frame(snp=lds_seg$SNP[lb4])

fwrite(lb1_snp, paste0(prefix,"_group1.txt"), row.names=F, quote=F, col.names=F)
fwrite(lb2_snp, paste0(prefix,"_group2.txt"), row.names=F, quote=F, col.names=F)
fwrite(lb3_snp, paste0(prefix,"_group3.txt"), row.names=F, quote=F, col.names=F)
fwrite(lb4_snp, paste0(prefix,"_group4.txt"), row.names=F, quote=F, col.names=F)