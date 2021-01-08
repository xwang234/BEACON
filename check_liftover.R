
library(data.table)
chr=1
outdir="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/"
#hg19=as.data.frame(fread(paste0(outdir,"chr22_all_noambiguous.map")))
hg38=as.data.frame(fread(paste0(outdir,"chr",chr,"_all_noambiguous_hg19tohg38_flip.bim")))
hg38$name=paste0(hg38$V1,":",hg38$V4,"_",hg38$V5,"_",hg38$V6)
#all(hg38$V2 %in% hg19$V2)
#traw38=as.data.frame(fread(paste0(outdir,"chr22_all_noambiguous_hg19tohg38_flip.traw")))

prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_noambiguous_chr"

gtexhg38=paste0(prefix,chr,".bim")
gtexhg38=as.data.frame(fread(gtexhg38))
gtexhg38$name=paste0(gtexhg38$V1,":",gtexhg38$V4,"_",gtexhg38$V5,"_",gtexhg38$V6)
trawgtex38=as.data.frame(fread(paste0(prefix,chr,".traw")))
rownames(trawgtex38)=paste0(trawgtex38$CHR,":",trawgtex38$POS,"_",trawgtex38$COUNTED,"_",trawgtex38$ALT)
trawgtex38=trawgtex38[,7:ncol(trawgtex38)]
tmp=unlist(strsplit(colnames(trawgtex38),"_"))
colnames(trawgtex38)=tmp[seq(1,length(tmp),2)]
#junction gtex
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")
idx=which(colnames(trawgtex38) %in% colnames(phenotype))
trawgtex38=trawgtex38[,idx]
#check missing
tmp=rowSums(trawgtex38)
sum(is.na(tmp))/length(tmp) #0.077
sum(is.na(trawgtex38))/nrow(trawgtex38)/ncol(trawgtex38) #0.0012
#check overlap
str11=paste0(hg38$V1,":",hg38$V4,"_",hg38$V5,"_",hg38$V6)
str12=paste0(hg38$V1,":",hg38$V4,"_",hg38$V6,"_",hg38$V5)
str2=rownames(trawgtex38)
sum(str2 %in% str11) #725896
sum(str2 %in% str12) #8953
sum(str2 %in% c(str11,str12))/length(str11) #0.921
sum(str2 %in% c(str11,str12))/length(str2) #0.258
idx=which(str2 %in% c(str11,str12))
str21=paste0(gtexhg38$V1[idx],":",gtexhg38$V4[idx],"_",gtexhg38$V5[idx],"_",gtexhg38$V6[idx])
str22=paste0(gtexhg38$V1[idx],":",gtexhg38$V4[idx],"_",gtexhg38$V6[idx],"_",gtexhg38$V5[idx])

BonnBE=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt"))
CambridgeBE=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Cambridge_autosomes.txt"))
CambridgeBE=CambridgeBE[CambridgeBE$CHR==chr,]
tmp=paste0(CambridgeBE$CHR,":",CambridgeBE$position,"_",CambridgeBE$non_effect_allele,"_",CambridgeBE$effect_allele)
sum(tmp %in% str21)
sum(tmp %in% str22)
