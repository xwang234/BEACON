
### analyze BCA GWAS output ###

rm(list=ls())
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output_pc10.RData")

mpoint <- ceiling(quantile(1:nrow(BCAout),0.95)/200)
pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(BCAout))

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,5],decreasing=T),5],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,7],decreasing=T),7],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)



plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,9],decreasing=T),9],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

snpfile <- read.table("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filteredMAF_20Feb2015.bim",header=F)

idx=snpfile$V1 %in% 1:23 #23:X
snpfile=snpfile[idx,]

snpfile[order(BCAout[,1]),][1:5,]
snpfile$V2[order(BCAout[,1])][1:5]
BCAout[,1][order(BCAout[,1])][1:5]

snpfile[order(BCAout[,3]),][1:5,]
snpfile$V2[order(BCAout[,3])][1:5]
BCAout[,3][order(BCAout[,3])][1:5]

snpfile[order(BCAout[,5]),][1:5,]
snpfile$V2[order(BCAout[,5])][1:5]
BCAout[,5][order(BCAout[,5])][1:5]


snpfile[order(BCAout[,7]),][1:5,]
BCAout[,7][order(BCAout[,7])][1:5]



nskip=16
eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca"
eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind"
#                         thesamples=geneexpsamplenames,nskip=16,opt=1)
# {
eigen <- read.table(eigfile,skip=nskip,stringsAsFactors = F)
eigsamples=read.table(eigsampfile,stringsAsFactors = F)
row.names(eigen) <- eigsamples$V2

sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
allsamples=intersect(row.names(eigen),sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
mean(sampletable$localid==row.names(eigen))
