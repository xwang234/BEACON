
#############################################
### Checking summary stat data from Puya  ###
#############################################

rm(list=ls())
summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"
library(data.table)
summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_autosomes.txt")
summarydat=as.data.frame(fread(summaryfile,header=T))
summarydat1=as.data.frame(fread(summaryfile1,header=T))
summaryfile2=paste0(summaryfolder,"EA_Bonn_autosomes.txt")
summaryfile3=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt")
summarydat2=as.data.frame(fread(summaryfile2,header=T))
summarydat3=as.data.frame(fread(summaryfile3,header=T))

summaryfile4=paste0(summaryfolder,"BE_Cambridge_autosomes.txt")
summarydat4=as.data.frame(fread(summaryfile4,header=T))


summaryfile5=paste0(summaryfolder,"EA_Cambridge_autosomes.txt")
summarydat5=as.data.frame(fread(summaryfile5,header=T))


#725330 19 rs10419226  0 18803172  A  B
#67473   2 kgp6039878  0   651290  B  A

tsnp1 <- "rs2903490"
tsnp1 <- "rs10419226"
tsnp1 <- "rs12420406"
sum(summarydat$SNP==tsnp1)
summarydat[summarydat$SNP==tsnp1,]
sum(summarydat1$rsid==tsnp1)
summarydat1[summarydat1$rsid==tsnp1,]
sum(summarydat2$SNP==tsnp1)
summarydat2[summarydat2$SNP==tsnp1,]
sum(summarydat3$SNP==tsnp1)
summarydat3[summarydat3$SNP==tsnp1,]
sum(summarydat4$SNP==tsnp1)
summarydat4[summarydat4$SNP==tsnp1,]
sum(summarydat5$SNP==tsnp1)
summarydat5[summarydat5$SNP==tsnp1,]


pp <- summarydat$P
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))

png("/fh/fast/dai_j/BEACON/Meta_summary_stat/qq_Bonn_BE.png")
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
title("Bonn: BE (n=1037) vs Control (n=3537)")
abline(0,1)
dev.off()

pp <- summarydat2$P
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
abline(0,1)

pp <- summarydat3$P
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
abline(0,1)

png("/fh/fast/dai_j/BEACON/Meta_summary_stat/qq_Oxford_BE.png")
pp <- summarydat1$pvalue
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
abline(0,1)
title("Oxford: BE (n=1851) vs Control (n=3496)")
dev.off()



png("/fh/fast/dai_j/BEACON/Meta_summary_stat/qq_Cambridge_BE.png")
pp <- summarydat4$P
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
abline(0,1)
title("Cambridge: BE (n=873) vs Control (n=3408)")
dev.off()


png("/fh/fast/dai_j/BEACON/Meta_summary_stat/qq_BEACON_BE.png")
pp <- summarydat5$P
nsnps <- length(pp)
mpoint <- ceiling(quantile(1:length(pp),0.995)/200)*200
pindex <- c(1,(1:ceiling(quantile(1:length(pp),0.995)/200))*200,(mpoint+1):length(pp))
plot((-log((nsnps:1)/nsnps,base=10))[pindex],(-log(pp[order(pp,decreasing=T)],base=10))[pindex],xlab="-log base 10 (expected null p-values)", ylab="-log base 10 (observed p-values)")
abline(0,1)
title("BEACON: BE (n=2413) vs Control (n=2185)")
dev.off()



