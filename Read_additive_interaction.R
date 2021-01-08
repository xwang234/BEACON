
rm(list=ls())

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/res_additive_test_genomewide.RData")
mpoint <- ceiling(quantile(1:nrow(EAres),0.95)/200)
pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(EAres))

plot(-log((nrow(BEres):1)/nrow(BEres),base=10)[pindex],-log(BEres[order(BEres[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BEres):1)/nrow(BEres),base=10)[pindex],-log(BEres[order(BEres[,2],decreasing=T),2],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BEres):1)/nrow(BEres),base=10)[pindex],-log(BEres[order(BEres[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(EAres):1)/nrow(EAres),base=10)[pindex],-log(EAres[order(EAres[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(EAres):1)/nrow(EAres),base=10)[pindex],-log(EAres[order(EAres[,2],decreasing=T),2],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(EAres):1)/nrow(EAres),base=10)[pindex],-log(EAres[order(EAres[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BEEAres):1)/nrow(BEEAres),base=10)[pindex],-log(BEEAres[order(BEEAres[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BEEAres):1)/nrow(BEEAres),base=10)[pindex],-log(BEEAres[order(BEEAres[,2],decreasing=T),2],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BEEAres):1)/nrow(BEEAres),base=10)[pindex],-log(BEEAres[order(BEEAres[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)







library(data.table)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output_pc10.RData")
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_GERD_interaction_output.RData")
#load(file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output.RData") 
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_cig_interaction_output.RData")
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_bmi_interaction_output.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_junction_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_stomach_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_blood_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_mucosa_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_muscularis_V7_eQTL_gene_pairs.RData")


lookuptable=fread("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")

#jointeqtl<- jointeqtl[jointeqtl$pval_nominal < 2e-11,]

lookuptable1 <- lookuptable[lookuptable$variant_id %in% jointeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% bloodeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% stomacheqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% mucosaeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% musculariseqtl$variant_id,]

idx <- match(musculariseqtl$variant_id,lookuptable1$variant_id)
idx <- match(mucosaeqtl$variant_id,lookuptable1$variant_id)
idx <- match(stomacheqtl$variant_id,lookuptable1$variant_id)
idx <- match(bloodeqtl$variant_id,lookuptable1$variant_id)
idx <- match(jointeqtl$variant_id,lookuptable1$variant_id)


lookuptable1 <- lookuptable1[idx,]

BCAout1 <- BEres[bim$V2 %in%lookuptable1$rs_id_dbSNP147_GRCh37p13, ]
BCAout1 <- EAres[bim$V2 %in%lookuptable1$rs_id_dbSNP147_GRCh37p13, ]
bim1 <- bim[bim$V2 %in%lookuptable1$rs_id_dbSNP147_GRCh37p13, ]


mpoint <- ceiling(quantile(1:nrow(BCAout1),0.95)/200)
pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(BCAout1))


plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,2],decreasing=T),2],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)
