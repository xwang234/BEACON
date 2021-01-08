#!/usr/bin/env Rscript

source("/fh/fast/stanford_j/Xiaoyu/QTL/code/functions.R")
getrsid=function(oldid=mqtl_23highrisk_cis$SNP)
{
  res=sapply(oldid,function(x){
    unlist(strsplit(x,":"))[1]
  })
  return(res)
}
#load("../result/Dong23SNPs.RData")
#500KB distance cutoff, based on 185 PRAD samples
mqtl_23highrisk_cis=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_23highrisk_cn_mutation_cis"))
qqplot(mqtl_23highrisk_cis$`p-value`)
idx=mqtl_23highrisk_cis$qvalue<0.05
unique(mqtl_23highrisk_cis$SNP[idx]) #"rs9918259"   "rs75783973"   "rs9257809" "rs10108511" "rs11789015" "rs199620551" "rs10423674"
unique(mqtl_23highrisk_cis$gene[idx]) #"cg20362242" "cg05482498" "cg13177375" "cg13835168" "cg21775007" "cg07076509" "cg12034943" "cg26584619"
idx=dong23snp$SNP %in% unique(getrsid(mqtl_23highrisk_cis$SNP[idx]))
View(dong23snp[idx,])

# mqtl_23highrisk_cis=readqtlres("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_23highrisk_cn_mutation_cis",
#                                snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_23highrisk_SNP_POS.txt",
#                                geposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_ME_rmprobes_POS.txt",fdrcutoff = 1)

#based on 89 EAC samples
mqtl_eac_23highrisk_cis=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_23highrisk_cn_mutation_cis"))
#include 21 SNPs, 8th and 10th were missing based on the cutoff 500KB
qqplot(mqtl_eac_23highrisk_cis$`p-value`)
idx=mqtl_eac_23highrisk_cis$qvalue<0.05
unique(mqtl_eac_23highrisk_cis$SNP[idx]) #"rs9257809"
unique(mqtl_eac_23highrisk_cis$gene[idx]) #"cg13177375"

mqtl_eac_23highrisk_cis=readqtlres("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_23highrisk_cn_mutation_cis",
                               snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_23highrisk_SNP_POS.txt",
                               geposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_ME_rmprobes_POS.txt",fdrcutoff = 1)

#use 250KB
mqtl_eac_23highrisk_cis1=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_23highrisk_cn_mutation_cis_250K"))
qqplot(mqtl_eac_23highrisk_cis1$`p-value`) 
idx=mqtl_eac_23highrisk_cis1$qvalue<0.05
unique(mqtl_eac_23highrisk_cis1$SNP[idx])
unique(mqtl_eac_23highrisk_cis1$gene[idx])

eqtl_eac_23highrisk_cis=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_EAC_highrisk_cn_mutation_cis"))
eqtl_eac_23highrisk_cis$SNP=getrsid(eqtl_eac_23highrisk_cis$SNP)
unique(eqtl_eac_23highrisk_cis$SNP[eqtl_eac_23highrisk_cis$qvalue<0.05]) #"rs9257809:29356331:A:G"
idx=which(eqtl_eac_23highrisk_cis$qvalue<0.05)
print(paste0(eqtl_eac_23highrisk_cis$SNP[idx],"-",eqtl_eac_23highrisk_cis$gene[idx],"-",eqtl_eac_23highrisk_cis$`p-value`[idx]))
View(eqtl_eac_23highrisk_cis[idx,])
dim(eqtl_eac_23highrisk_cis)
qqplot(eqtl_eac_23highrisk_cis$`p-value`)
length(unique(eqtl_eac_23highrisk_cis$SNP[eqtl_eac_23highrisk_cis$`p-value`<0.05]))

proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
eqtl_gtex_23highrisk_cis=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_gtex_23highrisk_cn_mutation_cis"))
eqtl_gtex_23highrisk_cis$SNP=getrsid(eqtl_gtex_23highrisk_cis$SNP)
eqtl_gtex_23highrisk_cis=eqtl_gtex_23highrisk_cis[eqtl_gtex_23highrisk_cis$gene %in% proteingenes,]
length(unique(eqtl_gtex_23highrisk_cis$SNP)) #23
eqtl_gtex_23highrisk_cis=eqtl_gtex_23highrisk_cis[order(eqtl_gtex_23highrisk_cis$`p-value`),]
unique(eqtl_gtex_23highrisk_cis$SNP[eqtl_gtex_23highrisk_cis$qvalue<0.05])
# "rs3072"     "rs7255"     "rs9257809"  "rs10108511"
idx=which(eqtl_gtex_23highrisk_cis$qvalue<0.05)
View(eqtl_gtex_23highrisk_cis[idx,])
print(paste0(eqtl_gtex_23highrisk_cis$SNP[idx],"-",eqtl_gtex_23highrisk_cis$gene[idx],"-",eqtl_gtex_23highrisk_cis$`p-value`[idx]))
dim(eqtl_gtex_23highrisk_cis)
qqplot(eqtl_gtex_23highrisk_cis$`p-value`)
length(unique(eqtl_gtex_23highrisk_cis$SNP[eqtl_gtex_23highrisk_cis$qvalue<0.05]))
length(unique(eqtl_gtex_23highrisk_cis$SNP[eqtl_gtex_23highrisk_cis$`p-value`<0.05])) #15

