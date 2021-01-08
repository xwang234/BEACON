#!/usr/bin/env bash
#used to do GWAS on imputed SNPs (selected by TCGA/GTEx models, model7 and model6)
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

prefix="TCGA"
prefix="GTEx"
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_${prefix}model/
infolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_${prefix}_PC4/
mkdir $outfolder

#remove multialleic
for i in {1..22}
do
  $plink --bfile ${infolder}chr${i}_select --recode --out ${outfolder}chr${i}_select
  $plink --file ${outfolder}chr${i}_select --bmerge ${infolder}chr${i}_select --merge-mode 6 --out duplicate
  $plink --bfile ${infolder}chr${i}_select --exclude ./duplicate.missnp --make-bed --out ${outfolder}chr${i}_select
done

rm ${outfolder}mergelist.txt
for i in {1..22}
do
echo ${outfolder}chr${i}_select  >>${outfolder}mergelist.txt
done

$plink --merge-list ${outfolder}mergelist.txt --make-bed --out ${outfolder}select


#for validate_gwas---
# models
#BE and control
#run perform_GWAS.R
$plink --bfile ${outfolder}select \
  --keep ${outfolder}imp_${prefix}_BE_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}imp_${prefix}_BE_CO
#run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}imp_${prefix}_BE_CO --covar  ${outfolder}imp_${prefix}_BE_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ${outfolder}imp_${prefix}_BE_CO

#EA and control
$plink --bfile ${outfolder}select \
  --keep ${outfolder}imp_${prefix}_EA_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}imp_${prefix}_EA_CO
#run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}imp_${prefix}_EA_CO --covar  ${outfolder}imp_${prefix}_EA_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ${outfolder}imp_${prefix}_EA_CO

#BEEA and control
$plink --bfile ${outfolder}select \
  --keep ${outfolder}imp_${prefix}_BEEA_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}imp_${prefix}_BEEA_CO
#run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}imp_${prefix}_BEEA_CO --covar  ${outfolder}imp_${prefix}_BEEA_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ${outfolder}imp_${prefix}_BEEA_CO
