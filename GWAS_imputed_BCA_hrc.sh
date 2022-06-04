#!/usr/bin/env bash
#used to do GWAS on imputed SNPs (selected by TCGA/GTEx models, model7 and model6)
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/

for chr in {1..22}
do
    $plink --bfile ${outfolder}chr${chr}_filter_hg19tohg38_flip --recode A-transpose --out ${outfolder}chr${chr}_filter_hg19tohg38_flip --memory 1200 &
done


# 
# 
# #remove multialleic
# for i in {1..22}
# do
#   $plink --bfile ${outfolder}chr${chr}_filter_hg19tohg38_flip --recode --out ${outfolder}chr${i}_select
#   $plink --file ${outfolder}chr${i}_select --bmerge ${infolder}chr${i}_select --merge-mode 6 --out duplicate
#   $plink --bfile ${infolder}chr${i}_select --exclude ./duplicate.missnp --make-bed --out ${outfolder}chr${i}_select
# done
# 
rm ${outfolder}mergelist.txt
for i in {1..22}
do
echo ${outfolder}chr${chr}_filter_hg19tohg38_flip  >>${outfolder}mergelist.txt
done
# 
$plink --merge-list ${outfolder}mergelist.txt --make-bed --out ${outfolder}select
sed -i -e 's/SEP//g' ${outfolder}select.fam
$plink --bfile ${outfolder}select --recode A-transpose --out ${outfolder}select
# 
# # models
# #BE and control
# #run perform_GWAS.R
$plink --bfile ${outfolder}select \
  --keep /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BE_CO_hrc_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}BE_CO
$plink --bfile ${outfolder}BE_CO --recode A-transpose --out ${outfolder}BE_CO
# #run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}BE_CO --covar /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BE_CO_hrc_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,pc5,pc6,sex \
  --logistic --beta --hide-covar --ci 0.95 --out ${outfolder}BE_CO
# 
# 
# #EA and control
$plink --bfile ${outfolder}select \
  --keep /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/EA_CO_hrc_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}EA_CO
# #run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}EA_CO --covar /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/EA_CO_hrc_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,pc5,pc6,sex \
  --logistic --beta --hide-covar --ci 0.95 --out ${outfolder}EA_CO
# 
# #run update_fam() in perform_GWAS.R
# $plink --bfile  ${outfolder}imp_${prefix}_EA_CO --covar  ${outfolder}imp_${prefix}_EA_selectedsamples_pheno_plink.txt \
#   --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ${outfolder}imp_${prefix}_EA_CO
# 
# #BEEA and control
$plink --bfile ${outfolder}select \
  --keep /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BEEA_CO_hrc_selectedsamples_plink.txt --make-bed \
  --out ${outfolder}BEEA_CO
# #run update_fam() in perform_GWAS.R
$plink --bfile  ${outfolder}BEEA_CO --covar  /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/BEEA_CO_hrc_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,pc5,pc6,sex --logistic --hide-covar --ci 0.95 --out ${outfolder}BEEA_CO

#get frequency
$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/select --freq --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/select
