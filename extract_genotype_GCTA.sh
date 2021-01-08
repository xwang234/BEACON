#!/usr/bin/env bash
#the code generates results based on genotyped data, use pheno,covar, outputs: _covar.hsq
#You can use this command to remove cryptic relatedness
#gcta64 --grm test --grm-cutoff 0.05 --make-grm --out test_rm05

#extract genotype data for GCTA
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result

cd $outfolder
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 --recode tab --snps-only --biallelic-only strict --maf 0.01 --hwe 0.05 --geno 0.01 --chr 1-22 --make-bed --out $outfolder/bca_filtered_04Apr2019

$plink --bfile $outfolder/bca_filtered_04Apr2019 --recode A-transpose --out $outfolder/bca_filtered_04Apr2019

gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.92.1beta6

prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BE --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE --keep $filelist --pca 20 --out GCTA_CO_BE --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_covar --thread-num 12

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA --keep $filelist --pca 20 --out GCTA_CO_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA --pheno $phenotype --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_covar --thread-num 12

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_BE_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA --keep $filelist --pca 20 --out GCTA_BE_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA --pheno $phenotype --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_covar --thread-num 12

#CO vs BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA --keep $filelist --pca 20 --out GCTA_CO_BEEA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar --thread-num 12


#subseting recurrent_HB_RF
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent1_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BE_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent1 --keep $filelist --pca 20 --out GCTA_CO_BE_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent1 --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent1 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_covar_recurrent1 --thread-num 12
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent0_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BE_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent0 --keep $filelist --pca 20 --out GCTA_CO_BE_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent0 --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent0 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_covar_recurrent0 --thread-num 12

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent1_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent1 --keep $filelist --pca 20 --out GCTA_CO_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent1 --pheno $phenotype --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent1 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_covar_recurrent1 --thread-num 12
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent0_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent0 --keep $filelist --pca 20 --out GCTA_CO_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent0 --pheno $phenotype --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent0 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_covar_recurrent0 --thread-num 12

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_recurrent1_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_BE_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent1 --keep $filelist --pca 20 --out GCTA_BE_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent1 --pheno $phenotype --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent1 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_covar_recurrent1 --thread-num 12
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_recurrent0_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_BE_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent0 --keep $filelist --pca 20 --out GCTA_BE_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent0 --pheno $phenotype --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent0 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_covar_recurrent0 --thread-num 12

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent1_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent1 --keep $filelist --pca 20 --out GCTA_CO_BEEA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent1 --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_recurrent1 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent1 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar_recurrent1 --thread-num 12
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent0_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent0 --keep $filelist --pca 20 --out GCTA_CO_BEEA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent0 --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_recurrent0 --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent0 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar_recurrent0 --thread-num 12

#Use recurrent as covariate
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --keep $filelist --pca 20 --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_covar_recurrent --thread-num 12

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --keep $filelist --pca 20 --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_covar_recurrent --thread-num 12

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --keep $filelist --pca 20 --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_covar_recurrent --thread-num 12

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --keep $filelist --pca 20 --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar_recurrent --thread-num 12




