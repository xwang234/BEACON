#!/usr/bin/env bash
#works on genotyped data, use pheno1,covar1,outputs _gentoyped.hsq

#You can use this command to remove cryptic relatedness
#gcta64 --grm test --grm-cutoff 0.05 --make-grm --out test_rm05

do_CGTA_gentoyped () {
  local prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019"
  local gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.92.1beta6
  local outprefix="$1"
  local filelist="$2"
  local phenotype="$3"
  local covar="$4"
  local qcovar="$5"
  local prevalence="$6"
  local nthread=16
  $gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out $outprefix --thread-num $nthread
  #$gctafolder/gcta64 --grm $outprefix --keep $filelist --pca 20 --out $outprefix --thread-num $nthread
  $gctafolder/gcta64 --grm $outprefix --pheno $phenotype --keep $filelist --reml --prevalence $prevalence --out $outprefix --thread-num $nthread
  $gctafolder/gcta64 --grm $outprefix --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence $prevalence --out $outprefix --thread-num $nthread
  
}

#extract genotype data for GCTA
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result

cd $outfolder
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 --recode tab --snps-only --biallelic-only strict --maf 0.01 --hwe 0.05 --geno 0.01 --chr 1-22 --make-bed --out $outfolder/bca_filtered_04Apr2019
$plink --bfile $outfolder/bca_filtered_04Apr2019 --recode A-transpose --out $outfolder/bca_filtered_04Apr2019

gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.92.1beta6

prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_genotyped $filelist $phenotype $covar $qcovar 0.016


#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_genotyped $filelist $phenotype $covar $qcovar 0.0025

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_genotyped $filelist $phenotype $covar $qcovar 0.01

#CO vs BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_genotyped $filelist $phenotype $covar $qcovar 0.016

#subseting recurrent_HB_RF
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_recurrent1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_recurrent0_genotyped $filelist $phenotype $covar $qcovar 0.016

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_recurrent1_genotyped $filelist $phenotype $covar $qcovar 0.0025

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_recurrent0_genotyped $filelist $phenotype $covar $qcovar 0.0025

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_recurrent1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_recurrent1_genotyped $filelist $phenotype $covar $qcovar 0.01

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_recurrent0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_recurrent0_genotyped $filelist $phenotype $covar $qcovar 0.01

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_recurrent1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_recurrent0_genotyped $filelist $phenotype $covar $qcovar 0.016


#subsetting smoking
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_smoke1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_smoke1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_smoke0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.smoke0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_smoke0_genotyped $filelist $phenotype $covar $qcovar 0.016

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_smoke1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke1.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_smoke1_genotyped $filelist $phenotype $covar $qcovar 0.0025

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_smoke0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.smoke0.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_smoke0_genotyped $filelist $phenotype $covar $qcovar 0.0025

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_smoke1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke1.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_smoke1_genotyped $filelist $phenotype $covar $qcovar 0.01

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_smoke0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.smoke0.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_smoke0_genotyped $filelist $phenotype $covar $qcovar 0.01

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_smoke1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_smoke1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_smoke0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_smoke0_genotyped $filelist $phenotype $covar $qcovar 0.016




#subsetting bmi
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_bmi1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_bmi1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_bmi0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.bmi0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BE_bmi0_genotyped $filelist $phenotype $covar $qcovar 0.016

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_bmi1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi1.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_bmi1_genotyped $filelist $phenotype $covar $qcovar 0.0025

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_bmi0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.bmi0.qcovar1"
do_CGTA_gentoyped GCTA_CO_EA_bmi0_genotyped $filelist $phenotype $covar $qcovar 0.0025

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_bmi1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi1.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_bmi1_genotyped $filelist $phenotype $covar $qcovar 0.01

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_bmi0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.bmi0.qcovar1"
do_CGTA_gentoyped GCTA_BE_EA_bmi0_genotyped $filelist $phenotype $covar $qcovar 0.01

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_bmi1_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_bmi1_genotyped $filelist $phenotype $covar $qcovar 0.016

filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_bmi0_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.qcovar1"
do_CGTA_gentoyped GCTA_CO_BEEA_bmi0_genotyped $filelist $phenotype $covar $qcovar 0.016



#Use recurrent as covariate
#CO VS BE
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar1"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --keep $filelist --pca 20 --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BE_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BE_covar1_recurrent --thread-num 12

#CO VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar1"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --keep $filelist --pca 20 --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_EA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.0025 --out GCTA_CO_EA_covar1_recurrent --thread-num 12

#BE VS EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_BE_EA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.BE.EA.qcovar1"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --keep $filelist --pca 20 --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_BE_EA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.01 --out GCTA_BE_EA_covar1_recurrent --thread-num 12

#CO VS BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --keep $filelist --pca 20 --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --pheno $phenotype --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_recurrent --thread-num 12
$gctafolder/gcta64 --grm GCTA_CO_BEEA_recurrent --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar1_recurrent --thread-num 12




