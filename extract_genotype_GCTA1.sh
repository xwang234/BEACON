#!/usr/bin/env bash
#the code is used to apply GCTA on imputed bca data, output imp_

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result
gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.92.1beta6
gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.93.0beta
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019
cd $outfolder

#first remove cryptic relatedness
#CO vs BE/EA
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar"
$gctafolder/gcta64 --bfile $prefix --autosome --maf 0.01 --keep $filelist --make-grm --out GCTA_CO_BEEA --thread-num 15
$gctafolder/gcta64 --grm GCTA_CO_BEEA --grm-cutoff 0.025 --make-grm --out GCTA_CO_BEEA_rm025 --thread-num 15
$gctafolder/gcta64 --grm GCTA_CO_BEEA --grm-cutoff 0.05 --make-grm --out GCTA_CO_BEEA_rm05 --thread-num 15
filelist=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt
phenotype=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1
covar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar1"
qcovar="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1"
$gctafolder/gcta64 --grm GCTA_CO_BEEA_rm05 --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence 0.016 --out GCTA_CO_BEEA_covar_rm05 --thread-num 12

#read genotype
filelist_plink=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1_plink.txt
for chr in {1..22}
do
  echo $chr
  vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz
  infofile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.info.gz
  infofile1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.info
  gunzip $infofile
  $plink --vcf $vcffile --const-fid 0 --make-bed --keep ${filelist_plink} --out s1
  #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp.multi-allelic
  $plink --bfile s1 --const-fid 0 --bmerge s1 --merge-mode 6
  $plink --bfile s1 --const-fid 0 --exclude plink.missnp --make-bed --out s2
  #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. (In fact, some tools, such as BEAGLE 4, require such duplicate variants to be merged/removed during preprocessing.) --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
  $plink --bfile s2 --const-fid 0 --list-duplicate-vars
  $plink --bfile s2 --const-fid 0 --exclude plink.dupvar --make-bed --out s3
  $plink --bfile s3 --const-fid 0 --mac 2 --maf 0.001 --biallelic-only strict --snps-only --hwe 0.000001 --geno 0.05 --keep ${filelist_plink} --qual-scores $infofile1 7 1 1 --qual-threshold 0.3 --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_gcta
  Rscript ../code/update_const_fid_fam.R /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_gcta.fam
done

#for chr in {1..22}
#do
#  echo $chr
#  $plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_gcta --exclude /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta.missnp --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_gcta 
#done


$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr1_gcta --merge-list ../result/BCA_GCTA_mergelist.txt --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta

do_impute_GCTA () {
  local prefix="$1"
  local filelist="$2"
  local phenotype="$3"
  local covar="$4"
  local qcovar="$5"
  local nthread="${6?40}"
 
  echo $prefix
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score-region 200 --keep $filelist --thread-num $nthread --out $prefix"
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score-region 200 --keep $filelist --thread-num $nthread --out $prefix
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score --ld-wind 200 --keep $filelist --thread-num $nthread --out $prefix
  Rscript ../code/stratify_GCTAsnps.R $outfolder/$prefix
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group1.txt --make-grm --thread-num $nthread --out ${prefix}_group1"
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group1.txt --make-grm --thread-num $nthread --out ${prefix}_group1
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group2.txt --make-grm --thread-num $nthread --out ${prefix}_group2
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group3.txt --make-grm --thread-num $nthread --out ${prefix}_group3"
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group3.txt --make-grm --thread-num $nthread --out ${prefix}_group3
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group4.txt --make-grm --thread-num $nthread --out ${prefix}_group4
  echo -e "${prefix}_group1\n${prefix}_group2\n${prefix}_group3\n${prefix}_group4" >${prefix}_mult_GRMs.txt
  echo "$gctafolder/gcta64 --reml --mgrm ${prefix}_mult_GRMs.txt --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml-no-constrain --thread-num $nthread --out $prefix"
  $gctafolder/gcta64 --reml --mgrm ${prefix}_mult_GRMs.txt --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --thread-num $nthread --reml-alg 1 --out $prefix
}

do_impute_GCTA imp_CO_BEEA /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1 40

do_impute_GCTA imp_CO_BEEA_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar1 40
do_impute_GCTA imp_CO_BEEA_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar1 40

do_impute_GCTA imp_CO_BEEA_bmi1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_bmi1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi1.qcovar1 40
do_impute_GCTA imp_CO_BEEA_bmi0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_bmi0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.bmi0.qcovar1 40

do_impute_GCTA imp_CO_BEEA_smoke1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_smoke1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke1.qcovar1 40
do_impute_GCTA imp_CO_BEEA_smoke0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_smoke0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.smoke0.qcovar1 40 

../code/CGTA_imputedata.sh imp_CO_BEEA /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_BEEA_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_BEEA_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_BE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_BE_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_BE_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.qcovar1 0.016

../code/CGTA_imputedata.sh imp_CO_EA /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar1 0.0025

../code/CGTA_imputedata.sh imp_CO_EA_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.qcovar1 0.0025

../code/CGTA_imputedata.sh imp_CO_EA_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.qcovar1 0.0025

#use 1grm
../code/CGTA_imputedata.sh imp1grm_CO_BEEA /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_BEEA_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_BEEA_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BEEA_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_BE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_BE_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_BE_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_BE_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.qcovar1 0.016

../code/CGTA_imputedata.sh imp1grm_CO_EA /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar1 0.0025

../code/CGTA_imputedata.sh imp1grm_CO_EA_recurrent1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent1_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.qcovar1 0.0025

../code/CGTA_imputedata.sh imp1grm_CO_EA_recurrent0 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/GCTA_CO_EA_recurrent0_individuals1.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.pheno1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.covar1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.qcovar1 0.0025

#$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score-region 200 --thread-num 40 --out test
#$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --extract snp_group1.txt --make-grm --thread-num 15 --out test_group1
#$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --extract snp_group2.txt --make-grm --thread-num 40 --out test_group2
#$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --extract snp_group3.txt --make-grm --thread-num 40 --out test_group3
#$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --extract snp_group4.txt --make-grm --thread-num 40 --out test_group4
#$gctafolder/gcta64 --reml --mgrm multi_GRMs.txt --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --prevalence 0.016 --reml-alg 1 --thread-num 40 --out test
