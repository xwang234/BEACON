#!/usr/bin/env bash
#SBATCH -t 4-5
#SBATCH --nodes=6
#SBATCH --mem=31G
#SBATCH --mem-per-cpu=28G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=xwang234@fhcrc.org

#http://gcta.freeforums.net/thread/243/variance-covaraince-matrix-positive-definite

outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result
gctafolder=/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.92.1beta6
cd $outfolder

do_impute_GCTA () {
  local prefix="$1"
  local filelist="$2"
  local phenotype="$3"
  local covar="$4"
  local qcovar="$5"
  local prevalence="$6"
  local nthread=16
 
  echo $prefix
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score-region 200 --keep $filelist --thread-num $nthread --out $prefix"
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --ld-score-region 200 --keep $filelist --thread-num $nthread --out $prefix
  #Rscript ../code/stratify_GCTAsnps.R $outfolder/$prefix
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group1.txt --make-grm --thread-num $nthread --out ${prefix}_group1"
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group1.txt --make-grm --thread-num $nthread --out ${prefix}_group1
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group2.txt --make-grm --thread-num $nthread --out ${prefix}_group2
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group3.txt --make-grm --thread-num $nthread --out ${prefix}_group3"
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group3.txt --make-grm --thread-num $nthread --out ${prefix}_group3
  #$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --keep $filelist --extract ${prefix}_group4.txt --make-grm --thread-num $nthread --out ${prefix}_group4
  echo -e "${prefix}_group1\n${prefix}_group2\n${prefix}_group3\n${prefix}_group4" >${prefix}_mult_GRMs.txt
  echo "$gctafolder/gcta64 --reml --mgrm ${prefix}_mult_GRMs.txt --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml-no-constrain --thread-num $nthread --out $prefix"
  $gctafolder/gcta64 --reml --mgrm ${prefix}_mult_GRMs.txt --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --thread-num $nthread --prevalence $prevalence --reml-alg 1 --out $prefix 
}

do_impute_GCTA_all () {
  local prefix="$1"
  local filelist="$2"
  local phenotype="$3"
  local covar="$4"
  local qcovar="$5"
  local prevalence="$6"
  local nthread=16
 
  echo $prefix
  echo "$gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --autosome --keep $filelist --make-grm --out $prefix --thread-num $nthread"
  $gctafolder/gcta64 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/gcta --autosome --keep $filelist --make-grm --out $prefix --thread-num $nthread
  echo "$gctafolder/gcta64 --grm $prefix --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence $prevalence --out $prefix --thread-num $nthread"
  $gctafolder/gcta64 --grm $prefix --pheno $phenotype --covar $covar --qcovar $qcovar --keep $filelist --reml --prevalence $prevalence --out $prefix --thread-num $nthread
}

do_impute_GCTA_all $1 $2 $3 $4 $5 $6
