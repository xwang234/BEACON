#!/usr/bin/env bash
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

prefix="$1" 
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/$prefix
if [[ ! -d "$outfolder" ]]; then mkdir "$outfolder"; fi

extract_vcf(){
  local chr="$1"
  echo $chr
  tabixreg=$outfolder/prediction_snps_tabix_chr${chr}.txt
  tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g_shapeit/chr${chr}.dose.vcf.gz >$outfolder/chr${chr}.select.dose.vcf
  xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g_shapeit/chr${chr}.dose.vcf.gz {} >> $outfolder/chr${chr}.select.dose.vcf
}

for chr in {1..22}; do extract_vcf $chr & done
wait
#read genotype
for chr in {1..22}
do
  echo $chr
  vcffile=$outfolder/chr${chr}.select.dose.vcf
  $plink --vcf $vcffile --const-fid 0 --snps-only --make-bed --out $outfolder/chr${chr}_select
  #cut -f 2 /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g_shapeit/chr${chr}_select.bim | sort | uniq -d > plink.dups
  #$plink --bfile $outfolder/chr${chr}_select --exclude plink.dups --make-bed --out $outfolder/chr${chr}_select
  $plink --bfile $outfolder/chr${chr}_select --recodeA --out $outfolder/chr${chr}_select
done

