#!/usr/bin/env bash
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

prefix="$1" 
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/$prefix
if [[ ! -d "$outfolder" ]]; then mkdir "$outfolder"; fi

chr=$2
start_loc=$3
end_loc=$4
echo $chr
tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz >$outfolder/dose.vcf
tabix /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz $chr:"$start_loc"-"$end_loc" >> $outfolder/dose.vcf


#read genotype
vcffile=$outfolder/dose.vcf
$plink --vcf $vcffile --const-fid 0 --snps-only --make-bed --out $outfolder/dose
  #cut -f 2 /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select.bim | sort | uniq -d > plink.dups
  #$plink --bfile $outfolder/chr${chr}_select --exclude plink.dups --make-bed --out $outfolder/chr${chr}_select
$plink --bfile $outfolder/dose --recodeA --out $outfolder/dose

