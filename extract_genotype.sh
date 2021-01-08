#!/usr/bin/env bash

#extract genotype data on 5snps of CAMBRIDGE samples
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/result

snps="rs11765529, rs7798662, rs4930068, rs2687201, rs10419226"

$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 --snps $snps --keep $outdir/cambridge_samples.txt --recode --tab --out $outdir/cambridge_genotype

$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 --snps $snps --recode --tab --out $outdir/cambridgeall_genotype

#extract genotype data on 23snps
snps=($(awk '{print $1}' ../result/Dong23SNPs_postion.txt))
chrs=($(awk '{print $2}' ../result/Dong23SNPs_postion.txt))
outfile=$outdir/Dong23SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
for ((i=0;i<${#snps[@]};i++))
do
  tabix /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_snp23_phase1/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
done


for chr in {1..22}
do
  echo $chr
  tabixreg=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/prediction_snps_tabix_chr${chr}.txt
  tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz >/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.select.dose.vcf
  xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz {} >> /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.select.dose.vcf
done

#read genotype
for chr in {1..22}
do
  echo $chr
  vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.select.dose.vcf
  $plink --vcf $vcffile --const-fid 0 --snps-only --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select --memory 64
  #cut -f 2 /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select.bim | sort | uniq -d > plink.dups
  #$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select --exclude plink.dups --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select
  $plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select --recodeA --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select --memory 640
done

snps="19:18817903, 19:18804294, 19:18803172"
$plink --vcf /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.dose.vcf.gz --const-fid 0 --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.dose
$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.dose -list-duplicate-vars 

$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.dose -exclude ./plink.dupvar --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.duprm.dose
$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.duprm.dose --snps $snps --recode A --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.3knownsnps

