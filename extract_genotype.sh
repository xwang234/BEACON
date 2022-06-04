#!/usr/bin/env bash

#extract genotype data on 5snps of CAMBRIDGE samples
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
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

#extract genotype data on 26snps based on new data (with extra controls)
ml tabix/0.2.6-GCCcore-8.3.0
#Beacon
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g
snps=($(awk '{print $1}' ../result/Dong26SNPs_postion.txt)) #Dong26SNPs_postion.txt has the same order as dong26snp (R data)
chrs=($(awk '{print $2}' ../result/Dong26SNPs_postion.txt))
for chr in {1..22}
do
  echo $chr
  prefix=$outdir/chr${chr}.dose
  tabix -p vcf -f $prefix.vcf.gz &
done

outfile=$outdir/Dong26SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
for ((i=0;i<${#snps[@]};i++))
do
  tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
done
#Cambridge
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g
for chr in {1..22}
do
  echo $chr
  prefix=$outdir/chr${chr}.dose
  tabix -p vcf -f $prefix.vcf.gz &
done

outfile=$outdir/Dong26SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
for ((i=0;i<${#snps[@]};i++))
do
  tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
done
#add "rs9918259"  "rs75783973" imputed from 1000 genome phase 1 (on chr5)
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000gphase1_maf001_snp
snps=($(awk '{print $1}' ../result/Dong26SNPs_postion.txt))
chrs=($(awk '{print $2}' ../result/Dong26SNPs_postion.txt))
chr=5
prefix=$outdir/chr${chr}.dose
tabix -p vcf -f $prefix.vcf.gz 

outfile=$outdir/Dong26SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
i=5
tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
i=6
tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile

#Cambridge
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000gphase1_maf001_snp
chr=5
prefix=$outdir/chr${chr}.dose
tabix -p vcf -f $prefix.vcf.gz 
outfile=$outdir/Dong26SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
i=5
tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
i=6
tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile

#extract dong26 from HRC
ml tabix/0.2.6-GCCcore-8.3.0
#Beacon
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc
snps=($(awk '{print $1}' ../result/Dong26SNPs_postion.txt)) #Dong26SNPs_postion.txt has the same order as dong26snp (R data)
chrs=($(awk '{print $2}' ../result/Dong26SNPs_postion.txt))
for chr in {1..22}
do
  echo $chr
  prefix=$outdir/chr${chr}.dose
  tabix -p vcf -f $prefix.vcf.gz &
done
wait
outfile=$outdir/Dong26SNPs_genotype.txt
#if [[ -f $outfile ]]; then rm $outfile; fi
for ((i=0;i<${#snps[@]};i++))
do
  #tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
  tmp=$(tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]})
  if [[ "$tmp" = "" ]]; then echo "$i+1";fi #check which snps are missing
done

#hrc_5missing=c(6,7,10,19,21)
# dong26snp[hrc_5missing,]
#            SNP Chr Position    pos38
# 6    rs9918259   5   663092   662977
# 7   rs75783973   5   668309   668194
# 10  rs76014404   6 62391538 61681633
# 19  rs66725070  15 58267416 57975218
# 21 rs199620551  19 18804294 18693484
#Cambridge
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_hrc
for chr in {1..22}
do
  echo $chr
  prefix=$outdir/chr${chr}.dose
  tabix -p vcf -f $prefix.vcf.gz &
done
wait
outfile=$outdir/Dong26SNPs_genotype.txt
if [[ -f $outfile ]]; then rm $outfile; fi
for ((i=0;i<${#snps[@]};i++))
do
  tabix $outdir/chr${chrs[$i]}.dose.vcf.gz ${snps[$i]} |tail -1 >>$outfile
done
