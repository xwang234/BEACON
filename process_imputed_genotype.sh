#!/usr/bin/env bash
#modified from process_bca_genotype.sh, to process imputed snps:cambridge_control and beacon_control
#first remove ambiguous SNPs, than hg19-hg38


#extract genotype data on all BCA samples,used to find common snps
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
#to use snpflip
ml Python
source /fh/fast/dai_j/CancerGenomics/Tools/snpflip/env/bin/activate
snpflip=/fh/fast/dai_j/CancerGenomics/Tools/snpflip/env/bin/snpflip
reference=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/reference/human_g1k_v37.fasta
hg38=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/hg38/hg38.fa



#run the following for each dataset

dataset="cambridgewtccc_1000g"
dataset="beacondbgapcontrol_1000g"
dataset="cambridgewtccc_qc_1000g"
dataset="beacondbgapcontrol_qc_1000g"

# #check missing
# checkmissing(){
#   local chr="$1"
#   echo $chr
#   vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}.dose.vcf.gz
#   $plink --vcf $vcffile  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 2560
# }
#filtering
filtering(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter
  vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}.dose.vcf.gz
  infofile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}.info.gz
  infofile1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}.info
  if [[ ! -f $infofile1 ]];then
     gunzip $infofile
  fi
  $plink --vcf $vcffile  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 2560
  Rscript /fh/fast/dai_j/BEACON/BEACON_GRANT/code/update_bimID.R tmp_s1_$chr.bim
  Rscript /fh/fast/dai_j/BEACON/BEACON_GRANT/code/update_infoID1.R $infofile1 tmp_s1_$chr.bim
  #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp.multi-allelic
  $plink --bfile tmp_s1_$chr  --bmerge tmp_s1_$chr --merge-mode 6 --out tmp_$chr   --memory 2560
  if [[ -f tmp_$chr.missnp ]];then
   $plink --bfile tmp_s1_$chr  --exclude tmp_$chr.missnp --make-bed --out tmp_s1_$chr --memory 2560
  fi
  #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. (In fact, some tools, such as BEAGLE 4, require such duplicate variants to be merged/removed during preprocessing.) --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
  $plink --bfile tmp_s1_$chr  --list-duplicate-vars --out tmp_$chr --memory 2560
  $plink --bfile tmp_s1_$chr  --exclude tmp_$chr.dupvar --make-bed --out tmp_s2_$chr --memory 2560
  $plink --bfile tmp_s2_$chr  --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --qual-scores $infofile1.newID 7 1 1 --qual-threshold 0.4 --make-bed --out $prefix --memory 2560
  #Rscript /fh/fast/dai_j/BEACON/BEACON_GRANT/code/update_const_fid_fam.R $prefix.fam
}
for chr in {1..22}
do
  filtering $chr &
done

#Two versions of following processes
#remove ambiguous version----------------------

for chr in {1..22}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter
  $snpflip -b $prefix.bim -f $reference -o snpfilp_output
  FILESIZE=$(stat -c%s snpfilp_output.reverse)
  #flip
  if [[ $FILESIZE -gt 0 ]];then
    $plink --bfile $prefix --flip snpflip_output.reverse --make-bed --out $prefix
  fi
  #remove variants in ambiguous strand
  FILESIZE1=$(stat -c%s snpfilp_output.ambiguous)
  if [[ $FILESIZE1 -gt 0 ]];then
    $plink --bfile $prefix --exclude snpfilp_output.ambiguous --make-bed --out ${prefix}_noambiguous
  fi
  #$plink --bfile ${prefix}_noambiguous --recode A-transpose --out ${prefix}_noambiguous
done

#to run lifover py
ml Python/2.7.15-foss-2018b
#convert hg19-hg38
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink #include liftover binary
chain=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/annotation/hg19ToHg38.over.chain.gz
liftover=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink/liftOverPlink.py
#https://github.com/sritchie73/liftOverPlink

#run extract_allbca_genotype.R to update variant ID (inlucding allels)
liftover(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous
  $plink --bfile $prefix --snps-only --biallelic-only --recode --maf 0.001 --out $prefix  --memory 2560 #biallelic only, maf 0.001
  python $liftover --map $prefix.map --out $outdir/chr${chr}_hg19tohg38_lifted --chain $chain
  python /fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink/rmBadLifts.py --map $outdir/chr${chr}_hg19tohg38_lifted.map --out $outdir/chr${chr}_hg19tohg38_good_lifted.map --log $outdir/chr${chr}_hg19tohg38_bad_lifted.dat
  cut -f 2 $outdir/chr${chr}_hg19tohg38_bad_lifted.dat > $outdir/chr${chr}_hg19tohg38_to_exclude.dat
  cut -f 4 $outdir/chr${chr}_hg19tohg38_lifted.bed.unlifted | sed "/^#/d" >> $outdir/chr${chr}_hg19tohg38_to_exclude.dat 
  # Note: this will clobber the lifted MAP file generated by `liftOverPlink`:
  $plink --file $prefix --recode --out $outdir/chr${chr}_hg19tohg38_lifted --exclude $outdir/chr${chr}_hg19tohg38_to_exclude.dat --memory 2560
  prefix1=${prefix}_hg19tohg38
  $plink --ped $outdir/chr${chr}_hg19tohg38_lifted.ped --map $outdir/chr${chr}_hg19tohg38_good_lifted.map --recode --out ${prefix1} --memory 2560
  
  $plink --file $prefix1 --chr $chr --snps-only --biallelic-only --make-bed --out ${prefix1} --memory 2560
  #$plink --bfile ${prefix1} --recode A-transpose --out ${prefix1} --memory 2560
  echo "chr${chr} done!"
}
for chr in {1..22}
do
  liftover $chr &
done

#snpflip
ml Python
snpflip(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38
  $snpflip -b $prefix.bim -f $hg38 -o snpfilp_output
  FILESIZE=$(stat -c%s snpfilp_output.reverse)
  #flip
  if [[ $FILESIZE -gt 0 ]];then
    $plink --bfile $prefix --flip snpfilp_output.reverse --make-bed --out ${prefix}_flip
  else
    cp $prefix.bed ${prefix}_flip.bed
    cp $prefix.bim ${prefix}_flip.bim
    cp $prefix.fam ${prefix}_flip.fam
  fi
  #remove variants in ambiguous strand
  FILESIZE1=$(stat -c%s snpfilp_output.ambiguous)
  if [[ $FILESIZE1 -gt 0 ]];then
    $plink --bfile ${prefix}_flip --exclude snpfilp_output.ambiguous --make-bed --out ${prefix}_flip
  fi
  
}
for chr in {1..22}
do
  snpflip $chr
  echo $chr
  prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip
  $plink --bfile ${prefix1} --recode A-transpose --out ${prefix1}
done

#convert to vcf for genotype extracting 
vcfsort=/fh/fast/dai_j/CancerGenomics/Tools/vcftools_0.1.12/bin/vcf-sort
bgzip=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/ensembl-vep/reference/htslib/bgzip

tovcf() {
  local chr="$1"
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip
   $plink --bfile $prefix --recode vcf --out $prefix 
  $vcfsort $prefix.vcf |$bgzip -c >$prefix.vcf.gz
}
for chr in {1..22}
do
  tovcf $chr 
done

#creat index file
ml tabix/0.2.6-GCCcore-8.3.0
for chr in {1..22}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip
  tabix -p vcf $prefix.vcf.gz
done

#keep ambiguous snp version---------------------------,final used

ml Python/2.7.15-foss-2018b
#convert hg19-hg38
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink #include liftover binary
chain=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/annotation/hg19ToHg38.over.chain.gz
liftover=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink/liftOverPlink.py
#https://github.com/sritchie73/liftOverPlink

liftover(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter
  $plink --bfile $prefix --snps-only --biallelic-only --recode --maf 0.001 --out $prefix  --memory 2560 #biallelic only, maf 0.001
  python $liftover --map $prefix.map --out $outdir/chr${chr}_hg19tohg38_lifted --chain $chain
  python /fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/liftOverPlink/rmBadLifts.py --map $outdir/chr${chr}_hg19tohg38_lifted.map --out $outdir/chr${chr}_hg19tohg38_good_lifted.map --log $outdir/chr${chr}_hg19tohg38_bad_lifted.dat
  cut -f 2 $outdir/chr${chr}_hg19tohg38_bad_lifted.dat > $outdir/chr${chr}_hg19tohg38_to_exclude.dat
  cut -f 4 $outdir/chr${chr}_hg19tohg38_lifted.bed.unlifted | sed "/^#/d" >> $outdir/chr${chr}_hg19tohg38_to_exclude.dat 
  # Note: this will clobber the lifted MAP file generated by `liftOverPlink`:
  $plink --file $prefix --recode --out $outdir/chr${chr}_filter_hg19tohg38_lifted --exclude $outdir/chr${chr}_hg19tohg38_to_exclude.dat --memory 2560
  prefix1=${prefix}_hg19tohg38
  $plink --ped $outdir/chr${chr}_filter_hg19tohg38_lifted.ped --map $outdir/chr${chr}_hg19tohg38_good_lifted.map --recode --out ${prefix1} --memory 2560
  
  $plink --file $prefix1 --chr $chr --snps-only --biallelic-only --make-bed --out ${prefix1} --memory 2560
  #$plink --bfile ${prefix1} --recode A-transpose --out ${prefix1} --memory 2560
  echo "chr${chr} done!"
}
for chr in {1..22}
do
  liftover $chr &
done

#snpflip
ml Python
snpflip(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38
  $snpflip -b $prefix.bim -f $hg38 -o ${outdir}/snpfilp_output_chr$chr
  FILESIZE=$(stat -c%s ${outdir}/snpfilp_output_chr$chr.reverse)
  #flip
  if [[ $FILESIZE -gt 0 ]];then
    $plink --bfile $prefix --flip ${outdir}/snpfilp_output_chr$chr.reverse --make-bed --out ${prefix}_flip  --memory 1000
  else
    cp $prefix.bed ${prefix}_flip.bed
    cp $prefix.bim ${prefix}_flip.bim
    cp $prefix.fam ${prefix}_flip.fam
  fi
  prefix1=${prefix}_flip
  $plink --bfile ${prefix1} --recode A-transpose --out ${prefix1} --memory 1000
}
for chr in {1..22}
do
  snpflip $chr &
  echo $chr
done

#convert to vcf for genotype extracting 
vcfsort=/fh/fast/dai_j/CancerGenomics/Tools/vcftools_0.1.12/bin/vcf-sort
bgzip=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/ensembl-vep/reference/htslib/bgzip
tovcf() {
  local chr="$1"
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip
   $plink --bfile $prefix --recode vcf --out $prefix 
  $vcfsort $prefix.vcf |$bgzip -c >$prefix.vcf.gz
}
for chr in {1..22}
do
  tovcf $chr 
done
#creat index file
ml tabix/0.2.6-GCCcore-8.3.0
for chr in {1..22}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip
  tabix -p vcf $prefix.vcf.gz
done

#merge beacon and cambridge-------------------
indir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation

#merge1
dataset1="cambridgewtccc_1000g"
dataset2="beacondbgapcontrol_1000g"
dataset3="merge_beacon_cambridge_1000g"

#merge2
dataset1="cambridgewtccc_qc_1000g"
dataset2="beacondbgapcontrol_qc_1000g"
dataset3="merge_beacon_cambridge_qc_1000g"
mergedataset() {
  local chr="$1"
  prefix1=$indir/$dataset1/chr${chr}_filter_hg19tohg38_flip
  prefix2=$indir/$dataset2/chr${chr}_filter_hg19tohg38_flip
  $$plink --bfile $prefix1  --bmerge $prefix2 --out $indir/$dataset3/chr${chr}_filter_hg19tohg38_flip --memory 2560
  $plink --bfile $indir/$dataset3/chr${chr}_filter_hg19tohg38_flip --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out $indir/$dataset3/chr${chr}_filter_hg19tohg38_flip --memory 2560
}
for chr in {1..22}
do
  mergedataset $chr &
done 


#convert to vcf for genotype extracting 
dataset="merge_beacon_cambridge_qc_1000g"
vcfsort=/fh/fast/dai_j/CancerGenomics/Tools/vcftools_0.1.12/bin/vcf-sort
bgzip=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/ensembl-vep/reference/htslib/bgzip
tovcf1() {
  local chr="$1"
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip
   $plink --bfile $prefix --recode vcf --out $prefix --memory 2560
}
for chr in {1..22}
do
  tovcf1 $chr &
done
wait
tovcf2() {
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip
  # $plink --bfile $prefix --recode vcf --out $prefix --memory 2560
  $vcfsort $prefix.vcf |$bgzip -c >$prefix.vcf.gz
}

for chr in {8..14}
do
  tovcf2 $chr 
done

for chr in {15..22}
do
  tovcf2 $chr 
done

#creat index file
ml tabix/0.2.6-GCCcore-8.3.0
for chr in {8..14}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip
  tabix -p vcf $prefix.vcf.gz
done