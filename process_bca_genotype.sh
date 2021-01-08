#!/usr/bin/env bash
#first remove ambiguous SNPs, than hg19-hg38

#ml Python #to use snpflip

#extract genotype data on all BCA samples,used to find common snps
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
snpflip=/fh/fast/dai_j/CancerGenomics/Tools/python3/snpflip/snpflip
reference=/fh/fast/stanford_j/Xiaoyu/Tools/reference/human_g1k_v37.fasta
hg38=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/hg38/hg38.fa

# for chr in {1..22}
# do
#   echo $chr
#   vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz 
#   #$plink --vcf $vcffile --const-fid 0 --snps-only --biallelic-only strict --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all
#   $snpflip -b /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all.bim -f $reference -o snpfilp_output
#   FILESIZE=$(stat -c%s snpfilp_output.reverse)
#   #flip
#   if [[ $FILESIZE -gt 0 ]];then
#     $plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all --flip snpflip_output.reverse --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all
#   fi
#   #remove variants in ambiguous strand
#   FILESIZE1=$(stat -c%s snpfilp_output.ambiguous)
#   if [[ $FILESIZE1 -gt 0 ]];then
#     $plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all --exclude snpfilp_output.ambiguous --make-bed --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all_noambiguous
#   fi
#   $plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all_noambiguous --recode A-transpose --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_all_noambiguous
# done

#filtering
filtering(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter
  vcffile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz
  infofile=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.info.gz
  infofile1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.info
  if [[ ! -f $infofile1 ]];then
    echo  gunzip $infofile
  fi
  $plink --vcf $vcffile --const-fid 0 --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out tmp_s1_$chr --memory 2560
  Rscript ../code/update_bimID.R tmp_s1_$chr.bim
  Rscript ../code/update_infoID.R $infofile1 tmp_s1_$chr.bim
  #If binary merging fails because at least one variant would have more than two alleles, a list of offending variant(s) will be written to plink.missnp.multi-allelic
  $plink --bfile tmp_s1_$chr --const-fid 0 --bmerge tmp_s1_$chr --merge-mode 6 --out tmp_$chr   --memory 2560
  if [[ -f tmp_$chr.missnp ]];then
   $plink --bfile tmp_s1_$chr --const-fid 0 --exclude tmp_$chr.missnp --make-bed --out tmp_s1_$chr --memory 2560
  fi
  #When multiple variants share the same bp coordinate and allele codes, it is likely that they are not actually distinct, and the duplicates should be merged or removed. (In fact, some tools, such as BEAGLE 4, require such duplicate variants to be merged/removed during preprocessing.) --list-duplicate-vars identifies them, and writes a report to plink.dupvar.
  $plink --bfile tmp_s1_$chr --const-fid 0 --list-duplicate-vars --out tmp_$chr --memory 2560
  $plink --bfile tmp_s1_$chr --const-fid 0 --exclude tmp_$chr.dupvar --make-bed --out tmp_s2_$chr --memory 2560
  $plink --bfile tmp_s2_$chr --const-fid 0 --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --qual-scores $infofile1.newID 7 1 1 --qual-threshold 0.4 --make-bed --out $prefix --memory 2560
  Rscript ../code/update_const_fid_fam.R $prefix.fam
}
for chr in {1..22}
do
  filtering $chr &
done


for chr in {1..22}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter
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

#ml python/2.7.15-foss-2018b
#convert hg19-hg38
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu/Tools/liftOverPlink #include liftover binary
chain=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/hg19ToHg38.over.chain.gz
liftover=/fh/fast/stanford_j/Xiaoyu/Tools/liftOverPlink/liftOverPlink.py
#https://github.com/sritchie73/liftOverPlink

#run extract_allbca_genotype.R to update variant ID (inlucding allels)
liftover(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous
  $plink --bfile $prefix --snps-only --biallelic-only --recode --maf 0.001 --out $prefix  --memory 2560 #biallelic only, maf 0.001
  python $liftover --map $prefix.map --out $outdir/chr${chr}_hg19tohg38_lifted --chain $chain
  python /fh/fast/stanford_j/Xiaoyu/Tools/liftOverPlink/rmBadLifts.py --map $outdir/chr${chr}_hg19tohg38_lifted.map --out $outdir/chr${chr}_hg19tohg38_good_lifted.map --log $outdir/chr${chr}_hg19tohg38_bad_lifted.dat
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
#ml Python
snpflip(){
  local chr="$1"
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous_hg19tohg38
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
  prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous_hg19tohg38_flip
  $plink --bfile ${prefix1} --recode A-transpose --out ${prefix1}
done

#convert to vcf for genotype extracting 
vcfsort=/fh/fast/dai_j/CancerGenomics/Tools/vcftools_0.1.12/bin/vcf-sort
tovcf() {
  local chr="$1"
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous_hg19tohg38_flip
   $plink --bfile $prefix --recode vcf --out $prefix 
  $vcfsort $prefix.vcf |bgzip -c >$prefix.vcf.gz
}
for chr in {1..22}
do
  tovcf $chr 
done
#creat index file
for chr in {1..22}
do
  echo $chr
  prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous_hg19tohg38_flip
  tabix -p vcf $prefix.vcf.gz
done


