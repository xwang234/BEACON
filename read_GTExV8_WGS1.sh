#!/usr/bin/env bash

#without removing ambiguous SNPs

#ml Python #to use snpflip
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
hg38=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/annotation/Homo_sapiens_assembly38.fasta #not working
hg38=/fh/fast/dai_j/CancerGenomics/Tools/database/reference/hg38/hg38.fa
snpflip=/fh/fast/dai_j/CancerGenomics/Tools/python3/snpflip/snpflip

##process on dbsnp151 hg1i and hg38, not used
#prefix=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180418 #hg38
#$plink --vcf $prefix.vcf.gz --allow-no-samples --snps-only --biallelic-only strict --make-bed --out $prefix 

#prefix=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180423
#$plink --vcf $prefix.vcf.gz --allow-no-samples --snps-only --biallelic-only strict --make-bed --out $prefix 

#prefix=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180418
#prefix1=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/dbsnp151_hg38/hg38
#for chr in {1..22}
#do
#  echo $chr
#  $plink --bfile $prefix --allow-no-samples --chr $chr --extract /fh/fast/stanford_j/Xiaoyu/Tools/annotation/dbsnp151_hg19_hg38_snpid.txt --make-bed --out ${prefix1}_${chr}
#done

#prefix=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_20180423
#prefix1=/fh/fast/stanford_j/Xiaoyu/Tools/annotation/dbsnp151_hg19/hg19
#for chr in {1..22}
#do
#  echo $chr
#  $plink --bfile $prefix --allow-no-samples --chr $chr --extract /fh/fast/stanford_j/Xiaoyu/Tools/annotation/dbsnp151_hg19_hg38_snpid.txt --make-bed --out ${prefix1}_${chr}
#done

wgsfile="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz"

prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze"
$plink --vcf $wgsfile --biallelic-only strict --out $prefix 


prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze"
prefix1="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter"
$plink --bfile $prefix --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.01 --snps-only --make-bed --out $prefix1
$plink --bfile $prefix1 --recode A-transpose --out $prefix1
for chr in {1..22}
do
  echo $chr
  $plink --bfile $prefix1 --chr $chr --biallelic-only strict --snps-only --make-bed --out ${prefix1}_chr${chr}
  $plink --bfile ${prefix1}_chr${chr} --recode A-transpose --out ${prefix1}_chr${chr}
done

#snpfilp
##snpflip(){
  
  prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter"
  sed 's/chr//g' $prefix.bim >$prefix.bim.copy #to apply snpflip, bim chr should not include "chr"
  #sed 's/chr//g' /fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.dict >/fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.dict.copy
  #sed 's/>chr//g' /fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.fasta >/fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.fasta.copy
  #sed 's/chr//g' /fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.fasta.fai >/fh/fast/stanford_j/Xiaoyu/Tools/annotation/Homo_sapiens_assembly38.fasta.fai.copy

  $snpflip -b $prefix.bim -f $hg38 -o snpfilp_output
  FILESIZE=$(stat -c%s snpfilp_output.reverse)
  #flip, actually no reverse strands were found, just copy
  if [[ $FILESIZE -gt 0 ]];then
    $plink --bfile $prefix --flip snpfilp_output.reverse --make-bed --out ${prefix}_flip
  else
    cp $prefix.bed ${prefix}_flip.bed
    cp $prefix.bim ${prefix}_flip.bim
    cp $prefix.fam ${prefix}_flip.fam
  fi
  #remove variants in ambiguous strand
  #FILESIZE1=$(stat -c%s snpfilp_output.ambiguous)
  #if [[ $FILESIZE1 -gt 0 ]];then
  #  $plink --bfile $prefix --exclude snpfilp_output.ambiguous --make-bed --out ${prefix}_noambiguous
  #fi
#}

plink1()
{
  local chr="$1"
  #echo $chr
  prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter_flip"
  $plink --bfile $prefix --chr $chr --biallelic-only strict --snps-only --mind 0.1 --geno 0.05 --hwe 0.00001 --maf 0.05 --make-bed --out ${prefix}_chr${chr} --memory 1000 #maf>0.05
  $plink --bfile ${prefix}_chr${chr} --recode A-transpose --out ${prefix}_chr${chr} --memory 1000
}

for chr in {1..22}
do
  echo $chr
  plink1 $chr &
  #$plink --bfile $prefix1 --chr $chr --biallelic-only strict --snps-only --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.05 --make-bed --out ${prefix1}_chr${chr} #maf>0.05
  #$plink --bfile ${prefix1}_chr${chr} --recode A-transpose --out ${prefix1}_chr${chr}
done

#eigenstrat
cd /fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin
prefix="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter"
prefix1=${prefix}maf05
$plink --bfile $prefix --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.01 --snps-only --make-bed --out $prefix1
$plink --bfile $prefix1 --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.01 --snps-only --recode --tab --out $prefix1
cp ${prefix1}.fam ${prefix1}.pedind
#gtex_gender.R includes genders for pedind file #not used
cp ${prefix1}.bim ${prefix1}.pedsnp
#sed 's/*/0/g' gtex.ped >gtex1.ped #remove * in the ped file
/fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin/convertf -p gtex.par.PED.EIGENSTRAT
perl /fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin/smartpca.perl -i GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.eigenstratgeno \
-a GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.snp \
-b GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.ind \
-k 15 -o GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.pca \
-e GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.eval \
-l GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.eigpca.log \
-p GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filtermaf05.plot -s 20

#work on new filters (maf=0.01,add insersion/deletion)
#maf0.01
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype
prefix="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze"
prefix1="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_maf001"
$plink --bfile $prefix --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out $prefix1


#$plink --bfile $prefix1 --recode A-transpose --out $prefix1
# for chr in {1..22}
# do
#   echo $chr
#   $plink --bfile $prefix1 --chr $chr --biallelic-only strict --make-bed --out ${prefix1}_chr${chr}
#   $plink --bfile ${prefix1}_chr${chr} --recode A-transpose --out ${prefix1}_chr${chr}
# done

#snpfilp
prefix="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_maf001"
sed 's/chr//g' $prefix.bim >$prefix.bim.copy #to apply snpflip, bim chr should not include "chr"
 
$snpflip -b $prefix.bim -f $hg38 -o snpfilp_output
FILESIZE=$(stat -c%s snpfilp_output.reverse)
#flip, actually no reverse strands were found, just copy
if [[ $FILESIZE -gt 0 ]];then
  $plink --bfile $prefix --flip snpfilp_output.reverse --make-bed --out ${prefix}_flip
else
  cp $prefix.bed ${prefix}_flip.bed
  cp $prefix.bim ${prefix}_flip.bim
  cp $prefix.fam ${prefix}_flip.fam
fi

plinkchr()
{
  local chr="$1"
  #echo $chr
  prefix="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_maf001_flip"
  $plink --bfile $prefix --chr $chr --biallelic-only strict --mind 0.1 --geno 0.05 --hwe 0.00001  --make-bed --out ${prefix}_chr${chr} --memory 2000 #maf>0.01
  $plink --bfile ${prefix}_chr${chr} --recode A-transpose --out ${prefix}_chr${chr} --memory 2000
  prefix1="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_maf005_flip"
  $plink --bfile $prefix --chr $chr --biallelic-only strict --maf 0.05 --mind 0.1 --geno 0.05 --hwe 0.00001  --make-bed --out ${prefix1}_chr${chr} --memory 2000 #maf>0.05
  $plink --bfile ${prefix1}_chr${chr} --recode A-transpose --out ${prefix1}_chr${chr} --memory 2000
}

for chr in {1..22}
do
  echo $chr
  plinkchr $chr &
  #$plink --bfile $prefix1 --chr $chr --biallelic-only strict --snps-only --mind 0.05 --geno 0.05 --hwe 0.00001 --maf 0.05 --make-bed --out ${prefix1}_chr${chr} #maf>0.05
  #$plink --bfile ${prefix1}_chr${chr} --recode A-transpose --out ${prefix1}_chr${chr}
done

#keep all snps maf=0
ml Python/3.7.4-GCCcore-8.3.0
source /fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1/bin/activate
snpflip=/fh/fast/dai_j/CancerGenomics/Tools/snpflip/env1/bin/snpflip
outdir=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype

#snpfilp
prefix="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze"
sed 's/chr//g' $prefix.bim >$prefix.bim.copy #to apply snpflip, bim chr should not include "chr"
 
$snpflip -b $prefix.bim -f $hg38 -o GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_snpflip_output
FILESIZE=$(stat -c%s GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_snpflip_output.reverse)
#flip, actually no reverse strands were found, just copy
if [[ $FILESIZE -gt 0 ]];then
  $plink --bfile $prefix --flip GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_snpflip_output.reverse --make-bed --out ${prefix}_flip
else
  cp $prefix.bed ${prefix}_flip.bed
  cp $prefix.bim ${prefix}_flip.bim
  cp $prefix.fam ${prefix}_flip.fam
fi

plinkchr()
{
  local chr="$1"
  #echo $chr
  prefix="$outdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_flip"
  $plink --bfile $prefix --chr $chr --mind 0.05 --geno 0.05 --hwe 0.00001  --make-bed --out ${prefix}_chr${chr} 
  $plink --bfile ${prefix}_chr${chr} --recode A-transpose --out ${prefix}_chr${chr} 
}

for chr in {1..22}
do
  echo $chr
  plinkchr $chr
done
