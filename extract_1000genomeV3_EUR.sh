
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/"

for chr in {1..22}
do
  echo $chr
  vcffile=${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  #$plink --vcf $vcffile --snps-only --make-bed --out ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --memory 6400
  #cut -f 2 /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_select.bim | sort | uniq -d > plink.dups
  #$plink --bfile ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --keep ${thousanddir}EURsamples.txt --make-bed --out ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR
  #$plink --bfile ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --recodeA --out ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --memory 6400
  $plink --bfile ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --recode A-transpose --out ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --memory 6400
  #$plink --bfile ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --recode --tab --out ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR
done

chr=5
tmp=($(awk '{ print $2 }' ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim))
for i in {1..100}
do
   echo ${tmp[$i]} >> ${thousanddir}test.snp
done

$plink --bfile ${thousanddir}ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR --extract ${thousanddir}test.snp --recode --tab --out ${thousanddir}test5
$plink --file ${thousanddir}test5 --make-bed --out ${thousanddir}test5
$plink --bfile ${thousanddir}test5 --recode A-transpose --out ${thousanddir}test5 --memory 6400
$plink --bfile ${thousanddir}test5 --recode A --out ${thousanddir}test5 --memory 6400
 