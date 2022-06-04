
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
cd $outfolder

$plink --bfile $outfolder/bca_filtered_30Nov2018 --recode tab --snps-only --biallelic-only strict --maf 0.05 --make-bed --out $outfolder/bca_filtered_10Jan2019
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin
cp ./bca_filtered_10Jan2019.fam ./bca_filtered_10Jan2019.pedind
cp ./bca_filtered_10Jan2019.bim ./bca_filtered_10Jan2019.pedsnp
/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin/convertf -p bca_filtered_10Jan2019.par.PED.EIGENSTRAT
perl /fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin/smartpca.perl -i bca_filtered_10Jan2019.eigenstratgeno -a bca_filtered_10Jan2019.snp -b bca_filtered_10Jan2019.ind -k 15 -o bca_filtered_10Jan2019.pca -e bca_filtered_10Jan2019.eval -l bca_filtered_10Jan2019.eigpca.log -p bca_filtered_10Jan2019.plot -s 6

#work on the genotyped data with new controls
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID
prefix2=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC
prefix3=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype
mergedataset() {
  local prefix1="$1"
  local prefix2="$2"
  local prefix3="$3"
  $plink --bfile $prefix1  --bmerge $prefix2 --out $prefix3 --memory 2560
  $plink --bfile $prefix3 --maf 0.05 --biallelic-only strict --snps-only --hwe 0.00001 --geno 0.05 --make-bed --out $prefix3 --memory 2560
  $plink --bfile $prefix3 --recode --out $prefix3 --memory 2560
}

#for combined
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype
#for individuals
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC
ml GCC/8.3.0
ml GSL/2.6-GCC-8.3.0
ml OpenBLAS/0.3.7-GCC-8.3.0
eigenstrat(){
  export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin
  #export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/EIG-7.2.1/bin
  prefix="$1"
  cp $prefix.fam $prefix.pedind
  cp $prefix.bim $prefix.pedsnp
  eigenstratconfig=$prefix.EIGENSTRAT
  if [[ -f $eigenstratconfig ]]; then rm $eigenstratconfig; fi
  echo "genotypename:    $prefix.ped" > $eigenstratconfig
  echo "snpname:         $prefix.pedsnp" >> $eigenstratconfig
  echo "indivname:       $prefix.pedind" >> $eigenstratconfig
  echo "outputformat:    EIGENSTRAT" >> $eigenstratconfig
  echo "genotypeoutname:    $prefix.eigenstratgeno" >> $eigenstratconfig
  echo "snpoutname:    $prefix.snp" >> $eigenstratconfig
  echo "indivoutname:    $prefix.ind" >> $eigenstratconfig
  echo "familynames:     NO" >> $eigenstratconfig
  echo "killr2:    YES" >> $eigenstratconfig  #prune 
  echo "r2thresh:    0.5" >> $eigenstratconfig
  if [[ ! -f $prefix.ped ]]; then $plink --bfile $prefix --recode --out $prefix ; fi
  echo "convertf----"
  /fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin/convertf -p $eigenstratconfig
  echo "smartpca-----"
  perl /fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/EIG-7.2.1/bin/smartpca.perl -i $prefix.eigenstratgeno -a $prefix.snp -b $prefix.ind -k 15 -o $prefix.pca -e prefix.eval -l $prefix.eigpca.log -p $prefix.plot -s 6
  #perl /fh/fast/dai_j/CancerGenomics/Tools/EIG-7.2.1/bin/smartpca.perl -i $prefix.eigenstratgeno -a $prefix.snp -b $prefix.ind -k 15 -o $prefix.pca -e prefix.eval -l $prefix.eigpca.log -p $prefix.plot -s 6
 
  echo "done"

}
