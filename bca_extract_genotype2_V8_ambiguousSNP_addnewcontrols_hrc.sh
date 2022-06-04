#!/usr/bin/env bash
#split each chr into 2
#use hg38 VCFs chr${chr}_filter_noambiguous_hg19tohg38_flip.vcf.gz
#to extract new beacon/cambridge SNPs

plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink

prefix="$1" 
#dataset="$2"
#use the merge data
dataset="merge_beacon_cambridge_hrc_maf001_snp"

outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/$prefix
if [[ ! -d "$outfolder" ]]; then mkdir "$outfolder"; fi

ml tabix/0.2.6-GCCcore-8.3.0
extract_vcf1(){
  local chr="$1"
  echo $chr
  if [[ -f $outfolder/chr${chr}_1.select.dose.vcf ]]; then rm $outfolder/chr${chr}_1.select.dose.vcf; fi 
  tabixreg=$outfolder/prediction_snps_tabix_chr${chr}_1.txt
  #tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip.vcf.gz >$outfolder/chr${chr}.select.dose.vcf
  #xargs -a $tabixreg -I {} echo {} |more                   
  #xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip.vcf.gz {} >> $outfolder/chr${chr}_1.select.dose.vcf
  xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip.vcf.gz {} >> $outfolder/chr${chr}_1.select.dose.vcf

}

extract_vcf2(){
  local chr="$1"
  echo $chr
  if [[ -f $outfolder/chr${chr}_2.select.dose.vcf ]]; then rm $outfolder/chr${chr}_2.select.dose.vcf; fi
  tabixreg=$outfolder/prediction_snps_tabix_chr${chr}_2.txt
  #tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_noambiguous_hg19tohg38_flip.vcf.gz >$outfolder/chr${chr}.select.dose.vcf
  #xargs -a $tabixreg -I {} echo {} |more
  xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip.vcf.gz {} >> $outfolder/chr${chr}_2.select.dose.vcf
}

plink_geno(){
  local chr="$1"
  #echo $chr
  vcffile=$outfolder/chr${chr}.test.select.dose.vcf
    #$plink --vcf $vcffile --const-fid 0 --snps-only --make-bed --out $outfolder/chr${chr}_select --memory 120
    $plink --vcf $vcffile --snps-only --make-bed --out $outfolder/chr${chr}_select --memory 120
    #cut -f 2 /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_select.bim | sort | uniq -d > plink.dups
    #$plink --bfile $outfolder/chr${chr}_select --exclude plink.dups --make-bed --out $outfolder/chr${chr}_select
    $plink --bfile $outfolder/chr${chr}_select --recodeA --out $outfolder/chr${chr}_select --memory 640
}
for chr in {1..22};
do
  if [[ -f "$outfolder/prediction_snps_tabix_chr${chr}_1.txt" ]];
  then
    echo $chr
    tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/"$dataset"/chr${chr}_filter_hg19tohg38_flip.vcf.gz >$outfolder/chr${chr}.head.select.dose.vcf
    extract_vcf1 $chr &
    extract_vcf2 $chr &
  fi
done
wait

#combine 2 pieces in each chr
for chr in {1..22};
do
  if [[ -f "$outfolder/prediction_snps_tabix_chr${chr}_1.txt" ]];
  then
    cat $outfolder/chr${chr}.head.select.dose.vcf $outfolder/chr${chr}_1.select.dose.vcf $outfolder/chr${chr}_2.select.dose.vcf >$outfolder/chr${chr}.test.select.dose.vcf
  fi
done

#generate plink form
for chr in {1..22}
do
  if [[ -f "$outfolder/prediction_snps_tabix_chr${chr}_1.txt" ]];
  then
    echo $chr
    plink_geno $chr &
  fi
done
wait

