#!/usr/bin/env bash
#only work on EA
#use Cambridge+BEACONcontrol as the validation, only work on BEEA
#BEACONcase+AMOS+other BEACONcontrol for discovery

#https://choishingwan.github.io/PRS-Tutorial/ldpred/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596916/
#https://biostat0903.github.io/DBSLMM/Scripts.html

plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink

#use imputed cambamos data
#extract cambamos data, target data---
#step1--------------
rm  ../result/BCA1000gnoambiguousmergelist.txt
for chr in {1..22}
do
echo /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous  >> ../result/BCA1000gnoambiguousmergelist.txt
done
$plink --merge-list ../result/BCA1000gnoambiguousmergelist.txt --make-bed \
--out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous

#do GWAS on Discovery
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous
$plink --bfile $prefix \
  --keep ../result/GWAS/Discovery_BEEA_CO_selectedsamples_plink.txt --make-bed \
  --out ../result/GWAS/Discovery_BEEA_CO_gwas
#update_fam() in R
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/Discovery_BEEA_CO_gwas
$plink --bfile  $prefix --covar  ../result/GWAS/Discovery_BEEA_CO_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --beta --hide-covar --ci 0.95 --out $prefix
  
#working on validation data
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous
$plink --bfile $prefix  \
 --keep ../result/validationBEEAsamples_plink.txt --make-bed --out $prefix1
 
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous
#QC
$plink \
    --bfile ${prefix} \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out ${prefix}_QC

$plink \
    --bfile ${prefix} \
    --keep ${prefix}_QC.fam \
    --extract ${prefix}_QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out ${prefix}_QC

$plink \
    --bfile ${prefix} \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.fam \
    --het \
    --out ${prefix}_QC

#Individuals that have a first or second degree relative in the sample (π>0.125) can be removed with the following command:
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous
$plink \
    --bfile $prefix \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.valid.sample \
    --rel-cutoff 0.125 \
    --out ${prefix}_QC
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous_QC
$plink --bfile $prefix  \
 --keep ../result/validationBEEAsamples_plink.txt --make-bed --out $prefix1


$plink \
    --bfile ${prefix1} \
    --keep ${prefix1}.valid.sample \
    --make-bed \
    --out ${prefix1}
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validationBEEA_filter_noambiguous_QC
$plink --bfile $prefix  \
 --keep ${prefix}.rel.id --make-bed --out $prefix
 


#step7---- 
#make BE/EA/BEEA files using PRS_LDpred.R
#work on base data using PRS_LDpred.R
#install ldpred
#ml Python
#python -m venv env
#source ldpredenv/bin/activate
#python -m pip install ldpred
#If you want to switch projects or otherwise leave your virtual environment, simply run:
#deactivate



#A1:effective allele,--only-hm3
run_ldpred () {
    meta="$1" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz
    target="$2" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
    N="$3" #sample size of base data
    echo $target
    if [[ -f $target.coord ]] ; then rm $target.coord; fi
    rm $target.weight*.txt
    #rm $target.score*.txt
    #rm $target.score*.txt.adj
    rm $target.ld*.gz
    echo "coord"
    ldpred coord \
     --rs MarkerName \
     --A1 Allele1 \
     --A2 Allele2 \
     --pos pos \
     --chr chr \
     --pval P \
     --only-hm3 \
     --eff Effect \
     --ssf-format CUSTOM \
     --eff_type LOGOR \
     --se StdErr \
     --N $N \
     --ssf $meta \
     --out $target.coord \
     --gf $target

    echo "gibbs"
    nsnp=$(wc -l $target.bim) #4937535 /3000=1645
    read -r -a tmp <<< "$nsnp"
    nldr=$((${tmp[0]}/3000))
    ldpred gibbs \
     --cf $target.coord \
     --ldr $nldr \
     --ldf $target.ld \
     --f 1 0.3 0.1 0.03 0.01 \
     --out $target.weight \
     --N $N
    
    #f is the probability that a marker is drawn from a Gaussian distribution, i.e., the fraction of causal markers.
    echo "score"
    ldpred score \
     --gf $target \
     --rf $target.weight \
     --cov-file $target.covariate \
     --f 1 0.3 0.1 0.03 0.01 \
     --out $target.score \
     --pf $target.pheno \
     --pf-format LSTANDARD 
}

#2646+3537+length(discoverysamples)=12169
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Discovery_Bonn_imp_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_BEEA 12169  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Validation_BEEA1.log  &



