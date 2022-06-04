#!/usr/bin/env bash
#add cambamos into the meta, only work on EA

#https://choishingwan.github.io/PRS-Tutorial/ldpred/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596916/
#https://biostat0903.github.io/DBSLMM/Scripts.html

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

#use imputed cambamos data
#extract cambamos data, target data---
#step1--------------
rm  ../result/BCA1000gnoambiguousmergelist.txt
for chr in {1..22}
do
echo /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous  >> ../result/BCA1000gnoambiguousmergelist.txt
done

#do GWAS on BEACON
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous
$plink --bfile $prefix \
  --keep ../result/GWAS/Beacon_EA_CO_selectedsamples_plink.txt --make-bed \
  --out ../result/GWAS/Beacon_EA_CO_gwas
#update_fam() in R
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/Beacon_EA_CO_gwas
$plink --bfile  $prefix --covar  ../result/GWAS/Beacon_EA_CO_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --beta --hide-covar --ci 0.95 --out $prefix
  

$plink --merge-list ../result/BCA1000gnoambiguousmergelist.txt --make-bed \
--out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous
$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous  \
--keep ../result/cambamossamples_plink.txt --make-bed \
--out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous

prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous
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
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous
$plink \
    --bfile $prefix \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.valid.sample \
    --rel-cutoff 0.125 \
    --out ${prefix}_QC
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous_QC
$plink --bfile $prefix  \
 --keep ../result/cambamossamples_plink.txt --make-bed --out $prefix1


$plink \
    --bfile ${prefix1} \
    --keep ${prefix1}.valid.sample \
    --make-bed \
    --out ${prefix1}
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/cambamos_filter_noambiguous_QC
$plink --bfile $prefix  \
 --keep ${prefix}.rel.id --make-bed --out $prefix
 


#step7---- 
#make BE/EA/BEEA files using PRS_LDpred.R
#work on base data using PRS_LDpred.R
#install ldpred
#ml Python
#python -m venv env
#source ./env/bin/activate
#python -m pip install ldpred
#If you want to switch projects or otherwise leave your virtual environment, simply run:
#deactivate



#A1:effective allele,--only-hm3
run_ldpred () {
    meta="$1" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz
    target="$2" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
    N="$3" #sample size of base data
    echo $target
    # if [[ -f $target.coord ]] ; then rm $target.coord; fi
    # rm $target.weight*.txt
    # #rm $target.score*.txt
    # #rm $target.score*.txt.adj
    # rm $target.ld*.gz
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


# run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE 14229  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_BE.log &
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA 8843 &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA.log & 
# run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA 11459 &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_BEEA.log &
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_imp_noAmos_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_noAmos 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_noAmos.log  &

#high quality SNPs
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.covariate
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_info_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info 8843 &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_info1.log & 


run_ldpred1 () {
    meta="$1" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz
    target="$2" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
    N="$3" #sample size of base data
    echo $target
    # if [[ -f $target.coord ]] ; then rm $target.coord; fi
    # rm $target.weight*.txt
    # #rm $target.score*.txt
    # #rm $target.score*.txt.adj
    # rm $target.ld*.gz
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
     --case-n case_N \
     --control-n control_N \
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
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn.covariate
run_ldpred1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_cn1.log  &
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info.covariate
run_ldpred1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_info_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_cn_info1.log  &



