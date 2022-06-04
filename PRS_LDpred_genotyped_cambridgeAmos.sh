#!/usr/bin/env bash
#Use genotyped data
#use Cambridge+AMOS as the validation, only work on EA

#https://choishingwan.github.io/PRS-Tutorial/ldpred/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596916/
#https://biostat0903.github.io/DBSLMM/Scripts.html

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

#extract cambamos(cambridge+amos) data, target data---
#step 1----
#remove ambiguous SNPs
#genotype A/B
#prefix=/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip
#ml Python #to use snpflip
snpflip=/fh/fast/dai_j/CancerGenomics/Tools/python3/snpflip/snpflip
reference=/fh/fast/stanford_j/Xiaoyu/Tools/reference/human_g1k_v37.fasta
$snpflip -b $prefix.bim -f $reference -o $prefix1
FILESIZE=$(stat -c%s ${prefix1}.reverse)
#flip
if [[ $FILESIZE -gt 0 ]];then
  $plink --bfile $prefix --flip ${prefix1}.reverse --make-bed --out ${prefix1}  --memory 1000
else
  cp $prefix.bed ${prefix1}.bed
  cp $prefix.bim ${prefix1}.bim
  cp $prefix.fam ${prefix1}.fam
fi
#remove variants on ambiguous strand
FILESIZE1=$(stat -c%s ${prefix1}.ambiguous)
if [[ $FILESIZE1 -gt 0 ]];then
  $plink --bfile $prefix1 --exclude ${prefix1}.ambiguous --make-bed --out ${prefix1}_noambiguous
fi

prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous
$plink --bfile $prefix  \
 --keep ../result/cambamossamples_plink.txt --make-bed --out $prefix1
 
#samples with extreme heterozygosity are typically removed prior to downstream analyses
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous
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
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous"
$plink \
    --bfile $prefix \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.valid.sample \
    --rel-cutoff 0.125 \
    --out ${prefix}_QC
    
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous"
prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous_QC"
$plink --bfile $prefix  \
 --keep ../result/cambamossamples_plink.txt --make-bed --out $prefix1

#step 3----      
$plink \
    --bfile ${prefix}_QC \
    --keep ${prefix}_QC.valid.sample \
    --make-bed \
    --out ${prefix1}
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/cambamos_filtered_30Nov2018_flip_noambiguous_QC
$plink --bfile $prefix  \
 --keep ${prefix}.rel.id --make-bed --out $prefix
 

 
#make BE/EA/BEEA files using PRS_LDpred.R
#work on base data using PRS_LDpred.R
#install ldpred
#ml Python/3.6.6-foss-2018b
python -m venv env
source ./env/bin/activate
python -m pip install --user ldpred
#If you want to switch projects or otherwise leave your virtual environment, simply run:
#deactivate
#If you want to re-enter the virtual environment just follow the same instructions above about activating a virtual environment. There’s no need to re-create the virtual environment.


# #ldpred=/fh/fast/dai_j/CancerGenomics/Tools/python3/lib/python3.6/site-packages/LDpred-1.0.11-py3.6.egg/ldpred/LDpred_fast.py
# #N numer of samples in base data
# ldpred coord \
# --rs MarkerName \
# --A1 Allele1 \
# --A2 Allele2 \
# --pos pos \
# --chr chr \
# --pval P \
# --eff Effect \
# --ssf-format CUSTOM \
# --N 14229 \
# --ssf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz \
# --out /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.coord \
# --gf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
# 
# nsnp=$(wc -l /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.bim) #4937535 /3000=1645
# read -r -a tmp <<< "$nsnp"
# nldr=$((${tmp[0]}/3000))
# ldpred gibbs \
#     --cf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.coord \
#     --ldr $nldr \
#     --ldf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.ld \
#     --out /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.weight \
#     --N 14229
#     
# ldpred score \
#     --gf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.ldpred \
#     --rf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.weight \
#     --out /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.score \
#     --pf /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE.pheno \
#     --pf-format LSTANDARD 

 
run_ldpred () {
    meta="$1" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz
    target="$2" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
    N="$3" #sample size of base data
    # echo $target
    # # echo "coord"
    # if [[ -f $target.coord ]] ; then rm $target.coord; fi
    # rm $target.weight*.txt
    # # rm $target.score*.txt
    # # rm $target.score*.txt.adj
    # # rm $target.ld*.gz
    # ldpred coord \
    #  --rs MarkerName \
    #  --A1 Allele1 \
    #  --A2 Allele2 \
    #  --pos pos \
    #  --chr chr \
    #  --pval P \
    #  --eff Effect \
    #  --eff_type LOGOR \
    #  --se StdErr \
    #  --ssf-format CUSTOM \
    #  --N $N \
    #  --ssf $meta \
    #  --out $target.coord \
    #  --gf $target
    # 
    # echo "gibbs"
    # nsnp=$(wc -l $target.bim) #4937535 /3000=1645
    # read -r -a tmp <<< "$nsnp"
    # nldr=$((${tmp[0]}/3000))
    # ldpred gibbs \
    #  --cf $target.coord \
    #  --ldr $nldr \
    #  --ldf $target.ld \
    #  --f 1 0.3 0.1 0.03 0.01 \
    #  --out $target.weight \
    #  --N $N
    
    #f is the probability that a marker is drawn from a Gaussian distribution, i.e., the fraction of causal markers.
    echo "score"
    ldpred score \
     --gf $target \
     --rf $target.weight \
     --cov-file $target.covariate \
     --f 1 0.3 0.1 0.03 0.01  \
     --out $target.score \
     --pf $target.pheno \
     --pf-format LSTANDARD 
}

run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_genotyped.log  &
# run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped 11459  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_BEEA_genotyped.log  &

#high quality SNPs
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.covariate
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_info_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_genotyped_info1.log  &

#add case_N,control_N
run_ldpred1 () {
    meta="$1" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz
    target="$2" #/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE
    N="$3" #sample size of base data
    # echo $target
    # # echo "coord"
    # if [[ -f $target.coord ]] ; then rm $target.coord; fi
    # rm $target.weight*.txt
    # # rm $target.score*.txt
    # # rm $target.score*.txt.adj
    # # rm $target.ld*.gz
    # ldpred coord \
    #  --rs MarkerName \
    #  --A1 Allele1 \
    #  --A2 Allele2 \
    #  --pos pos \
    #  --chr chr \
    #  --pval P \
    #  --eff Effect \
    #  --eff_type LOGOR \
    #  --se StdErr \
    #  --ssf-format CUSTOM \
    #  --N $N \
    #  --case-n case_N \
    #  --control-n control_N \
    #  --ssf $meta \
    #  --out $target.coord \
    #  --gf $target
    # 
    # echo "gibbs"
    # nsnp=$(wc -l $target.bim) #4937535 /3000=1645
    # read -r -a tmp <<< "$nsnp"
    # nldr=$((${tmp[0]}/3000))
    # ldpred gibbs \
    #  --cf $target.coord \
    #  --ldr $nldr \
    #  --ldf $target.ld \
    #  --f 1 0.3 0.1 0.03 0.01 \
    #  --out $target.weight \
    #  --N $N
    
    #f is the probability that a marker is drawn from a Gaussian distribution, i.e., the fraction of causal markers.
    echo "score"
    ldpred score \
     --gf $target \
     --rf $target.weight \
     --cov-file $target.covariate \
     --f 1 0.3 0.1 0.03 0.01  \
     --out $target.score \
     --pf $target.pheno \
     --pf-format LSTANDARD 
}

cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn.covariate
run_ldpred1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_genotyped_cn1.log  &
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.bim /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info.bim
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.bed /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.fam /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info.fam
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.pheno /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info.pheno
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info.covariate /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info.covariate
run_ldpred1 /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_Bonn_info_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info 8843  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_CambAmos_EA_genotyped_cn_info1.log  &

