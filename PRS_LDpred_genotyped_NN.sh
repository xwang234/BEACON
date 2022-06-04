#!/usr/bin/env bash
#Use genotyped data
#to compare with NN only work on EA

#https://choishingwan.github.io/PRS-Tutorial/ldpred/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4596916/
#https://biostat0903.github.io/DBSLMM/Scripts.html

plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink

#working on validation data
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped
$plink --bfile $prefix  \
 --keep ../result/PRSDNN_test_plinksamples.txt --make-bed --out $prefix1
 
#samples with extreme heterozygosity are typically removed prior to downstream analyses
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped
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
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped"
$plink \
    --bfile $prefix \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.valid.sample \
    --rel-cutoff 0.125 \
    --out ${prefix}_QC
    
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped"
prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped_QC"
$plink --bfile $prefix  \
 --keep ../result/PRSDNN_test_plinksamples.txt --make-bed --out $prefix1

#step 3----      
$plink \
    --bfile ${prefix}_QC \
    --keep ${prefix}_QC.valid.sample \
    --make-bed \
    --out ${prefix1}
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtest_genotyped_QC
$plink --bfile $prefix  \
 --keep ${prefix}.rel.id --make-bed --out $prefix
 

 
#make BE/EA/BEEA files using PRS_LDpred.R
#work on base data using PRS_LDpred.R
#install ldpred
ml Python/3.7.4-GCCcore-8.3.0
python -m venv ldpredenv1.0.9
source ./ldpredenv1.0.9/bin/activate
python -m pip install ldpred==1.0.9
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
    echo $target
    echo "coord"
    # if [[ -f $target.coord ]] ; then rm $target.coord; fi
    # rm $target.weight*.txt
    # #rm $target.score*.txt
    # rm $target.score*.txt.adj
    # rm $target.ld*.gz
    # ldpred coord \
    #  --rs MarkerName \
    #  --A1 Allele2 \
    #  --A2 Allele1 \
    #  --pos pos \
    #  --chr chr \
    #  --pval P \
    #  --eff Effect \
    #  --eff_type OR \
    #  --se StdErr \
    #  --ssf-format CUSTOM \
    #  --N $N \
    #  --ssf $meta \
    #  --out $target.coord \
    #  --gf $target

    echo "gibbs"
    nsnp=$(wc -l $target.bim) #4937535 /3000=1645
    read -r -a tmp <<< "$nsnp"
    nldr=$((${tmp[0]}/3000))
    #f is the probability that a marker is drawn from a Gaussian distribution, i.e., the fraction of causal markers.
    #--f 1 0.3 0.1 0.03 0.01 0.003 0.001  \
    ldpred gibbs \
     --cf $target.coord \
     --ldr $nldr \
     --ldf $target.ld \
     --f 1 0.3 0.1 0.03 0.05 0.01 0.003 0.001\
     --out $target.weight \
     --N $N
    
    echo "score"
    ldpred score \
     --gf $target \
     --rf $target.weight \
     --cov-file $target.covariate \
     --f 1 0.3 0.1 0.03 0.05 0.01 0.003 0.001 \
     --out $target.score \
     --pf $target.pheno \
     --pf-format LSTANDARD 
}
wc -l ../result/PRSDNN_train_plinksamples.txt #4578
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtrain.assoc.logistic.metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/DNNvalidition_EA_genotyped 4578  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_DNNvalidition_EA_genotyped.log  &


