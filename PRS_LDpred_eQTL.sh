#!/usr/bin/env bash

plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

#step1--------------
rm  ../result/BCA1000gnoambiguousmergelist.txt
for chr in {1..22}
do
echo /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}_filter_noambiguous  >> ../result/BCA1000gnoambiguousmergelist.txt
done

$plink --merge-list ../result/BCA1000gnoambiguousmergelist.txt --make-bed \
--out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous

$plink --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous  \
--keep ../result/beaconsamples_plink.txt --make-bed \
--out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous

#step2
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous
$plink --bfile $prefix --extract ../result/beacon_eQTLkeepSNPs.txt --make-bed --out ${prefix}_eQTL

prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL
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

#step4
#Individuals that have a first or second degree relative in the sample (π>0.125) can be removed with the following command:
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL
$plink \
    --bfile $prefix \
    --extract ${prefix}_QC.prune.in \
    --keep ${prefix}_QC.valid.sample \
    --rel-cutoff 0.125 \
    --out ${prefix}_QC
  
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL
prefix1=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL_QC
$plink --bfile $prefix  \
 --keep ../result/beaconsamples_plink.txt --make-bed --out $prefix1


$plink \
    --bfile ${prefix1} \
    --keep ${prefix1}.valid.sample \
    --make-bed \
    --out ${prefix1}
    
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL_QC
$plink --bfile $prefix  \
 --keep ${prefix}.rel.id --make-bed --out $prefix

#last step
#A1:effective allele
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
    # echo "coord"
    # ldpred coord \
    #  --rs MarkerName \
    #  --A1 Allele1 \
    #  --A2 Allele2 \
    #  --pos pos \
    #  --chr chr \
    #  --pval P \
    #  --eff Effect \
    #  --ssf-format CUSTOM \
    #  --eff_type LOGOR \
    #  --se StdErr \
    #  --N $N \
    #  --ssf $meta \
    #  --out $target.coord \
    #  --gf $target

    echo "gibbs"
    nsnp=$(wc -l $target.bim) #4937535 /3000=1645
    read -r -a tmp <<< "$nsnp"
    nldr=$((${tmp[0]}/3000))
    ldpred gibbs \
     --cf $target.coord \
     --ldr $nldr \
     --ldf $target.ld \
     --f 0.01 0.003 0.001 0.0003 0.0001 \
     --out $target.weight \
     --N $N
    
    #f is the probability that a marker is drawn from a Gaussian distribution, i.e., the fraction of causal markers.
    echo "score"
    ldpred score \
     --gf $target \
     --rf $target.weight \
     --cov-file $target.covariate \
     --f 0.01 0.003 0.001 0.0003 0.0001 \
     --out $target.score \
     --pf $target.pheno \
     --pf-format LSTANDARD 
}


run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_eQTL 14229  &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_BE_eQTL.log &
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_eQTL 9549 &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_EA_eQTL.log & 
run_ldpred /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt.gz /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL 11459 &> /fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_Beacon_BEEA_eQTL.log &

