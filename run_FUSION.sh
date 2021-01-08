#!/usr/bin/env bash

#Add the bundled GCTA binary gcta_nr_robust to path (coded by Po-Ru Loh for robust non-linear optimization)
export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master
#Download and install PLINK2, add plink to path
export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/plink2

fusiondir=/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master
#an example
Rscript $fusiondir/FUSION.assoc_test.R \
--sumstats ../result/fusion/sumstats/PGC2.SCZ.sumstats \
--weights ../result/fusion/WEIGHTS/GTEx.Whole_Blood.pos \
--weights_dir ../result/fusion/WEIGHTS/ \
--ref_ld_chr /fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/LDREF/1000G.EUR. \
--chr 22 \
--out ../result/fusion/result/PGC2.SCZ.22.dat

#to compute weights for GTEX V8 data------
#https://github.com/gusevlab/fusion_twas/blob/master/examples/GTEX_v7.sh

#process genotype, get HM3 genotype data
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.bed \
  /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid.bed
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.bim \
  /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid.bim  
cp /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.fam \
  /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid.fam

plink2 --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_rsid \
  --extract /fh/fast/dai_j/BEACON/BEACON_GRANT/data/HM3/snplist.txt --make-bed \
  --out /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_HM3
  
cd /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx
GCTA="gcta_nr_robust"
PLINK="plink2 --allow-no-sex"
GEMMA="/fh/fast/dai_j/CancerGenomics/Tools/GEMMA/gemma-0.98.1-linux-static"
FUSION="/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master"
#for bslmm,this is a workaround because GEMMA requires results to go into an output subdirectory 
ln -s ./ output
RED='\033[0;31m'
NC='\033[0m' # No Color

PRE="Esophagus_Gastroesophageal_Junction"
PRE="Stomach"
PRE="Whole_Blood"
PRE="Adipose_Visceral_Omentum"


NR=10

PRE=$1
NR=${SLURM_ARRAY_TASK_ID}

#prepare folders and files ---

mkdir tmp/$PRE.$NR

rm -f HSQ/$PRE.$NR.hsq HSQ/TSPEC.$PRE.$NR.hsq
zcat expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz |head -n 1 |tr '\t' '\n' >expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz.HEADER


COVAR="expression_covariates/$PRE.v8.EUR.covariates.txt"
#generate plink format covariate file

ncov=$(head -n 1 $COVAR | awk '{print NF}')
COVAR1="expression_covariates/$PRE.v8.EUR.covariates1.txt"
COVAR2="expression_covariates/$PRE.v8.EUR.covariates_plink.txt"
for ((i=1;i<=$ncov;i++));do cut -f"$i" $COVAR | paste -s;done >$COVAR1
cut -f1 $COVAR1 |paste - $COVAR1 >$COVAR2
rm -f $COVAR1

#main loop on all the genes----
#the position file
posfile=WEIGHTS/$PRE.pos
if [[ -f $posfile ]];then rm $posfile; fi
#zcat expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz | tail -n+2 | awk -v i=$NR 'NR > (i-1)*100 && NR <= i*100' |  while read PARAM; do
#put results in one folder
ngene=0
zcat expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz | tail -n+2 |  while read PARAM; do
ngene=$(( ngene + 1 ))
tmp=$(( $ngene % 200 ))
if [[ $tmp -eq 0 ]];then echo -e "${RED}Working on $ngene genes...${NC}";fi
#echo $PARAM >test.txt
#PARAM=`cat test.txt`

GNAME=`echo $PARAM | awk '{ print $4 }'`
CHR=`echo $PARAM | awk '{ print $1 }'`
P0=`echo $PARAM | awk '{ p=$2 - 500e3; if(p<0) p=0; print p; }'`
P1=`echo $PARAM | awk '{ print $3 + 500e3 }'`

OUT="tmp/$PRE.$NR/$PRE.$GNAME"

echo $GNAME $CHR $P0 $P1
echo $PRE/$PRE.$GNAME.wgt.RDat $CHR $P0 $P1 >> $posfile
echo $PARAM | tr ' ' '\n' | paste expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz.HEADER expression_matrices/$PRE.v8.EUR.normalized_expression.bed.gz.HEADER - | tail -n+5 > $OUT.pheno
rm -f $OUT.bed
$PLINK --silent --bfile /fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_HM3 --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out $OUT --pheno $OUT.pheno --keep $OUT.pheno --maf 0.0001
if [ ! -f $OUT.bed ]; then
continue
fi

if [[ ! -d "WEIGHTS/$PRE" ]];then mkdir WEIGHTS/$PRE;fi

FINAL_OUT="WEIGHTS/$PRE/$PRE.$GNAME"
Rscript $FUSION/FUSION.compute_weights.R \
--bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA \
--models blup,lasso,top1,enet --covar $COVAR2 --hsq_p 1.0 --PATH_plink /fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

if [ ! -f $FINAL_OUT.wgt.RDat ]; then continue; fi

cat $FINAL_OUT.hsq | awk -vw="$PRE $w" '{ print w,$0 }' >> HSQ/$PRE.$NR.hsq
rm -f $FINAL_OUT.hsq

# # Get best eQTL for this gene
# META=`cat GTEx.v7.meta.id | awk -vg=$GNAME '$3 == g { print $2 }'`
# if [ "$META" == "" ]; then
# continue
# fi
# cat GTEx.v7.meta.covar | awk -vm=$META '$1 == m' | tr '\t' '\n' | paste GTEx.v7.meta.covar.HEADER - | tail -n+2 | sort -k1,1 \
# | join -1 1 -2 1 - <(cut -f2- $COVAR | sort -k1,1) | awk '$2 != "NA" { print $1,$0 }' | tr -s ' ' '\t' > $OUT.meta.covar
# 
# FINAL_OUT="WEIGHTS/$PRE.TSPEC/$PRE.$GNAME"
# 
# Rscript $FUSION/FUSION.compute_weights.R \
# --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet --covar $OUT.meta.covar --resid
# cat $FINAL_OUT.hsq | awk -vw="$PRE $w" '{ print w,$0 }' >> HSQ/TSPEC.$PRE.$NR.hsq
# rm -f $FINAL_OUT.hsq

#rm $OUT.*
done

#rm -fr tmp/$PRE.$NR
#touch HSQ/$PRE.$NR.done

# #change the sumstat file, too slow
#  ml Anaconda3/2019.07
#  cd /fh/fast/dai_j/CancerGenomics/Tools/ldsc
#  source activate ldsc
#  
#  cd -
# transform_sumstats () {
#   local Ncas="$1"
#   local Ncon="$2"
#   local input="$3"
#   local out="$4"
#   /fh/fast/dai_j/CancerGenomics/Tools/ldsc/munge_sumstats.py \
# --out $out \
# --merge-alleles /fh/fast/dai_j/CancerGenomics/Tools/ldsc/w_hm3.snplist \
# --N-cas $Ncas \
# --N-con $Ncon \
# --sumstats $out 
# }
# 
# transform_sumstats 3925 2185 /fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt ../result/BEEA_BEACON_autosomes.sumstats
# /fh/fast/dai_j/CancerGenomics/Tools/ldsc/munge_sumstats.py \
# --out ../result/BEEA_BEACON_autosomes_N_hm3.sumstat \
# --merge-alleles /fh/fast/dai_j/CancerGenomics/Tools/ldsc/w_hm3.snplist \
# --N 6110 \
# --sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N_hm3.txt

#create sumstats using run_FUSION.py

#Run assoc
fussion_assoc(){
  local sumstats="$1"
  local weights="$2"
  local chr="$3"
  local outprefix="$4"
  Rscript $fusiondir/FUSION.assoc_test.R \
--sumstats $sumstats \
--weights $weights \
--weights_dir /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/ \
--ref_ld_chr /fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/LDREF/1000G.EUR. \
--chr $chr \
--out $outprefix.$chr.dat
}
# for chr in {1..21}
# do
#   fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N_hm3.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/GTEx.Whole_Blood.pos $chr /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BEACON_BEEA_Blood &
# done
# for chr in {1..22}
# do
#   fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N_hm3.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Gastroesophageal_Junction.P01.pos $chr /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BEACON_BEEA_Junction &
# done

for chr in {1..22}
do
  fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Gastroesophageal_Junction.P01.pos $chr /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction &
done

run_fussion_assoc(){
  local sumstats="$1"
  local posfile="$2"
  local outprefix="$3"
  declare -a pids
  for chr in {1..22}
  do
    echo $chr
    fussion_assoc $sumstats $posfile $chr $outprefix &
    pids[${chr}]=$!
  done
  for pid in ${pids[*]}; do
    wait $pid
  done
  echo "DONE"
}

run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Gastroesophageal_Junction.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Junction
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/GTEx.Whole_Blood.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Blood
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Stomach.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Stomach
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Mucosa.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Mucosa
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Muscularis.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Muscularis
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BEEA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Adipose_Visceral_Omentum.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BEEA_Adiplose

run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Gastroesophageal_Junction.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Junction
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/GTEx.Whole_Blood.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Blood
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Stomach.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Stomach
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Mucosa.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Mucosa
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Muscularis.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Muscularis
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_EA.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Adipose_Visceral_Omentum.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_EA_Adiplose

run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Gastroesophageal_Junction.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Junction
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/GTEx.Whole_Blood.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Blood
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Stomach.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Stomach
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Mucosa.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Mucosa
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Esophagus_Muscularis.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Muscularis
run_fussion_assoc /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_BE.sumstats /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/WEIGHTS/Adipose_Visceral_Omentum.P01.pos /fh/fast/dai_j/BEACON/BEACON_GRANT/result/fusion/result/BCA_BE_Adiplose

