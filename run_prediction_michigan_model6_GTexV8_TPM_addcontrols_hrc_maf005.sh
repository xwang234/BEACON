#!/usr/bin/env bash

prefix="GTExV8mucosadata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_mucosa_HRC_maf005 #cv 10 times
outprefix="dist500K_GTEx_mucosa_HRC_MAF005" #cv 100 times

prefix="GTExV8junctiondata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_junction_HRC_maf005
outprefix=dist500K_GTEx_junction_HRC_MAF005

prefix="GTExV8muscularisdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_muscularis_HRC_maf005
outprefix=dist500K_GTEx_muscularis_HRC_MAF005

prefix="GTExV8stomachdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_stomach_HRC_maf005
outprefix=dist500K_GTEx_stomach_HRC_MAF005

prefix="GTExV8blooddata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_blood_HRC_maf005
outprefix=dist500K_GTEx_blood_HRC_MAF005

prefix="GTExV8adiposedata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
#outprefix=dist500K_GTEx_adipose_HRC_maf005
outprefix=dist500K_GTEx_adipose_HRC_MAF005

rdata=$prefix.RData
if [[ -f ../result/$rdata ]]
then
  salloc -t 6-1 --constraint=gizmok -n 71 mpirun -n 1 Rscript ./prediction_michigan_model6_GTexV8_TPM_addcontrols_hrc_maf005.R $rdata $outprefix &
  echo $rdata
fi
  
#rdata=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8mucosadata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
#rdata=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8muscularisdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")


prefix="GTExV8mucosadata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix="dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor" #cv 100 times

prefix="GTExV8junctiondata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_junction_HRC_MAF005_rmhighcor

prefix="GTExV8muscularisdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor

prefix="GTExV8stomachdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_stomach_HRC_MAF005_rmhighcor

prefix="GTExV8blooddata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_blood_HRC_MAF005_rmhighcor

prefix="GTExV8adiposedata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_adipose_HRC_MAF005_rmhighcor

rdata=$prefix.RData
if [[ -f ../result/$rdata ]]
then
  salloc -t 6-1 --constraint=gizmok -n 81 mpirun -n 1 Rscript ./prediction_michigan_model6_GTexV8_TPM_addcontrols_rmhighcor_hrc_maf005.R $rdata $outprefix &
  echo $rdata
fi

#include all covar----
prefix="GTExV8mucosadata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix="dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor_allcovar" #cv 100 times

prefix="GTExV8junctiondata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar

prefix="GTExV8muscularisdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor_allcovar

prefix="GTExV8stomachdata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_stomach_HRC_MAF005_rmhighcor_allcovar

prefix="GTExV8blooddata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_blood_HRC_MAF005_rmhighcor_allcovar

prefix="GTExV8adiposedata_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction"
outprefix=dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar

rdata=$prefix.RData
if [[ -f ../result/$rdata ]]
then
  salloc -t 6-1 --constraint=gizmok -n 101 mpirun -n 1 Rscript ./prediction_michigan_model6_GTexV8_TPM_addcontrols_rmhighcor_hrc_allcovr_maf005.R $rdata $outprefix &
  echo $rdata
fi
