#!/usr/bin/env bash
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink

/fh/fast/dai_j/CancerGenomics/Tools/GenGen-1.0.1/convert_bim_allele.pl \
  /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015.bim \
  /fh/fast/dai_j/CancerGenomics/Tools/GenGen-1.0.1/gengenlib/hh550_610.snptable \
  -intype  ilmnab -outtype dbsnp -outfile ../result/GWAS/test.dbsnp.bim
#BE and control
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 \
  --keep ../result/GWAS/BE_CO_selectedsamples_plink.txt --make-bed \
  --out ../result/GWAS/BE_CO_19Feb2015
#run update_fam() in perform_GWAS.R
  
$plink --bfile  ../result/GWAS/BE_CO_19Feb2015 --covar  ../result/GWAS/BE_CO_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ../result/GWAS/BE_CO_19Feb2015





#EA and control
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 \
  --keep ../result/GWAS/EA_CO_selectedsamples_plink.txt --make-bed \
  --out ../result/GWAS/EA_CO_19Feb2015
#run update_fam() in perform_GWAS.R
  
$plink --bfile  ../result/GWAS/EA_CO_19Feb2015 --covar  ../result/GWAS/EA_CO_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ../result/GWAS/EA_CO_19Feb2015


#BEEA and control
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 \
  --keep ../result/GWAS/BEEA_CO_selectedsamples_plink.txt --make-bed \
  --out ../result/GWAS/BEEA_CO_19Feb2015
#run update_fam() in perform_GWAS.R
$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 \
  --keep ../result/GWAS/BEEA_CO_selectedsamples_plink.txt --recode \
  --out ../result/GWAS/BEEA_CO_19Feb2015
  
$plink --bfile  ../result/GWAS/BEEA_CO_19Feb2015 --covar  ../result/GWAS/BEEA_CO_selectedsamples_pheno_plink.txt \
  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ../result/GWAS/BEEA_CO_19Feb2015
  
#output
#   .assoc.linear, .assoc.logistic (multi-covariate association analysis report)
# Produced by --linear/--logistic.
# 
# A text file with a header line, and T lines per variant typically with the following nine fields (where T is normally the number of terms, but the 'genotypic' and 'hethom' modifiers and the --tests flag can change this):
# 
# CHR	Chromosome code. Not present with 'no-snp' modifier.
# SNP	Variant identifier. Not present with 'no-snp'.
# BP	Base-pair coordinate. Not present with 'no-snp'.
# A1	Allele 1 (usually minor). Not present with 'no-snp'.
# TEST	Test identifier
# NMISS	Number of observations (nonmissing genotype, phenotype, and covariates)
# 'BETA'/'OR'	Regression coefficient (--linear, "--logistic beta") or odds ratio (--logistic without 'beta')
# STAT	T-statistic
# P	Asymptotic p-value for t-statistic
# If --ci 0.xy has also been specified, the following three fields are inserted before 'STAT':
# 
# SE	Standard error of beta (log-odds) estimate
# Lxy	Bottom of xy% symmetric approx. confidence interval
# Hxy	Top of xy% approx. confidence interval

awk '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }' ../result/GWAS/BEEA_CO_19Feb2015.ped >test1.ped
#compute LD score
for chr in {1..22}
do 
  echo $chr
  $plink --bfile ../result/GWAS/BE_CO_19Feb2015 --chr $chr --make-bed --out ../result/GWAS/BE_w_ld_chr/$chr
  $plink --bfile ../result/GWAS/BE_w_ld_chr/$chr --cm-map /fh/fast/dai_j/CancerGenomics/Tools/database/1000genome/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt $chr \
  --make-bed --out ../result/GWAS/BE_w_ld_chr/$chr
done
for chr in {1..22}
do
  echo $chr
  /fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
	--bfile ../result/GWAS/BE_w_ld_chr/$chr --l2  --ld-wind-cm 1 \
	--out ../result/GWAS/BE_w_ld_chr/$chr
done

for chr in {1..22}
do 
  echo $chr
  $plink --bfile ../result/GWAS/EA_CO_19Feb2015 --chr $chr --make-bed --out ../result/GWAS/EA_w_ld_chr/$chr
  $plink --bfile ../result/GWAS/EA_w_ld_chr/$chr --cm-map /fh/fast/dai_j/CancerGenomics/Tools/database/1000genome/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt $chr \
  --make-bed --out ../result/GWAS/EA_w_ld_chr/$chr
done
for chr in {1..22}
do
  echo $chr
  /fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
	--bfile ../result/GWAS/EA_w_ld_chr/$chr --l2  --ld-wind-cm 1 \
	--out ../result/GWAS/EA_w_ld_chr/$chr
done

#compute h2
gzip ../result/GWAS/BEEA.sumstats
ml anaconda2
/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
--h2 ../result/GWAS/BEEA.sumstats.gz \
--ref-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/ \
--w-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/ \
--out ../result/GWAS/BEEA_h2_1000g

#gzip ../result/GWAS/BE.sumstats
/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
--h2 ../result/GWAS/BE.sumstats.gz \
--ref-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/  \
--w-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/  \
--samp-prev 0.6008042 \
--pop-prev 0.016 \
--out ../result/GWAS/BE_h2_1000g
/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
--h2 ../result/GWAS/BE.sumstats.gz \
--ref-ld-chr ../result/GWAS/BE_w_ld_chr/ \
--w-ld-chr ../result/GWAS/BE_w_ld_chr/ \
--samp-prev 0.6008042 \
--pop-prev 0.016 \
--out ../result/GWAS/BE_h2

gzip ../result/GWAS/BE_hm3.sumstats
/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py \
--h2 ../result/GWAS/BE_hm3.sumstats.gz \
--ref-ld-chr ../result/GWAS/BE_w_ld_chr/ \
--w-ld-chr ../result/GWAS/BE_w_ld_chr/ \
--samp-prev 0.6008042 \
--pop-prev 0.016 \
--out ../result/GWAS/BE_h2_hm3

#gzip ../result/GWAS/EA.sumstats

$plink --bfile /fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015 \
  --geno 0.03 --maf 0.01 --mind 0.03 --hwe 0.0001 --make-bed \
  --out ../result/bca_filtered_puya

$plink --bfile  ../result/bca_filtered_puya --genome --indep-pairwise 50 5 0.2 --min 0.2 --out ../result/bca_filtered_puya_ibd
$plink --bfile  ../result/bca_filtered_puya --remove ../result/bca_filtered_puya_ibd.genome --make-bed --out ../result/bca_filtered_puya1


