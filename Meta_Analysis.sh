#METAL commands used for meta analysis
#METAL program
#/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#GWAS data were processed in PRS_LDpred.R

#To add sample size
#sed 's/$/\tgreen/' ./test.txt
#sed '0,/Apple/{s/Apple/Banana/}' input_filename
#Bonn:1037 BE, 1609 EA, 2646 BEEA, 3537 controls
#Cambridge:873 BE, 995 EA, 1868 BEEA, 3408 controls
#Oxford:1851 BE, 3496 controls
#Beacon:2185 controls, 1512 EA, 2185 controls 3925 EA/BE
#Meta BEACON 2406 BE, 1508 EA, 3914 EA/BE, 6718 controls
add_samplesize()
{
  local N="$1"
  local infile="$2"
  local outfile="$3"
  sed "s/$/\ $N/" $infile > $outfile
  sed -i "0,/$N/{s/$N/N/}" $outfile
}
add_samplesize 4574 /fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_N.txt
add_samplesize 4281 /fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Cambridge_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt
add_samplesize 9124 /fh/fast/dai_j/BEACON/Meta_summary_stat/BE_BEACON_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_BEACON_autosomes_N.txt
add_samplesize 5374 /fh/fast/dai_j/BEACON/Meta_summary_stat/BE_oxford_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_oxford_autosomes_N.txt
add_samplesize 5146 /fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Bonn_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_N.txt
add_samplesize 4403 /fh/fast/dai_j/BEACON/Meta_summary_stat/EA_Cambridge_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt
add_samplesize 8226 /fh/fast/dai_j/BEACON/Meta_summary_stat/EA_BEACON_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_BEACON_autosomes_N.txt
add_samplesize 6183 /fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt
add_samplesize 5276 /fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt
add_samplesize 10632 /fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt


# #BE studies-----, use p-value
# #load the first input file
# MARKER SNP
# ALLELE non_effect_allele effect_allele
# EFFECT BETA
# PVALUE P
# WEIGHT N
# PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_N.txt
# 
# #load the second input file
# MARKER rsid
# ALLELE non-effect-allele effect-allele
# EFFECT beta
# PVALUE pvalue
# WEIGHT N
# PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_oxford_autosomes_N.txt
# 
# #load the third input file
# MARKER SNP
# ALLELE non_effect_allele effect_allele
# EFFECT beta
# PVALUE P
# WEIGHT N
# PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt
# 
# OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_METAANALYSIS_BE_ .tbl
# MINWEIGHT 10000
# ANALYZE 

#Use All SNPs from studies------
#BE studies-----,use standard error
#load the first input file
MARKER SNP
ALLELE non_effect_allele effect_allele
PVALUE P
EFFECT BETA
STDERRLABEL SE
SCHEME STDERR
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_N.txt

#load the second input file
MARKER rsid
ALLELE non-effect-allele effect-allele
PVALUE pvalue
EFFECT beta
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_oxford_autosomes_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_METAANALYSIS_BE_ .tbl
MINWEIGHT 10000
ANALYZE 

QUIT

#BEEA studies---, standard error
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Cambridge_METAANALYSIS_BEEA_ .tbl
MINWEIGHT 10000
ANALYZE 

QUIT

#EA studies---, standard error
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Cambridge_METAANALYSIS_EA_ .tbl
MINWEIGHT 6000
ANALYZE 

QUIT

#Use common SNPs from studies------
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#BE studies-----,use standard error
#load the first input file
MARKER SNP
ALLELE non_effect_allele effect_allele
PVALUE P
EFFECT BETA
STDERRLABEL SE
SCHEME STDERR
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_comsnp_N.txt

#load the second input file
MARKER rsid
ALLELE non-effect-allele effect-allele
PVALUE pvalue
EFFECT beta
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_oxford_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_METAANALYSIS_BE_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#BEEA studies---, standard error
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_comsnp_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Cambridge_METAANALYSIS_BEEA_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_comsnp_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Cambridge_METAANALYSIS_EA_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_BEACON_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_EA_comsnp .tbl

ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#BE studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_BEACON_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BE_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#BEEA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BEEA_comsnp .tbl
ANALYZE 

QUIT

#change strand of BEACON
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#BEEA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE effect_allele non_effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BEEA_comsnp_new .tbl
ANALYZE 

QUIT


#use BEACON +BONN
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_METAANALYSIS_EA_comsnp .tbl
ANALYZE 

QUIT

#use BEACON +BONN, high quality SNP
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_BB_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_BB_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_METAANALYSIS_EA_info_comsnp .tbl
ANALYZE 

QUIT

#use Cambridge +BONN, high quality SNP
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_BC_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_BC_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Cambridge_METAANALYSIS_EA_info_comsnp .tbl
ANALYZE 

QUIT

#use BEACON +BONN, geneotyped, beacon without amos
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_genotyped_BB_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_genotyped_BB_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_METAANALYSIS_EA_genotypedNoAmos_comsnp .tbl
ANALYZE 

QUIT

#use BEACON +BONN, imputed, beacon without amos
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_imp_BB_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_imp_BB_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_METAANALYSIS_EA_impNoAmos_comsnp .tbl
ANALYZE 

QUIT

#use Discovery +BONN
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_genotyped_BD_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_genotyped_BD_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Discovery_Bonn_METAANALYSIS_EA_genotyped_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_genotyped_BD_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Discovery_genotyped_BD_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Discovery_Bonn_METAANALYSIS_BEEA_genotyped_comsnp .tbl
ANALYZE 

QUIT

#use Discovery +BONN
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_BD_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_BD_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Discovery_Bonn_METAANALYSIS_EA_comsnp .tbl
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
#EA studies---, standard error, add Beacon
#load the secnod input file
MARKER SNP
ALLELE non_effect_allele effect_allele
SCHEME STDERR
STDERRLABEL SE
EFFECT BETA
PVALUE P
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_BD_autosomes_comsnp_N.txt

#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Discovery_BD_autosomes_comsnp_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Discovery_Bonn_METAANALYSIS_BEEA_comsnp .tbl
ANALYZE 

QUIT


#BEACON+CAMBRIDGE
/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_BEACON_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_EA_ .tbl
MINWEIGHT 10000
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_BEACON_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BE_ .tbl
MINWEIGHT 10000
ANALYZE 

QUIT

/fh/fast/dai_j/CancerGenomics/Tools/METAL/build/bin/metal
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt
#load the third input file
MARKER SNP
ALLELE non_effect_allele effect_allele
EFFECT beta
PVALUE P
STDERRLABEL se
WEIGHT N
PROCESS /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt

OUTFILE /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEACON_Cambridge_METAANALYSIS_BEEA_ .tbl
MINWEIGHT 12000
ANALYZE 

QUIT
