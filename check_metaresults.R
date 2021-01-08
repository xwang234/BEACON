#!/usr/bin/env Rscript
library(data.table)
BEEAmeta=as.data.frame(fread("../result/Bonn_Cambridge_METAANALYSIS_BEEA_comsnp1.tbl"))
BEEAmeta[which.min(BEEAmeta$`P-value`),]
snps=c("rs7255","rs2464469","rs17451754","rs17749155","rs10108511")
for (i in 1:length(snps))
{
  idx=which(BEEAmeta$MarkerName==snps[i])
  if (length(idx)>0)
  {
    print(BEEAmeta[idx,])
  }
}
# MarkerName Allele1 Allele2  Effect StdErr   P-value Direction
# 6426356     rs7255       t       c -0.1225 0.0298 4.038e-05        --
#   MarkerName Allele1 Allele2 Effect StdErr P-value Direction
# 7869655  rs2464469       a       g  0.085 0.0292 0.00357        ++
#   MarkerName Allele1 Allele2 Effect StdErr   P-value Direction
# 1012609 rs17451754       a       g 0.1986 0.0418 2.054e-06        ++
#   MarkerName Allele1 Allele2  Effect StdErr  P-value Direction
# 10488654 rs17749155       a       g -0.1358 0.0423 0.001312        --
#   MarkerName Allele1 Allele2  Effect StdErr   P-value Direction
# 10524281 rs10108511       t       c -0.1112 0.0286 0.0001021        --

BEEAmeta=as.data.frame(fread(input="../result/BEEA_Bonn_Cambridge_metastat.txt.gz"))

BEEABonn=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt"))
BEEACambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt"))
for (i in 1:length(snps))
{
  idx=which(BEEABonn$SNP==snps[i])
  if (length(idx)>0)
  {
    print(BEEABonn[idx,c("SNP","non_effect_allele","effect_allele","OR","P","BETA","SE")])
  }
}
# SNP non_effect_allele effect_allele       OR          P     BETA        SE
# 834542 rs7255                 T             C 0.893153 0.00156999 -0.13035 0.0412311
# SNP non_effect_allele effect_allele       OR          P      BETA        SE
# 8185061 rs2464469                 G             A 0.865184 0.00631877 -0.107429 0.0393399
# SNP non_effect_allele effect_allele       OR           P      BETA        SE
# 4797236 rs17451754                 G             A 0.794338 2.94989e-05 -0.228502 0.0547009
# SNP non_effect_allele effect_allele      OR         P     BETA        SE
# 5004555 rs17749155                 G             A 1.13041 0.0196747 0.138512 0.0593834
# SNP non_effect_allele effect_allele      OR          P      BETA        SE
# 5010918 rs10108511                 T             C 0.88556 0.00603387 -0.104454 0.0380396

for (i in 1:length(snps))
{
  idx=which(BEEACambridge$SNP==snps[i])
  if (length(idx)>0)
  {
    print(BEEACambridge[idx,c("SNP","non_effect_allele","effect_allele","P","beta","se")])
  }
}
# SNP          P      beta        se
# 1262699 rs7255 0.00843291 -0.113909 0.0432573
# SNP        P       beta        se
# 10838875 rs2464469 0.185203 -0.0575936 0.0434628
# SNP         P     beta        se
# 5727883 rs17451754 0.0153613 -0.15656 0.0649404
# SNP         P    beta        se
# 6801549 rs17749155 0.0274039 0.13301 0.0601567
# SNP          P      beta        se
# 7080064 rs10108511 0.00571086 -0.120082 0.0434696

#Allele1 is non-effect allele

BeaconBEEA=as.data.frame(fread("../result/Beacon_BEEA.bim"))

