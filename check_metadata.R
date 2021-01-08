
library(data.table)
BonnBE=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Bonn_autosomes.txt"))
quantile(BonnBE$freq_non_effect_allele_controls)
# 0%       25%       50%       75%      100% 
# 0.0010001 0.6297340 0.8799011 0.9761530 0.9989999 
quantile(BonnBE$freq_non_effect_allele_cases)
# 0%       25%       50%       75%      100% 
# 0.0000251 0.6298894 0.8797990 0.9760573 0.9999995
quantile(BonnBE$freq_effect_allele_controls)
# 0%        25%        50%        75%       100% 
# 0.00100010 0.02384697 0.12009900 0.37026600 0.99900000 
quantile(BonnBE$info)
# 0%      25%      50%      75%     100% 
# 0.400001 0.710738 0.919692 0.979656 1.000000

BonnBEEA=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Bonn_autosomes.txt"))
sum(BonnBE$SNP %in% BonnBEEA$SNP)
CambridgeBE=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BE_Cambridge_autosomes.txt"))
CambridgeBEEA=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt"))
quantile(CambridgeBE$all_maf)
# 0%        25%        50%        75%       100% 
# 0.00100001 0.00614612 0.04669250 0.21074800 0.50000000 
quantile(CambridgeBE$info)
# 0%      25%      50%      75%     100% 
# 0.400000 0.780014 0.962504 0.994745 1.000000
quantile(CambridgeBE$cases_maf)
# 0%        25%        50%        75%       100% 
# 0.00000000 0.00614204 0.04650170 0.21080000 0.54435100 
quantile(CambridgeBEEA$all_maf)
# 0%        25%        50%        75%       100% 
# 0.00100001 0.00615275 0.04669620 0.21077400 0.50000000
chr=1
bimdat=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr",chr,"_all_noambiguous.bim")))
