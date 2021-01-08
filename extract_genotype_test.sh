
#!/usr/bin/env bash
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC
chr=19
tabixreg=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC/prediction_snps_tabix_test.txt
tabix -H /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz >$outfolder/chr${chr}.test.dose.vcf
xargs -a $tabixreg -I {} tabix -f /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr${chr}.dose.vcf.gz {} >> $outfolder/chr${chr}.test.dose.vcf
$plink --vcf $outfolder/chr${chr}.test.dose.vcf --const-fid 0 --snps-only --make-bed --out $outfolder/chr${chr}_test
$plink --bfile $outfolder/chr${chr}_test --recodeA --out $outfolder/chr${chr}_test
