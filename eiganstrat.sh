
outfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result
plink=/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink
cd $outfolder

$plink --bfile $outfolder/bca_filtered_30Nov2018 --recode tab --snps-only --biallelic-only strict --maf 0.05 --make-bed --out $outfolder/bca_filtered_10Jan2019
export PATH=$PATH:/fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin
cp ./bca_filtered_10Jan2019.fam ./bca_filtered_10Jan2019.pedind
cp ./bca_filtered_10Jan2019.bim ./bca_filtered_10Jan2019.pedsnp
/fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin/convertf -p bca_filtered_10Jan2019.par.PED.EIGENSTRAT
perl /fh/fast/stanford_j/Xiaoyu/Tools/EIG-7.2.1/bin/smartpca.perl -i bca_filtered_10Jan2019.eigenstratgeno -a bca_filtered_10Jan2019.snp -b bca_filtered_10Jan2019.ind -k 15 -o bca_filtered_10Jan2019.pca -e bca_filtered_10Jan2019.eval -l bca_filtered_10Jan2019.eigpca.log -p bca_filtered_10Jan2019.plot -s 6

