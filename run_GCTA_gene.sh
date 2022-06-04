organs=("junction" "muscularis" "adipose" "blood" "mucosa" "stomach")
resultfolder=/fh/fast/dai_j/BEACON/BEACON_GRANT/result
# for organ in ${organs[@]}
# do
#  echo $organ
#  ./GCTA_gene.R $organ "mgrm" > $resultfolder/GCTA_${organ}_mgrm_peer.txt 2>&1 &
# done
# 
# for organ in ${organs[@]}
# do
#  echo $organ
#  ./GCTA_gene.R $organ "1grm" > $resultfolder/GCTA_${organ}_1grm_peer.txt 2>&1 &
# done

#MPI
for organ in ${organs[@]}
do
 echo $organ
 salloc -t 2-1 --constraint=gizmok -n 45 mpirun -n 1 Rscript ./GCTA_gene.R $organ 1grm > $resultfolder/GCTA_${organ}_1grm.txt 2>&1 &
done