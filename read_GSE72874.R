library(data.table)
GE=data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GSE72874/GSE72874-GPL10558_series_matrix.txt",skip=62,header = T))
rownames(GE)=GE$ID_REF
GE=GE[,-1]
dim(GE)
#[1] 47208    70
ME=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GSE72874/GSE72874-GPL13534_series_matrix.txt",skip=63,sep="\t",header = T))
rownames(ME)=ME$ID_REF
ME=ME[,-1]
