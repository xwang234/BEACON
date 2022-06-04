
library(data.table)
ped=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.ped")))
bim=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.bim")))
traw=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.traw")))
raw=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.raw")))
map=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.map")))
fam=as.data.frame(fread(paste0("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/test5.fam")))

tmp1=raw[,7:ncol(raw)]
tmp2=ped[,7:ncol(ped)]
for (j in 1:nrow(bim))
{
  idx1=which(tmp2[,j]==paste0(bim$V6[j]," ",bim$V6[j]))
  tmp2[idx1,j]=0
  idx1=which(tmp2[,j]==paste0(bim$V6[j]," ",bim$V5[j]))
  tmp2[idx1,j]=1
  idx1=which(tmp2[,j]==paste0(bim$V5[j]," ",bim$V6[j]))
  tmp2[idx1,j]=1
  idx1=which(tmp2[,j]==paste0(bim$V5[j]," ",bim$V5[j]))
  tmp2[idx1,j]=2
}
all(tmp1==tmp2)