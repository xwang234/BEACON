chr=unlist(strsplit(selsnps[1],":"))[1]
#read data from 1000 genome v3 EUR
thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
#refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
refbimstr1=paste0(refbim$V1,":",refbim$V4,"_",refbim$V5,"_",refbim$V6)
refbimstr2=paste0(refbim$V1,":",refbim$V4,"_",refbim$V6,"_",refbim$V5)
refbimstr0=paste0(refbim$V1,":",refbim$V4)

tmp=intersect(selsnpstr0,BE_Bonnsummarydatstr)
missingsnps=selsnpstr0[which(!selsnpstr0 %in% BE_Bonnsummarydatstr)]
avai_refsnps=refbimstr0[which(refbimstr0 %in% missingsnps)]
refreplacesnps=data.frame(snpstr0=avai_refsnps,replacestr0=NA,cor=NA)
refreplacesnps=data.frame(snpstr0=selsnpstr0,replacestr0=NA,cor=NA,missing=!selsnpstr0 %in% BE_Bonnsummarydatstr)
for (k in 1:nrow(refreplacesnps))
{
  idx=which(refbimstr0==refreplacesnps$snpstr0[k])
  if (length(idx)>0)
  {
    tmp1=cor(t(refraw[,7:ncol(refraw)]),unlist(refraw[idx,7:ncol(refraw)]))
    
    tmp2=order(tmp1[,1],decreasing = T)
    
    tmp2=tmp2[tmp2!=idx]
    idx1=which(refbimstr0 %in% BE_Bonnsummarydatstr)
    tmp3=tmp2[tmp2 %in% idx1]
    
    refreplacesnps$replacestr0[k]=refbimstr0[tmp3[1]]
    refreplacesnps$cor[k]=tmp1[tmp3[1],1]
  }
}

sum(refreplacesnps$replacestr0 %in% BE_Bonnsummarydatstr)
refreplacesnps$inBonn=refreplacesnps$replacestr0 %in% BE_Bonnsummarydatstr

#use the closest snps in BonnBE
refreplacesnps1=data.frame(snpstr0=missingsnps,replacestr0=NA,dist=NA)
for (k in 1:nrow(refreplacesnps1))
{
  tmp1=as.integer(unlist(strsplit(refreplacesnps1$snpstr0[k],":"))[2])
  idx=which(BE_Bonnsummarydat$chr==chr)
  idx1=which.min(abs(BE_Bonnsummarydat$position[idx]-tmp1))
  refreplacesnps1$replacestr0[k]=BE_Bonnsummarydatstr[idx[idx1]]
  refreplacesnps1$dist[k]=min(abs(BE_Bonnsummarydat$position[idx]-tmp1))
  
}

tmp=intersect(selsnpstr0,BE_Bonnsummarydatstr)
#tmp=c(tmp,refreplacesnps$replacestr0[refreplacesnps$inBonn & refreplacesnps$cor>=0.7])
tmp=c(tmp,refreplacesnps1$replacestr0[refreplacesnps1$dist<=50])
idx1=match(tmp,BE_Bonnsummarydatstr) #use id from summary data
BonnBEsnpstr=BE_Bonnsummarydatstr1[idx1]

Z1=Z2=NULL
tmp1=intersect(BonnBEsnpstr,refbimstr1)
tmp2=intersect(BonnBEsnpstr,refbimstr2)
idx1=match(tmp1,refbimstr1)
if (length(tmp1)>0)
{
  Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
  rownames(Z1)=tmp1
}
idx2=match(tmp2,refbimstr2)
if (length(tmp2)>0)
{
  Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
  #change id
  rownames(Z2)=tmp2
  Z2=2-Z2
}
Z3=as.data.frame(t(rbind(Z1,Z2)))
tmp1=intersect(selsnpstr0,BE_Bonnsummarydatstr)
tmp=colnames(Z3)
if (length(tmp)<length(tmp1)) warning(paste0(genename,",BE:there are ",length(tmp1)-length(tmp)," not found in 1000G ref."))

if (length(tmp)>1)
{
  valsnps1 <- colnames(Z3)
  selsnpstr1 <- colnames(Z3)
  val1 <- BE_Bonnsummarydat[match(tmp,BE_Bonnsummarydatstr1),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z3[,match(selsnpstr1,colnames(Z3))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[3,1] <- length(tmp)
  result[3,2] <- mean(val1$P<0.05)
  result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
}else if (length(tmp)==1)
{
  result[3,3]=BE_Bonnsummarydat$P[which(BE_Bonnsummarydatstr==tmp)]
}

rownames(result)[3]="BE_Bonn"
result[3,]
