
rm(list=ls())
if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data//bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  for (i in 1:ncol(sampletable))
  {
    idx=which(sampletable[,i]==-9)
    sampletable[idx,i]=NA
  }
}

fam=read.table("../result/bca_filtered_04Apr2019.fam",stringsAsFactors = F)
sum(sampletable$phenoBE_bc==2 & sampletable$phenoEA_bc==2,na.rm=T) #0

all(sampletable$localid %in% fam$V2) #T
sum(fam$V2 %in% sampletable$localid) #8208
pheno=data.frame(famid=fam$V1,subjid=fam$V2,pheno=-9,stringsAsFactors = F)
idx1=which(pheno$subjid %in% sampletable$localid[sampletable$phenoBE_bc==1 | sampletable$phenoEA_bc==1]) #control
idx2=which(pheno$subjid %in% sampletable$localid[sampletable$phenoBE_bc==2]) #BE case
idx3=which(pheno$subjid %in% sampletable$localid[sampletable$phenoEA_bc==2]) #EA case
length(intersect(idx1,idx2))
length(intersect(idx3,idx2))
length(intersect(idx1,idx3))
pheno$pheno[idx1]=0
pheno$pheno[c(idx2,idx3)]=1
table(pheno$pheno)
# -9    0    1 
# 1244 2185 5810 
write.table(pheno,file="../result/bca_filtered_04Apr2019.GCTA.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx2),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BE.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx3),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.EA.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx2,idx3),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.pheno",sep=" ",col.names = F,row.names = F,quote=F)
pheno$EA=-9
pheno$EA[c(idx1,idx2)]=0
pheno$EA[idx3]=1
write.table(pheno[c(idx2,idx3),c(1,2,4)],file="../result/bca_filtered_04Apr2019.GCTA.BE.EA.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx2),1:2],file="../result/GCTA_CO_BE_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx3),1:2],file="../result/GCTA_CO_EA_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx2,idx3),1:2],file="../result/GCTA_BE_EA_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[c(idx1,idx2,idx3),1:2],file="../result/GCTA_CO_BEEA_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
#subseting recurrent_HB_RF
idx=match(pheno$subjid,sampletable$localid)
pheno$recurrent_HB_RF=NA
pheno$recurrent_HB_RF=sampletable$recurrent_HB_RF[idx]
idx4=which(pheno$recurrent_HB_RF==1)
write.table(pheno[intersect(c(idx1,idx2),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx3),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx2,idx3),idx4),c(1,2,4)],file="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2,idx3),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2),idx4),1:2],file="../result/GCTA_CO_BE_recurrent1_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx3),idx4),1:2],file="../result/GCTA_CO_EA_recurrent1_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx2,idx3),idx4),1:2],file="../result/GCTA_BE_EA_recurrent1_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2,idx3),idx4),1:2],file="../result/GCTA_CO_BEEA_recurrent1_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
idx4=which(pheno$recurrent_HB_RF==0)
write.table(pheno[intersect(c(idx1,idx2),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx3),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx2,idx3),idx4),c(1,2,4)],file="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2,idx3),idx4),1:3],file="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.pheno",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2),idx4),1:2],file="../result/GCTA_CO_BE_recurrent0_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx3),idx4),1:2],file="../result/GCTA_CO_EA_recurrent0_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx2,idx3),idx4),1:2],file="../result/GCTA_BE_EA_recurrent0_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)
write.table(pheno[intersect(c(idx1,idx2,idx3),idx4),1:2],file="../result/GCTA_CO_BEEA_recurrent0_individuals.txt",sep=" ",col.names = F,row.names = F,quote=F)

#generate the cov input,sex,site
generate_cov=function(idx=c(idx1,idx2),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.covar")
{
  res=pheno[idx,1:2]
  idxsampletable=match(pheno$subjid[idx],sampletable$localid)
  if (sum(is.na(idxsampletable))>0) warning("some ids are not mapped")
  res$sex=sampletable$sex[idxsampletable]
  res$sex[res$sex==1]="Male"
  res$sex[res$sex==2]="Female"
  # res$site=sampletable$site[idxsampletable]
  # if (sum(is.na(res$site))>0) warning("some sites are missing")
  # res$site=paste0("Site",res$site)
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}
generate_cov()
generate_cov(idx=c(idx1,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.covar")
generate_cov(idx=c(idx2,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.covar")
generate_cov(idx=c(idx1,idx2,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.covar")
#subseting recurrent_HB_RF
idx4=which(pheno$recurrent_HB_RF==1)
generate_cov(idx=intersect(c(idx1,idx2),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.covar")
generate_cov(idx=intersect(c(idx1,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.covar")
generate_cov(idx=intersect(c(idx2,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.covar")
generate_cov(idx=intersect(c(idx1,idx2,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.covar")
idx4=which(pheno$recurrent_HB_RF==0)
generate_cov(idx=intersect(c(idx1,idx2),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.covar")
generate_cov(idx=intersect(c(idx1,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.covar")
generate_cov(idx=intersect(c(idx2,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.covar")
generate_cov(idx=intersect(c(idx1,idx2,idx3),idx4),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.covar")

#use recurrent as covariate
generate_cov1=function(idx=c(idx1,idx2),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent.covar")
{
  res=pheno[idx,1:2]
  idxsampletable=match(pheno$subjid[idx],sampletable$localid)
  if (sum(is.na(idxsampletable))>0) warning("some ids are not mapped")
  res$sex=sampletable$sex[idxsampletable]
  res$sex[res$sex==1]="Male"
  res$sex[res$sex==2]="Female"
  res$recurrent_HB_RF=sampletable$recurrent_HB_RF[idxsampletable]
  res$recurrent_HB_RF[res$recurrent_HB_RF==0]="recurrent_HB_RF0"
  res$recurrent_HB_RF[res$recurrent_HB_RF==1]="recurrent_HB_RF1"
  # res$site=sampletable$site[idxsampletable]
  # if (sum(is.na(res$site))>0) warning("some sites are missing")
  # res$site=paste0("Site",res$site)
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}
generate_cov1()
generate_cov1(idx=c(idx1,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent.covar")
generate_cov1(idx=c(idx2,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent.covar")
generate_cov1(idx=c(idx1,idx2,idx3),outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent.covar")

generate_qcov=function(idx=c(idx1,idx2),pcafile="../result/GCTA_CO_BE.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar")
{
  res=pheno[idx,1:2]
  pcadat=read.table(pcafile,stringsAsFactors = F)
  idxa=match(pheno$subjid[idx],pcadat$V2)
  pcadat=pcadat[idxa,]
  res=cbind(res,pcadat[,3:ncol(pcadat)])
  idxsampletable=match(pheno$subjid[idx],sampletable$localid)
  res$age=sampletable$age[idxsampletable]
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}

generate_qcov(idx=c(idx1,idx2),pcafile="../result/GCTA_CO_BE.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar")
generate_qcov(idx=c(idx1,idx3),pcafile="../result/GCTA_CO_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.qcovar")
generate_qcov(idx=c(idx2,idx3),pcafile="../result/GCTA_BE_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.qcovar")
generate_qcov(idx=c(idx1,idx2,idx3),pcafile="../result/GCTA_CO_BEEA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.qcovar")

#subseting recurrent_HB_RF
idx4=which(pheno$recurrent_HB_RF==1)
generate_qcov(idx=intersect(c(idx1,idx2),idx4),pcafile="../result/GCTA_CO_BE.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent1.qcovar")
generate_qcov(idx=intersect(c(idx1,idx3),idx4),pcafile="../result/GCTA_CO_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent1.qcovar")
generate_qcov(idx=intersect(c(idx2,idx3),idx4),pcafile="../result/GCTA_BE_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent1.qcovar")
generate_qcov(idx=intersect(c(idx1,idx2,idx3),idx4),pcafile="../result/GCTA_CO_BEEA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent1.qcovar")
idx4=which(pheno$recurrent_HB_RF==0)
generate_qcov(idx=intersect(c(idx1,idx2),idx4),pcafile="../result/GCTA_CO_BE.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.recurrent0.qcovar")
generate_qcov(idx=intersect(c(idx1,idx3),idx4),pcafile="../result/GCTA_CO_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.EA.recurrent0.qcovar")
generate_qcov(idx=intersect(c(idx2,idx3),idx4),pcafile="../result/GCTA_BE_EA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.BE.EA.recurrent0.qcovar")
generate_qcov(idx=intersect(c(idx1,idx2,idx3),idx4),pcafile="../result/GCTA_CO_BEEA.eigenvec",outfile="../result/bca_filtered_04Apr2019.GCTA.CO.BEEA.recurrent0.qcovar")


check_covariate=function(y=c(rep(0,length(idx1)),rep(1,length(idx2))),covfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.covar",qcovfile="../result/bca_filtered_04Apr2019.GCTA.CO.BE.qcovar")
{
  covar=read.table(covfile,stringsAsFactors = F)
  qcovar=read.table(qcovfile,stringsAsFactors = F)
  tmp=cbind(covar[,3:ncol(covar)],qcovar[,3:ncol(qcovar)])
  colnames(tmp)=c("Sex","Site",paste0("PC",1:20),"Age")
  fit=glm(y~.,data=tmp)
  tmp1=model.matrix(~.,data=tmp)
  
}

readhsq=function(hsqfile="../result/GCTA_CO_BE_covar.hsq")
{
  tmp=read.table(hsqfile,sep="\t",fill=T,stringsAsFactors = F,header = T)
  
  return(data.frame(H=tmp$Variance[c(4,7)],SE=tmp$SE[c(4,7)],pvalue=tmp$Variance[12],stringsAsFactors = F))
}
readhsq()
# H       SE pvalue
# 1 0.564408 0.062474      0
# 2 0.364057 0.040298      0
readhsq(hsqfile="../result/GCTA_CO_EA_covar.hsq")
# H       SE pvalue
# 1 0.618551 0.072198      0
# 2 0.256594 0.029950      0
readhsq(hsqfile="../result/GCTA_BE_EA_covar.hsq")
# H       SE    pvalue
# 1 0.156066 0.059618 0.0034029
# 2 0.087784 0.033534 0.0034029
readhsq(hsqfile="../result/GCTA_CO_BEEA_covar.hsq")
# H       SE pvalue
# 1 0.462828 0.043705      0
# 2 0.359769 0.033973      0
readhsq(hsqfile="../result/GCTA_CO_BE_covar_recurrent1.hsq")
# H       SE   pvalue
#1 0.325742 0.218424 0.061372
#2 0.286741 0.192272 0.061372
## 1 0.378872 0.195269 0.023097
readhsq(hsqfile="../result/GCTA_CO_BE_covar_recurrent0.hsq")
# H       SE     pvalue
# 1 0.531605 0.137236 3.9806e-05
# 2 0.337186 0.087046 3.9806e-05
readhsq(hsqfile="../result/GCTA_CO_EA_covar_recurrent1.hsq")
# H       SE     pvalue
# 1 0.884951 0.278931 0.00077176
# 2 0.442589 0.139501 0.00077176
readhsq(hsqfile="../result/GCTA_CO_EA_covar_recurrent0.hsq")
# H       SE     pvalue
# 1 0.60361 0.138062 2.3584e-06
# 2 0.26003 0.059476 2.3584e-06
readhsq(hsqfile="../result/GCTA_BE_EA_covar_recurrent1.hsq")
# H       SE  pvalue
# 1 0.044969 0.09361 0.31198
readhsq(hsqfile="../result/GCTA_BE_EA_covar_recurrent0.hsq")
# H       SE    pvalue
# 1 0.22888 0.095639 0.0081457
readhsq(hsqfile="../result/GCTA_CO_BEEA_covar_recurrent1.hsq")
# H       SE    pvalue
# 1 0.495421 0.178913 0.0025483
readhsq(hsqfile="../result/GCTA_CO_BEEA_covar_recurrent0.hsq")
# H       SE     pvalue
# 1 0.334252 0.062528 1.2708e-08
readhsq(hsqfile="../result/GCTA_CO_BE_covar_recurrent.hsq")
# H       SE     pvalue
# 1 0.269713 0.053639 1.1938e-07
readhsq(hsqfile="../result/GCTA_CO_EA_covar_recurrent.hsq")
# H       SE     pvalue
# 1 0.214692 0.038196 1.2717e-09
readhsq(hsqfile="../result/GCTA_BE_EA_covar_recurrent.hsq")
# H       SE    pvalue
# 1 0.107966 0.046573 0.0079086
readhsq(hsqfile="../result/GCTA_CO_BEEA_covar_recurrent.hsq")
# H       SE   pvalue
# 1 0.272575 0.042966 1.01e-11