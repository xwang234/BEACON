#!/usr/bin/env Rscript
#check gwas results

#first do gwas in terminal
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
#generate covariate files
generate_cov=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC",
                        nskip=16)
{
  eigfile=paste0(prefix,".pca")
  eigsampfile=paste0(prefix,".pedind")
  eigsamples=read.table(eigsampfile,stringsAsFactors = F)
  covariate=eigsamples[,c(1,2,6,5)]
  colnames(covariate)=c("FID","IID","affected","sex")
  covariate$IID=gsub("SEP","",covariate$IID)
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  colnames(tmp)=paste0("pc",1:ncol(tmp))
  rownames(tmp)=eigsamples$V2
  tmp=tmp[,1:4]
  covariate=cbind(covariate,tmp)
  covariatefile=paste0(prefix,".covariate")
  write.table(covariate,file=covariatefile,quote=F,row.names = F)
}

generate_cov()
generate_cov(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID")


# plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
#genotyped
# do_gwas(){
#   local prefix="$1"
#   $plink --bfile  $prefix --covar  $prefix.covariate \
#   --covar-name pc1,pc2,pc3,pc4,sex --logistic --hide-covar --ci 0.95 --out $prefix
# }

prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g"
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g"
merge_imp()
{
  local prefix="$1"
  file1=$prefix/filter_hg19tohg38_flip_allchrs.txt
  rm $file1
  for chr in {1..22}
  do
  echo $prefix/chr${chr}_filter_hg19tohg38_flip  >> $file1
  done
  $plink --merge-list $file1 --make-bed \
  --out $prefix/filter_hg19tohg38_flip

  
}
#imputed
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/filter_hg19tohg38_flip"
prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC"

prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/filter_hg19tohg38_flip"
prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID"

update_fam=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/filter_hg19tohg38_flip",opt="BEEA")
{
  if (!file.exists(paste0(prefix,".fam0"))) file.copy(paste0(prefix,".fam"),paste0(prefix,".fam0"))
  fam=read.table(paste0(prefix,".fam"))
  fam$V6=1
  fam$V2=gsub("SEP","",fam$V2)
  if (opt=="BEEA")
  {
    comsamples=intersect(sampletable$localid,fam$V2)
    idx1=match(comsamples,sampletable$localid)
    idx2=match(comsamples,fam$V2)
    fam$V6[idx2]=sampletable$phenoEABE_bca[idx1]
    write.table(fam,file=paste0(prefix,".fam"),quote=F,row.names = F,col.names = F)
  }
}
update_fam()
update_fam(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/filter_hg19tohg38_flip")
do_gwas(){
  local prefix="$1" #genotype
  local prefix1="$2" #covariate
  $plink --bfile  $prefix --covar  $prefix1.covariate --allow-no-sex \
  --covar-name pc1,pc2,pc3,pc4,sex --logistic --hide-covar --ci 0.95 --out $prefix
}

tmp=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.fam")
colnames(tmp)[1:2]=c("fam","localid")
tmp1=merge(tmp,sampletable,by=c("fam","localid"))
all(tmp1$V6==tmp1$phenoEABE_bca) #T:BEEA
#compare results
compare_gwas=function(file1="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/filter_hg19tohg38_flip.assoc.logistic",
                      file2="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt")
{
  tmp1=as.data.frame(fread(file1))
  tmp2=as.data.frame(fread(file2))
  # tmp3=sapply(1:nrow(tmp2),function(x){nchar(tmp2$non_effect_allele[x])})
  # tmp4=sapply(1:nrow(tmp2),function(x){nchar(tmp2$effect_allele[x])})
  # sum(tmp3>1 |tmp4>1)/nrow(tmp2) #0.06648272
  # tmp1_snp=sapply(1:nrow(tmp1),function(x){
  #   unlist(strsplit(tmp1$SNP[x],"_"))[1]
  # })
  
  # tmp2_snp=paste0(tmp2$CHR,":",tmp2$position)
  # com_snp=intersect(tmp1_snp,tmp2_snp)
  # idx1=match(com_snp,tmp1_snp)
  # idx2=match(com_snp,tmp2_snp)
  # plot(-log10(tmp1$P[idx1]),-log10(tmp2$P[idx2]))
  
  tmp1_snp=tmp1$SNP
  tmp2_snp1=paste0(tmp2$CHR,":",tmp2$position,"_",tmp2$non_effect_allele,"_",tmp2$effect_allele)
  tmp2_snp2=paste0(tmp2$CHR,":",tmp2$position,"_",tmp2$effect_allele,"_",tmp2$non_effect_allele)
  
  com_snp1=intersect(tmp1_snp,tmp2_snp1)
  idx1=match(com_snp1,tmp1_snp)
  idx2=match(com_snp1,tmp2_snp1)
  com_snp2=intersect(tmp1_snp,tmp2_snp2)
  idx3=match(com_snp2,tmp1_snp)
  idx4=match(com_snp2,tmp2_snp2)
  plot(-log10(tmp1$P[c(idx1,idx3)]),-log10(tmp2$P[c(idx2,idx4)]),xlab="New",ylab="Old")
  abline(0,1,col="red")
  cor1=round(cor(-log10(tmp1$P[c(idx1,idx3)]),-log10(tmp2$P[c(idx2,idx4)])),2)
  #text(1.5,5.5,paste0("correlation=",cor1),cex=1.3)
  text(1.5,6.5,paste0("correlation=",cor1),cex=1.3)
}

compare_gwas(file1="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/filter_hg19tohg38_flip.assoc.logistic",
             file2="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt")

qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log10)",
       ylab="Observed p-value (log10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

qqplot_gwas=function(file1="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/filter_hg19tohg38_flip.assoc.logistic",
                      file2="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_Cambridge_autosomes.txt")
{
  tmp1=as.data.frame(fread(file1))
  tmp2=as.data.frame(fread(file2))
  # tmp1_snp=sapply(1:nrow(tmp1),function(x){
  #   unlist(strsplit(tmp1$SNP[x],"_"))[1]
  # })
  
  # tmp2_snp=paste0(tmp2$CHR,":",tmp2$position)
  # com_snp=intersect(tmp1_snp,tmp2_snp)
  # idx1=match(com_snp,tmp1_snp)
  # idx2=match(com_snp,tmp2_snp)
  # plot(-log10(tmp1$P[idx1]),-log10(tmp2$P[idx2]))
  
  tmp1_snp=tmp1$SNP
  tmp2_snp1=paste0(tmp2$CHR,":",tmp2$position,"_",tmp2$non_effect_allele,"_",tmp2$effect_allele)
  tmp2_snp2=paste0(tmp2$CHR,":",tmp2$position,"_",tmp2$effect_allele,"_",tmp2$non_effect_allele)
  
  com_snp1=intersect(tmp1_snp,tmp2_snp1)
  idx1=match(com_snp1,tmp1_snp)
  idx2=match(com_snp1,tmp2_snp1)
  com_snp2=intersect(tmp1_snp,tmp2_snp2)
  idx3=match(com_snp2,tmp1_snp)
  idx4=match(com_snp2,tmp2_snp2)
  par(mar=c(6,6,2,1))
  qqplot(tmp1$P[c(idx1,idx3)])
  qqplot(tmp2$P[c(idx2,idx4)])
  
}
qqplot_gwas(file1="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/filter_hg19tohg38_flip.assoc.logistic",
             file2="/fh/fast/dai_j/BEACON/Meta_summary_stat/BEEA_BEACON_autosomes.txt")

#compare controls using PCA
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

#add covariate table (pc1-pc4)
readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/merge_beacon_cambridge_genotype.pedind",
                        nskip=16)
{
  eigsamples=read.table(eigsampfile,stringsAsFactors = F)
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  colnames(tmp)=paste0("pc",1:ncol(tmp))
  rownames(tmp)=eigsamples$V2
  tmp$sex="M"
  tmp$sex[eigsamples$V5==2]="F"
  return(tmp)
}
tmp=read.table(eigfile,stringsAsFactors = F)
tmp=tmp[2:16,1]
pcavarexplained=tmp/sum(tmp)
pcavarexplained
# [1] 0.31023245 0.09634963 0.08770073 0.06452556 0.06091837 0.05633626 0.04568182 0.04027799 0.03601621 0.03516664 0.03463740 0.03385747 0.03313324
# [14] 0.03272935 0.03243687
plot(pcavarexplained*100,ylab="Variance explained (%)",cex.axis=1.3,cex.lab=1.3,type="l",xlab="PC index")
points(pcavarexplained*100)
covariatetable=readeigenstrat()
rownames(covariatetable)=gsub("SEP","",rownames(covariatetable))
#add case/control
tmp=covariatetable
tmp$phenoBE_bca=tmp$phenoEA_bca=tmp$phenoEABE_bca=1
comsamples=intersect(sampletable$localid,rownames(tmp))
idx1=match(comsamples,rownames(tmp))
idx2=match(comsamples,sampletable$localid)
tmp$phenoBE_bca[idx1]=sampletable$phenoBE_bca[idx2]
tmp$phenoEA_bca[idx1]=sampletable$phenoEA_bca[idx2]
tmp$phenoEABE_bca[idx1]=sampletable$phenoEABE_bca[idx2]
tmp$phenoBE_bca[tmp$phenoBE_bca==-9]=NA
tmp$phenoEA_bca[tmp$phenoEA_bca==-9]=NA
tmp$phenoEABE_bca[tmp$phenoEABE_bca==-9]=NA
covariatetable=tmp
covariatetable$sex=factor(covariatetable$sex)

amos=read.table("../data/AdditionalControlPuya/GENEVA_Melanoma_controls_caucasians_cleaned2.fam")
sum(rownames(covariatetable) %in% amos$V2) #1038
beacon=read.table("../data/AdditionalControlPuya/BEACON_CRIC_GENEVAMel_PD.fam")
beacon$V2=gsub("_","",beacon$V2)
cambridge=read.table("../data/AdditionalControlPuya/cambridge_cases_WTCCC_controls.fam")
cambridge$V2=gsub("_","",cambridge$V2)
covariatetable$type=NA
covariatetable$type[rownames(covariatetable) %in% sampletable$localid]="Old beacon cases"
covariatetable$type[which(covariatetable$type=="Old beacon cases" & covariatetable$phenoEABE_bca==1)]="Old beacon controls"
covariatetable$type[(rownames(covariatetable) %in% cambridge$V2) & (is.na(covariatetable$type))]="New WTCCC2 controls"
covariatetable$type[rownames(covariatetable) %in% beacon$V2 & is.na(covariatetable$type)]="New dbGap controls"


#covariatetable$type[which(rownames(covariatetable) %in% amos$V2)]="New AMOS controls"

#new dbgap controls
# cric=read.table("../data/AdditionalControlPuya/CRIC_Non-hispanic_white_cleaned2_sexchecked.fam")
# cric$V2=gsub("_","",cric$V2)
# geneva=read.table("../data/AdditionalControlPuya/GENEVA_Melanoma_controls_caucasians_cleaned2.fam")
# geneva$V2=gsub("_","",geneva$V2)
# pd=read.table("../data/AdditionalControlPuya/PD_ENV_controls_caucasians_cleaned2_sexchecked.fam")
# pd$V2=gsub("_","",pd$V2)
# 
# covariatetable$type=NA
# covariatetable$type[rownames(covariatetable) %in% cric$V2]="dbGap controls1"
# covariatetable$type[rownames(covariatetable) %in% geneva$V2]="dbGap controls2"
# covariatetable$type[rownames(covariatetable) %in% pd$V2]="dbGap controls3"

allcolors=c("blue","red","green","black","brown4","darkseagreen","darkseagreen1")
allpch=c(15,16,17,18,19)

library(scales)
plot.pca=function(dat=covariatetable[covariatetable$type!= "Old beacon cases",],yinc=1.2,pc1=1,pc2=2,main="",xlim=NULL,ylim=NULL)
{
  dat=dat[!is.na(dat$type),]
  types=factor(dat$type,levels = unique(dat$type))
 
  if (is.null(xlim))
  {
    xmin=min(dat[,pc1])
    xmax=max(dat[,pc1])
    ymin=min(dat[,pc2])
    ymax=max(dat[,pc2])*yinc
    pch=allpch[types]
    par(mar=c(6,6,2,1))
    plot(dat[,pc1],dat[,pc2],col=alpha(allcolors[types],0.05),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
         ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("PC",pc1),ylab=paste0("PC",pc2),main=main,bty='l',cex=1.3)
    #text(x=dat[,1],y=dat[,2],rownames(dat),col=colors[clinicaltable$`Platium response`],cex=1)
    
  }else
  {
    plot(dat[,pc1],dat[,pc2],col=alpha(allcolors[types],0.05),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=xlim,
         ylim=ylim,xlab=paste0("PC",pc1),ylab=paste0("PC",pc2),main=main,bty='l',cex=1.3)
  }
    legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=allpch,cex=1.3,bty = "n",ncol=2)
}

plot.pca()
plot.pca(pc2=3)
plot.pca(pc1=2,pc2=3)
plot.pca(dat=covariatetable[!covariatetable$type %in% c("Old beacon cases","New WTCCC2 controls"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("Old beacon controls"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("New WTCCC2 controls"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("New dbGap controls"),])
covariatetable$type[which(covariatetable$phenoEABE_bca==2 & rownames(covariatetable) %in% cambridge$V2)]="Cambridge cases"
covariatetable$type[which(covariatetable$phenoEABE_bca==2 & rownames(covariatetable) %in% beacon$V2)]="BEACON cases"

#covariatetable$type[which(covariatetable$type %in% c("Old beacon controls","New dbGap controls"))]="BEACON controls"
plot.pca(dat=covariatetable[covariatetable$type %in% c("New dbGap controls","BEACON cases","Old beacon controls"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("New dbGap controls","BEACON cases","Old beacon controls"),],pc2=3)

plot.pca(dat=covariatetable[covariatetable$type %in% c("BEACON controls","BEACON cases"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("BEACON controls","BEACON cases"),],pc2=3)
plot.pca(dat=covariatetable[covariatetable$type %in% c("New WTCCC2 controls","Cambridge cases"),])
plot.pca(dat=covariatetable[covariatetable$type %in% c("New WTCCC2 controls","Cambridge cases"),],pc2=3)

allcolors[1]="red"
plot.pca(dat=covariatetable[covariatetable$type %in% c("New WTCCC2 controls"),],pc2=3,xlim=c(-0.03,0.01),ylim=c(-0.04,0.04))
allcolors[1]="blue"
plot.pca(dat=covariatetable[covariatetable$type %in% c("Cambridge cases"),],pc2=3,xlim=c(-0.03,0.01),ylim=c(-0.04,0.04))


# #work on PCs generated by cambridge/beacon individually
# #Cambridge
# covariatetable=readeigenstrat(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC.pca",
#                               eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/Cambridge_WTCCC.pedind")
# rownames(covariatetable)=gsub("SEP","",rownames(covariatetable))
# tmp=covariatetable
# tmp$phenoBE_bca=tmp$phenoEA_bca=tmp$phenoEABE_bca=1
# comsamples=intersect(sampletable$localid,rownames(tmp))
# idx1=match(comsamples,rownames(tmp))
# idx2=match(comsamples,sampletable$localid)
# tmp$phenoBE_bca[idx1]=sampletable$phenoBE_bca[idx2]
# tmp$phenoEA_bca[idx1]=sampletable$phenoEA_bca[idx2]
# tmp$phenoEABE_bca[idx1]=sampletable$phenoEABE_bca[idx2]
# tmp$phenoBE_bca[tmp$phenoBE_bca==-9]=NA
# tmp$phenoEA_bca[tmp$phenoEA_bca==-9]=NA
# tmp$phenoEABE_bca[tmp$phenoEABE_bca==-9]=NA
# covariatetable=tmp
# covariatetable$sex=factor(covariatetable$sex)
# covariatetable$type="Cambridge cases"
# covariatetable$type[which(covariatetable$phenoEABE_bca==1)]="Cambridge new controls"
# plot.pca()
# 
# #BEACON
# covariatetable=readeigenstrat(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.pca",
#                               eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.pedind")
# rownames(covariatetable)=gsub("SEP","",rownames(covariatetable))
# tmp=covariatetable
# tmp$phenoBE_bca=tmp$phenoEA_bca=tmp$phenoEABE_bca=1
# comsamples=intersect(sampletable$localid,rownames(tmp))
# idx1=match(comsamples,rownames(tmp))
# idx2=match(comsamples,sampletable$localid)
# tmp$phenoBE_bca[idx1]=sampletable$phenoBE_bca[idx2]
# tmp$phenoEA_bca[idx1]=sampletable$phenoEA_bca[idx2]
# tmp$phenoEABE_bca[idx1]=sampletable$phenoEABE_bca[idx2]
# tmp$phenoBE_bca[tmp$phenoBE_bca==-9]=NA
# tmp$phenoEA_bca[tmp$phenoEA_bca==-9]=NA
# tmp$phenoEABE_bca[tmp$phenoEABE_bca==-9]=NA
# covariatetable=tmp
# covariatetable$sex=factor(covariatetable$sex)
# covariatetable$type="Old beacon cases"
# covariatetable$type[which(covariatetable$phenoEABE_bca==1 & rownames(covariatetable) %in% sampletable$localid)]="Old beacon controls"
# covariatetable$type[which(covariatetable$phenoEABE_bca==1 & ! rownames(covariatetable) %in% sampletable$localid)]="New dbGap controls"
# 
# plot.pca(dat=covariatetable)

#work on old PCs
dat=sampletable[,c("ev1_bca","ev2_bca","ev3_bca","ev4_bca","site","phenoEABE_bca")]
dat$type=NA
dat$type[which(dat$phenoEABE_bca==1)]="Old beacon controls"
plot.pca(dat=dat)
plot.pca(dat=dat,pc2=3)
