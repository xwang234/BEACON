data_vr<-read.table('/fh/fast/dai_j/BEACON/PLINKinputCombo_b_20Feb2014.txt',header=TRUE,sep="\t")

readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  eigsamples=read.table(eigsampfile,stringsAsFactors = F)
  eigsamples=eigsamples$V2
  if (opt==1)
  {
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=eigsamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}

eigenstratmatrix=readeigenstrat(opt=2)
eigenstratmatrix=t(eigenstratmatrix)
sum(data_vr$localid %in% rownames(eigenstratmatrix))
idx=which(!data_vr$localid %in% rownames(eigenstratmatrix))
table(data_vr$site[idx])
# 15  17  18  23  27 
# 1   1   1   1 227 

availsamples=data_vr$localid[data_vr$localid %in% rownames(eigenstratmatrix)]
idx=match(availsamples,data_vr$localid)
dat=data_vr[idx,]
idx=dat$ev1_b==-9 |dat$ev2_b==-9|dat$ev3_b==-9|dat$ev4_b==-9
availsamples=dat$localid[!idx]
idx=match(availsamples,data_vr$localid)
dat=data_vr[idx,]
idx=match(availsamples,rownames(eigenstratmatrix))
eigenstrat=eigenstratmatrix[idx,]
all(rownames(eigenstrat)==dat$localid)
cor(eigenstrat[,1],dat$ev1_b) #[1] 0.9202314
cor(eigenstrat[,4],dat$ev2_b)
cor(eigenstrat[,4],dat$ev3_b)



load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Dong23SNPs.RData")



PC <- PC[row.names(PC)%in% row.names(dat),]

PC <- PC[match(row.names(dat),row.names(PC)),]



dat <- cbind(dat,PC)



out <- matrix(0,23,3)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[10,4]
  
  out[i,2] <- summary(fit2)$coef[10,4]
  
  out[i,3] <- summary(fit3)$coef[10,4]
  
} 

dat1=dat
for (i in 1:ncol(dat1)) 
{
  idx=which(dat1[,i]==-9)
  if (length(idx)>0) dat1[idx,i]=NA
}

#only consider BE
out1=rep(NA,23)
idx=complete.cases(dat1[,c("phenoBE_bc","recurrent_HB_RF","sex","age")])
dat2=dat1[idx,]
for (i in 1:23)
{
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat2[,i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat2)
  out1[i] <- summary(fit1)$coef[10,4]
}

#for the data from paper

library('aod')
library('MASS')
library('gtools')
library('npmlreg')
library('rootSolve')
library('sandwich')
library('xlsx')
library('ncdf')
library('quadprog')
#
nc<-open.ncdf("/fh/fast/dai_j/BEACON/beagessAll.geno.subj.plink.filtered.v6_vH.nc")
snp_anot<-read.csv("/fh/fast/dai_j/BEACON/SNP_annotation.csv",header=TRUE,sep=",")
data_id<-read.table('/fh/fast/dai_j/BEACON/beagessAll.geno.subj.plink.filtered.v6_vH.fam')
data_vr<-read.table('/fh/fast/dai_j/BEACON/PLINKinputCombo_b_20Feb2014.txt',header=TRUE,sep="\t")
filter<-read.table('/fh/fast/dai_j/BEACON/SNP_qual_maf_filter_eabe_vs_co_extract.txt')
# 
data_nv<-data.frame(read.table("/fh/fast/dai_j/BEACON/eabe-newvar-29Apr2014.txt",header=TRUE))

snps<-subset(snp_anot,(snp_anot$rsID)%in%(filter$V1))
snp<-snps[snps$chromosome<=23,]
pos<-snp$position
sampleID<-order(data_id$V2)
data_vr$sampleID<-sampleID
data_nv$sampleID<-sampleID
#
pdata0<-data_vr[order(data_vr$sampleID,decreasing=FALSE),]
pdata1<-data_nv[order(data_nv$sampleID,decreasing=FALSE),]
pdata<-pdata0
pdata$bmi_cat3<-pdata1$bmi_cat3
pdata$hb_rf_comb<-pdata1$hb_rf_comb
#
pdata$hbrf<-pdata0$recurrent_HB_RF
#pdata$hbrf<-pdata$hb_rf_comb
#
g<-get.var.ncdf(nc,"chromosome",start=c(1),count=c(-1))
SNP.pos<-get.var.ncdf(nc,"position",start=c(1),count=c(-1))
sampleNames<-get.var.ncdf(nc,"sampleID",start=c(1),count=c(-1))
#
chrom<-as.numeric(as.vector(snp$chromosome))
snptotal<-length(pos)
#table(pdata$hb_rf_comb)
rsid<-snp$rsID
#
cc_name<-c("phenoBE_b","phenoEA_b","phenoEABE_b")
varname<-c("bmi_cat3","cig_smk_ever","recurrent_HB_RF")
rsids<-c("rs2687201","rs10419226","rs11789015","rs3072","rs2701108","rs9936833","rs9257809")
nvar<-length(varname)

l=2
genx<-get.var.ncdf(nc,"genotype",start=c(min(which(SNP.pos==pos[rsid==rsids[l]])),1),count=c(1,-1))
#
xx<-genx[!is.na(genx)&(genx>=0)]
all0<-as.numeric(xx==0)
all1<-as.numeric(xx==1)
all2<-as.numeric(xx==2)
GFs1<-as.vector(apply(cbind(all0,all1,all2),2,mean))
#
if(GFs1[1]<GFs1[3]){genx<-2-genx}
#
xx<-genx[!is.na(genx)&(genx>=0)]
all0<-as.numeric(xx==0)
all1<-as.numeric(xx==1)
all2<-as.numeric(xx==2)
GFs<-as.vector(apply(cbind(all0,all1,all2),2,mean))
#
AF1<-GFs[1]+0.5*GFs[2]
AF2<-GFs[3]+0.5*GFs[2]
AFs<-c(AF1,AF2)
#
MAF<-min(AFs)
#
#table(genx)
#
#*********************************************************************************************************
#
gdata <- cbind(sampleNames,genx)
gdata <- gdata[sampleNames %in% pdata$sampleID,]
gdata <- data.frame(gdata)
names(gdata) <- c("sampleID","gene")
adat <- merge(pdata,gdata,by="sampleID",all.x=T)
data<-adat[,order(names(adat))]
data$be<-as.numeric(data$phenoBE_b==2)
data$ea<-as.numeric(data$phenoEA_b==2)
data$eabe<-as.numeric(data$phenoEABE_b==2)
#
australia<-c(17,25)
europe<-c(11,13,16)
america<-c(12,14,15,18,19,20,21,22,23,27)
#
region<-rep(0,length(data[,1]))
region[which(data$site  %in% australia)]<-rep(1,length(which(data$site  %in% australia)))
region[which(data$site  %in% europe)]<-rep(2,length(which(data$site  %in% europe)))
region[which(data$site  %in% america)]<-rep(3,length(which(data$site  %in% america)))
#
data$region<-region
#


k<-1
data$cc<-as.vector(data[,which(names(data)==cc_name[k])])
data$cstatus<-as.numeric(data$cc==2)
data$bmi<-data$bmi_cat3
data$smk<-data$cig_smk_ever
data$hbrf<-data$recurrent_HB_RF
data1=data
data1=data[data$gene>=0,]
for (i in 1:ncol(data1))
{
  idx=which(data1[,i]==-9)
  if (length(idx)>0) data1[idx,i]=NA
}
idx=complete.cases(data1[,c("age","sex","ev1_b","ev2_b","ev3_b","ev4_b","hbrf","cc")])
fit2=glm(I(cc==2)~as.numeric(gene)*hbrf+factor(sex)+age + ev1_b + ev2_b + ev3_b + ev4_b,family=binomial,data=data1)
