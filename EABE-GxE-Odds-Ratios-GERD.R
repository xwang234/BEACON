##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#                       EABE: BMI, SMOKING and HBRF joint-interaction with gene
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
rm(list=ls())
#
#
#
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
#names(data_nv)
#
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
#
#
#********************************************************************************************************************
#
#
#
#********************************************************************************************************************
#
#
#
#names(pdata)
#
cc_name<-c("phenoBE_b","phenoEA_b","phenoEABE_b")
varname<-c("bmi_cat3","cig_smk_ever","recurrent_HB_RF")
rsids<-c("rs2687201","rs10419226","rs11789015","rs3072","rs2701108","rs9936833","rs9257809")
nvar<-length(varname)
#
#
#********************************************************************************************************************
#
#   Site missing data percentage
#
#********************************************************************************************************************
#
#
site_missing<-sapply(1:nvar,function(i){
    #i<-1
    sdata<-pdata
    sdata$expo<-as.vector(sdata[,which(names(sdata)==varname[i])]) 
    sdata$miss<-as.numeric(sdata$expo==-9 | as.numeric(is.na(sdata$expo))==1) 
    tab1<-table(sdata$miss,sdata$site)
    round(100*c(tab1[2,]/apply(tab1,2,sum),mean(sdata$miss)),1)
  })     
#
with(pdata[pdata$phenoBE_b!=-9,],table(phenoBE_b,site))
with(pdata[pdata$phenoEA_b!=-9,],table(phenoEA_b,site))
with(pdata[pdata$phenoEA_b!=-9,],table(phenoEABE_b,site))
#
#
#********************************************************************************************************************
#
#
#********************************************************************************************************************
#
expo_summary<-function(i){
  #
  sdata<-pdata
  sdata$expo<-as.vector(sdata[,which(names(sdata)==varname[i])]) 
  #i<-1;k<-1
  #
  summar<-sapply(1:2,function(k){
  #
  #
  sdata$cc<-as.vector(sdata[,which(names(sdata)==cc_name[k])])
  sdata$cstatus<-as.numeric(sdata$cc==2)  
  sdatas<-sdata[as.numeric(sdata$cc!=-9 & sdata$ev1_b!=-9 & sdata$ev2_b!=-9 & sdata$ev3_b!=-9 & sdata$ev4_b!=-9 & 
                             sdata$expo!=-9 & as.numeric(is.na(sdata$expo))==0)==1,]  
  # 
  tab<-table(sdatas$expo,sdatas$cstatus)
  fi<-glm(cstatus~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(expo),family=binomial,data=sdatas)
  ri<-summary(fi)$coefficients
  vi<-vcovHC(fi,type="HC0")
  #ri
  be_exp<-as.vector(ri[-c(1:7),1])
  se_exp<-as.vector(ri[-c(1:7),2])
  pb_exp<-as.vector(ri[-c(1:7),4])
  or_exp<-exp(be_exp)
  or_exp_lb<-exp(be_exp-1.96*se_exp)
  or_exp_ub<-exp(be_exp+1.96*se_exp)
  ORs<-cbind(or_exp, or_exp_lb,or_exp_ub)
  Tabs<-cbind(tab[,2],tab[,1])
  Est<-data.frame(be_exp,se_exp,pb_exp)
  colnames(ORs)<-c("OR","OR_lb","OR_ub")
  colnames(Tabs)<-c("case","control")
  list(Tabs=Tabs,Est=Est,ORs=ORs)
  #
  })
  summar
}
#
#
#********************************************************************************************************************
#
i<-1
#
Risksummary<-expo_summary(i)
BE<-Risksummary[,1]
EA<-Risksummary[,2]
#
BE$Tabs
BE$ORs
#
EA$Tabs
EA$ORs
#
#
#
#********************************************************************************************************************
#
#
fit_be<-glm(I(phenoBE_b==2)~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(bmi_cat3),family=binomial,data=pdata[pdata$phenoBE_b!=-9 & 
                         pdata$ev1_b!=-9 & pdata$ev2_b!=-9 & pdata$ev3_b!=-9 & pdata$ev4_b!=-9 & pdata$bmi_cat3!=-9, ])        
fit_ea<-glm(I(phenoEA_b==2)~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(bmi_cat3),family=binomial,data=pdata[pdata$phenoEA_b!=-9 & 
                         pdata$ev1_b!=-9 & pdata$ev2_b!=-9 & pdata$ev3_b!=-9 & pdata$ev4_b!=-9 & pdata$bmi_cat3!=-9, ])   
summary(fit_be)$coefficients
summary(fit_ea)$coefficients
#
#
#
#
#********************************************************************************************************************
#
#
#********************************************************************************************************************
#
#
snp_maf<-get(load("EABE-MAF.Rdata"))
maf_s<-snp_maf[snp_maf[,6]>=0.25,c(1,4:6)]
len_maf<-length(maf_s[,1])
#
#
SNP_sample<-function(l){
  #
  set.seed(1000+l)
  #
  pvals<-rep(0,14)
  snpmaf0<-maf_s[sample(1:len_maf,7,replace=F),]
  #
  it<-0
  snpmaf<-snpmaf0
  #
  while(min(pvals<=0.05)){
  #
  it<-it+1
  set.seed(5000+it+l)  # s3
  #set.seed(5000+it+l)  #s
  snpmaf<-maf_s[sample(1:len_maf,7,replace=F),]
  pval<-matrix(0,ncol=2,nrow=7)
  pos_maf<-snpmaf[,3]
  #
  for(r in 1:7){
  #r<-1
  genx<-get.var.ncdf(nc,"genotype",start=c(min(which(SNP.pos==pos_maf[r])),1),count=c(1,-1))
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
  cc_name<-c("phenoBE_b","phenoEA_b")
  #
  pval[r,]<-as.vector(sapply(1:2,function(k){
    #k<-1
    data$cc<-as.vector(data[,which(names(data)==cc_name[k])])
    data$cstatus<-as.numeric(data$cc==2)
      #
      obs<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                      & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                      & data$cc!=-9)  
      #
      datas<-data[obs==1,]
      x<-datas$gene
      y<-datas$cstatus
      x1<-as.numeric(x>0)
      sex<-as.numeric(datas$sex==2)
      datas$y<-y
      datas$x<-x
      fit_mar<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+x,family=binomial,data=datas)
      res_mar<-summary(fit_mar)$coefficients
      res_mar[8,4]
  })
  )
  }   
  #
  pvals<-as.vector(pval) 
  }
  #it;pvals
  #
  snpmaf
  }
#
#
#
#********************************************************************************************************************
#
#
rsids_sample<-SNP_sample(10)
rsid_sample<-as.vector(rsids_sample$rsID)
#
#
rsid1<-c("rs2687201","rs10419226","rs11789015","rs3072","rs2701108","rs9936833","rs9257809")
rsids<-rsid1
#rsids<-rsid_sample
#
#
#********************************************************************************************************************
#
#
Table4<-function(l){
  #l<-1
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
  #
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
  cc_name<-c("phenoBE_b","phenoEA_b")
  varname<-c("bmi_cat3","cig_smk_ever","recurrent_HB_RF")
  #
  k<-1;i<-3
  #k<-2;i<-1
  data$cc<-as.vector(data[,which(names(data)==cc_name[k])])
  data$cstatus<-as.numeric(data$cc==2)
  #
  data$expo<-as.vector(data[,which(names(data)==varname[i])])
  obs<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                  & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                  & data$cc!=-9 & data$expo!=-9 & !is.na(data$expo))  
  #
  obs_mar<-as.numeric(!is.na(data$gene) & data$gene>=0 & data$gene!=3 & data$cc!=-9 & data$expo!=-9 & !is.na(data$expo)) 
  #
  datas<-data[obs==1,]
  datas1<-data[obs_mar==1,]
  datas1$cat<-with(datas1,ifelse(gene==0 & expo==0,1,ifelse(gene>0 & expo==0,2,ifelse(gene==0 & expo==1,3,4))))
  tab<-cbind(table(with(datas1[datas1$cstatus==1,],ifelse(gene==0 & expo==0,1,ifelse(gene>0 & expo==0,2,ifelse(gene==0 & expo==1,3,4))))),
             table(with(datas1[datas1$cstatus==0,],ifelse(gene==0 & expo==0,1,ifelse(gene>0 & expo==0,2,ifelse(gene==0 & expo==1,3,4))))))
  #
  #
  datas$GGERD<-with(datas,ifelse(gene==0 & expo==0,1,ifelse(gene>0 & expo==0,2,ifelse(gene==0 & expo==1,3,4))))
  datas$gene1<-as.numeric(datas$gene>0)
  #
  fit_t4<-glm(cstatus~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(GGERD),family=binomial,data=datas)
  re_t4<-summary(fit_t4)$coefficients
  ll4<-length(re_t4[,1])
  par_exp<-c((ll4-2):ll4)
  be_exp<-re_t4[par_exp,1]
  se_exp<-re_t4[par_exp,2]
  or_exp<-exp(be_exp)
  or_exp_lb<-exp(be_exp-1.96*se_exp)
  or_exp_ub<-exp(be_exp+1.96*se_exp)
  OR_GERD<-data.frame(or_exp,or_exp_lb,or_exp_ub)
  #
  ffint<-summary(glm(cstatus~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+gene1+expo+gene1*expo,family=binomial,data=datas))
  rrint<-ffint$coefficients
  #
  par_exp<-c(8:10)
  vvint<-vcov(ffint)
  #vvint<-vcovHC(ffint,type="HC0")
  #
  be_exp<- rrint[par_exp,1]
  va_exp<-vvint[par_exp,par_exp]
  #
  g0<-c(1,0,0)
  g1<-c(1,0,1)
  g2<-c(0,1,0)
  g3<-c(0,1,1)
  #
  be_exp0<-be_exp%*%g0
  be_exp1<-be_exp%*%g1
  be_exp2<-be_exp%*%g2
  be_exp3<-be_exp%*%g3
  #
  se_exp0<-sqrt(t(g0)%*%va_exp%*%g0)
  se_exp1<-sqrt(t(g1)%*%va_exp%*%g1)
  se_exp2<-sqrt(t(g2)%*%va_exp%*%g2)
  se_exp3<-sqrt(t(g3)%*%va_exp%*%g3)
  #
  or_exp0<-exp(be_exp0)
  or_exp1<-exp(be_exp1)
  or_exp2<-exp(be_exp2)
  or_exp3<-exp(be_exp3)
  #
  or_exp0_lb<-exp(be_exp0-1.96*se_exp0)
  or_exp0_ub<-exp(be_exp0+1.96*se_exp0)
  or_exp1_lb<-exp(be_exp1-1.96*se_exp1)
  or_exp1_ub<-exp(be_exp1+1.96*se_exp1)
  or_exp2_lb<-exp(be_exp2-1.96*se_exp2)
  or_exp2_ub<-exp(be_exp2+1.96*se_exp2)
  or_exp3_lb<-exp(be_exp3-1.96*se_exp3)
  or_exp3_ub<-exp(be_exp3+1.96*se_exp3)
  #
  gor_exp<-rbind(c(or_exp0,or_exp0_lb,or_exp0_ub),
                 c(or_exp1,or_exp1_lb,or_exp1_ub),
                 c(or_exp2,or_exp2_lb,or_exp2_ub),
                 c(or_exp3,or_exp3_lb,or_exp3_ub))
  #
  list(Tab_BE=tab,OR_GERD=OR_GERD, GG_OR=gor_exp)
  #
  
}
Table4(1)
  #
  #
  #
  #*********************************************************************************************************
  #
  # 
  #
  #*********************************************************************************************************
  #
  #
  cc_name<-c("phenoBE_b","phenoEA_b","phenoEABE_b")
  varname<-c("bmi_cat3","cig_smk_ever","recurrent_HB_RF")
  rsids<-c("rs2687201","rs10419226","rs11789015","rs3072","rs2701108","rs9936833","rs9257809")
  #rsids<-rsid_sample
  nvar<-length(varname)
  #names(pdata);table(pdata$phenoEABE_b)
  #
  #
  #********************************************************************************************************************
  #
  #
Tables<-function(l){
#
#
#l<-1
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
cc_name1<-c("phenoBE_b","phenoEA_b")
varname1<-c("bmi_recent_healthy","bmi_cat3","cig_smk_ever","heartburn","heartburn_freq4",
           "recurrent_heartburn","reflux","reflux_freq4","recurrent_reflux","recurrent_HB_RF")
#
data_all<-data
#
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#
#
#
#
#data<-data_all
#data<-data_all[which(data_all$site  %in% america),]
data<-data_all[which(data_all$site  %in% europe),]
#data<-data_all[which(data_all$site  %in% australia),]
#
#
#
#
#*********************************************************************************************************
#
sites<-unique(data$site)
sites<-sort(sites)
data1<-data
#
site_mr<-sapply(1:nvar,function(i){
  #
  #i<-1
  data1$expo<-as.vector(data1[,which(names(data1)==varname[i])])
  data1$missing<-as.numeric(data1$expo==-9|as.numeric(is.na(data1$expo))==1)
  #
  mrates<-rep(0,length(sites))
  for(s in 1:length(sites)){
    data_site<-data1[data1$site==sites[s],]
    
    mrates[s]<-mean(data_site$missing)
  }
  mrates})
#
site_mrate<-data.frame(100*round(site_mr,3))
sitename<-c("England-Sheffield","Kaiser","Sweden-Karolinska","Mayo","EGA-WA","Ireland-FINBAR","Australia-Queensland",
            "Toronto","UNC","WA reflux","Canada-Nova Scotia","EGA-NJ","USC Keck","Australia-wide","WA Reid")
#
#write.csv(site_mrate,file="EABE-GxE-exposure-missing rate by site-01-26-2014.csv")
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#

  Table1<-function(k){
  #k<-1
  data$cc<-as.vector(data[,which(names(data)==cc_name[k])])
  data$cstatus<-as.numeric(data$cc==2)
  
  oddratio<-matrix(0,ncol=13,nrow=nvar)
  #
  for(i in 1:nvar){
    #
    data$expo<-as.vector(data[,which(names(data)==varname[i])])
    obs<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                    & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                    & data$cc!=-9 & data$expo!=-9 & !is.na(data$expo))  
    #
    obs1<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                     & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                     & data$cc!=-9) 
    #
    #*********************************************************************************************************
    #
    #
    #table(datas$education12)
    #australia<-c(17,25)
    #europe<-c(11,13,16)
    #america<-c(12,14,15,18,19,20,21,22,23,27)
    #
    datas<-data[obs==1,]
    datas1<-data[obs1==1,]
    mrate<-length(which(datas1$expo==-9|as.numeric(is.na(datas1$expo))==1))/length(datas1[,1])
    #
    x<-datas$gene
    y<-datas$cstatus
    x1<-as.numeric(x>0)
                   sex<-as.numeric(datas$sex==2)
                   datas$y<-y
                   datas$x<-x
                   datas$x1<-x1
                   n<-length(datas[,1])
                   #
                   #
                   #*********************************************************************************************************
                   #
                   #  EABE: Odds ratios SNP: rs2687201
                   #
                   #*********************************************************************************************************
                   #
                   #
                   fit_int2<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(x)+expo+factor(x)*expo,family=binomial,data=datas)
                   #
                   re_int<-summary(fit_int2)$coefficients
                   #va_int<-vcov(fit_int)
                   va_int<-vcovHC(fit_int2,type="HC0")
                   #
                   par_exp<-c(10:12)
                   be_exp<-re_int[par_exp,1]
                   va_exp<-va_int[par_exp,par_exp]
                   #
                   g0<-c(1,0,0)
                   g1<-c(1,1,0)
                   g2<-c(1,0,1)
                   #
                   be_exp0<-be_exp%*%g0
                   be_exp1<-be_exp%*%g1
                   be_exp2<-be_exp%*%g2
                   #
                   se_exp0<-sqrt(t(g0)%*%va_exp%*%g0)
                   se_exp1<-sqrt(t(g1)%*%va_exp%*%g1)
                   se_exp2<-sqrt(t(g2)%*%va_exp%*%g2)
                   #
                   or_exp0<-exp(be_exp0)
                   or_exp1<-exp(be_exp1)
                   or_exp2<-exp(be_exp2)
                   #
                   or_exp0_lb<-exp(be_exp0-1.96*se_exp0)
                   or_exp0_ub<-exp(be_exp0+1.96*se_exp0)
                   or_exp1_lb<-exp(be_exp1-1.96*se_exp1)
                   or_exp1_ub<-exp(be_exp1+1.96*se_exp1)
                   or_exp2_lb<-exp(be_exp2-1.96*se_exp2)
                   or_exp2_ub<-exp(be_exp2+1.96*se_exp2)
                   #
                   or_exp<-c(or_exp0,or_exp0_lb,or_exp0_ub,or_exp1,or_exp1_lb,or_exp1_ub,or_exp2,or_exp2_lb,or_exp2_ub)
                   prob_chi2<-(wald.test(b=coef(fit_int2),Sigma=vcovHC(fit_int2,type="HC0"),Terms=par_exp[-1])$result)$chi2
                   #
                   fit_int3<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+x+expo+x*expo,family=binomial,data=datas)
                   ptrend<-(summary(fit_int3)$coefficients)[10,4]
                   #
                   fit_int4<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(x)+expo+factor(x1)*expo,family=binomial,data=datas)
                   res_int4<-summary(fit_int4)$coefficients
                   p_int4<-res_int4[11,4]
                   #
                  
    #
    design1<-model.matrix(fit_int3)
    design2<-model.matrix(fit_int4)
    l2<-length(design2[1,])
    des1<-design1
    des2<-design2[,-c(l2-1)]
    design<-cbind(des1,des2)
    beta<-rbind(summary(fit_int3)$coefficients,summary(fit_int4)$coefficients)[,1]
    #
    fscore<-function(b){
      #
      expb1<-exp(as.vector(des1%*%b[1:10]))
      expb2<-exp(as.vector(des2%*%b[11:21]))
      prob1<-expb1/(1+expb1)
      prob2<-expb2/(1+expb2)
      cbind(des1*(y-prob1),des2*(y-prob2))
    }
    #
    score<-function(b){apply(fscore(b),2,sum)}
    fisher<-function(b){-gradient(score,b,centered =FALSE,pert=1e-8)}
    vars<-function(b){(n-1)*var(fscore(b))}
    #
    sand<-function(b){
      GG1<-solve(fisher(b))
      BB1<-vars(b)
      VA1<-GG1%*%BB1%*%t(GG1)
      VA1
    }
    sandvar<-sand(beta)
    par_int<-c(10,21)
    vcovm<-sandvar[par_int,par_int]
    corr_int<-cov2cor(vcovm)[1,2]
    #
    fit_int5<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(x)+expo+x*expo,family=binomial,data=datas)
    res_int5<-summary(fit_int5)$coefficients
    p_int5<-res_int5[11,4]
    #
    
    oddratio[i,]<-c(round(100*mrate,2),round(or_exp,3),round(p_int4,4),round(ptrend,4),corr_int)
    
  }
  #
  rownames(oddratio)<-varname
  colnames(oddratio)<-c("MR","OR0","OR0_lb","OR0_ub","OR1","OR1_lb","OR1_ub","OR2","OR2_lb","OR2_ub","Prob","P_trend","correlation")
  oddratio
  }
#
#names(datas)
#table(datas$education12)
#
#*********************************************************************************************************
#  
#
#*********************************************************************************************************
#  
#
  Table2<-function(k){
  #k<-1
  data$cc<-as.vector(data[,which(names(data)==cc_name[k])])
  data$cstatus<-as.numeric(data$cc==2)
  data$bmi<-data$bmi_cat3
  data$smk<-data$cig_smk_ever
  data$hbrf<-data$recurrent_HB_RF
  #
  obs<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                  & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                  & data$cc!=-9 & data$hbrf!=-9 & !is.na(data$hbrf) & data$bmi!=-9 & !is.na(data$bmi)
                  & data$smk!=-9 & !is.na(data$smk))  
  #
  obs1<-as.numeric(!is.na(data$ev1_b) & !is.na(data$ev2_b) & !is.na(data$ev3_b) & !is.na(data$ev4_b) & !is.na(data$gene) 
                   & data$ev1_b!=-9 & data$ev2_b!=-9 & data$ev3_b!=-9 & data$ev4_b!=-9 & data$gene>=0 & data$gene!=3 
                   & data$cc!=-9)
  #
  #*********************************************************************************************************
  #
  #
  #
  datas<-data[obs==1,]
  datas1<-data[obs1==1,]
  mrate<-length(which(datas1$bmi==-9|as.numeric(is.na(datas1$bmi))==1|
                        datas1$smk==-9|as.numeric(is.na(datas1$smk))==1|
                        datas1$hbrf==-9|as.numeric(is.na(datas1$hbrf))==1))/length(datas1[,1])
  #
  x<-datas$gene
  y<-datas$cstatus
  x1<-as.numeric(x>0)
  sex<-as.numeric(datas$sex==2)
  datas$y<-y
  datas$x<-x
  datas$x1<-x1
  n<-length(datas[,1])
  #
  #
  #*********************************************************************************************************
  #
  fit_int2<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(x)+bmi+smk+hbrf+
                  factor(x)*bmi+factor(x)*smk+factor(x)*hbrf,family=binomial,data=datas)
  re_int<-summary(fit_int2)$coefficients
  #
  #table(datas$hbrf)
  #
  va_int<-vcovHC(fit_int2,type="HC0")
  #
  par_bmi<-c(10,13,14)
  par_smk<-c(11,15,16)
  par_hbx<-c(12,17,18)
  #
  be_bmi<-re_int[par_bmi,1]
  va_bmi<-va_int[par_bmi,par_bmi]
  #
  be_smk<-re_int[par_smk,1]
  va_smk<-va_int[par_smk,par_smk]
  #
  be_hbx<-re_int[par_hbx,1]
  va_hbx<-va_int[par_hbx,par_hbx]
  #
  g0<-c(1,0,0)
  g1<-c(1,1,0)
  g2<-c(1,0,1)
  #
  be_bmi0<-be_bmi%*%g0
  be_bmi1<-be_bmi%*%g1
  be_bmi2<-be_bmi%*%g2
  #
  be_smk0<-be_smk%*%g0
  be_smk1<-be_smk%*%g1
  be_smk2<-be_smk%*%g2
  #
  be_hbx0<-be_hbx%*%g0
  be_hbx1<-be_hbx%*%g1
  be_hbx2<-be_hbx%*%g2
  #
  se_bmi0<-sqrt(t(g0)%*%va_bmi%*%g0)
  se_bmi1<-sqrt(t(g1)%*%va_bmi%*%g1)
  se_bmi2<-sqrt(t(g2)%*%va_bmi%*%g2)
  #
  se_smk0<-sqrt(t(g0)%*%va_smk%*%g0)
  se_smk1<-sqrt(t(g1)%*%va_smk%*%g1)
  se_smk2<-sqrt(t(g2)%*%va_smk%*%g2)
  #
  se_hbx0<-sqrt(t(g0)%*%va_hbx%*%g0)
  se_hbx1<-sqrt(t(g1)%*%va_hbx%*%g1)
  se_hbx2<-sqrt(t(g2)%*%va_hbx%*%g2)
  #
  or_bmi0<-exp(be_bmi0)
  or_bmi1<-exp(be_bmi1)
  or_bmi2<-exp(be_bmi2)
  #
  or_smk0<-exp(be_smk0)
  or_smk1<-exp(be_smk1)
  or_smk2<-exp(be_smk2)
  #
  or_hbx0<-exp(be_hbx0)
  or_hbx1<-exp(be_hbx1)
  or_hbx2<-exp(be_hbx2)
  #
  or_bmi0_lb<-exp(be_bmi0-1.96*se_bmi0)
  or_bmi0_ub<-exp(be_bmi0+1.96*se_bmi0)
  or_bmi1_lb<-exp(be_bmi1-1.96*se_bmi1)
  or_bmi1_ub<-exp(be_bmi1+1.96*se_bmi1)
  or_bmi2_lb<-exp(be_bmi2-1.96*se_bmi2)
  or_bmi2_ub<-exp(be_bmi2+1.96*se_bmi2)
  #
  or_smk0_lb<-exp(be_smk0-1.96*se_smk0)
  or_smk0_ub<-exp(be_smk0+1.96*se_smk0)
  or_smk1_lb<-exp(be_smk1-1.96*se_smk1)
  or_smk1_ub<-exp(be_smk1+1.96*se_smk1)
  or_smk2_lb<-exp(be_smk2-1.96*se_smk2)
  or_smk2_ub<-exp(be_smk2+1.96*se_smk2)
  #
  or_hbx0_lb<-exp(be_hbx0-1.96*se_hbx0)
  or_hbx0_ub<-exp(be_hbx0+1.96*se_hbx0)
  or_hbx1_lb<-exp(be_hbx1-1.96*se_hbx1)
  or_hbx1_ub<-exp(be_hbx1+1.96*se_hbx1)
  or_hbx2_lb<-exp(be_hbx2-1.96*se_hbx2)
  or_hbx2_ub<-exp(be_hbx2+1.96*se_hbx2)
  #
  pbmi_chi2<-round(((wald.test(b=coef(fit_int2),Sigma=vcovHC(fit_int2,type="HC0"),Terms=par_bmi[-1])$result)$chi2)[3],4)
  psmk_chi2<-round(((wald.test(b=coef(fit_int2),Sigma=vcovHC(fit_int2,type="HC0"),Terms=par_smk[-1])$result)$chi2)[3],4)
  phbx_chi2<-round(((wald.test(b=coef(fit_int2),Sigma=vcovHC(fit_int2,type="HC0"),Terms=par_hbx[-1])$result)$chi2)[3],4)
  #
  fit_int5<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+factor(x)+bmi+smk+hbrf+
                  factor(x1)*bmi+factor(x1)*smk+factor(x1)*hbrf,family=binomial,data=datas)
  res_int5<-summary(fit_int5)$coefficients
  p_int5<-res_int5[13:15,4]
  pbmi_int<- p_int5[1]
  psmk_int<- p_int5[2]
  phbx_int<- p_int5[3]
  #
  or_bmi<-c(round(100*mrate,2),round(c(or_bmi0,or_bmi0_lb,or_bmi0_ub,or_bmi1,or_bmi1_lb,or_bmi1_ub,or_bmi2,or_bmi2_lb,or_bmi2_ub),3),pbmi_int)
  or_smk<-c(round(100*mrate,2),round(c(or_smk0,or_smk0_lb,or_smk0_ub,or_smk1,or_smk1_lb,or_smk1_ub,or_smk2,or_smk2_lb,or_smk2_ub),3),psmk_int)  
  or_hbx<-c(round(100*mrate,2),round(c(or_hbx0,or_hbx0_lb,or_hbx0_ub,or_hbx1,or_hbx1_lb,or_hbx1_ub,or_hbx2,or_hbx2_lb,or_hbx2_ub),3),phbx_int)  
  #
  fit_int3<-glm(y~age+factor(sex)+ev1_b+ev2_b+ev3_b+ev4_b+x+bmi+smk+hbrf+x*bmi+x*smk+x*hbrf,family=binomial,data=datas)
  ptrend<-round((summary(fit_int3)$coefficients)[12:14,4],4)
  
  OR<-cbind(rbind(or_bmi,or_smk,or_hbx),ptrend)
  rownames(OR)<-varname
  colnames(OR)<-c("MR","OR0","OR0_lb","OR0_ub","OR1","OR1_lb","OR1_ub","OR2","OR2_lb","OR2_ub","Prob","Ptrend")
  OR
}
#
#*********************************************************************************************************
#  
#
#*********************************************************************************************************
#  
#
BE_Table1<-Table1(1)
EA_Table1<-Table1(2)
BE_Table2<-Table2(1)
EA_Table2<-Table2(2)
EABE_Table1<-Table1(3)
EABE_Table2<-Table2(3)
#
list(BE_Table1=BE_Table1,EA_Table1=EA_Table1,EABE_Table1=EABE_Table1,
     BE_Table2=BE_Table2,EA_Table2=EA_Table2,EABE_Table2=EABE_Table2)
}
#
#
#*********************************************************************************************************
#  
#
#*********************************************************************************************************
#  
#
SNP1<-Tables(1)
SNP2<-Tables(2)
SNP3<-Tables(3)
SNP4<-Tables(4)
SNP5<-Tables(5)
SNP6<-Tables(6)
SNP7<-Tables(7)
#
#   Table 1
BE_Table1_SNP1<-SNP1$BE_Table1
BE_Table1_SNP2<-SNP2$BE_Table1
BE_Table1_SNP3<-SNP3$BE_Table1
BE_Table1_SNP4<-SNP4$BE_Table1
BE_Table1_SNP5<-SNP5$BE_Table1
BE_Table1_SNP6<-SNP6$BE_Table1
BE_Table1_SNP7<-SNP7$BE_Table1
#
EA_Table1_SNP1<-SNP1$EA_Table1
EA_Table1_SNP2<-SNP2$EA_Table1
EA_Table1_SNP3<-SNP3$EA_Table1
EA_Table1_SNP4<-SNP4$EA_Table1
EA_Table1_SNP5<-SNP5$EA_Table1
EA_Table1_SNP6<-SNP6$EA_Table1
EA_Table1_SNP7<-SNP7$EA_Table1
#
EABE_Table1_SNP1<-SNP1$EABE_Table1
EABE_Table1_SNP2<-SNP2$EABE_Table1
EABE_Table1_SNP3<-SNP3$EABE_Table1
EABE_Table1_SNP4<-SNP4$EABE_Table1
EABE_Table1_SNP5<-SNP5$EABE_Table1
EABE_Table1_SNP6<-SNP6$EABE_Table1
EABE_Table1_SNP7<-SNP7$EABE_Table1
#
#   Table 2
BE_Table2_SNP1<-SNP1$BE_Table2
BE_Table2_SNP2<-SNP2$BE_Table2
BE_Table2_SNP3<-SNP3$BE_Table2
BE_Table2_SNP4<-SNP4$BE_Table2
BE_Table2_SNP5<-SNP5$BE_Table2
BE_Table2_SNP6<-SNP6$BE_Table2
BE_Table2_SNP7<-SNP7$BE_Table2
#
EA_Table2_SNP1<-SNP1$EA_Table2
EA_Table2_SNP2<-SNP2$EA_Table2
EA_Table2_SNP3<-SNP3$EA_Table2
EA_Table2_SNP4<-SNP4$EA_Table2
EA_Table2_SNP5<-SNP5$EA_Table2
EA_Table2_SNP6<-SNP6$EA_Table2
EA_Table2_SNP7<-SNP7$EA_Table2
#
EABE_Table2_SNP1<-SNP1$EABE_Table2
EABE_Table2_SNP2<-SNP2$EABE_Table2
EABE_Table2_SNP3<-SNP3$EABE_Table2
EABE_Table2_SNP4<-SNP4$EABE_Table2
EABE_Table2_SNP5<-SNP5$EABE_Table2
EABE_Table2_SNP6<-SNP6$EABE_Table2
EABE_Table2_SNP7<-SNP7$EABE_Table2
#
#
#*********************************************************************************************************
# 
Prob_BE_Table1<-cbind(BE_Table1_SNP1[,11],BE_Table1_SNP2[,11],BE_Table1_SNP3[,11],
                      BE_Table1_SNP4[,11],BE_Table1_SNP5[,11],BE_Table1_SNP6[,11],BE_Table1_SNP7[,11])
Ptrd_BE_Table1<-cbind(BE_Table1_SNP1[,12],BE_Table1_SNP2[,12],BE_Table1_SNP3[,12],
                      BE_Table1_SNP4[,12],BE_Table1_SNP5[,12],BE_Table1_SNP6[,12],BE_Table1_SNP7[,12])
# 
Prob_EA_Table1<-cbind(EA_Table1_SNP1[,11],EA_Table1_SNP2[,11],EA_Table1_SNP3[,11],
                      EA_Table1_SNP4[,11],EA_Table1_SNP5[,11],EA_Table1_SNP6[,11],EA_Table1_SNP7[,11])
Ptrd_EA_Table1<-cbind(EA_Table1_SNP1[,12],EA_Table1_SNP2[,12],EA_Table1_SNP3[,12],
                      EA_Table1_SNP4[,12],EA_Table1_SNP5[,12],EA_Table1_SNP6[,12],EA_Table1_SNP7[,12])
#

corr_BE<-cbind(BE_Table1_SNP1[,13],BE_Table1_SNP2[,13],BE_Table1_SNP3[,13],BE_Table1_SNP4[,13],
               BE_Table1_SNP5[,13],BE_Table1_SNP6[,13],BE_Table1_SNP7[,13])
corr_EA<-cbind(EA_Table1_SNP1[,13],EA_Table1_SNP2[,13],EA_Table1_SNP3[,13],EA_Table1_SNP4[,13],
               EA_Table1_SNP5[,13],EA_Table1_SNP6[,13],EA_Table1_SNP7[,13])
colnames(corr_BE)<-rsids
colnames(corr_EA)<-rsids
#
#write.csv(corr_BE,file="EABE-GxE-BEACON-interaction-correlation-BE-04-21-2015.csv")
#write.csv(corr_EA,file="EABE-GxE-BEACON-interaction-correlation-EA-04-21-2015.csv")
#
#*********************************************************************************************************
# 
#
ntest<-21
SNP_s<-rbind(rsids,rsids,rsids)
Prob_Table1<-c(as.vector(Prob_BE_Table1),as.vector(Ptrd_BE_Table1),as.vector(Prob_EA_Table1),as.vector(Ptrd_EA_Table1))
SNP_Table1<-c(as.vector(SNP_s),as.vector(SNP_s),as.vector(SNP_s),as.vector(SNP_s))
#
CC_Table1<-c(rep("BE",42),rep("EA",42))
Prob_Type<-c(rep("Prob",21),rep("Ptrd",21),rep("Prob",21),rep("Ptrd",21))
Exposure<-as.vector(cbind(varname,varname,varname,varname,varname,varname,varname))
Exposures<-c(rep(Exposure,4))
#
#
Prob_BE_Table2<-cbind(BE_Table2_SNP1[,11],BE_Table2_SNP2[,11],BE_Table2_SNP3[,11],
                      BE_Table2_SNP4[,11],BE_Table2_SNP5[,11],BE_Table2_SNP6[,11],BE_Table2_SNP7[,11])
Ptrd_BE_Table2<-cbind(BE_Table2_SNP1[,12],BE_Table2_SNP2[,12],BE_Table2_SNP3[,12],
                      BE_Table2_SNP4[,12],BE_Table2_SNP5[,12],BE_Table2_SNP6[,12],BE_Table2_SNP7[,12])
#
Prob_EA_Table2<-cbind(EA_Table2_SNP1[,11],EA_Table2_SNP2[,11],EA_Table2_SNP3[,11],
                      EA_Table2_SNP4[,11],EA_Table2_SNP5[,11],EA_Table2_SNP6[,11],EA_Table2_SNP7[,11])
Ptrd_EA_Table2<-cbind(EA_Table2_SNP1[,12],EA_Table2_SNP2[,12],EA_Table2_SNP3[,12],
                      EA_Table2_SNP4[,12],EA_Table2_SNP5[,12],EA_Table2_SNP6[,12],EA_Table2_SNP7[,12])
#
#
colnames(Prob_BE_Table1)<-rsids
colnames(Ptrd_BE_Table1)<-rsids
colnames(Prob_EA_Table1)<-rsids
colnames(Ptrd_EA_Table1)<-rsids
# 
colnames(Prob_BE_Table2)<-rsids
colnames(Ptrd_BE_Table2)<-rsids
colnames(Prob_EA_Table2)<-rsids
colnames(Ptrd_EA_Table2)<-rsids
#
#*********************************************************************************************************
# 
#
#*********************************************************************************************************
# 
#
#
fdr<-function(pCrude){
#
tt<-length(pCrude)
pBonferroni<-pmin(tt*pCrude,1)
pFDR<-rep(0,tt);for(i in 1:tt){
  #i<-64;tt*pCrude[i]/sum(as.numeric(pCrude<=pCrude[i]))
  pFDR[i]<-min(tt*pCrude[i]/sum(as.numeric(pCrude<=pCrude[i])),1)}
col<-1:tt
pCrude_ord<-pCrude[order(pCrude,decreasing=FALSE)]

FDR_val<-cbind(pCrude,pBonferroni,pFDR)
FDR_sig<-cbind(col,pCrude[order(pCrude,decreasing=FALSE)],pFDR[order(pCrude,decreasing=FALSE)],col*0.05/tt)
list(FDR_val=FDR_val,FDR_sig=FDR_sig)
}
#
#
#
#*********************************************************************************************************
# 
#   Table 1
#
#*********************************************************************************************************
#
#
SNP_s<-rbind(rsids,rsids,rsids)
Prob_Table1<-c(as.vector(Prob_BE_Table1),as.vector(Prob_EA_Table1))
SNP_Table1<-c(as.vector(SNP_s),as.vector(SNP_s))
#
CC_Table1<-c(rep("BE",21),rep("EA",21))
Prob_Type<-c(rep("Prob",21),rep("Prob",21))
Exposure<-as.vector(cbind(varname,varname,varname,varname,varname,varname,varname))
Exposures<-c(rep(Exposure,2))
#
#
FDR_prob<-fdr(Prob_Table1)$FDR_val
FDR_Table1_prob<-data.frame(CC_Table1,Exposures,SNP_Table1,Prob_Type,Prob_Table1,FDR_prob)
#
FDR_Table1_prob[SNP_Table1=="rs2687201",]
FDR_Table1_prob[SNP_Table1=="rs10419226",]
#
fdr(Prob_Table1)$FDR_sig
fdr(Ptrd_Table1)$FDR_sig
#
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#

SNP_s<-rbind(rsids,rsids,rsids)
Ptrd_Table1<-c(as.vector(Ptrd_BE_Table1),as.vector(Ptrd_EA_Table1))
SNP_Table1<-c(as.vector(SNP_s),as.vector(SNP_s))
#
CC_Table1<-c(rep("BE",21),rep("EA",21))
Prob_Type<-c(rep("Ptrd",21),rep("Ptrd",21))
Exposure<-as.vector(cbind(varname,varname,varname,varname,varname,varname,varname))
Exposures<-c(rep(Exposure,2))
#
#
FDR_ptrd<-fdr(Ptrd_Table1)$FDR_val
FDR_Table1_ptrd<-data.frame(CC_Table1,Exposures,SNP_Table1,Prob_Type,Ptrd_Table1,FDR_ptrd)
#
FDR_Table1_ptrd[SNP_Table1=="rs2687201",]
FDR_Table1_ptrd[SNP_Table1=="rs10419226",]
#
#
library('stats')
library("fdrtool")
library("qvalue")
#
fdr_prob_pa<-p.adjust(Prob_Table1, "fdr")
fdr_ptrd_pa<-p.adjust(Ptrd_Table1, "fdr")
#
fdr_prob_ft<-fdrtool(Prob_Table1, statistic="Prob_Table1")$qval
fdr_ptrd_ft<-fdrtool(Ptrd_Table1, statistic="Ptrd_Table1")$qval
#
#
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#
snp1_ea<-as.numeric(CC_Table1=="EA" & SNP_Table1=="rs2687201")
snp1_be<-as.numeric(CC_Table1=="BE" & SNP_Table1=="rs2687201"  & Exposure=="recurrent_HB_RF")
snp2_ea<-as.numeric(CC_Table1=="EA" & SNP_Table1=="rs10419226" & Exposure=="recurrent_HB_RF")
#
ord_prob<-order(Prob_Table1,decreasing=F)
snp1_be_prob<-snp1_be[ord_prob]
snp1_ea_prob<-snp1_ea[ord_prob]
snp2_ea_prob<-snp2_ea[ord_prob]
#
ord_ptrd<-order(Ptrd_Table1,decreasing=F)
snp1_be_ptrd<-snp1_be[ord_ptrd]
snp1_ea_ptrd<-snp1_ea[ord_ptrd]
snp2_ea_ptrd<-snp2_ea[ord_ptrd]
#
#
Prob_1<-Prob_Table1[ord_prob]
Ptrd_1<-Ptrd_Table1[ord_ptrd]
#
#
#
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#
#
# QQ-plot
#
ggd.qqplot_prob = function(pvector, main=NULL, ...) {
  #
  o = -log10(sort(pvector,decreasing=F))
  #
  xe<-rep(0,length(o))
  for(i in 1:length(o)){
      ue<-rep(0,1000)
      for(j in 1:1000){
      ue[j]<-sort(runif(length(o)),decreasing=F)[i]}
    xe[i]<-mean(ue)
  }
  #
  e<--log10(xe)
  #e = -log10(1:length(o)/length(o))
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,2), ylim=c(0,3.8))
  lines(e,e,col="blue")
  points(e[snp1_be_prob==1],o[snp1_be_prob==1],col="red",pch=19,cex=1)
  #points(e[snp1_ea==1],o[snp1_ea==1],col="brown",pch=19,cex=1)
  points(e[snp2_ea_prob==1],o[snp2_ea_prob==1],col="brown",pch=19,cex=1)
  legend(0,3.8,c("Interaction between rs2687201 and weekly heartburn or reflux in BE", 
  #              "Interaction between rs2687201 and risk factors in EA", 
                 "Interaction between rs10419226 and weekly heartburn or reflux in EA", 
                 "Other interactions with risk factors"), col=c("red","brown","black"),pch=19,cex=1)
}
#cyan4
#
ggd.qqplot_ptrd = function(pvector, main=NULL, ...) {
  #
  o = -log10(sort(pvector,decreasing=F))
  #
  xe<-rep(0,length(o))
  for(i in 1:length(o)){
    ue<-rep(0,1000)
    for(j in 1:1000){
      ue[j]<-sort(runif(length(o)),decreasing=F)[i]}
    xe[i]<-mean(ue)
  }
  #
  e<--log10(xe)
  #e = -log10(1:length(o)/length(o))
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,2), ylim=c(0,3.8))
  lines(e,e,col="blue")
  points(e[snp1_be_ptrd==1],o[snp1_be_ptrd==1],col="red",pch=19,cex=1)
  #points(e[snp1_ea==1],o[snp1_ea==1],col="brown",pch=19,cex=1)
  points(e[snp2_ea_ptrd==1],o[snp2_ea_ptrd==1],col="brown",pch=19,cex=1)
  legend(0,3.8,c("Interaction between rs2687201 and weekly heartburn or reflux in BE", 
                 #              "Interaction between rs2687201 and risk factors in EA", 
                 "Interaction between rs10419226 and weekly heartburn or reflux in EA", 
                 "Other interactions with risk factors"), col=c("red","brown","black"),pch=19,cex=1)
}
#cyan4
#tiff("EABE-GxE-qqplot-04-15-2015.tif", width=8,height=10,units="in",compression='lzw',res=600)
#ggd.qqplot(Prob_Table1, " " )
#
#dev.off()
#
#
#******************************************************************************************************
#
#
tiff("EABE-GxE-qqplot-7snps-42 tests-04-17-2015.tif", width=16,height=10,units="in",compression='lzw',res=600)
par(mfrow=c(1,2))
#ggd.qqplot0(Prob_Table1, "(a)" )
#ggd.qqplot0(Ptrd_Table1, "(b)" )
ggd.qqplot_prob(Prob_1, "(a)" )
ggd.qqplot_ptrd(Ptrd_1, "(b)" )
#
dev.off()
#
#
#
#******************************************************************************************************
#
#
#******************************************************************************************************
#
#
#prob_random<-Prob_Table1
#
ggd.qqplot0 = function(pvector, main=NULL, ...) {
  #
  o = -log10(sort(pvector,decreasing=F))
  #
  xe<-rep(0,length(o))
  for(i in 1:length(o)){
    ue<-rep(0,1000)
    for(j in 1:1000){
      ue[j]<-sort(runif(length(o)),decreasing=F)[i]}
    xe[i]<-mean(ue)
  }
  #
  e<--log10(xe)
  #e = -log10(1:length(o)/length(o))
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,2), ylim=c(0,3.8))
  lines(e,e,col="blue")
}
#cyan4
#
#******************************************************************************************************
#
#tiff("EABE-GxE-qqplot-random-a-04-17-2015.tif", width=16,height=10,units="in",compression='lzw',res=600)
#ggd.qqplot0(Prob_Table1, "(a)" )
#
#dev.off()
#
#******************************************************************************************************
#
#
#******************************************************************************************************
#
#
#
tiff("GEinteraction-BEACON-qqplot-7reported-random-SNPs-42Tests-04-17-2015.tif", 
                        width=16,height=10,units="in",compression='lzw',res=600)
par(mfrow=c(2,2))
ggd.qqplot_prob(Prob_1,  "(a)" )
ggd.qqplot0(Prob_Table1, "(b)" )
ggd.qqplot_ptrd(Ptrd_1,  "(c)" )
ggd.qqplot0(Ptrd_Table1, "(d)" )
#
dev.off()
#
#
#
#******************************************************************************************************
#
#
#******************************************************************************************************
#
#























k<-5
#
fdr_prb_be_table1<-fdr(Prob_BE_Table1[k,])
fdr_trd_be_table1<-fdr(Ptrd_BE_Table1[k,])
#
fdr_prb_ea_table1<-fdr(Prob_EA_Table1[k,])
fdr_trd_ea_table1<-fdr(Ptrd_EA_Table1[k,])
#
prb_be_Table1<-c(Prob_BE_Table1[1,],Prob_BE_Table1[2,],Prob_BE_Table1[3,],
                 Prob_BE_Table1[4,],Prob_BE_Table1[5,])
trd_be_Table1<-c(Ptrd_BE_Table1[1,],Ptrd_BE_Table1[2,],Ptrd_BE_Table1[3,],
                  Ptrd_BE_Table1[4,],Ptrd_BE_Table1[5,])
#
fdr(prb_be_Table1)
#length(prb_be_Table1)

fdr_trd_be_table1<-fdr(Ptrd_BE_Table1[k,])

#
#*********************************************************************************************************
# 
#  Table 1: BE
fdr_prb_be_table1$FDR_val[1,c(1,3)]
fdr_trd_be_table1$FDR_val[1,c(1,3)]
#
#   Table 1: EA
fdr_prb_ea_table1$FDR_val[1,c(1,3)]
fdr_trd_ea_table1$FDR_val[1,c(1,3)]
#
#p.adjust(Prob_BE_Table1[k,],method="fdr",n=6)
#p.adjust(Ptrd_BE_Table1[k,],method="fdr",n=6)
#
#
#
#
#
#*********************************************************************************************************
# 
#   Table 2
#
#*********************************************************************************************************
# 
# 
k<-3
fdr_prb_be_table2<-fdr(Prob_BE_Table2[k,])
fdr_trd_be_table2<-fdr(Ptrd_BE_Table2[k,])
#
fdr_prb_ea_table2<-fdr(Prob_EA_Table2[k,])
fdr_trd_ea_table2<-fdr(Ptrd_EA_Table2[k,])
#
#*********************************************************************************************************
# 
#  Table 2: BE
fdr_prb_be_table2$FDR_val[1,c(1,3)]
fdr_trd_be_table2$FDR_val[1,c(1,3)]
#
#   Table 2: EA
fdr_prb_ea_table2$FDR_val[1,c(1,3)]
fdr_trd_ea_table2$FDR_val[1,c(1,3)]
#
#
#*********************************************************************************************************
#
#
#*********************************************************************************************************
#
#
#
#
#
l<-1;
SNP<-Tables(l)
BE1<-SNP$BE_Table1
EA1<-SNP$EA_Table1
BE2<-SNP$BE_Table2
EA2<-SNP$EA_Table2
#
write.csv(BE1,file="EABE-OR-BE-Table1-SNP1-02-06-2015.csv")
write.csv(EA1,file="EABE-OR-EA-Table1-SNP1-02-06-2015.csv")
write.csv(BE2,file="EABE-OR-BE-Table2-SNP1-02-06-2015.csv")
write.csv(EA2,file="EABE-OR-EA-Table2-SNP1-02-06-2015.csv")
#
#*********************************************************************************************************
#  
#
#*********************************************************************************************************
#  
