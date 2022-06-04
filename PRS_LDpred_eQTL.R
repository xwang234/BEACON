#!/usr/bin/env Rscript
# generate data for LDpred,for imputed beacon data 
library(data.table)
#step0 meta analysis
BEmeta=as.data.frame(fread("../result/Bonn_Oxford_Cambridge_METAANALYSIS_BE_comsnp1.tbl"))
colnames(BEmeta)[which(colnames(BEmeta)=="P-value")]="P"
BECambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_comsnp_N.txt"))
idx=match(BEmeta$MarkerName,BECambridge$SNP)
BEmeta$chr=BECambridge$CHR[idx]
BEmeta$pos=BECambridge$position[idx]
write.table(BEmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt",col.names = T,
            row.names=F,quote=F,sep="\t")
cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt"
system(cmd)

EAmeta=as.data.frame(fread("../result/Bonn_Cambridge_METAANALYSIS_EA_comsnp1.tbl"))
colnames(EAmeta)[which(colnames(EAmeta)=="P-value")]="P"
EACambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_comsnp_N.txt"))
idx=match(EAmeta$MarkerName,EACambridge$SNP)
EAmeta$chr=EACambridge$CHR[idx]
EAmeta$pos=EACambridge$position[idx]
write.table(EAmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_Cambridge_metastat.txt",col.names = T,
            row.names=F,quote=F,sep="\t")
cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_Cambridge_metastat.txt"
system(cmd)

BEEAmeta=as.data.frame(fread("../result/Bonn_Cambridge_METAANALYSIS_BEEA_comsnp1.tbl"))
colnames(BEEAmeta)[which(colnames(BEEAmeta)=="P-value")]="P"
BEEACambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_comsnp_N.txt"))
idx=match(BEEAmeta$MarkerName,BEEACambridge$SNP)
BEEAmeta$chr=BEEACambridge$CHR[idx]
BEEAmeta$pos=BEEACambridge$position[idx]
write.table(BEEAmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt",col.names = T,
            row.names=F,quote=F,sep="\t")
system("rm /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt.gz")
cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt"
system(cmd)


#step2-------
#remove Cambridge from BCA
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
#site 30:Cambridge, 55:AMOS
idx=which(!sampletable$site %in% c(30,55))
beaconsamples=sampletable[idx,1:2]
table(sampletable$phenoBE_bca[idx])
# -9    1    2 
# 1683 2182 2413 
table(sampletable$phenoEA_bca[idx])
# -9    1    2 
# 2581 2185 1512 
table(sampletable$phenoEABE_bca[idx])
# -9    1    2 
# 174 2183 3921 
write.table(beaconsamples,file="../result/beaconsamples_plink.txt",row.names = F,col.names = F,sep="\t",quote=F)

#step 1---get eQTLs
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_junction=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_blood=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_stomach=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
alleqtl=unique(c(snp_junction,snp_blood,snp_stomach))
alleqtl=alleqtl[!is.na(alleqtl)]
#overlift hg38-->hg19
#hg19lookuptable=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt"))
hg38lookuptable=as.data.frame(fread(input="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/V8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz"))
hg38lookuptable0=hg38lookuptable
hg38lookuptable=hg38lookuptable[!is.na(hg38lookuptable$variant_id_b37),]
tmp=gsub("chr","",hg38lookuptable$chr)
hg38name=paste0(tmp,":",hg38lookuptable$variant_pos,"_",hg38lookuptable$ref,"_",hg38lookuptable$alt)
hg38name1=paste0(tmp,":",hg38lookuptable$variant_pos,"_",hg38lookuptable$alt,"_",hg38lookuptable$ref)
sum(alleqtl %in% hg38name) #220587
sum(alleqtl %in% hg38name1) #737458
tmp1=unlist(unlist(strsplit(hg38lookuptable$variant_id_b37,"_")))
hg19chr=tmp1[seq(1,length(tmp1),5)]
hg19pos=tmp1[seq(2,length(tmp1),5)]
sum(hg38lookuptable$variant_pos==hg19pos) #392935,0.85% of all
hg19ref=tmp1[seq(3,length(tmp1),5)]
hg19alt=tmp1[seq(4,length(tmp1),5)]
hg19name=paste0(hg19chr,":",hg19pos,"_",hg19ref,"_",hg19alt)
hg19posname=paste0(hg19chr,":",hg19pos)
hg19name1=paste0(hg19chr,":",hg19pos,"_",hg19alt,"_",hg19ref)
idx1=which(hg38name %in% alleqtl)
idx2=which(hg38name1 %in% alleqtl)
idx=c(idx1,idx2)
idx=idx[order(idx)]
hg19eqtl=hg19name[idx]
hg19eqtl1=hg19name1[idx]
hg19eqtlposname=hg19posname[idx]
beaconhg19=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous.bim"))
sum(beaconhg19$V2 %in% hg19eqtl) #192514
sum(beaconhg19$V2 %in% hg19eqtl1) #628083
beaconhg19posname=paste0(beaconhg19$V1,":",beaconhg19$V4)
sum(beaconhg19posname %in% hg19eqtlposname) #820599
idx1=which((beaconhg19posname %in% hg19eqtlposname))
tmp=data.frame(snp=beaconhg19$V2[idx1])
fwrite(tmp,file="../result/beacon_eQTLkeepSNPs.txt",row.names = F,col.names = F,quote=F)

#step3
keep_nohetsamples=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL_QC")
{
  dat <- read.table(paste0(prefix,".het"), header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  write.table(valid[,c(1,2)], paste0(prefix,".valid.sample"), quote=F, row.names=F) 
}

#step5----
generate_plink=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_eQTL_QC",
                        opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_eQTL",metadat=BEmeta)
{
  #
  bim <- read.table(paste0(prefix,".bim"),stringsAsFactors = F)
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
  metadat$Allele1=toupper(metadat$Allele1) #Allele1 is non-effect
  metadat$Allele2=toupper(metadat$Allele2)
  metadat$A2=metadat$Allele1
  metadat$A1=metadat$Allele2
  colnames(metadat)[which(colnames(metadat)=="MarkerName")]="SNP"
  colnames(metadat)[which(colnames(metadat)=="chr")]="CHR"
  colnames(metadat)[which(colnames(metadat)=="pos")]="BP"
  #update bim SNP ID based on metadat
  bim_posname1=paste0(bim$CHR,":",bim$BP,"_",bim$B.A1,"_",bim$B.A2)
  bim_posname2=paste0(bim$CHR,":",bim$BP,"_",bim$B.A2,"_",bim$B.A1)
  meta_posname=paste0(metadat$CHR,":",metadat$BP,"_",metadat$Allele1,"_",metadat$Allele2)
  idx1=which(bim_posname1 %in% meta_posname)
  idx2=which(bim_posname2 %in% meta_posname)
  tmp=bim_posname1[idx1]
  idx3=match(tmp,bim_posname1)
  idx4=match(tmp,meta_posname)
  bim$SNP[idx3]=metadat$SNP[idx4]
  tmp=bim_posname2[idx2]
  idx3=match(tmp,bim_posname2)
  idx4=match(tmp,meta_posname)
  bim$SNP[idx3]=metadat$SNP[idx4]
  print(sum(bim$SNP %in% metadat$SNP))#number of overlapped SNPs
  write.table(bim,file=paste0(prefix1,".bim"),col.names = F,row.names = F,sep=" ",quote=F)
  # #Identify SNPs that require strand flipping
  # info <- merge(bim, metadat, by = c("SNP", "CHR", "BP"))
  # # Function for finding the complementary allele
  # complement <- function(x) {
  #   switch (
  #     x,
  #     "A" = "T",
  #     "C" = "G",
  #     "T" = "A",
  #     "G" = "C",
  #     return(NA)
  #   )
  # }
  # # Get SNPs that have the same alleles across base and target
  # info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
  # # Identify SNPs that are complementary between base and target
  # info$C.A1 <- sapply(info$B.A1, complement)
  # info$C.A2 <- sapply(info$B.A2, complement)
  # info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
  # # Update the complementary alleles in the bim file
  # # This allow us to match the allele in subsequent analysis
  # complement.snps <- bim$SNP %in% info.complement$SNP
  # bim[complement.snps,]$B.A1 <-
  #   sapply(bim[complement.snps,]$B.A1, complement)
  # bim[complement.snps,]$B.A2 <-
  #   sapply(bim[complement.snps,]$B.A2, complement)
  # 
  # 
  # # Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)
  # info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
  # # Update the recode SNPs
  # recode.snps <- bim$SNP %in% info.recode$SNP
  # tmp <- bim[recode.snps,]$B.A1
  # bim[recode.snps,]$B.A1 <- bim[recode.snps,]$B.A2
  # bim[recode.snps,]$B.A2 <- tmp
  # 
  # # identify SNPs that need recoding & complement
  # info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
  # # Update the recode + strand flip SNPs
  # com.snps <- bim$SNP %in% info.crecode$SNP
  # tmp <- bim[com.snps,]$B.A1
  # bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
  # bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
  # 
  # # Output updated bim file
  # write.table(
  #   bim,
  #   paste0(prefix1,"_QC.adj.bim"),
  #   quote = F,
  #   row.names = F,
  #   col.names = F,
  #   sep="\t"
  # )
  # cmd=paste0("cp ",prefix,".bim ",prefix1,".bim.bk")
  # system(cmd,wait=T)
  # cmd=paste0("cp ",prefix1,"_QC.adj.bim ",prefix1,".bim")
  # system(cmd,wait=T)
  
  for (ext in c(".bed",".fam"))
  {
    #file.copy(paste0(prefix,ext),paste0(prefix1,ext))
    cmd=paste0("cp ",paste0(prefix,ext)," ",paste0(prefix1,ext))
    system(cmd,wait = T)
  }
  fam=read.table(paste0(prefix1,".fam"))
  idx=match(fam$V2,sampletable$localid)
  if (opt=="BE")
  {
    fam$V6=sampletable$phenoBE_bca[idx]
  }
  if (opt=="EA")
  {
    fam$V6=sampletable$phenoEA_bca[idx]
  }
  if (opt=="BEEA")
  {
    fam$V6=sampletable$phenoEABE_bca[idx]
  }
  write.table(fam,file=paste0(prefix1,".fam"),row.names = F,col.names = F,sep=" ",quote=F)
  #generate phenotype file
  tmp=fam[,c(1,2,6)]
  write.table(tmp,file=paste0(prefix1,".pheno"),row.names = F,col.names = F,sep=" ",quote=F)
}
generate_plink(opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_eQTL")
generate_plink(opt="EA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_eQTL",metadat=EAmeta)
generate_plink(opt="BEEA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL",metadat=BEEAmeta)


generate_covariate=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_eQTL")
{
  fam=read.table(paste0(prefix,".fam"),stringsAsFactors = F)
  idx=match(fam$V2,sampletable$localid)
  Covariate=sampletable[idx,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
  colnames(Covariate)[which(colnames(Covariate)=="ev1_bca")]="pc1"
  colnames(Covariate)[which(colnames(Covariate)=="ev2_bca")]="pc2"
  colnames(Covariate)[which(colnames(Covariate)=="ev3_bca")]="pc3"
  colnames(Covariate)[which(colnames(Covariate)=="ev4_bca")]="pc4"
  rownames(Covariate)=sampletable$localid[idx]
  Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)
  Covariate=data.frame(IID=rownames(Covariate),Covariate)
  write.table(Covariate,file=paste0(prefix,".covariate"),col.names = T,row.names = F,sep=" ",quote=F)
}
generate_covariate()
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_eQTL")
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL")




#after runing PRS_LDpred.sh, read results:
#library(DescTools)
#p.threshold <- c(1,0.3,0.1)
p.threshold <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001)
library(pROC)
plot_ROC=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_genotyped",opt=1)
{
  phenotype <- read.table(paste0(prefix,".pheno"), header=F)
  colnames(phenotype)=c("FID","IID","case")
  if (opt==1)
  {
    for(i in p.threshold)
    {
      # Go through each p-value threshold .score_LDpred_p1.0000e-01.txt
      tmp=formatC(i, format = "e", digits =4)
      LDpredfile=paste0(prefix,".score_LDpred_p",tmp,".txt")
      if (file.exists(LDpredfile))
      {
        prs <- read.table(LDpredfile, header=T,sep=",",stringsAsFactors = F)
        pheno.prs <- merge(phenotype, prs, by=c("IID"))
        pheno.prs$true_phens[pheno.prs$true_phens==-9]=NA
        pheno.prs$case[pheno.prs$case==-9]=NA
        testroc<- roc(pheno.prs$true_phens, pheno.prs$PRS)
        print(paste0(i,":auc=",testroc$auc))
        #plot(testroc, print.auc=TRUE,main=tmp)
      }
    }
  }
  
  if (opt==2) #LDpred-inf
  {
    LDpredfile=paste0(prefix,".score_LDpred-inf.txt")
    prs <- read.table(LDpredfile, header=T,sep=",",stringsAsFactors = F)
    pheno.prs <- merge(phenotype, prs, by=c("IID"))
    pheno.prs$true_phens[pheno.prs$true_phens==-9]=NA
    pheno.prs$case[pheno.prs$case==-9]=NA
    testroc<- roc(pheno.prs$true_phens, pheno.prs$PRS)
    plot(testroc, print.auc=TRUE,main="inf")
  }
  
}
par(mfrow=c(2,2))
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL",opt=2)

require(MASS)
#PRS high, risk high, for genotypted
PRS_odds=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_eQTL",p=0.1,opt=1)
{
  phenotype <- read.table(paste0(prefix,".pheno"), header=F)
  colnames(phenotype)=c("FID","IID","case")
  if (opt==1)
  {
    tmp=formatC(p, format = "e", digits =4)
    LDpredfile=paste0(prefix,".score_LDpred_p",tmp,".txt.adj")
  }else
  {
    LDpredfile=paste0(prefix,".score_LDpred-inf.txt")
  }
  
  prs <- read.table(LDpredfile, header=T,sep=",")
  #prs$case=prs$true_phens
  pheno.prs <- merge(phenotype, prs, by=c("IID"))
  #quantiles:
  qt=quantile(pheno.prs$PRS)
  qtidx=cut(pheno.prs$PRS,qt)
  pheno.prs$qtidx=qtidx
  pheno.prs$group=as.numeric(qtidx)
  pheno.prs$true_phens[pheno.prs$true_phens==-9]=NA
  pheno.prs$case[pheno.prs$case==-9]=NA
  if (sum(pheno.prs$case==2,na.rm = T)>0)
  {
    pheno.prs$case[which(pheno.prs$case==1)]=0
    pheno.prs$case[which(pheno.prs$case==2)]=1
  }
  tmp=pheno.prs[,c("IID","PRS")]
  write.table(tmp,file=paste0(prefix,".PRS.txt"),sep="\t",row.names = F,quote=F)
  for (i in 2:4)
  {
    print(paste0(i,"th quantile:"))
    dat=pheno.prs[pheno.prs$group %in% c(1,i),]
    dat$group=factor(dat$group,levels = c(1,i))
    fit=glm(I(dat$case==1)~factor(dat$group),family = binomial)
    print(exp(cbind(coef(fit), confint(fit)))[2,])
  }
}

#BE---
PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL",opt=2)
