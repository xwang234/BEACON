#!/usr/bin/env Rscript
# generate data for LDpred,for imputed validation data 
#use Cambridge+BEACONcontrol as the validation, only work on EA
#BEACONcase+AMOS+other BEACONcontrol for discovery

library(data.table)
#generate disc/valid samples
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
#site 30:Cambridge, 55:AMOS
idx=which(sampletable$site<30)
table(sampletable$phenoEA_bca[idx])
# -9    1    2 
# 2566 2185 1512
idx=which(sampletable$site==30)
table(sampletable$phenoEA_bca[idx])
# -9    2 
# 936 1003 
idx=which(sampletable$site==55)
table(sampletable$phenoEA_bca[idx])
# 1 
# 1022 
table(sampletable$site,sampletable$phenoEA_bca)
#      -9    1    2
# 11  174    0  102
# 12  303  214    0
# 13    2  116   63
# 14  822    0  502
# 15    0   88   42
# 16  202  218  193
# 17  332  323    0
# 18   36  259  247
# 19  101    0    0
# 20  160  167    0
# 21  122   92   54
# 22    2   26   14
# 23    8  437   59
# 25    0  245  236
# 27  302    0    0
# 30  936    0 1003
# 55    0 1022    0
# NA   15    0    0
#sites: control will be picked to combine with Cambridge to form validationset
controlsites=c(12,17,20,23)
table(sampletable$phenoEA_bca[sampletable$site %in% controlsites])
# -9    1    2 
# 803 1141   59 
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
validation_cases=sampletable$localid[sampletable$site %in% 30 & !is.na(sampletable$phenoEA_bca)]
validation_controls=sampletable$localid[which(sampletable$site %in% controlsites & sampletable$phenoEA_bca==1)]
validationsamples=c(validation_cases,validation_controls)
idx=match(validationsamples,sampletable$localid)
discoverysamples=sampletable$localid[!is.na(sampletable$phenoEA_bca)]
discoverysamples=discoverysamples[! discoverysamples %in% validationsamples]
idx=match(discoverysamples,sampletable$localid)
table(sampletable$phenoEA_bca[idx])
# 1    2 
# 2066 1512 
idx=which(sampletable$localid %in% discoverysamples)
sampletable1=sampletable[idx,]
discsamples=sampletable[idx,1:2]
table(sampletable$phenoBE_bca[idx])
# 1 
# 2065 
table(sampletable$phenoEA_bca[idx])
# 1    2 
# 2066 1512 
#write.table(discsamples,file="../result/discsamples_plink.txt",row.names = F,col.names = F,sep="\t",quote=F)

EAsamples=sampletable1$localid[which(sampletable1$phenoEA_bc==2)]
COsamples=sampletable1$localid[which(sampletable1$phenoEA_bc==1)]
fam=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/bca_filter_noambiguous.fam",header = F,stringsAsFactors = F)
COsamples=COsamples[COsamples %in% fam$V2]
EAsamples=EAsamples[EAsamples %in% fam$V2]
filesforGWAS=function(casesamples=EAsamples,controlsamples=COsamples,prefix="Discovery_EA_CO")
{
  idx1=match(casesamples,fam$V2)
  tmp1=data.frame(FID=fam$V1[idx1],IID=fam$V2[idx1],affected=2,stringsAsFactors = F)
  idx2=match(controlsamples,fam$V2)
  tmp2=data.frame(FID=fam$V1[idx2],IID=fam$V2[idx2],affected=1,stringsAsFactors = F)
  tmp=rbind(tmp1,tmp2)
  idx=match(tmp$IID,sampletable1$localid)
  tmp$pc1=sampletable1$ev1_bca[idx]
  tmp$pc2=sampletable1$ev2_bca[idx]
  tmp$pc3=sampletable1$ev3_bca[idx]
  tmp$pc4=sampletable1$ev4_bca[idx]
  tmp$age=sampletable1$age[idx]
  tmp$sex=sampletable1$sex[idx]
  
  #for validate_twas, ind ID changed
  tmp1=tmp
  tmp1$IID=paste0(tmp$FID,"_",tmp$IID)
  tmp1$FID=0
  write.table(tmp1[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink1.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp[,c(1,2)],file=paste0("../result/GWAS/",prefix,"_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp,file=paste0("../result/GWAS/",prefix,"_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
  write.table(tmp1,file=paste0("../result/GWAS/",prefix,"_selectedsamples_pheno_plink1.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
filesforGWAS()
#check files:
tmp=read.table("../result/GWAS/Discovery_EA_CO_gwas.fam")
tmp1=read.table("../result/GWAS/Discovery_EA_CO_selectedsamples_pheno_plink.txt",header = T)

all(tmp$V2 %in% tmp1$IID) #T
udpate_fam=function(famfile="../result/GWAS/Discovery_EA_CO_gwas.fam",phenofile="../result/GWAS/Discovery_EA_CO_selectedsamples_pheno_plink.txt")
{
  fam=read.table(famfile,stringsAsFactors = F)
  pheno=read.table(phenofile,header=T,stringsAsFactors = F)
  idx=match(fam$V2,pheno$IID)
  fam$V5=pheno$sex[idx]
  fam$V6=pheno$affected[idx]
  write.table(fam,file=famfile,col.names=F,row.names = F,quote=F)
}
udpate_fam()
#after running gwas, update the result

#after running gwas, update the result
update_gwas=function(gwasfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/Discovery_EA_CO_gwas.assoc.logistic",
                     bimfile="../result/GWAS/Discovery_EA_CO_gwas.bim")
{
  gwas=data.frame(fread(gwasfile,sep=" "))
  #check gwas
  # tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Beacon_autosomes_comsnp_N.txt"))
  # tmp1=merge(gwas,tmp,by="SNP")
  
  bim=data.frame(fread(bimfile))
  all(gwas$SNP==bim$V2)
  all(gwas$A1 %in% bim$V5)
  idx=match(bim$V2,gwas$SNP)
  gwas$A2=bim$V6[idx]
  colnames(gwas)[colnames(gwas)=="A1"]="effect_allele"
  colnames(gwas)[colnames(gwas)=="A2"]="non_effect_allele"
  colnames(gwas)[colnames(gwas)=="BETA"]="beta"
  colnames(gwas)[colnames(gwas)=="SE"]="se"
  colnames(gwas)[colnames(gwas)=="BP"]="position"
  gwas$N=gwas$NMISS
  gwas1=as.data.frame(fread("/fh/fast/dai_j/BEACON/Meta_summary_stat/EA_BEACON_autosomes.txt"))
  posname1_gwas=paste0(gwas$CHR,"_",gwas$position,"_",gwas$effect_allele,"_",gwas$non_effect_allele)
  posname2_gwas=paste0(gwas$CHR,"_",gwas$position,"_",gwas$non_effect_allele,"_",gwas$effect_allele)
  posname1_gwas1=paste0(gwas1$CHR,"_",gwas1$position,"_",gwas1$effect_allele,"_",gwas1$non_effect_allele)
  tmp=intersect(posname1_gwas, posname1_gwas1)
  idx1=match(tmp,posname1_gwas)
  idx2=match(tmp,posname1_gwas1)
  gwas$SNP[idx1]=gwas1$SNP[idx2]
  tmp=intersect(posname2_gwas, posname1_gwas1)
  idx1=match(tmp,posname2_gwas)
  idx2=match(tmp,posname1_gwas1)
  gwas$SNP[idx1]=gwas1$SNP[idx2]
  
  
  idx=is.na(gwas$P)
  gwas=gwas[!idx,]
  fwrite(gwas,file=gwasfile,quote=F,sep=" ")
}

#run Meta_Analysis.R and Meta_Analysis.sh to get meta results





#step0 meta analysis
# BEmeta=as.data.frame(fread("../result/Bonn_Oxford_Cambridge_METAANALYSIS_BE_comsnp1.tbl"))
# colnames(BEmeta)[which(colnames(BEmeta)=="P-value")]="P"
# BECambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_comsnp_N.txt"))
# idx=match(BEmeta$MarkerName,BECambridge$SNP)
# BEmeta$chr=BECambridge$CHR[idx]
# BEmeta$pos=BECambridge$position[idx]
# write.table(BEmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt",col.names = T,
#             row.names=F,quote=F,sep="\t")
# cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/Bonn_Oxford_Cambridge_metastat.txt"
# system(cmd)
# 
EAmeta=as.data.frame(fread("../result/Discovery_Bonn_METAANALYSIS_EA_comsnp1.tbl"))
colnames(EAmeta)[which(colnames(EAmeta)=="P-value")]="P"
EABeacon=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_BD_autosomes_comsnp_N.txt"))
idx=match(EAmeta$MarkerName,EABeacon$SNP)
EAmeta$chr=EABeacon$CHR[idx]
EAmeta$pos=EABeacon$position[idx]
write.table(EAmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_Bonn_imp_metastat.txt",col.names = T,
            row.names=F,quote=F,sep="\t")
cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_Bonn_imp_metastat.txt"
system(cmd)
# 
# BEEAmeta=as.data.frame(fread("../result/Bonn_Cambridge_METAANALYSIS_BEEA_comsnp1.tbl"))
# colnames(BEEAmeta)[which(colnames(BEEAmeta)=="P-value")]="P"
# BEEACambridge=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_comsnp_N.txt"))
# idx=match(BEEAmeta$MarkerName,BEEACambridge$SNP)
# BEEAmeta$chr=BEEACambridge$CHR[idx]
# BEEAmeta$pos=BEEACambridge$position[idx]
# write.table(BEEAmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt",col.names = T,
#             row.names=F,quote=F,sep="\t")
# system("rm /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt.gz")
# cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_Cambridge_metastat.txt"
# system(cmd)

#step2-------
keep_nohetsamples=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validation_filter_noambiguous_QC")
{
  dat <- read.table(paste0(prefix,".het"), header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  write.table(valid[,c(1,2)], paste0(prefix,".valid.sample"), quote=F, row.names=F) 
}


generate_plink=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/validation_filter_noambiguous_QC",
                        opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE",metadat=BEmeta)
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
#generate_plink(opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE")
generate_plink(opt="EA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA",metadat=EAmeta)
#generate_plink(opt="BEEA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA",metadat=BEEAmeta)



generate_covariate=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE")
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
#generate_covariate()
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA")
#generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA")




#after runing PRS_LDpred.sh, read results:
#library(DescTools)
#p.threshold <- c(1,0.3,0.1)
p.threshold <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001)
# Read in the covariates
read_LDpred=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA")
{
  # Read in the phenotype file 
  phenotype <- read.table(paste0(prefix,".pheno"), header=F)
  colnames(phenotype)=c("FID","IID","case")
  
  
  covariate <- read.table(paste0(prefix,".covariate"), header=T,stringsAsFactors = F)
  # Now merge the files
  pheno <- merge(phenotype, covariate, by=c("IID"))
  # We can then calculate the null model (model with PRS) using a linear regression 
  
  null.model <- glm(I(case==1)~., data=pheno[,!colnames(pheno)%in%c("FID","IID")],family = "binomial")
  null.model1 <- glm(I(case==1)~1, data=pheno[,!colnames(pheno)%in%c("FID","IID")],family = "binomial")
  # And the R2 of the null model is 
  null.r2 <- 1-logLik(null.model)/logLik(null.model1)
  prs.result <- NULL
  for(i in p.threshold){
    # Go through each p-value threshold .score_LDpred_p1.0000e-01.txt
    tmp=formatC(i, format = "e", digits =4)
    LDpredfile=paste0(prefix,".score_LDpred_p",tmp,".txt")
    if (file.exists(LDpredfile))
    {
      #prs <- read.table(paste0(prefix,".score_LDpred_p",tmp,".txt"), header=T,sep=",")
      prs <- read.table(paste0(prefix,".score_LDpred_p",tmp,".txt.adj"), header=T,sep=",")
      # Merge the prs with the phenotype matrix
      # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
      # relevant columns
      pheno.prs <- merge(pheno, prs, by=c("IID"))
      # Now perform a linear regression on Height with PRS and the covariates
      # ignoring the FID and IID from our model
      
      model <- glm(I(case==1)~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID","true_phens","cov_prs")],family = "binomial")
      model1 <- glm(I(case==1)~1, data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID","true_phens","cov_prs")],family = "binomial")
      # model R2 is obtained as 
      model.r2 <- 1-logLik(model)/logLik(model1)
      # R2 of PRS is simply calculated as the model R2 minus the null R2
      prs.r2 <- model.r2-null.r2
      #prs.r2=1-logLik(model)/logLik(null.model)
      # We can also obtain the coeffcient and p-value of association of PRS as follow
      prs.coef <- summary(model)$coeff["PRS",]
      #prs.coef <- summary(model)$coeff["cov_prs",]
      prs.beta <- as.numeric(prs.coef[1])
      prs.se <- as.numeric(prs.coef[2])
      prs.p <- as.numeric(prs.coef[4])
      # We can then store the results
      prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
    }
  }
  # Best result is:
  prs.result[which.max(prs.result$R2),]
  return(prs.result)
}

plot_LDpred=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA")
{
  prs.result=read_LDpred(prefix=prefix)
  #pdf(paste0(prefix,".LDpred.pdf"),width=12,height=8)
  png(paste0(prefix,".LDpred.png"), height=10, width=10, res=300, unit="in")
  # First, obtain the colorings based on the p-value
  col <- suppressWarnings(colorRampPalette(c("dodgerblue", "firebrick")))
  # We want the color gradient to match the ranking of p-values
  prs.result <- prs.result[order(-log10(prs.result$P)),]
  prs.result$color <-  col(nrow(prs.result))
  prs.result <- prs.result[order(prs.result$Threshold),]
  # generate a pretty format for p-value output
  prs.result$print.p <- round(prs.result$P, digits = 3)
  prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0 ] <-
    format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0 ], digits = 2)
  prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
  # Generate the axis labels
  xlab <- expression(italic(P) - value ~ threshold ~ (italic(P)[T]))
  ylab <- expression(paste("PRS model fit:  ", R ^ 2))
  # Setup the drawing area
  layout(t(1:2), widths=c(8.8,1.2))
  par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
       oma=c(0,0.5,0,0),
       mar=c(4,6,0.5,0.5))
  # Plotting the bars
  b<- barplot(height=prs.result$R2, 
              col=prs.result$color, 
              border=NA, 
              ylim=c(0, max(prs.result$R2)*1.25), 
              axes = F, ann=F)
  # Plot the axis labels and axis ticks
  odd <- seq(0,nrow(prs.result)+1,2)
  even <- seq(1,nrow(prs.result),2)
  axis(side=1, at=b[odd], labels=prs.result$Threshold[odd], lwd=2)
  axis(side=1, at=b[even], labels=prs.result$Threshold[even],lwd=2)
  axis(side=1, at=c(0,b[1],2*b[length(b)]-b[length(b)-1]), labels=c("","",""), lwd=2, lwd.tick=0)
  # Write the p-value on top of each bar
  text( parse(text=paste(
    prs.result$print.p)), 
    x = b+0.1, 
    y =  prs.result$R2+ (max(prs.result$R2)*1.05-max(prs.result$R2)), 
    srt = 45)
  # Now plot the axis lines
  box(bty='L', lwd=2)
  axis(2,las=2, lwd=2)
  # Plot the axis titles
  title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
  title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
  # Generate plot area for the legend
  par(cex.lab=1.5, cex.axis=1.25, font.lab=2, 
      mar=c(20,0,20,4))
  prs.result <- prs.result[order(-log10(prs.result$P)),]
  image(1, -log10(prs.result$P), t(seq_along(-log10(prs.result$P))), col=prs.result$color, axes=F,ann=F)
  axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
  # plot legend title
  title(bquote(atop(-log[10] ~ model, italic(P) - value), ), 
        line=2, cex=1.5, font=2, adj=0)
  dev.off()
  return(prs.result)
}

BE_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE")
BEEA_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA")
EA_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA")