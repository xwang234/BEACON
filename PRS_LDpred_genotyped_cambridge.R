#!/usr/bin/env Rscript
# generate data for LDpred,for genotyped validation data 
#use Cambridge+BEACONcontrol as the validation, only work on EA
#BEACONcase+AMOS+other BEACONcontrol for discovery

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
write.table(discsamples,file="../result/discsamples_plink.txt",row.names = F,col.names = F,sep="\t",quote=F)

EAsamples=sampletable1$localid[which(sampletable1$phenoEA_bc==2)]
COsamples=sampletable1$localid[which(sampletable1$phenoEA_bc==1)]
fam=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous.fam")
COsamples=COsamples[COsamples %in% fam$V2]
filesforGWAS=function(casesamples=EAsamples,controlsamples=COsamples,prefix="Discovery_genotyped_EA_CO")
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
tmp=read.table("../result/GWAS/Discovery_genotyped_EA_CO_gwas.fam")
tmp1=read.table("../result/GWAS/Discovery_genotyped_EA_CO_selectedsamples_pheno_plink.txt",header = T)

all(tmp$V2 %in% tmp1$IID) #T
udpate_fam=function(famfile="../result/GWAS/Discovery_genotyped_EA_CO_gwas.fam",phenofile="../result/GWAS/Discovery_genotyped_EA_CO_selectedsamples_pheno_plink.txt")
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

update_gwas=function(gwasfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/Discovery_genotyped_EA_CO_gwas.assoc.logistic",
                     bimfile="../result/GWAS/Discovery_genotyped_EA_CO_gwas.bim")
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
  idx=is.na(gwas$P)
  gwas=gwas[!idx,]
  fwrite(gwas,file=gwasfile,quote=F,sep=" ")
}

#run Meta_Analysis.R and Meta_Analysis.sh to get meta results


library(data.table)

EAmeta=as.data.frame(fread("../result/Discovery_Bonn_METAANALYSIS_EA_genotyped_comsnp1.tbl"))
colnames(EAmeta)[which(colnames(EAmeta)=="P-value")]="P"
EABeacon=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_genotyped_BD_autosomes_comsnp_N.txt"))
idx=match(EAmeta$MarkerName,EABeacon$SNP)
EAmeta$chr=EABeacon$CHR[idx]
EAmeta$pos=EABeacon$position[idx]
write.table(EAmeta,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_Bonn_genotyped_metastat.txt",col.names = T,
            row.names=F,quote=F,sep="\t")
cmd="gzip /fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Discovery_Bonn_genotyped_metastat.txt"
system(cmd)

#work on validation------

#Pick validation samples
library(readxl)
#step 2---
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
idx=which(sampletable$localid %in% validationsamples)
sampletable1=sampletable[idx,]
valisamples=sampletable[idx,1:2]
table(sampletable$phenoBE_bca[idx])
# -9    1 
# 1005 1139 
table(sampletable$phenoEA_bca[idx])
# 1    2 
# 1141 1003 

write.table(valisamples,file="../result/validationsamples_plink.txt",row.names = F,col.names = F,sep="\t",quote=F)


keep_nohetsamples=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/validation_filtered_30Nov2018_flip_noambiguous_QC")
{
  dat <- read.table(paste0(prefix,".het"), header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  write.table(valid[,c(1,2)], paste0(prefix,".valid.sample"), quote=F, row.names=F) 
}

#step 4----
#generate validation BE/EA/BEEA plink files
generate_plink=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/validation_filtered_30Nov2018_flip_noambiguous_QC",
                        opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_BE_genotyped",metadat=BEmeta)
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
#generate_plink(opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_genotyped")
generate_plink(opt="EA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA_genotyped",metadat=EAmeta)
#generate_plink(opt="BEEA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped",metadat=BEEAmeta)


generate_covariate=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_BE_genotyped")
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
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA_genotyped")
#generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped")




#after runing PRS_LDpred.sh, read results:
#library(DescTools)
p.threshold <- c(1,0.3,0.1)
p.threshold <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001)
# Read in the covariates
# read_LDpred=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped")
# {
#   # Read in the phenotype file 
#   phenotype <- read.table(paste0(prefix,".pheno"), header=F)
#   colnames(phenotype)=c("FID","IID","case")
#   
#   
#   covariate <- read.table(paste0(prefix,".covariate"), header=T,stringsAsFactors = F)
#   # Now merge the files
#   pheno <- merge(phenotype, covariate, by=c("IID"))
#   # We can then calculate the null model (model with PRS) using a linear regression 
#   
#   null.model <- glm(case~., data=pheno[,!colnames(pheno)%in%c("FID","IID")],family = "binomial")
#   null.model1 <- glm(case~1, data=pheno[,!colnames(pheno)%in%c("FID","IID")],family = "binomial")
#   # And the R2 of the null model is 
#   null.r2 <- 1-logLik(null.model)/logLik(null.model1)
#   prs.result <- NULL
#   for(i in p.threshold){
#     # Go through each p-value threshold .score_LDpred_p1.0000e-01.txt
#     tmp=formatC(i, format = "e", digits =4)
#     LDpredfile=paste0(prefix,".score_LDpred_p",tmp,".txt")
#     if (file.exists(LDpredfile))
#     {
#       #prs <- read.table(paste0(prefix,".score_LDpred_p",tmp,".txt"), header=T,sep=",")
#       prs <- read.table(paste0(prefix,".score_LDpred_p",tmp,".txt.adj"), header=T,sep=",")
#       # Merge the prs with the phenotype matrix
#       # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
#       # relevant columns
#       pheno.prs <- merge(pheno, prs, by=c("IID"))
#       # Now perform a linear regression on Height with PRS and the covariates
#       # ignoring the FID and IID from our model
#       
#       model <- glm(case~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID","true_phens","cov_prs")],family = "binomial")
#       model1 <- glm(case~1, data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID","true_phens","cov_prs")],family = "binomial")
#       # model R2 is obtained as 
#       model.r2 <- 1-logLik(model)/logLik(model1)
#       # R2 of PRS is simply calculated as the model R2 minus the null R2
#       prs.r2 <- model.r2-null.r2
#       #prs.r2=1-logLik(model)/logLik(null.model)
#       # We can also obtain the coeffcient and p-value of association of PRS as follow
#       prs.coef <- summary(model)$coeff["PRS",]
#       #prs.coef <- summary(model)$coeff["cov_prs",]
#       prs.beta <- as.numeric(prs.coef[1])
#       prs.se <- as.numeric(prs.coef[2])
#       prs.p <- as.numeric(prs.coef[4])
#       # We can then store the results
#       prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
#     }
#   }
#   # Best result is:
#   prs.result[which.max(prs.result$R2),]
#   return(prs.result)
# }
  
# plot_LDpred=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped")
# {
#   prs.result=read_LDpred(prefix=prefix)
#   #pdf(paste0(prefix,".LDpred.pdf"),width=12,height=8)
#   png(paste0(prefix,".LDpred.png"), height=10, width=10, res=300, unit="in")
#   # First, obtain the colorings based on the p-value
#   col <- suppressWarnings(colorRampPalette(c("dodgerblue", "firebrick")))
#   # We want the color gradient to match the ranking of p-values
#   prs.result <- prs.result[order(-log10(prs.result$P)),]
#   prs.result$color <-  col(nrow(prs.result))
#   prs.result <- prs.result[order(prs.result$Threshold),]
#   # generate a pretty format for p-value output
#   prs.result$print.p <- round(prs.result$P, digits = 3)
#   prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0 ] <-
#     format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0 ], digits = 2)
#   prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
#   # Generate the axis labels
#   xlab <- expression(italic(P) - value ~ threshold ~ (italic(P)[T]))
#   ylab <- expression(paste("PRS model fit:  ", R ^ 2))
#   # Setup the drawing area
#   layout(t(1:2), widths=c(8.8,1.2))
#   par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
#        oma=c(0,0.5,0,0),
#        mar=c(4,6,0.5,0.5))
#   # Plotting the bars
#   b<- barplot(height=prs.result$R2, 
#               col=prs.result$color, 
#               border=NA, 
#               ylim=c(0, max(prs.result$R2)*1.25), 
#               axes = F, ann=F)
#   # Plot the axis labels and axis ticks
#   odd <- seq(0,nrow(prs.result)+1,2)
#   even <- seq(1,nrow(prs.result),2)
#   axis(side=1, at=b[odd], labels=prs.result$Threshold[odd], lwd=2)
#   axis(side=1, at=b[even], labels=prs.result$Threshold[even],lwd=2)
#   axis(side=1, at=c(0,b[1],2*b[length(b)]-b[length(b)-1]), labels=c("","",""), lwd=2, lwd.tick=0)
#   # Write the p-value on top of each bar
#   text( parse(text=paste(
#     prs.result$print.p)), 
#     x = b+0.1, 
#     y =  prs.result$R2+ (max(prs.result$R2)*1.05-max(prs.result$R2)), 
#     srt = 45)
#   # Now plot the axis lines
#   box(bty='L', lwd=2)
#   axis(2,las=2, lwd=2)
#   # Plot the axis titles
#   title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
#   title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
#   # Generate plot area for the legend
#   par(cex.lab=1.5, cex.axis=1.25, font.lab=2, 
#       mar=c(20,0,20,4))
#   prs.result <- prs.result[order(-log10(prs.result$P)),]
#   image(1, -log10(prs.result$P), t(seq_along(-log10(prs.result$P))), col=prs.result$color, axes=F,ann=F)
#   axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
#   # plot legend title
#   title(bquote(atop(-log[10] ~ model, italic(P) - value), ), 
#         line=2, cex=1.5, font=2, adj=0)
#   dev.off()
#   return(prs.result)
# }
# 
# BE_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_genotyped")
# BEEA_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_genotyped")
# EA_LDpred=plot_LDpred(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped")
# save(BE_LDpred,EA_LDpred,BEEA_LDpred,file="../result/LDpred_geneotyped_result.RData")

library(pROC)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
plot_ROC=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE_genotyped",opt=1,opt1="Beacon")
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
        idx=match(pheno.prs$IID,sampletable$localid)
        sampletable1=sampletable[idx,]
        sampletable1$prs=pheno.prs$PRS
        sampletable1$case=pheno.prs$case
        if (opt1=="Beacon") #validation dataset
        {
          fit1 <- glm(I(case==2)~prs,family=binomial,data=sampletable1,y=T)
          roc1=roc(fit1$y,fit1$fitted.values,quiet = T)
          #remove bmi_recent_healthy
          dat=sampletable1[,c("case","age","sex","recurrent_HB_RF","cig_smk_ever","nsaid_ever","prs","site")]
          idx=complete.cases(dat)
          if (sum(idx)>0)
          {
            fit3 <- glm(I(case==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+prs,family=binomial,data=sampletable1,y=T)
            roc3=roc(fit3$y,fit3$fitted.values,quiet = T)
            print(paste0("p=",i,":auc_prs=",round(roc1$auc,3)," auc_env=",round(roc3$auc,3)))
          }else
          {
            print(paste0("p=",i,":auc_prs=",round(roc1$auc,3)))
          }
          
        }else
        {
          fit1 <- glm(I(case==2)~prs,family=binomial,data=sampletable1[sampletable1$site>=30,],y=T)
          #fit3 <- glm(I(case==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+prs,family=binomial,data=sampletable1[sampletable1$site>=30,],y=T)
          roc1=roc(fit1$y,fit1$fitted.values,quiet = T)
          #roc3=roc(fit3$y,fit3$fitted.values,quiet = T)
          print(paste0("p=",i,":auc_prs=",round(roc1$auc,3)))
        }
        
        #plot(roc1, print.auc=TRUE,main=tmp)
      }
    }
  }
  
  if (opt==2) #LDpred-inf
  {
    LDpredfile=paste0(prefix,".score_LDpred-inf.txt")
    prs <- read.table(LDpredfile, header=T,sep=",",stringsAsFactors = F)
    pheno.prs <- merge(phenotype, prs, by=c("IID"))
    testroc<- roc(pheno.prs$true_phens, pheno.prs$PRS)
    plot(testroc, print.auc=TRUE,main="inf")
  }
  
}
#par(mfrow=c(2,2))

plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped")
# [1] "p=1:auc_prs=0.562 auc_env=0.755"
# [1] "p=0.3:auc_prs=0.563 auc_env=0.755"
# [1] "p=0.1:auc_prs=0.563 auc_env=0.755"
# [1] "p=0.03:auc_prs=0.566 auc_env=0.755"
# [1] "p=0.01:auc_prs=0.565 auc_env=0.754"
# [1] "p=0.003:auc_prs=0.525 auc_env=0.744"
# [1] "p=0.001:auc_prs=0.509 auc_env=0.743"
# [1] "p=3e-04:auc_prs=0.508 auc_env=0.743"
# [1] "p=1e-04:auc_prs=0.529 auc_env=0.744"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped_cn")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped_info")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped_cn_info")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA")
# [1] "p=1:auc_prs=0.565 auc_env=0.753"
# [1] "p=0.3:auc_prs=0.565 auc_env=0.753"
# [1] "p=0.1:auc_prs=0.565 auc_env=0.753"
# [1] "p=0.03:auc_prs=0.566 auc_env=0.753"
# [1] "p=0.01:auc_prs=0.568 auc_env=0.754"
# [1] "p=0.003:auc_prs=0.527 auc_env=0.744"
# [1] "p=0.001:auc_prs=0.552 auc_env=0.75"
# [1] "p=3e-04:auc_prs=0.542 auc_env=0.747"
# [1] "p=1e-04:auc_prs=0.538 auc_env=0.746"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_cn") #add case_N,control_N
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_info")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_cn_info")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_noAmos",opt1="Cambridge")
# [1] "p=1:auc_prs=0.597"
# [1] "p=0.3:auc_prs=0.598"
# [1] "p=0.1:auc_prs=0.6"
# [1] "p=0.03:auc_prs=0.601"
# [1] "p=0.01:auc_prs=0.591"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_noAmos",opt1="Cambridge")
# [1] "p=1:auc_prs=0.601"
# [1] "p=0.3:auc_prs=0.601"
# [1] "p=0.1:auc_prs=0.602"
# [1] "p=0.03:auc_prs=0.603"
# [1] "p=0.01:auc_prs=0.606"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA_genotyped")
# [1] "p=1:auc_prs=0.603 auc_env=0.813"
# [1] "p=0.3:auc_prs=0.603 auc_env=0.813"
# [1] "p=0.1:auc_prs=0.604 auc_env=0.813"
# [1] "p=0.03:auc_prs=0.603 auc_env=0.813"
# [1] "p=0.01:auc_prs=0.596 auc_env=0.812"
# [1] "p=0.003:auc_prs=0.576 auc_env=0.808"
# [1] "p=0.001:auc_prs=0.562 auc_env=0.805"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA")
# [1] "p=1:auc_prs=0.602 auc_env=0.813"
# [1] "p=0.3:auc_prs=0.602 auc_env=0.813"
# [1] "p=0.1:auc_prs=0.603 auc_env=0.813"
# [1] "p=0.03:auc_prs=0.604 auc_env=0.814"
# [1] "p=0.01:auc_prs=0.607 auc_env=0.815"
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_N5000",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_N10000",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_N15000",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_info",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_cn_info",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_noAmos",opt1="Cambridge")


plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_info",opt1="Cambridge")
plot_ROC(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_cn_info",opt1="Cambridge")

#calculate ordds ratios from PRS
require(MASS)
#PRS high, risk high, for genotypted
PRS_odds=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE",p=0.1)
{
  phenotype <- read.table(paste0(prefix,".pheno"), header=F)
  colnames(phenotype)=c("FID","IID","case")
  
  tmp=formatC(p, format = "e", digits =4)
  LDpredfile=paste0(prefix,".score_LDpred_p",tmp,".txt.adj")
  if (file.exists(LDpredfile))
  {
    prs <- read.table(LDpredfile, header=T,sep=",")
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
    #return(tmp)
  }
}

#EA---

PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA",p=0.01)
# [1] "2th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.306476 1.082987 1.576619 
# [1] "3th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.435085 1.189547 1.732120 
# [1] "4th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.999902 1.658746 2.413559 
PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_genotyped_noAmos",p=0.03)
# [1] "2th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.451905 1.134177 1.860319 
# [1] "3th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.988629 1.550527 2.554916 
# [1] "4th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   2.498519 1.944306 3.218334 
PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/CambAmos_EA_noAmos",p=0.03)

PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA_genotyped",p=0.1)

PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA",p=0.01)
# [1] "2th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   1.344442 1.049057 1.724538 
# [1] "3th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   2.152021 1.682578 2.758004 
# [1] "4th quantile:"
# Waiting for profiling to be done...
# 2.5 %   97.5 % 
#   2.653871 2.072655 3.406600 
PRS_odds(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA_eQTL",p=0.1)
BEprs=read.table("../result/Beacon_BE_genotyped.PRS.txt",header = T,stringsAsFactors = F)
EAprs=read.table("../result/Beacon_EA_genotyped.PRS.txt",header = T,stringsAsFactors = F)
BEEAprs=read.table("../result/Beacon_BEEA_genotyped.PRS.txt",header = T,stringsAsFactors = F)
save(BEprs,EAprs,BEEAprs,file="../result/LDpred_genotyped_PRS.RData")

BEprs=read.table("../result/Beacon_BE.PRS.txt",header = T,stringsAsFactors = F)
EAprs=read.table("../result/Beacon_EA.PRS.txt",header = T,stringsAsFactors = F)
BEEAprs=read.table("../result/Beacon_BEEA.PRS.txt",header = T,stringsAsFactors = F)
#save(BEprs,EAprs,BEEAprs,file="../result/LDpred_PRS.RData")

BEprs=read.table("../result/Beacon_BE_eQTL.PRS.txt",header = T,stringsAsFactors = F)
EAprs=read.table("../result/Beacon_EA_eQTL.PRS.txt",header = T,stringsAsFactors = F)
BEEAprs=read.table("../result/Beacon_BEEA_eQTL.PRS.txt",header = T,stringsAsFactors = F)
save(BEprs,EAprs,BEEAprs,file="../result/LDpred_eQTL_PRS.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_genotyped_PRS.RData")
plotroc=function(EAfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA.PRS.txt")
{
  library(readxl)
  sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
  sampletable <- data.frame(sampletable)
  for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
  
  EAprs=read.table(EAfile,header=T,stringsAsFactors = F)
  names(EAprs) <- c("localid","EA.prs")
  
  sampletable1 <- merge(sampletable,EAprs,by="localid")
  
  #png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_ROC.png")
  fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable1,y=T,x=T)
  #summary(fit1)
  fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable1,y=T)
  roc2 <- roc(fit2$y~fit2$fitted.values)
  print(roc2)
  #summary(fit2)
  library(pROC)
  roc1 <- roc(fit1$y~fit1$fitted.values)
  #remove "bmi_recent_healthy"
  dat=sampletable1[,c("phenoEA_bca","age","sex","recurrent_HB_RF","bmi_recent_healthy","cig_smk_ever","nsaid_ever","EA.prs")]
  idx=complete.cases(dat)
  if (sum(idx)>0)
  {
    fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable1,y=T)
    #summary(fit3)
    
    plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.2,cex.lab=1.2)
    roc2 <- roc(fit3$y~fit3$fitted.values)
    plot(roc2,ylim=c(0,1),print.auc=F,col=4,add=T)
    text(0.4,0.5,paste0("PRS: AUC=",round(roc1$auc,2)),col=3,cex=1.2)
    text(0.7,0.98,paste0("PRS+Environment: AUC=",round(roc2$auc,2)),col=4,cex=1.2)
  }else
  {
    plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.2,cex.lab=1.2)
    text(0.4,0.5,paste0("PRS: AUC=",round(roc1$auc,3)),col=3,cex=1.2)
  }
  
  #dev.off()
}
plotroc()
