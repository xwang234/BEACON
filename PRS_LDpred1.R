#!/usr/bin/env Rscript
# generate data for LDpred,for imputed beacon data 
#add beacon into the meta, only work on EA
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
keep_nohetsamples=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_QC")
{
  dat <- read.table(paste0(prefix,".het"), header=T) # Read in the EUR.het file, specify it has header
  m <- mean(dat$F) # Calculate the mean  
  s <- sd(dat$F) # Calculate the SD
  valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
  write.table(valid[,c(1,2)], paste0(prefix,".valid.sample"), quote=F, row.names=F) 
}


generate_plink=function(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/beacon_filter_noambiguous_QC",
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
generate_plink(opt="BE",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BE")
generate_plink(opt="EA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA1",metadat=EAmeta)
generate_plink(opt="BEEA",prefix1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA",metadat=BEEAmeta)



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
generate_covariate()
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA1")
generate_covariate(prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_BEEA")




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