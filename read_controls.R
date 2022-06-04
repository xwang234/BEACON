#this code is used to read new controls

# The uploaded data includes genotype files as described below (no imputed files):
#   
#   
#   Plink files for BEACON + dbgap controls  (cleaned version):
#   
#   BEACON_CRIC_GENEVAMel_PD.fam
# 
# BEACON_CRIC_GENEVAMel_PD.bim
# 
# BEACON_CRIC_GENEVAMel_PD.bed
# 
# 
# Plink files for Cambridge + WTCCC controls (cleaned version):
#   
#   cambridge_cases_WTCCC_controls.fam
# 
# cambridge_cases_WTCCC_controls.bed
# 
# cambridge_cases_WTCCC_controls.bim
# 
# 
# 
# Please note that the above Cambridge files include 5142 WTCCC controls,  but only 3408 were used with Cambridge cases in our lancet Oncol paper (the rest were used with Oxford data). To avoid sample overlap with the Oxford data, please use  phenotypes from the following file (exclude those coded -9 in the BE, EA, and BEEA fields in the file below):
#   
#   cambridge_cases_WTCCC_controls.sample
# 
# 
# To know which individuals are dbgap controls, they are listed in the following files for your reference:
#   
#   CRIC_Non-hispanic_white_cleaned2_sexchecked.fam
# 
# GENEVA_Melanoma_controls_caucasians_cleaned2.fam
# 
# PD_ENV_controls_caucasians_cleaned2_sexchecked.fam


sampletable=readxl::read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
# 11=England-Sheffield
# 12=Kaiser
# 13=Sweden - Karolinska
# 14=Mayo
# 15=EGA-WA
# 16=Ireland - FINBAR
# 17=Australia - Queensland
# 18=Toronto
# 19=UNC
# 20=WA reflux
# 21=Canada - Nova Scotia
# 22=EGA-NJ
# 23=USC Keck
# 25=Australia-wide
# 27=WA Reid
# 30=Cambridge
# 55=AMOS
table(sampletable$site)
# 11   12   13   14   15   16   17   18   19   20   21   22   23   25   27   30   55   NA 
# 276  517  181 1324  130  613  655  542  101  327  268   42  504  481  302 1939 1022   15

for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  sampletable[idx,i]=NA
}

BEACON_CRIC=read.table("../data/AdditionalControlPuya/BEACON_CRIC_GENEVAMel_PD.fam")
sum(BEACON_CRIC$V2 %in% sampletable$localid) #6091
tmp=sampletable$localid[sampletable$localid %in% BEACON_CRIC$V2]
idx=match(tmp,sampletable$localid)
table(sampletable$phenoEABE_bca[idx])
table(sampletable$phenoEA_bca[idx])
table(sampletable$phenoBE_bca[idx])
sum(!BEACON_CRIC$V2 %in% sampletable$localid) #4541 New
BEACON_CRIC_snp=read.table("../data/AdditionalControlPuya/BEACON_CRIC_GENEVAMel_PD.bim")
dim(BEACON_CRIC_snp) #[1] 719534      6
idx1=sampletable$localid %in% BEACON_CRIC$V2
table(sampletable$site[idx1])
# 11   12   13   14   15   16   17   18   19   20   21   22   23   25   27 
# 269  454  179 1312  129  607  646  502  100  323  260   40  495  481  294
BEACON_CRIC_new=BEACON_CRIC$V2[!BEACON_CRIC$V2 %in% sampletable$localid]
#some sample ID contain underscore, need to be replaced to run plink
# --update-ids expects input with the following four fields:
#   Old family ID
# Old within-family ID
# New family ID
# New within-family ID
sum(grepl("_",BEACON_CRIC$V1)) #0
sum(grepl("_",BEACON_CRIC$V2)) #1522
sum(grepl("SEP",BEACON_CRIC$V2))#0

BEACON_CRIC1=BEACON_CRIC
BEACON_CRIC1$V2=gsub("_","SEP",BEACON_CRIC$V2)
BEACON_CRIC_new1=BEACON_CRIC1$V2[!BEACON_CRIC$V2 %in% sampletable$localid]
write.table(BEACON_CRIC1,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.fam",col.names = F,row.names = F,sep="\t",quote=F)
file.copy("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_CRIC_GENEVAMel_PD.bed","/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.bed")
file.copy("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_CRIC_GENEVAMel_PD.bim","/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.bim")

Cambridge_WTCCC=read.table("../data/AdditionalControlPuya/cambridge_cases_WTCCC_controls.fam")
sum(Cambridge_WTCCC$V2 %in% sampletable$localid) #1885
sum(!Cambridge_WTCCC$V2 %in% sampletable$localid) #5142 WTCCC new
sum(Cambridge_WTCCC$V2 %in% BEACON_CRIC$V2) #0
Cambridge=Cambridge_WTCCC$V2[Cambridge_WTCCC$V2 %in% sampletable$localid]
Cambridge_WTCCC_snp=read.table("../data/AdditionalControlPuya/cambridge_cases_WTCCC_controls.bim")
dim(Cambridge_WTCCC_snp) #[1] 482715      6
idx2=sampletable$localid %in% Cambridge_WTCCC$V2
table(sampletable$site[idx2])
# 30 
# 1885 
Cambridge_WTCCC_pheno=read.table("../data/AdditionalControlPuya/cambridge_cases_WTCCC_controls_phenotypes")
WTCCC_new=Cambridge_WTCCC$V2[!Cambridge_WTCCC$V2 %in% sampletable$localid]
idx=match(WTCCC_new,Cambridge_WTCCC_pheno$V2)
table(Cambridge_WTCCC_pheno$V25[idx])
# -9    0 
# 1734 3408 
table(Cambridge_WTCCC_pheno$V26[idx])
table(Cambridge_WTCCC_pheno$V27[idx])
sum(grepl("SEP",Cambridge_WTCCC$V2)) #0
sum(grepl("_",Cambridge_WTCCC$V2)) #5142
sum(grepl("_",Cambridge_WTCCC$V1)) #0
Cambridge_WTCCC1=Cambridge_WTCCC
Cambridge_WTCCC1$V2=gsub("_","SEP",Cambridge_WTCCC1$V2)
write.table(Cambridge_WTCCC1,file="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/Cambridge_WTCCC_newID.fam",col.names = F,row.names = F,sep="\t",quote=F)
file.copy("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/cambridge_cases_WTCCC_controls.bed","/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/Cambridge_WTCCC_newID.bed")
file.copy("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/cambridge_cases_WTCCC_controls.bim","/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/Cambridge_WTCCC_newID.bim")


dbgap_control1=read.table("../data/AdditionalControlPuya/CRIC_Non-hispanic_white_cleaned2_sexchecked.fam")
dim(dbgap_control1) # 1544    6
sum(dbgap_control1$V2 %in% BEACON_CRIC$V2) #1522
dbgap_control2=read.table("../data/AdditionalControlPuya/GENEVA_Melanoma_controls_caucasians_cleaned2.fam")
dim(dbgap_control2) # 1047    6
sum(dbgap_control2$V2 %in% BEACON_CRIC$V2) #1038

dbgap_control3=read.table("../data/AdditionalControlPuya/PD_ENV_controls_caucasians_cleaned2_sexchecked.fam")
dim(dbgap_control3) # 1986    6
sum(dbgap_control3$V2 %in% BEACON_CRIC$V2) #1981

dbgap_control=c(dbgap_control1$V2,dbgap_control2$V2,dbgap_control3$V2)
length(dbgap_control) #4577
sum(dbgap_control %in% BEACON_CRIC_new) #4541
all(BEACON_CRIC_new %in% dbgap_control) #T
all(dbgap_control %in% BEACON_CRIC$V2) #F

#From Puya lancet supplement
# To increase the statistical power for detecting risk loci, an
# additional 4,541 unscreened controls of European ancestry genotyped on the Illumina HumanOmni1-Quad array
# were obtained from dbGaP and merged with the BEACON data.
# dbgap_newcontrol=BEACON_CRIC_new1
# idx=match(dbgap_newcontrol,BEACON_CRIC1$V2)
# tmp=BEACON_CRIC1[idx,]
# write.table(BEACON_CRIC1[idx,1:2],file="../result/BEACON_dbgapcontrol_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote=F)
# 
# tmp1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/BEACON_dbgapcontrol.fam")
# idx=tmp$V2 %in% tmp1$V2
# tmp2=tmp$V2[!tmp$V2 %in% tmp1$V2]
# idx=match(tmp2,tmp$V2)
# # Following data cleaning and
# # merging, 873 BE cases, 995 EA cases (1,868 BE/EA combined), and 3,408 unscreened controls, all of European
# # descent, were available for this study (see Appendix table 1)
WTCCC_newcontrol=Cambridge_WTCCC_pheno$V2[which(Cambridge_WTCCC_pheno$V25==0)]
Cambridge_WTCCCnew=c(Cambridge,WTCCC_newcontrol)
idx=match(Cambridge_WTCCCnew,Cambridge_WTCCC$V2)
write.table(Cambridge_WTCCC1[idx,1:2],file="../result/Cambridge_WTCCC_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote=F)

#check number of samples in each chr data (generated by convert_ped_2vcf.sh),quality control
number_samples=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf/",
                        prefix="BEACON_dbgapcontrol_chr")
{
  for (i in 1:22)
  {
    fam=read.table(paste0(outfolder,prefix,i,".fam"))
    if (i==1)
    {
      sample=fam$V2
    }else
    {
      sample=intersect(sample,fam$V2)
    }
    
    bim=read.table(paste0(outfolder,prefix,i,".bim"))
    print(paste0(i,":",nrow(fam),":",nrow(bim)))
    
  }
  return(sample)
}
qc_beaconsamples=number_samples()
idx=match(qc_beaconsamples,BEACON_CRIC1$V2)
write.table(BEACON_CRIC1[idx,1:2],file="../result/BEACON_dbgapcontrol_QCplinksamples.txt",col.names = F,row.names = F,sep="\t",quote=F)

qc_wtcccsamples=number_samples(prefix="Cambridge_WTCCCcontrol_chr")
idx=match(qc_wtcccsamples,Cambridge_WTCCC1$V2)
write.table(Cambridge_WTCCC1[idx,1:2],file="../result/Cambridge_WTCCC_QCplinksamples.txt",col.names = F,row.names = F,sep="\t",quote=F)

outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/imputation_vcf"
prefix="BEACON_dbgapcontrol_newID"
tmp=paste0(outfolder,"/",prefix,"_chr",chr,"_missing.imiss")
tmp1=read.table(tmp,header=T)
