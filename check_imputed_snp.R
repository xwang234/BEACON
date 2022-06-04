
#
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdata_ambiguous_TPM_for_prediction.RData")
snppos_stomach=snppos
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8muscularisdata_ambiguous_TPM_for_prediction.RData")
snppos_muscularis=snppos

#readsnp
readsnp=function(dataset="cambridgewtccc_1000g")
{
  snp=snppos=NULL
  for (chr in 1:22)
  {
    cat(chr,'..')
    
    bim1file=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter_flip_chr",chr,".bim")
    bim1=fread(bim1file)
    bim1=as.data.frame(bim1)
    bim2file=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/",dataset,"/chr",chr,"_filter_hg19tohg38_flip.bim")
    bim2=fread(bim2file)
    bim2=as.data.frame(bim2)
    str11=paste0(bim1$V4,"_",bim1$V5,"_",bim1$V6)
    str12=paste0(bim1$V4,"_",bim1$V6,"_",bim1$V5)
    str2=paste0(bim2$V4,"_",bim2$V5,"_",bim2$V6)
    idx=str11 %in% str2 | str12 %in% str2
    traw1file=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_filter_flip_chr",chr,".traw.gz")
    traw1=fread(input=traw1file)
    traw1=as.data.frame(traw1)
    tmp=paste0(traw1$CHR,":",traw1$POS,"_",traw1$COUNTED,"_",traw1$ALT)
    rownames(traw1)=tmp
    traw1=traw1[,7:ncol(traw1)]
    tmp=unlist(strsplit(colnames(traw1),"_"))
    colnames(traw1)=tmp[seq(1,length(tmp),2)]
    idx1=match(colnames(phenotype),colnames(traw1))
    if (sum(is.na(idx1))>0) warning(paste0("some RNAseq samples are missing in genotype data,",chr))
    traw1=traw1[,idx1]
    snp_=traw1[idx,]
    snppos_=data.frame(chr=bim1$V1[idx],pos=bim1$V4[idx],minor=bim1$V5[idx],major=bim1$V6[idx],stringsAsFactors = F)
    rownames(snppos_)=paste0(chr,":",snppos_$pos,"_",snppos_$minor,"_",snppos_$major)
    snp=rbind.data.frame(snp,snp_)
    snppos=rbind.data.frame(snppos,snppos_)
  }
  
  k <- which(is.na(snp), arr.ind=TRUE)
  length(k)/nrow(snp)/ncol(snp)
  snp[k] <- rowMeans(snp, na.rm=TRUE)[k[,1]]
  all(rownames(snppos)==rownames(snp))
  if (any(colnames(snp)!=colnames(phenotype))) warning("samples name not match!")
  all(colnames(snp)==rownames(covariate))
  return(list(snp=snp,snppos=snppos))
}
cambridge1000g=readsnp()
cambridge1000g_qc=readsnp(dataset="cambridgewtccc_qc_1000g")
beacon1000g=readsnp(dataset="beacondbgapcontrol_1000g")
beacon1000g_qc=readsnp(dataset="beacondbgapcontrol_qc_1000g")
merge1000g=readsnp(dataset="merge_beacon_cambridge_1000g")
merge1000g_qc=readsnp(dataset="merge_beacon_cambridge_qc_1000g")

overlapsnps=function(snppos2=cambridge1000g$snppos,snppos1=snppos_muscularis)
{
  str11=paste0(snppos1$chr,":",snppos1$pos,"_",snppos1$minor,"_",snppos1$major)
  str12=paste0(snppos1$chr,":",snppos1$pos,"_",snppos1$major,"_",snppos1$minor)
  str2=paste0(snppos2$chr,":",snppos2$pos,"_",snppos2$minor,"_",snppos2$major)
  idx=str11 %in% str2 | str12 %in% str2
  print(paste0("snp1 has ",length(str11)," snps"))
  print(paste0("snp2 has ",length(str2)," snps"))
  print(paste0(sum(idx)," snps are overlap(",round(sum(idx)/length(idx)*100),"%)"))
}

overlapsnps()
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 5206460 snps"
# [1] "5115603 snps are overlap(98%)"
overlapsnps(snppos2=cambridge1000g_qc$snppos)
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 5250578 snps"
# [1] "5159493 snps are overlap(99%)"
overlapsnps(snppos2=beacon1000g$snppos)
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 4765617 snps"
# [1] "4678457 snps are overlap(90%)"
overlapsnps(snppos2=beacon1000g_qc$snppos)
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 5293138 snps"
# [1] "5198990 snps are overlap(100%)"
overlapsnps(snppos1=cambridge1000g$snppos,snppos2 = beacon1000g$snppos)
# [1] "snp1 has 5206460 snps"
# [1] "snp2 has 4765617 snps"
# [1] "4665331 snps are overlap(90%)"
overlapsnps(snppos1=cambridge1000g_qc$snppos,snppos2 = beacon1000g_qc$snppos)
# [1] "snp1 has 5250578 snps"
# [1] "snp2 has 5293138 snps"
# [1] "5231570 snps are overlap(100%)"

overlapsnps(snppos2=merge1000g$snppos)
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 4634237 snps"
# [1] "4560873 snps are overlap(87%)"
overlapsnps(snppos2=merge1000g_qc$snppos)
# [1] "snp1 has 5218673 snps"
# [1] "snp2 has 5196891 snps"
# [1] "5117484 snps are overlap(98%)"

overlapmodelsnps=function(prefix="dist500K_GTEx_muscularis_June11",snppos=cambridge1000g$snppos,opt=1)
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  str1=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed = T)))
  str1=str1[!is.na(str1)]
  str21=paste0(snppos$chr,":",snppos$pos,"_",snppos$minor,"_",snppos$major)
  str22=paste0(snppos$chr,":",snppos$pos,"_",snppos$major,"_",snppos$minor)
  idx=str1 %in% str21 | str1 %in% str22
  print(paste0("model has ",length(str1)," snps"))
  if (opt==1) print(paste0("snp2 has ",length(str21)," snps"))
  print(paste0(sum(idx)," snps are overlap (",round(sum(idx)/length(idx)*100),"%)"))
}

overlapallmodelsnps=function(snppos=cambridge1000g$snppos)
{
  tmp=c("dist500K_GTEx_June11","dist500K_GTEx_stomach_June11","dist500K_GTEx_muscularis_June11","dist500K_GTEx_mucosa_June11","dist500K_GTEx_blood_June11","dist500K_GTEx_adipose_June11")
  tmp1=c("Junction","Stomach","Muscularis","Mucosa","Blood","Adipose")
  for (i in 1:length(tmp))
  {
    print(paste0(tmp1[i],"------"))
    if (i==1)
    {
      overlapmodelsnps(prefix=tmp[i],snppos=snppos)
    }else
    {
      overlapmodelsnps(prefix=tmp[i],snppos=snppos,opt=2)
    }
    
  }
}
overlapallmodelsnps()
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 5206460 snps"
# [1] "466267 snps are overlap (97%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "465899 snps are overlap (97%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "488260 snps are overlap (97%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "496257 snps are overlap (97%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "459000 snps are overlap (98%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "500091 snps are overlap (97%)"

overlapallmodelsnps(snppos=cambridge1000g_qc$snppos)
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 5250578 snps"
# [1] "469042 snps are overlap (98%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "468794 snps are overlap (98%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "491229 snps are overlap (98%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "499035 snps are overlap (98%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "461229 snps are overlap (98%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "503166 snps are overlap (98%)"
overlapallmodelsnps(snppos=beacon1000g$snppos)
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 4765617 snps"
# [1] "434507 snps are overlap (91%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "433890 snps are overlap (91%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "455925 snps are overlap (91%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "462349 snps are overlap (91%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "426752 snps are overlap (91%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "465817 snps are overlap (91%)"
overlapallmodelsnps(snppos=beacon1000g_qc$snppos)
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 5293138 snps"
# [1] "475139 snps are overlap (99%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "474881 snps are overlap (99%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "497655 snps are overlap (99%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "505444 snps are overlap (99%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "467135 snps are overlap (99%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "509695 snps are overlap (99%)"
overlapallmodelsnps(snppos=merge1000g$snppos)
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 4634237 snps"
# [1] "422591 snps are overlap (88%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "421827 snps are overlap (88%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "443328 snps are overlap (88%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "449891 snps are overlap (88%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "415609 snps are overlap (88%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "452961 snps are overlap (88%)"
overlapallmodelsnps(snppos=merge1000g_qc$snppos)
# [1] "Junction------"
# [1] "model has 478676 snps"
# [1] "snp2 has 5196891 snps"
# [1] "465167 snps are overlap (97%)"
# [1] "Stomach------"
# [1] "model has 478391 snps"
# [1] "464919 snps are overlap (97%)"
# [1] "Muscularis------"
# [1] "model has 501373 snps"
# [1] "487205 snps are overlap (97%)"
# [1] "Mucosa------"
# [1] "model has 509246 snps"
# [1] "494902 snps are overlap (97%)"
# [1] "Blood------"
# [1] "model has 470515 snps"
# [1] "457446 snps are overlap (97%)"
# [1] "Adipose------"
# [1] "model has 513441 snps"
# [1] "499059 snps are overlap (97%)"

#check samples
BEACON_CRIC1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/AdditionalControlPuya/BEACON_dbgapcontrol_newID.fam")
BEACONsamples=gsub("SEP","_",BEACON_CRIC1$V2)
imputedBEACONchr17=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_1000g/chr17_filter_hg19tohg38_flip.fam")
chr17imputedBEACONsamples=gsub("SEP","",imputedBEACONchr17$V2)
imputedBEACON=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_1000g/chr1_filter_hg19tohg38_flip.fam")
imputedBEACONsamples=gsub("SEP","",imputedBEACON$V2)
imputedqcBEACON=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/chr1_filter_hg19tohg38_flip.fam")
imputedqcBEACONsamples=gsub("SEP","",imputedqcBEACON$V2)

BEACONqc_CRIC1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/chr17_filter_hg19tohg38_flip.fam")
Cambridge_WTCCC1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Cambridge_WTCCC_plinksamples.txt")
Cambridgesamples=gsub("SEP","",Cambridge_WTCCC1$V2)
imputedCambridge=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_1000g/chr17_filter_hg19tohg38_flip.fam")
imputedCambridgesamples=gsub("SEP","",imputedCambridge$V2)
imputedqcCambridge=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/chr17_filter_hg19tohg38_flip.fam")
imputedqcCambridgesamples=gsub("SEP","",imputedqcCambridge$V2)

nrow(BEACON_CRIC1)+nrow(Cambridge_WTCCC1) #15925
mergesample1000g=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/","merge_beacon_cambridge_1000g","/chr",1,"_filter_hg19tohg38_flip.fam"))
mergesample=gsub("SEP","",mergesample1000g$V2)
mergesample1000g1=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/","merge_beacon_cambridge_1000g","/chr",17,"_filter_hg19tohg38_flip.fam"))
chr17mergesample=gsub("SEP","",mergesample1000g1$V2)
mergesampleqc1000g=read.table(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/","merge_beacon_cambridge_qc_1000g","/chr",1,"_filter_hg19tohg38_flip.fam"))
mergeqcsample=gsub("SEP","",mergesampleqc1000g$V2)

#Beacon 2,413 BE cases, 1,512 EA cases and 2,185 controls +4541 controls
#Cambridge  873 BE cases, 995 EA cases (1,868 BE/EA combined), and 3,408 unscreened controls

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
sampletable_beacon=sampletable$localid[!sampletable$site %in% c(30,55)]
sampletable_cambridge=sampletable$localid[sampletable$site %in% c(30)]
sampletable_amos=sampletable$localid[sampletable$site %in% c(55)]
idx=match(sampletable_beacon,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #2413 BE
table(sampletable$phenoEA_bca[idx]) #1512 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2185 controls

idx=match(sampletable_cambridge,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #882 BE
table(sampletable$phenoEA_bca[idx]) #1003 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #0 controls

BEACONsamples1=BEACONsamples[BEACONsamples %in% sampletable$localid]
BEACONsamples2=BEACONsamples[!BEACONsamples %in% sampletable$localid] #4541 new controls
idx=match(BEACONsamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #2406 BE
table(sampletable$phenoEA_bca[idx]) #1508 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2177 controls

chr17imputedBEACONsamples1=chr17imputedBEACONsamples[chr17imputedBEACONsamples %in% sampletable$localid]
chr17imputedBEACONsamples2=chr17imputedBEACONsamples[!chr17imputedBEACONsamples %in% sampletable$localid] #4539 new controls (miss 2)
idx=match(chr17imputedBEACONsamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #2405 BE (miss 1)
table(sampletable$phenoEA_bca[idx]) #1508 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2176 controls (miss 1)

imputedBEACONsamples1=imputedBEACONsamples[imputedBEACONsamples %in% sampletable$localid]
imputedBEACONsamples2=imputedBEACONsamples[!imputedBEACONsamples %in% sampletable$localid] #4541 new controls
idx=match(imputedBEACONsamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #2406 BE
table(sampletable$phenoEA_bca[idx]) #1508 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2177

imputedqcBEACONsamples1=imputedqcBEACONsamples[imputedqcBEACONsamples %in% sampletable$localid]
imputedqcBEACONsamples2=imputedqcBEACONsamples[!imputedqcBEACONsamples %in% sampletable$localid] #4530 new controls
idx=match(imputedqcBEACONsamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #2401 BE
table(sampletable$phenoEA_bca[idx]) #1507 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2171

Cambridgesamples1=Cambridgesamples[Cambridgesamples %in% sampletable$localid]
sum(sampletable_cambridge %in% Cambridgesamples1)
Cambridgesamples2=Cambridgesamples[!Cambridgesamples %in% sampletable$localid] #3408 new controls
idx=match(Cambridgesamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #882 BE
table(sampletable$phenoEA_bca[idx]) #1003 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #0 controls

imputedCambridgesamples1=imputedCambridgesamples[imputedCambridgesamples %in% sampletable$localid]
imputedCambridgesamples2=imputedCambridgesamples[!imputedCambridgesamples %in% sampletable$localid] #3408 new controls
idx=match(imputedCambridgesamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #882 BE
table(sampletable$phenoEA_bca[idx]) #1003 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #0

imputedqcCambridgesamples1=imputedqcCambridgesamples[imputedqcCambridgesamples %in% sampletable$localid]
imputedqcCambridgesamples2=imputedqcCambridgesamples[!imputedqcCambridgesamples %in% sampletable$localid] #3408 new controls
idx=match(imputedqcCambridgesamples1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #881 BE
table(sampletable$phenoEA_bca[idx]) #1003 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #0

mergesample1=mergesample[mergesample %in% sampletable$localid]
idx=match(mergesample1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #3288 BE
table(sampletable$phenoEA_bca[idx]) #2511 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2177
mergesample2=mergesample[!mergesample %in% sampletable$localid] #7949=4541+3408

chr17mergesample1=chr17mergesample[chr17mergesample %in% sampletable$localid]
idx=match(chr17mergesample1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #3287 BE (miss 1)
table(sampletable$phenoEA_bca[idx]) #2511 EA
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2176 (miss 1)
chr17mergesample2=chr17mergesample[!chr17mergesample %in% sampletable$localid] #7947 (miss 2)

mergeqcsample1=mergeqcsample[mergeqcsample %in% sampletable$localid]
idx=match(mergeqcsample1,sampletable$localid)
table(sampletable$phenoBE_bca[idx]) #3282 BE (miss 6)
table(sampletable$phenoEA_bca[idx]) #2510 EA (miss 1)
sum(sampletable$phenoEA_bca[idx]==1 | sampletable$phenoBE_bca[idx]==1) #2171 (miss 6)
mergeqcsample2=mergeqcsample[!mergeqcsample %in% sampletable$localid] #7938 (miss 11)


#compare genotypes used for TWAS----------------------
prefix1="dist500K_GTEx_mucosa_Jan19"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_addcontrols_for_prediction.RData")
snp1=snp
snppos1=snppos

prefix2="dist500K_GTEx_mucosa_June11"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
snp2=snp
snppos2=snppos

sum(rownames(snppos1) %in% rownames(snppos2))
comsnps=intersect(rownames(snppos1),rownames(snppos2))
tmp1=snp1[match(comsnps,rownames(snppos1)),]
tmp2=snp2[match(comsnps,rownames(snppos2)),]
sum(tmp1!=tmp2) #0 all GTex genotypes are the same, no need to check this

load(paste0(outfolder1,"/bca_extractgenotype.RData"))
bcagenotype1=bcagenotype
load(paste0(outfolder2,"/bca_extractgenotype.RData"))
bcagenotype2=bcagenotype
tmp=sapply(1:ncol(bcagenotype2),function(x){
  tmp=unlist(strsplit(colnames(bcagenotype2)[[x]],"_"))
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(bcagenotype2)=tmp
sum(rownames(bcagenotype1) %in% rownames(bcagenotype2))
comsnps=intersect(rownames(bcagenotype1),rownames(bcagenotype2))
comsamples=intersect(colnames(bcagenotype1),colnames(bcagenotype2))
tmp1=bcagenotype1[match(comsnps,rownames(bcagenotype1)),match(comsamples,colnames(bcagenotype1))]
tmp2=bcagenotype2[match(comsnps,rownames(bcagenotype2)),match(comsamples,colnames(bcagenotype2))]
#sum(tmp1!=tmp2)
samplecor=sapply(1:ncol(tmp1),function(x){
  cor(tmp1[,x],tmp2[,x],use="complete")
})
quantile(samplecor)
# 0%       25%       50%       75%      100% 
# 0.9618792 0.9733894 0.9782746 0.9790821 0.9870781 
tmp=table(tmp1[,1],tmp2[,1])
sum(diag(tmp))/sum(tmp)
sampleacc=sapply(1:ncol(tmp1),function(x){
  if (length(unique(tmp1[,x]))==3 & length(unique(tmp2[,x]))==3)
  {
    tmp=table(tmp1[,x],tmp2[,x])
    tmp1=sum(diag(tmp))/sum(tmp)
  }else
  {
    tmp1=NA
  }
  tmp1
})
quantile(sampleacc)

#compare gene models
outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix1)
outfolder2=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix2)
load(paste0(outfolder1,"/preidiction_michigan_model.RData"))
res_min1=res_min
load(paste0(outfolder2,"/preidiction_michigan_model.RData"))
res_min2=res_min

comp_res_min=data.frame(numtotalsnp1=res_min1$numtotalsnp,numtotalsnp2=res_min2$numtotalsnp,numselectedsnp1=res_min1$numselectedsnp,
                        numselectedsnp2=res_min2$numselectedsnp,numcomsnp=0)
for (i in 1:nrow(res_min1))
{
  if (i %%2000==0) cat(i,'..')
  tmp=intersect(unlist(strsplit(res_min1$selectedsnps[i],"|",fixed = T)),unlist(strsplit(res_min2$selectedsnps[i],"|",fixed = T)))
  comp_res_min$numcomsnp[i]=length(tmp[!is.na(tmp)])
}
idx=which(comp_res_min$numselectedsnp1!=0 & comp_res_min$numselectedsnp2!=0) #22649
quantile(comp_res_min$numcomsnp[idx]/comp_res_min$numselectedsnp2[idx])
# 0%       25%       50%       75%      100% 
# 0.0000000 0.8214286 0.9047619 0.9642857 1.0000000
idx1=idx[which(comp_res_min$numcomsnp[idx]==0)] #16
quantile(comp_res_min$numselectedsnp1[idx1])

#compare predicted geneexp
load(paste0(outfolder1,"/bca_predict_geneexp.RData"))
predict_min1=predict_min[,3:ncol(predict_min)]
load(paste0(outfolder2,"/bca_predict_geneexp.RData"))
predict_min2=predict_min[,3:ncol(predict_min)]
geneexpsamplenames=strsplit(colnames(predict_min)[3:ncol(predict_min)],"_") #use localid
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_min2)=geneexpsamplenames
comsamples=intersect(colnames(predict_min1),colnames(predict_min2))
comgenes=intersect(rownames(predict_min1),rownames(predict_min2))
tmp1=predict_min1[match(comgenes,rownames(predict_min1)),match(comsamples,colnames(predict_min1))]
tmp2=predict_min2[match(comgenes,rownames(predict_min2)),match(comsamples,colnames(predict_min2))]
samplecor=sapply(1:ncol(tmp1),function(x){
  cor(tmp1[,x],tmp2[,x],use="complete")
})
quantile(samplecor)
# 0%       25%       50%       75%      100% 
# 0.8346116 0.8484896 0.8519404 0.8545010 0.8661318 

#"missing" samples
tmp=sampletable$localid[!sampletable$localid %in% colnames(predict_min1)]
idx=match(tmp,sampletable$localid)
table(sampletable$site[idx]) #AMOS
idx=which(sampletable$site==55)

melanoma=read.table("../data/AdditionalControlPuya/GENEVA_Melanoma_controls_caucasians_cleaned2.fam")
sum(melanoma$V2 %in% tmp) #0
sum(melanoma$V2 %in% sampletable$localid) #0

#check correlation between pcs
pca1=covariatetable[,c('pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','pc9','pc10')]
pca2=sampletable[,c('ev1_bca','ev2_bca','ev3_bca','ev4_bca')]
rownames(pca2)=sampletable$localid
colnames(pca2)=c('pc1','pc2','pc3','pc4')
comsamples=intersect(rownames(pca1),rownames(pca2))
idx1=match(comsamples,rownames(pca1))
idx2=match(comsamples,rownames(pca2))
tmp1=pca1[idx1,]
tmp2=pca2[idx2,]
tmp=cor(tmp1,tmp2)
#          pc1         pc2           pc3          pc4
# pc1  0.988635945 0.012985708  3.152872e-05 -0.006644249
# pc2 -0.003826767 0.147137512 -9.031539e-02 -0.001317920
# pc3 -0.005128448 0.942083511  4.495601e-02  0.102220004
# pc4 -0.032793925 0.001991751 -4.752410e-01  0.033535805
