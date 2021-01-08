#!/usr/bin/env Rscript
library(data.table)
plink="/fh/fast/stanford_j/Xiaoyu/Tools/plink-1.07-x86_64/plink1.9/plink"
sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
fam=read.table("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015.fam",stringsAsFactors = F)

BEsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$phenoBE_bc==2)]
EAsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$phenoEA_bc==2)]
COsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & (sampletable$phenoBE_bc==1 | sampletable$phenoEA_bc==1))]

GERDsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$recurrent_HB_RF==1)]
NGERDsamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$recurrent_HB_RF==0)]
Smokersamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$cig_smk_ever==1)]
NSmokersamples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$cig_smk_ever==0)]
BMI1samples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$bmi_recent_healthy<=30)]
BMI0samples=sampletable$localid[which(sampletable$localid %in% fam$V2 & sampletable$bmi_recent_healthy>30)]

outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDSC/"
filesforGWAS=function(casesamples=BEsamples,controlsamples=COsamples,pcafile="../result/bca_filtered_10Jan2019.pca.evec",prefix="BE_CO")
{
  idx1=match(casesamples,fam$V2)
  tmp1=data.frame(FID=fam$V1[idx1],IID=fam$V2[idx1],affected=2,stringsAsFactors = F)
  idx2=match(controlsamples,fam$V2)
  tmp2=data.frame(FID=fam$V1[idx2],IID=fam$V2[idx2],affected=1,stringsAsFactors = F)
  tmp=rbind(tmp1,tmp2)
  pcadat=read.table(pcafile,stringsAsFactors = F)
  tmp=tmp[tmp$IID %in% pcadat$V1,]
  idx=match(tmp$IID,pcadat$V1)
  tmp$pc1=pcadat$V2[idx]
  tmp$pc2=pcadat$V3[idx]
  tmp$pc3=pcadat$V4[idx]
  tmp$pc4=pcadat$V5[idx]
  idx=match(tmp$IID,sampletable$localid)
  tmp$age=sampletable$age[idx]
  tmp$sex=sampletable$sex[idx]  
  write.table(tmp[,c(1,2)],file=paste0(outfolder,prefix,"_selectedsamples_plink.txt"),row.names = F,col.names = F,sep="\t",quote=F)
  write.table(tmp,file=paste0(outfolder,prefix,"_selectedsamples_pheno_plink.txt"),row.names = F,col.names = T,sep="\t",quote=F)
}
update_fam=function(prefix="BEEA")
{
  famfile=paste0(outfolder,prefix,".fam")
  phenofile=paste0(outfolder,prefix,"_selectedsamples_pheno_plink.txt")
  fam=read.table(famfile,stringsAsFactors = F)
  pheno=read.table(phenofile,header=T,stringsAsFactors = F)
  idx=match(fam$V2,pheno$IID)
  fam$V6=pheno$affected[idx]
  write.table(fam,file=famfile,col.names=F,row.names = F,quote=F)
}
create_sumstats=function(prefix="BEEA",opt="nohm3")
{
  gwasfile=paste0(outfolder,prefix,".assoc.logistic")
  bimfile=paste0(outfolder,prefix,".bim")
  gwasdat=as.data.frame(fread(gwasfile))
  bimdat=as.data.frame(fread(bimfile))
  idx=which(bimdat$V1 %in% c(1:23))
  gwasdat=gwasdat[idx,]
  bimdat=bimdat[idx,]
  all(gwasdat$BP==bimdat$V4)
  if (opt=="hm3")
  {
    hm3=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/w_hm3.snplist"))
    idx=which(bimdat$V2 %in% hm3$SNP)
    prefix=paste0(prefix,"_hm3")
    gwasdat=gwasdat[idx,]
    bimdat=bimdat[idx,]
  }
  
  res=data.frame(SNP=bimdat$V2,hg19chr=bimdat$V1,bp=bimdat$V4,A1=bimdat$V5,A2=bimdat$V6,
                 or=gwasdat$OR,se=gwasdat$SE,pval=gwasdat$P,Z=gwasdat$STAT,N=gwasdat$NMISS,stringsAsFactors = F)
  write.table(res,file=paste0(outfolder,prefix,".sumstats"),col.names = T,row.names = F,sep=" ",quote=F)
}
readh2=function(prefix="BEEA",opt="1000g")
{
  tmp=read.table(paste0(outfolder,prefix,"_h2_",opt,".log"),stringsAsFactors =F,fill=T)
  h2=as.numeric(tmp$V5[27])
  se=tmp$V6[27]
  se=gsub("(","",se,fixed = T)
  se=as.numeric(gsub(")","",se,fixed = T))
  return(list(h2=h2,se=se))
}
#remember to run module load anaconda2 before start rstudio
compute_ldsc=function(plinkprefix="/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bca_filtered_19Feb2015",
                      casesamples=c(BEsamples,EAsamples),controlsamples=COsamples,prefix="BEEA",prevalence=0.016,
                      opt="1000g")
{
  print(prefix)
  print("Generate plink files")
  filesforGWAS(casesamples = casesamples,controlsamples = controlsamples,prefix=prefix)
  cmd=paste0(plink," --bfile ",plinkprefix," --keep ",paste0(outfolder,prefix,"_selectedsamples_plink.txt")," --make-bed --out ",paste0(outfolder,prefix))
  system(cmd)
  update_fam(prefix=prefix)
  print("GWAS")
  cmd=paste0(plink," --bfile  ",paste0(outfolder,prefix)," --covar  ",paste0(outfolder,prefix,"_selectedsamples_pheno_plink.txt"),"  --covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out ",paste0(outfolder,prefix))
  system(cmd)
  tmp=read.table(paste0(outfolder,prefix,".log"),fill=T,stringsAsFactors = F,sep="\n")
  idx=which(grepl("are cases and",tmp[,1]))
  tmp=tmp[idx,1]
  tmp=unlist(strsplit(tmp," "))
  idx=which(grepl("cases",tmp))
  ncase=as.integer(tmp[idx-2])
  idx=which(grepl("controls",tmp))
  ncontrol=as.integer(tmp[idx-2])
  samp_prev=as.character(ncase/(ncase+ncontrol))
  create_sumstats(prefix=prefix)
  if (file.exists(paste0(outfolder,prefix,".sumstats.gz"))) file.remove(paste0(outfolder,prefix,".sumstats.gz"))
  cmd=paste0("gzip ",paste0(outfolder,prefix,".sumstats"))
  system(cmd)
  print("ldsc")
  if (opt=="1000g")
  {
    cmd1=paste0("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py ",
               "--h2 ",paste0(outfolder,prefix,".sumstats.gz"),
               " --ref-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/ ",
               "--w-ld-chr /fh/fast/dai_j/CancerGenomics/Tools/ldsc/eur_w_ld_chr/  ",
               "--samp-prev ",samp_prev,
               " --pop-prev ",prevalence,
               " --out ",paste0(outfolder,prefix,"_h2_",opt))
  }
  if (opt=="computeld")
  {
    ldfolder=paste0(outfolder,prefix,"w_ld_chr/")
    if (! dir.exists(ldfolder)) dir.create(ldfolder)
    for (chr in 1:22)
    {
      cmd=paste0(plink," --bfile ",paste0(outfolder,prefix)," --chr ",chr," --make-bed --out ",ldfolder,chr)
      system(cmd)
      cmd=paste0(plink," --bfile ",ldfolder,chr," --cm-map /fh/fast/dai_j/CancerGenomics/Tools/database/1000genome/1000GP_Phase3/genetic_map_chr",chr,"_combined_b37.txt ",chr," --make-bed --out ",ldfolder,chr)
      system(cmd)
      cmd=paste0("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py --bfile ",ldfolder,chr," --l2  --ld-wind-cm 1 --out ",ldfolder,chr)
      system(cmd)
    }
    cmd1=paste0("/fh/fast/dai_j/CancerGenomics/Tools/ldsc/ldsc.py ",
               "--h2 ",paste0(outfolder,prefix,".sumstats.gz"),
               " --ref-ld-chr ",ldfolder,
               " --w-ld-chr   ",ldfolder,
               " --samp-prev ",samp_prev,
               " --pop-prev ",prevalence,
               " --out ",paste0(outfolder,prefix,"_h2_",opt))
  }
  system(cmd1)
  h2res=readh2(prefix=prefix,opt=opt)
  res=data.frame(prefix=prefix,ncase=ncase,ncontrol=ncontrol,h2=h2res$h2,se=h2res$se,stringsAsFactors = F)
  return(res)
}
BEEA_res=compute_ldsc()
BEEA_recurrent0_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),NGERDsamples),
                                 controlsamples = intersect(COsamples,NGERDsamples),prefix="BEEA_recurrent0")  
BEEA_recurrent1_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),GERDsamples),
                                 controlsamples = intersect(COsamples,GERDsamples),prefix="BEEA_recurrent1") 
BEEA_bmi0_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),BMI0samples),
                                 controlsamples = intersect(COsamples,BMI0samples),prefix="BEEA_bmi0")
BEEA_bmi1_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),BMI1samples),
                           controlsamples = intersect(COsamples,BMI1samples),prefix="BEEA_bmi1")
BEEA_smoke0_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),NSmokersamples),
                           controlsamples = intersect(COsamples,NSmokersamples),prefix="BEEA_smoke0")
BEEA_smoke1_res=compute_ldsc(casesamples = intersect(c(BEsamples,EAsamples),Smokersamples),
                             controlsamples = intersect(COsamples,Smokersamples),prefix="BEEA_smoke1")
ldsc_res=rbind(BEEA_res,BEEA_recurrent0_res,BEEA_recurrent1_res,BEEA_bmi0_res,BEEA_bmi1_res,BEEA_smoke0_res,
               BEEA_smoke1_res)
write.csv(ldsc_res,file="../result/ldsc_res.csv",row.names=F)
BE_puya_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
                    casesamples=BEsamples,controlsamples=COsamples,prefix="BE_puya",prevalence=0.016,
                    opt="computeld")
EA_puya_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
                         casesamples=EAsamples,controlsamples=COsamples,prefix="EA_puya",prevalence=0.0025,
                         opt="computeld")
BEEA_puya_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
                         casesamples=c(BEsamples,EAsamples),controlsamples=COsamples,prefix="BEEA_puya",prevalence=0.016,
                         opt="computeld")
BEEA_puya_recurrent0_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),NGERDsamples),
                                 controlsamples = intersect(COsamples,NGERDsamples),prefix="BEEA_puya_recurrent0",opt="computeld")  
BEEA_puya_recurrent1_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),GERDsamples),
                                 controlsamples = intersect(COsamples,GERDsamples),prefix="BEEA_puya_recurrent1",opt="computeld") 
BEEA_puya_bmi0_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),BMI0samples),
                           controlsamples = intersect(COsamples,BMI0samples),prefix="BEEA_puya_bmi0",opt="computeld")
BEEA_puya_bmi1_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),BMI1samples),
                           controlsamples = intersect(COsamples,BMI1samples),prefix="BEEA_puya_bmi1",opt="computeld")
BEEA_puya_smoke0_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),NSmokersamples),
                             controlsamples = intersect(COsamples,NSmokersamples),prefix="BEEA_puya_smoke0",opt="computeld")
BEEA_puya_smoke1_res=compute_ldsc(plinkprefix="../result/bca_filtered_puya1",
  casesamples = intersect(c(BEsamples,EAsamples),Smokersamples),
                             controlsamples = intersect(COsamples,Smokersamples),prefix="BEEA_puya_smoke1",opt="computeld")
ldsc_puya_res=rbind(BE_puya_res,EA_puya_res,BEEA_puya_res,BEEA_puya_recurrent0_res,BEEA_puya_recurrent1_res,BEEA_puya_bmi0_res,
                    BEEA_puya_bmi1_res,BEEA_puya_smoke0_res,BEEA_puya_smoke1_res)
write.csv(ldsc_puya_res,file="../result/ldsc_puya_res.csv",row.names=F)
