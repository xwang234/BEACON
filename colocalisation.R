#!/usr/bin/env Rscript
#https://cran.r-project.org/web/packages/coloc/vignettes/a06_SuSiE.html
#https://cran.r-project.org/web/packages/coloc/vignettes/a02_data.html

#load summary data
library(data.table)
library(GenomicRanges)
#summary stat from validation
summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"
BE_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_Bonn_autosomes.txt"),header=T))
BE_Oxfordsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_oxford_autosomes.txt"),header=T))

EA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"EA_Bonn_autosomes.txt"),header=T))
BEEA_Bonnsummarydat=as.data.frame(fread(paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"),header=T))

tmp=intersect(BE_Bonnsummarydat$SNP,EA_Bonnsummarydat$SNP)

tmp=intersect("SNP",colnames(BE_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="pvalue")]="P"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="beta")]="BETA"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Bonnsummarydat)[which(colnames(BE_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

tmp=intersect("SNP",colnames(BE_Oxfordsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="rsid")]="SNP"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="pvalue")]="P"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="beta")]="BETA"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="non-effect-allele")]="non_effect_allele"
}



tmp=intersect("SNP",colnames(EA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="pvalue")]="P"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="beta")]="BETA"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(EA_Bonnsummarydat)[which(colnames(EA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}


tmp=intersect("SNP",colnames(BEEA_Bonnsummarydat))
if (length(tmp)==0)
{
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="rsid")]="SNP"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="pvalue")]="P"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="beta")]="BETA"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="effect-allele")]="effect_allele"
  colnames(BEEA_Bonnsummarydat)[which(colnames(BEEA_Bonnsummarydat)=="non-effect-allele")]="non_effect_allele"
}

add_chr_tosummarydat=function(dat1=BE_Bonnsummarydat)
{
  dat2=BE_Oxfordsummarydat
  if (!"chr" %in% colnames(dat1))
  {
    dat1$chr=NA
    tmp=intersect(dat1$SNP,dat2$SNP)
    idx1=match(tmp,dat1$SNP)
    idx2=match(tmp,dat2$SNP)
    dat1$chr[idx1]=dat2$chr[idx2]
    # idx=which(is.na(dat1$chr))
    # tmp=intersect(dat1$SNP[idx],dbsnp$V2)
    # idx1=match(tmp,dat1$SNP)
    # idx2=match(tmp,dbsnp$V2)
    # dat1$chr[idx1]=dbsnp$V1[idx2]
    print(paste0(sum(is.na(dat1$chr))/nrow(dat1),"of all snps are missing chr"))
  }
  dat1$snpname=dat1$snpname1=NA
  dat1$snpname=paste0(dat1$chr,":",dat1$position,"_",dat1$effect_allele,"_",dat1$non_effect_allele)
  dat1$snpname1=paste0(dat1$chr,":",dat1$position,"_",dat1$non_effect_allele,"_",dat1$effect_allele)
  return(dat1)
}
BE_Bonnsummarydat=add_chr_tosummarydat(dat1=BE_Bonnsummarydat)
EA_Bonnsummarydat=add_chr_tosummarydat(dat1=EA_Bonnsummarydat)
BEEA_Bonnsummarydat=add_chr_tosummarydat(dat1=BEEA_Bonnsummarydat)
BE_Oxfordsummarydat=add_chr_tosummarydat(dat1=BE_Oxfordsummarydat)

#change cooridates of hg38 to hg19. Summarydata use hg19, BCA use hg38
hg38tohg19=function(snpnames=rownames(snp)[1:10])
{
  library(rtracklayer)
  library(GenomicRanges)
  chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
  tmp=unlist(strsplit(snpnames,":"))
  chr0=tmp[1]
  chr=paste0("chr",tmp[1])
  tmp=tmp[seq(2,length(tmp),2)]
  tmp=unlist(strsplit(tmp,"_"))
  pos=as.integer(tmp[seq(1,length(tmp),3)])
  alt=tmp[seq(2,length(tmp),3)]
  ref=tmp[seq(3,length(tmp),3)]
  gr_dat=GRanges(seqnames = chr,ranges=IRanges(start=pos,width=1))
  tmp=liftOver(gr_dat,chain)
  newsnpnames=newpos=rep(NA,length(tmp))
  for (i in 1:length(tmp))
  {
    tmp1=unlist(tmp[i])
    if (length(tmp1)==0)
    {
      warning(paste0(snpnames[i]," not transformed"))
    }else
    {
      if (length(tmp1)==1)
      {
        newpos[i]=start(tmp1)
      }else
      {
        warning(paste0(snpnames[i]," transformed to ",length(tmp1)," snps"))
        newpos[i]=start(tmp1)[1]
      }
    }
  }
  newsnpnames=paste0(chr0,":",newpos,"_",alt,"_",ref)
  res=data.frame(snphg38=snpnames,snphg19=newsnpnames,stringsAsFactors = F)
  return(res)
}

library("biomaRt")
#mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
#snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp") #hg38
snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp",host="uswest.ensembl.org")

#faster than using biomart to find snpnames, but can't find all the results
find_snpname=function(selsnps,summarydat=BE_Bonnsummarydat)
{
  selsnps_hg19=hg38tohg19(snpnames = selsnps)
  selsnps_hg19$selsnps=selsnps
  selsnps_hg19$snp=NA
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname1)
  idx1=match(tmp,selsnps_hg19$snphg19)
  idx2=match(tmp,summarydat$snpname1)
  selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
  return(selsnps_hg19)
}



library(coloc)

#the table including interesting genes
#need trait and orgon info for each gene
table3=read.csv("../result/TWAS_table3.csv")
genes=c("EXOC3",
        "ZNF641",
        "HSP90AA1"
)
table3=table3[table3$gene %in% genes,]

#GTEx and Bonn, use selected SNPs
colcal=function(i=1)
{
  gene=table3$gene[i]
  trait=table3$Trait[i]
  organ=table3$organ[i]
  print(gene)
  
  #first work on GTEx
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
  idx=which(rownames(res_min)==gene)
  selsnps=unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T))
  
  coeffs=as.numeric(unlist(strsplit(res_min$selectedsnps_coeff[idx],"|",fixed=T)))
  GTEx_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
  idx2=which(rownames(phenotype)==gene)
  Y=unlist(phenotype[idx2,])
  for (j in 1:length(selsnps))
  {
    idx1=which(rownames(snp)==selsnps[j])
    dat=cbind.data.frame(snp=t(snp[idx1,]),covariate)
    colnames(dat)[1]="snp"
    
    fit=lm(Y~.,data=dat)
    thecoeff=summary(fit)$coefficients
    GTEx_beta$beta[j]=thecoeff[2,1]
    GTEx_beta$varbeta[j]=thecoeff[2,2]^2 #variance
    GTEx_beta$pvalue[j]=thecoeff[2,4]
  }
  
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  selsnpinfo=tmp
  selsnpinfo$hg19_position=NA
  for (j in 1:nrow(selsnpinfo))
  {
    tmp=unlist(strsplit(selsnpinfo$snphg19[j],"_"))
    selsnpinfo$hg19_position[j]=as.numeric(unlist(strsplit(tmp[1],":"))[2])
  }
  GTExdat=list()
  GTExdat$beta=GTEx_beta$beta
  names(GTExdat$beta)=selsnpinfo$snp
  GTExdat$varbeta=GTEx_beta$varbeta
  names(GTExdat$varbeta)=selsnpinfo$snp
  GTExdat$N=length(Y)
  GTExdat$sdY=sd(Y)
  GTExdat$type="quant"
  GTExdat$snp=selsnpinfo$snp
  GTExdat$postion=selsnpinfo$hg19_position
  
  
  #validation GWAS data
  type=table3$Trait[i]
  Bonn_N=NA
  if (type=="BE")
  {
    sumdat=BE_Bonnsummarydat
    Bonn_N=1037+3537
  }
  if (type=="EA")
  {
    sumdat=EA_Bonnsummarydat
    Bonn_N=1609+3537
  }
  if (type=="BEEA")
  {
    sumdat=BEEA_Bonnsummarydat
    Bonn_N=1037+1609+3537
  }
  idx1=which(selsnpinfo$snphg19 %in% sumdat$snpname)
  idx2=which(selsnpinfo$snphg19 %in% sumdat$snpname1)
  Bonn_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
  for (j in 1:length(selsnps))
  {
    idx1=which(sumdat$snpname==selsnpinfo$snphg19[j])
    if (length(idx1)>0)
    {
     Bonn_beta$beta[j]=sumdat$BETA[idx1]
     Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
     Bonn_beta$pvalue[j]=sumdat$P[idx1]
    }
    idx1=which(sumdat$snpname1==selsnpinfo$snphg19[j])
    if (length(idx1)>0)
    {
      Bonn_beta$beta[j]=-sumdat$BETA[idx1]
      Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
      Bonn_beta$pvalue[j]=sumdat$P[idx1]
    }
  }
  
  Bonndat=list()
  Bonndat$beta=Bonn_beta$beta
  names(Bonndat$beta)=selsnpinfo$snp
  Bonndat$varbeta=Bonn_beta$varbeta
  names(Bonndat$varbeta)=selsnpinfo$snp
  Bonndat$N=Bonn_N
  Bonndat$type="cc"
  Bonndat$snp=selsnpinfo$snp
  Bonndat$postion=selsnpinfo$hg19_position
  
  #compute LD
  chr=unlist(strsplit(selsnps[1],":"))[1]
  #read data from 1000 genome v3 EUR
  thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
  refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
  #refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
  refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
  refbimstr1=paste0(refbim$V1,":",refbim$V4,"_",refbim$V5,"_",refbim$V6)
  refbimstr2=paste0(refbim$V1,":",refbim$V4,"_",refbim$V6,"_",refbim$V5)
  
  
  
  Z1=Z2=NULL
  tmp1=intersect(selsnpinfo$snphg19,refbimstr1)
  tmp2=intersect(selsnpinfo$snphg19,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  tmp=selsnpinfo$snphg19[selsnpinfo$snphg19 %in% colnames(Z3)]
  if (length(tmp)<length(selsnps)) print(paste0(gene,":",length(selsnps)-length(tmp)," snps are missing in 1000genome control"))
  Z3=Z3[,match(tmp,colnames(Z3))]
  LD=cor(Z3)
  idx=match(colnames(Z3),selsnpinfo$snphg19)
  rownames(LD)=colnames(LD)=selsnpinfo$snp[idx]
  
  tmp=Bonndat$snp[!is.na(Bonndat$beta)]
  tmp1=intersect(tmp,rownames(LD))
  idx=match(tmp1,Bonndat$snp)
  Bonndat$beta=Bonndat$beta[idx]
  Bonndat$varbeta=Bonndat$varbeta[idx]
  Bonndat$snp=Bonndat$snp[idx]
  Bonndat$postion=Bonndat$postion[idx]
  GTExdat$LD=LD
  Bonndat$LD=LD[idx,idx]
  Bonndat$LD=LD[idx,idx]
  
  idx=match(tmp1,GTExdat$snp)
  GTExdat$beta=GTExdat$beta[idx]
  GTExdat$varbeta=GTExdat$varbeta[idx]
  GTExdat$snp=GTExdat$snp[idx]
  GTExdat$postion=GTExdat$postion[idx]
  GTExdat$LD=LD[idx,idx]
  
  check_dataset(GTExdat)
  check_dataset(Bonndat)
  my.res <- coloc.abf(dataset1=GTExdat, dataset2=Bonndat)
  print(my.res)
  pdf(paste0("../result/coloc_",gene,"_abf.pdf"),width=12)
  sensitivity(my.res,"H4 > 0.8")
  dev.off()
  
  S3=runsusie(GTExdat)
  S4=runsusie(Bonndat)
  susie.res=NULL
  if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(S3,S4)
    print(susie.res$summary)
  }
  if (!is.null(susie.res$summary))
  {
    if(requireNamespace("susieR",quietly=TRUE) & nrow(susie.res$summary)>1) {
      pdf(paste0("../result/coloc_",gene,"_susie.pdf"),width=12)
      sensitivity(susie.res,"H4 > 0.8",row=1,dataset1=GTExdat,dataset2=Bonndat)
      sensitivity(susie.res,"H4 > 0.8",row=2,dataset1=GTExdat,dataset2=Bonndat)
      dev.off()
    }
  }
  
  return(list(GTExdat=GTExdat,Bonndat=Bonndat,GTEx_beta=GTEx_beta,Bonn_beta=Bonn_beta))
}
res1=colcal(1)
res2=colcal(2)
res3=colcal(3)
save(res1,res2,res3,file="../result/colocalisationRES.RData")
# #example:
# library(coloc)
# data(coloc_test_data)
# attach(coloc_test_data) ## datasets D1, D2, D3 and D4
# S3=runsusie(D3)
# S4=runsusie(D4)
# if(requireNamespace("susieR",quietly=TRUE)) {
#   susie.res=coloc.susie(S3,S4)
#   print(susie.res$summary)
# }
# if(requireNamespace("susieR",quietly=TRUE)) {
#   sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=D3,dataset2=D4)
#   sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=D3,dataset2=D4)
# }

#GTEx and BCA, use all common SNPs
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  if (length(idx)>0)
    sampletable[idx,i]=NA
}

#BCA covariate table, principle components
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
  tmp1=cbind(eigsamples,tmp)
  return(tmp1)
}

covariatetable=readeigenstrat()
rownames(covariatetable)=covariatetable$V2=gsub("SEP","",rownames(covariatetable))
#all(colnames(predict_min)[3:ncol(predict_min)] %in% rownames(covariatetable))
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
genes=table3$gene
traits=table3$Trait
organs=table3$organ
#to check 2 genes above and under ZNF641


colcal1=function(i=1)
{
  gene=genes[i]
  trait=traits[i]
  organ=organs[i]
  print(gene)
  
  #first work on GTEx
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1))
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2))
  idx=which(rownames(phenotype)==gene)
  tmp=distance(gr_snp,gr_pos[idx])
  idx=which(tmp<5e5)
  selsnps=rownames(snp)[idx]
  
  GTEx_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
  idx2=which(rownames(phenotype)==gene)
  Y=unlist(phenotype[idx2,])
  for (j in 1:length(selsnps))
  {
    idx1=which(rownames(snp)==selsnps[j])
    dat=cbind.data.frame(snp=t(snp[idx1,]),covariate)
    colnames(dat)[1]="snp"
    
    fit=lm(Y~.,data=dat)
    thecoeff=summary(fit)$coefficients
    GTEx_beta$beta[j]=thecoeff[2,1]
    GTEx_beta$varbeta[j]=thecoeff[2,2]^2 #variance
    GTEx_beta$pvalue[j]=thecoeff[2,4]
  }
  
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      #sometimes biomart has connection problems
      tmp2=tryCatch(
        {
          getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        },
        error=function(e)
        {
          return(F)
        }
      )
      
      #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
      if (class(tmp2)[1]=="data.frame")
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  selsnpinfo=tmp
  selsnpinfo$hg19_position=NA
  selsnpinfo$snphg19_=NA
  for (j in 1:nrow(selsnpinfo))
  {
    tmp=unlist(strsplit(selsnpinfo$snphg19[j],"_"))
    selsnpinfo$hg19_position[j]=as.numeric(unlist(strsplit(tmp[1],":"))[2])
    tmp=unlist(strsplit(selsnpinfo$snphg19[j],"_"))
    selsnpinfo$snphg19_[j]=paste(tmp[c(1,3,2)],collapse = "_")
  }
  GTExdat=list()
  GTExdat$beta=GTEx_beta$beta
  names(GTExdat$beta)=selsnpinfo$snp
  GTExdat$varbeta=GTEx_beta$varbeta
  names(GTExdat$varbeta)=selsnpinfo$snp
  GTExdat$N=length(Y)
  GTExdat$sdY=sd(Y)
  GTExdat$type="quant"
  GTExdat$snp=selsnpinfo$snp
  GTExdat$postion=selsnpinfo$hg19_position
  
  #BCA data
  chr=unlist(strsplit(selsnps[1],":"))[1]
  BCAsnp=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/chr",chr,"_filter_hg19tohg38_flip.traw")))
  rownames(BCAsnp)=BCAsnp$SNP
  colnames(BCAsnp)=gsub("SEP","",colnames(BCAsnp))
  BCAsnp=BCAsnp[,7:ncol(BCAsnp)]
  if (!all(colnames(BCAsnp) %in% paste0(covariatetable$V1,"_",covariatetable$V2))& ncol(BCAsnp)!=nrow(covariatetable)) print("some samples are missing")
  mycovariatetable=covariatetable[match(colnames(BCAsnp),paste0(covariatetable$V1,"_",covariatetable$V2)),]
  tmp1=rownames(BCAsnp)[rownames(BCAsnp) %in% selsnpinfo$snphg19]
  tmp2=rownames(BCAsnp)[rownames(BCAsnp) %in% selsnpinfo$snphg19_]
  if (length(tmp1)+length(tmp2)<length(selsnps)) print("some snps are missing in BCA")
  if (length(tmp2)>0)
  {
    dat1=BCAsnp[match(tmp1,rownames(BCAsnp)),]
    dat2=BCAsnp[match(tmp2,rownames(BCAsnp)),]
    dat2=2-dat2
    idx=match(rownames(dat2),selsnpinfo$snphg19_)
    rownames(dat2)=selsnpinfo$snphg19[idx]
    dat3=rbind(dat1,dat2)
    dat3=dat3[match(selsnpinfo$snphg19,rownames(dat3)),]
  }else
  {
    dat3=BCAsnp[match(tmp1,rownames(BCAsnp)),]
  }
  
  type=traits[i]
  if (type=="BE")
  {
    Y_BCA=mycovariatetable$phenoBE_bca
    N_BCA=2401+881+6701+3408
  }
  if (type=="EA")
  {
    Y_BCA=mycovariatetable$phenoEA_bca
    N_BCA=1507+1003+6701+3408
  }
  if (type=="BEEA")
  {
    Y_BCA=mycovariatetable$phenoEABE_bca
    N_BCA=2401+881+1507+1003+6701+3408
  }
  BCA_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
  BCA_beta$freq_non_effect_allele_cases=BCA_beta$freq_non_effect_allele_controls=NA
  for (j in 1:nrow(dat3))
  {
    dat=data.frame(snp=dat3[j,])
    dat=cbind(t(dat),mycovariatetable[,c("pc1","pc2","pc3","pc4","pc5","pc6")])
    colnames(dat)[1]="snp"
    fit=glm(I(Y_BCA==2)~.,data=dat,family="binomial")
    tmp=summary(fit)$coefficients
    BCA_beta$beta[j]=tmp[2,1]
    BCA_beta$varbeta[j]=tmp[2,2]^2
    BCA_beta$pvalue[j]=tmp[2,4]
    tmp=unlist(dat3[j,which(Y_BCA==1)])
    BCA_beta$freq_non_effect_allele_controls[j]=(sum(tmp==0)*2+sum(tmp==1))/length(tmp)/2
    tmp=unlist(dat3[j,which(Y_BCA==2)])
    BCA_beta$freq_non_effect_allele_cases[j]=(sum(tmp==0)*2+sum(tmp==1))/length(tmp)/2
    
  }
  
  BCAdat=list()
  BCAdat$beta=BCA_beta$beta
  names(BCAdat$beta)=selsnpinfo$snp
  BCAdat$varbeta=BCA_beta$varbeta
  names(BCAdat$varbeta)=selsnpinfo$snp
  BCAdat$N=N_BCA
  BCAdat$type="cc"
  BCAdat$snp=selsnpinfo$snp
  BCAdat$postion=selsnpinfo$hg19_position
  
  #compute LD
  chr=unlist(strsplit(selsnps[1],":"))[1]
  #read data from 1000 genome v3 EUR
  thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
  refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
  #refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
  refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
  #refbimstr=paste0(refbim$V1,":",refbim$V4)
  refbimstr1=paste0(refbim$V1,":",refbim$V4,"_",refbim$V5,"_",refbim$V6)
  refbimstr2=paste0(refbim$V1,":",refbim$V4,"_",refbim$V6,"_",refbim$V5)
  
  Z1=Z2=NULL
  tmp1=intersect(selsnpinfo$snphg19,refbimstr1)
  tmp2=intersect(selsnpinfo$snphg19,refbimstr2)
  idx1=match(tmp1,refbimstr1)
  if (length(tmp1)>0)
  {
    Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    rownames(Z1)=tmp1
  }
  idx2=match(tmp2,refbimstr2)
  if (length(tmp2)>0)
  {
    Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #change id
    rownames(Z2)=tmp2
    Z2=2-Z2
  }
  Z3=as.data.frame(t(rbind(Z1,Z2)))
  tmp=selsnpinfo$snphg19[selsnpinfo$snphg19 %in% colnames(Z3)]
  if (length(tmp)<length(selsnps)) print(paste0(gene,":",length(selsnps)-length(tmp)," snps are missing in 1000genome control"))
  Z3=Z3[,match(tmp,colnames(Z3))]
  tmp=rowMeans(Z3)
  Z3=Z3[!is.na(tmp),]
  LD=cor(Z3)
  idx=match(colnames(Z3),selsnpinfo$snphg19)
  rownames(LD)=colnames(LD)=selsnpinfo$snp[idx]
  
  tmp=selsnpinfo$snp[idx]
  idx=match(tmp,BCAdat$snp)
  BCAdat$beta=BCAdat$beta[idx]
  BCAdat$varbeta=BCAdat$varbeta[idx]
  BCAdat$snp=BCAdat$snp[idx]
  BCAdat$postion=BCAdat$postion[idx]
  BCAdat$LD=LD
  
  idx=match(tmp,GTExdat$snp)
  GTExdat$beta=GTExdat$beta[idx]
  GTExdat$varbeta=GTExdat$varbeta[idx]
  GTExdat$snp=GTExdat$snp[idx]
  GTExdat$postion=GTExdat$postion[idx]
  GTExdat$LD=LD
  
  check_dataset(GTExdat)
  check_dataset(BCAdat)
  my.res <- coloc.abf(dataset1=GTExdat, dataset2=BCAdat)
  print(my.res)
  pdf(paste0("../result/coloc1_",gene,"_abf.pdf"),width=12)
  sensitivity(my.res,"H4 > 0.8")
  dev.off()
  
  S3=runsusie(GTExdat)
  S4=runsusie(BCAdat)
  susie.res=NULL
  if(requireNamespace("susieR",quietly=TRUE)) {
    susie.res=coloc.susie(S3,S4)
    print(susie.res$summary)
  }
  if (!is.null(susie.res$summary))
  {
    if(requireNamespace("susieR",quietly=TRUE) & nrow(susie.res$summary)>1) {
      pdf(paste0("../result/coloc1_",gene,"_susie.pdf"),width=12)
      sensitivity(susie.res,"H4 > 0.8",row=1,dataset1=GTExdat,dataset2=Bonndat)
      sensitivity(susie.res,"H4 > 0.8",row=2,dataset1=GTExdat,dataset2=Bonndat)
      dev.off()
    }
  }
  
  
  #validation Bonn data
  type=traits[i]
  Bonn_N=NA
  if (type=="BE")
  {
    sumdat=BE_Bonnsummarydat
    Bonn_N=1037+3537
  }
  if (type=="EA")
  {
    sumdat=EA_Bonnsummarydat
    Bonn_N=1609+3537
  }
  if (type=="BEEA")
  {
    sumdat=BEEA_Bonnsummarydat
    Bonn_N=1037+1609+3537
  }
  idx1=which(selsnpinfo$snphg19 %in% sumdat$snpname)
  idx2=which(selsnpinfo$snphg19 %in% sumdat$snpname1)
  Bonn_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
  Bonn_beta$freq_non_effect_allele_cases=Bonn_beta$freq_non_effect_allele_controls=NA
  for (j in 1:length(selsnps))
  {
    idx1=which(sumdat$snpname==selsnpinfo$snphg19[j])
    if (length(idx1)>0)
    {
      Bonn_beta$beta[j]=sumdat$BETA[idx1]
      Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
      Bonn_beta$pvalue[j]=sumdat$P[idx1]
      Bonn_beta$freq_non_effect_allele_cases[j]=sumdat$freq_non_effect_allele_cases[idx1]
      Bonn_beta$freq_non_effect_allele_controls[j]=sumdat$freq_non_effect_allele_controls[idx1]
    }
    idx1=which(sumdat$snpname1==selsnpinfo$snphg19[j])
    if (length(idx1)>0)
    {
      Bonn_beta$beta[j]=-sumdat$BETA[idx1]
      Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
      Bonn_beta$pvalue[j]=sumdat$P[idx1]
      Bonn_beta$freq_non_effect_allele_cases[j]=1-sumdat$freq_non_effect_allele_cases[idx1]
      Bonn_beta$freq_non_effect_allele_controls[j]=1-sumdat$freq_non_effect_allele_controls[idx1]
    }
  }
  
  Bonndat=list()
  Bonndat$beta=Bonn_beta$beta
  names(Bonndat$beta)=selsnpinfo$snp
  Bonndat$varbeta=Bonn_beta$varbeta
  names(Bonndat$varbeta)=selsnpinfo$snp
  Bonndat$N=Bonn_N
  Bonndat$type="cc"
  Bonndat$snp=selsnpinfo$snp
  Bonndat$postion=selsnpinfo$hg19_position
  
  
  allbeta=data.frame(snp=selsnpinfo$snp,chr=chr,poistion=selsnpinfo$hg19_position,effect_allele=NA,non_effect_allel=NA,Beta_eQTL=GTEx_beta$beta,Pvalue_eQTL=GTEx_beta$pvalue,
                     Beta_BCA=BCA_beta$beta,Pvalue_BCA=BCA_beta$pvalue,Beta_Bonn=Bonn_beta$beta,Pvalue_Bonn=Bonn_beta$pvalue,selected=F,
                     freq_non_effect_allele_cases_Bonn=Bonn_beta$freq_non_effect_allele_cases,freq_non_effect_allele_controls_Bonn=Bonn_beta$freq_non_effect_allele_controls,
                     freq_non_effect_allele_cases_BCA=BCA_beta$freq_non_effect_allele_cases,freq_non_effect_allele_controls_BCA=BCA_beta$freq_non_effect_allele_controls)
  for (j in 1:nrow(allbeta))
  {
    tmp=unlist(strsplit(selsnpinfo$selsnps[j],"_"))
    allbeta$effect_allele[j]=tmp[2]
    allbeta$non_effect_allele[j]=tmp[3]
  }
  idx=which(rownames(res_min)==gene)
  selsnps1=unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T))
  idx=match(selsnps1,selsnps)
  allbeta$selected[idx]=T
  
  return(list(GTExdat=GTExdat,Bonndat=Bonndat,BCAdat=BCAdat,GTEx_beta=GTEx_beta,Bonn_beta=Bonn_beta,BCA_beta=BCA_beta,allbeta=allbeta,selsnpinfo=selsnpinfo))
}

res1=colcal1(1)
write.csv(res1$allbeta,file="../result/EXOC3_allsnps.csv",row.names = F)
res2=colcal1(2)
write.csv(res2$allbeta,file="../result/ZNF641_allsnps.csv",row.names = F)
res3=colcal1(3)
write.csv(res3$allbeta,file="../result/HSP90AA1_allsnps.csv",row.names = F)
save(res1,res2,res3,file="../result/colocalisationRES1.RData")

par(mar=c(6,6,2,1))
allbeta1=res1$allbeta

par(mfrow=c(1,2))
plot(allbeta1$freq_non_effect_allele_controls_BCA,allbeta1$freq_non_effect_allele_controls_Bonn,xlab="Frequency of non_effect_allele_control (BCA)",ylab="Frequency of non_effect_allele_control (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta1$freq_non_effect_allele_cases_BCA,allbeta1$freq_non_effect_allele_cases_Bonn,xlab="Frequency of non_effect_allele_case (BCA)",ylab="Frequency of non_effect_allele_case (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
# idx=which(allbeta1$Beta_BCA*allbeta1$Beta_Bonn<0)
# points(allbeta1$freq_non_effect_allele_controls_BCA[idx],allbeta1$freq_non_effect_allele_controls_Bonn[idx],col="blue")
# View(allbeta1[idx,])
par(mfrow=c(1,2))
plot(allbeta1$Beta_BCA,allbeta1$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta1$freq_non_effect_allele_cases_BCA-allbeta1$freq_non_effect_allele_controls_BCA,allbeta1$freq_non_effect_allele_cases_Bonn-allbeta1$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

plot(allbeta1$Beta_eQTL,allbeta1$Beta_BCA,xlab="Beta (eQTL)",ylab="Beta (BCA)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta1$Beta_eQTL,allbeta1$Beta_Bonn,xlab="Beta (eQTL)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

# tmp1=allbeta1$Beta_BCA
# tmp2=allbeta1$Beta_Bonn
# tmp1[idx]=-tmp1[idx]
# plot(tmp1,tmp2)

allbeta2=res2$allbeta
plot(allbeta2$freq_non_effect_allele_controls_BCA,allbeta2$freq_non_effect_allele_controls_Bonn,xlab="Frequency of non_effect_allele_control (BCA)",ylab="Frequency of non_effect_allele_control (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta2$freq_non_effect_allele_cases_BCA,allbeta2$freq_non_effect_allele_cases_Bonn,xlab="Frequency of non_effect_allele_case (BCA)",ylab="Frequency of non_effect_allele_case (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
# idx=which(allbeta2$Beta_BCA*allbeta2$Beta_Bonn<0)
# points(allbeta2$freq_non_effect_allele_controls_BCA[idx],allbeta2$freq_non_effect_allele_controls_Bonn[idx],col="blue")
# View(allbeta2[idx,])

plot(allbeta2$Beta_BCA,allbeta2$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta2$freq_non_effect_allele_cases_BCA-allbeta2$freq_non_effect_allele_controls_BCA,allbeta2$freq_non_effect_allele_cases_Bonn-allbeta2$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

plot(allbeta2$Beta_eQTL,allbeta2$Beta_BCA,xlab="Beta (eQTL)",ylab="Beta (BCA)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta2$Beta_eQTL,allbeta2$Beta_Bonn,xlab="Beta (eQTL)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

#check a snp j=801, 12:48186927_T_G 12:48580710_T_G 12:48186927_T_G rs6580657      48580710 12:48580710_G_T
thebim=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/chr12_filter_hg19tohg38_flip.bim"))
idx=which(thebim$V4==48186927)
thebim[idx,]
# V1              V2 V3       V4 V5 V6
# 94766 12 12:48580710_T_G  0 48186927  T  G
thetraw=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/chr12_filter_hg19tohg38_flip.traw"))
tmp=unlist(thetraw[idx,7:ncol(thetraw)])
(sum(tmp==0)*2+sum(tmp==1))/length(tmp)/2 #0.86 of G
theinfo=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc/chr12.info"))
idx=which(grepl("12:48580710",theinfo$SNP))
theinfo[idx,]
#          SNP            REF(0) ALT(1) ALT_Frq   MAF AvgCall     Rsq  Genotyped LooRsq EmpR EmpRsq Dose0 Dose1
# 662540 12:48580710:T:G      T      G 0.85761 0.14239 0.99754 0.98023   Imputed      -    -      -     -     -
cmd=paste0("echo '12:48580710-48580710' >../result/text.txt")
system(cmd)
# tabix  /fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc/chr12.dose.vcf.gz "12:48580710-48580710" > ../result/test.txt
test=as.data.frame(fread("../result/test.txt",sep="\t"))
test[1,1:9]
# V1       V2              V3 V4 V5 V6   V7                                        V8           V9
# 1 12 48580710 12:48580710:T:G  T  G  . PASS AF=0.85761;MAF=0.14239;R2=0.98023;IMPUTED GT:DS:HDS:GP
idx1=which(grepl("0|1",test[1,10:ncol(test)],fixed = T)) #VCF dose count on ALT
idx2=which(grepl("1|0",test[1,10:ncol(test)],fixed = T))
idx0=which(grepl("0|0",test[1,10:ncol(test)],fixed = T))
idx3=which(grepl("1|1",test[1,10:ncol(test)],fixed = T))
(length(idx1)+length(idx2)+length(idx0)*2)/(length(idx0)+length(idx1)+length(idx2)+length(idx3))/2 # 0.1400226 of T
res2$selsnpinfo[801,]
allbeta2[801,]

allbeta3=res3$allbeta
plot(allbeta3$freq_non_effect_allele_controls_BCA,allbeta3$freq_non_effect_allele_controls_Bonn,xlab="Frequency of non_effect_allele_control (BCA)",ylab="Frequency of non_effect_allele_control (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta3$freq_non_effect_allele_cases_BCA,allbeta3$freq_non_effect_allele_cases_Bonn,xlab="Frequency of non_effect_allele_case (BCA)",ylab="Frequency of non_effect_allele_case (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
# idx=which(allbeta3$Beta_BCA*allbeta3$Beta_Bonn<0)
# points(allbeta3$freq_non_effect_allele_controls_BCA[idx],allbeta3$freq_non_effect_allele_controls_Bonn[idx],col="blue")
# View(allbeta3[idx,])

plot(allbeta3$Beta_BCA,allbeta3$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta3$freq_non_effect_allele_cases_BCA-allbeta3$freq_non_effect_allele_controls_BCA,allbeta3$freq_non_effect_allele_cases_Bonn-allbeta3$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

plot(allbeta3$Beta_eQTL,allbeta3$Beta_BCA,xlab="Beta (eQTL)",ylab="Beta (BCA)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta3$Beta_eQTL,allbeta3$Beta_Bonn,xlab="Beta (eQTL)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)


load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
idx=order(gtexv8_ge_anno$Chromosome,gtexv8_ge_anno$start)
gtexv8_ge_anno=gtexv8_ge_anno[idx,]
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
idx=which(proteingenes=="ZNF641")
proteingenes[(idx-4):(idx+4)]
#"ASB8"     "CCDC184"  "OR10AD1"  "H1FNT"    "ZNF641"   "ANP32D"   "C12orf54" "OR8S1"    "LALBA"  
genes=c("ASB8","CCDC184","OR10AD1","H1FNT","ANP32D","C12orf54","OR8S1","LALBA")
organ="junction"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
sum(genes %in% rownames(res_min))
genes[genes %in% rownames(res_min)]
# "ASB8"     "CCDC184"  "OR10AD1"  "H1FNT"    "C12orf54" "OR8S1"
genes=genes[genes %in% rownames(res_min)]
organs=rep("junction",length(genes))
traits=rep("BEEA",length(genes))
res1_=colcal1(1)
res2_=colcal1(2)
res3_=colcal1(3)
res4_=colcal1(4)
res5_=colcal1(5)
res6_=colcal1(6)
save(res1_,res2_,res3_,res4_,res5_,res6_,file="../result/colocalisationRES2.RData")

allbeta=res4_$allbeta
plot(allbeta$freq_non_effect_allele_controls_BCA,allbeta$freq_non_effect_allele_controls_Bonn,xlab="Frequency of non_effect_allele_control (BCA)",ylab="Frequency of non_effect_allele_control (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta$freq_non_effect_allele_cases_BCA,allbeta$freq_non_effect_allele_cases_Bonn,xlab="Frequency of non_effect_allele_case (BCA)",ylab="Frequency of non_effect_allele_case (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)


plot(allbeta$Beta_BCA,allbeta$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta$freq_non_effect_allele_cases_BCA-allbeta$freq_non_effect_allele_controls_BCA,allbeta$freq_non_effect_allele_cases_Bonn-allbeta$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

allbeta=res5_$allbeta
plot(allbeta$Beta_BCA,allbeta$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta$freq_non_effect_allele_cases_BCA-allbeta$freq_non_effect_allele_controls_BCA,allbeta$freq_non_effect_allele_cases_Bonn-allbeta$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)

#check beta between BCA and Bonn
organ="junction"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1))
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2))
check_beta=function(genes,traits)
{
  oldchr=0
  allbeta0=NULL
  for (i in 1:length(genes))
  {
    gene=genes[i]
    trait=traits[i]
    #organ=organs[i]
    print(gene)
    
    # #first work on GTEx
    # outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
    # load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
    # load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
    # gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1))
    # gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2))
    idx=which(rownames(phenotype)==gene)
    tmp=distance(gr_snp,gr_pos[idx])
    idx=which(tmp<5e5)
    selsnps=rownames(snp)[idx]
    
    # GTEx_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
    # idx2=which(rownames(phenotype)==gene)
    # Y=unlist(phenotype[idx2,])
    # for (j in 1:length(selsnps))
    # {
    #   idx1=which(rownames(snp)==selsnps[j])
    #   dat=cbind.data.frame(snp=t(snp[idx1,]),covariate)
    #   colnames(dat)[1]="snp"
    #   
    #   fit=lm(Y~.,data=dat)
    #   thecoeff=summary(fit)$coefficients
    #   GTEx_beta$beta[j]=thecoeff[2,1]
    #   GTEx_beta$varbeta[j]=thecoeff[2,2]^2 #variance
    #   GTEx_beta$pvalue[j]=thecoeff[2,4]
    # }
    
    tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
    idx=which(is.na(tmp$snp))
    if (length(idx)>0) #use biomart, slower method
    {
      pos=allele1=allele2=rep(NA,length(selsnps))
      for (j in idx)
      {
        #cat(j,"\n")
        tmp1=unlist(strsplit(tmp$selsnps[j],":"))
        chr=tmp1[1]
        tmp1=unlist(strsplit(tmp1[2],"_"))
        allalleles=tmp1[2:3]
        allele1[j]=tmp1[2] #minor allele
        allele2[j]=tmp1[3]
        pos[j]=as.numeric(tmp1[1])
        #find rsid based on position and allles
        #tmp=listAttributes(snpmart)
        #sometimes biomart has connection problems
        tmp2=tryCatch(
          {
            getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
          },
          error=function(e)
          {
            return(F)
          }
        )
        
        #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
        if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
        if (class(tmp2)[1]=="data.frame")
        {
          for (k in 1:nrow(tmp2))
          {
            myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
            if (all(allalleles %in% myalleles))
            {
              tmp$snp[j]=tmp2$refsnp_id[k]
              break
            }
          }
        }
      }
    }
    selsnpinfo=tmp
    selsnpinfo$hg19_position=NA
    selsnpinfo$snphg19_=NA
    for (j in 1:nrow(selsnpinfo))
    {
      tmp=unlist(strsplit(selsnpinfo$snphg19[j],"_"))
      selsnpinfo$hg19_position[j]=as.numeric(unlist(strsplit(tmp[1],":"))[2])
      tmp=unlist(strsplit(selsnpinfo$snphg19[j],"_"))
      selsnpinfo$snphg19_[j]=paste(tmp[c(1,3,2)],collapse = "_")
    }
    # GTExdat=list()
    # GTExdat$beta=GTEx_beta$beta
    # names(GTExdat$beta)=selsnpinfo$snp
    # GTExdat$varbeta=GTEx_beta$varbeta
    # names(GTExdat$varbeta)=selsnpinfo$snp
    # GTExdat$N=length(Y)
    # GTExdat$sdY=sd(Y)
    # GTExdat$type="quant"
    # GTExdat$snp=selsnpinfo$snp
    # GTExdat$postion=selsnpinfo$hg19_position
    
    #BCA data
    cat("BCA..")
    chr=unlist(strsplit(selsnps[1],":"))[1]
    if (oldchr!=chr)
    {
      BCAsnp=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/merge_beacon_cambridge_hrc_maf005_snp/chr",chr,"_filter_hg19tohg38_flip.traw")))
      rownames(BCAsnp)=BCAsnp$SNP
      colnames(BCAsnp)=gsub("SEP","",colnames(BCAsnp))
      BCAsnp=BCAsnp[,7:ncol(BCAsnp)]
      oldchr=chr
    }
    
    if (!all(colnames(BCAsnp) %in% paste0(covariatetable$V1,"_",covariatetable$V2))& ncol(BCAsnp)!=nrow(covariatetable)) print("some samples are missing")
    mycovariatetable=covariatetable[match(colnames(BCAsnp),paste0(covariatetable$V1,"_",covariatetable$V2)),]
    tmp1=rownames(BCAsnp)[rownames(BCAsnp) %in% selsnpinfo$snphg19]
    tmp2=rownames(BCAsnp)[rownames(BCAsnp) %in% selsnpinfo$snphg19_]
    selsnpinfo$BCAflip=F
    selsnpinfo$BCAflip[selsnpinfo$snphg19_ %in% rownames(BCAsnp)]=T
    if (length(tmp1)+length(tmp2)<length(selsnps)) print("some snps are missing in BCA")
    if (length(tmp2)>0)
    {
      dat1=BCAsnp[match(tmp1,rownames(BCAsnp)),]
      dat2=BCAsnp[match(tmp2,rownames(BCAsnp)),]
      dat2=2-dat2
      idx=match(rownames(dat2),selsnpinfo$snphg19_)
      rownames(dat2)=selsnpinfo$snphg19[idx]
      dat3=rbind(dat1,dat2)
      dat3=dat3[match(selsnpinfo$snphg19,rownames(dat3)),]
    }else
    {
      dat3=BCAsnp[match(tmp1,rownames(BCAsnp)),]
    }
    
    type=traits[i]
    if (type=="BE")
    {
      Y_BCA=mycovariatetable$phenoBE_bca
      N_BCA=2401+881+6701+3408
    }
    if (type=="EA")
    {
      Y_BCA=mycovariatetable$phenoEA_bca
      N_BCA=1507+1003+6701+3408
    }
    if (type=="BEEA")
    {
      Y_BCA=mycovariatetable$phenoEABE_bca
      N_BCA=2401+881+1507+1003+6701+3408
    }
    BCA_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
    BCA_beta$freq_non_effect_allele_cases=BCA_beta$freq_non_effect_allele_controls=NA
    for (j in 1:nrow(dat3))
    {
      dat=data.frame(snp=dat3[j,])
      dat=cbind(t(dat),mycovariatetable[,c("pc1","pc2","pc3","pc4","pc5","pc6")])
      colnames(dat)[1]="snp"
      fit=glm(I(Y_BCA==2)~.,data=dat,family="binomial")
      tmp=summary(fit)$coefficients
      BCA_beta$beta[j]=tmp[2,1]
      BCA_beta$varbeta[j]=tmp[2,2]^2
      BCA_beta$pvalue[j]=tmp[2,4]
      tmp=unlist(dat3[j,which(Y_BCA==1)])
      BCA_beta$freq_non_effect_allele_controls[j]=(sum(tmp==0)*2+sum(tmp==1))/length(tmp)/2
      tmp=unlist(dat3[j,which(Y_BCA==2)])
      BCA_beta$freq_non_effect_allele_cases[j]=(sum(tmp==0)*2+sum(tmp==1))/length(tmp)/2
      
    }
    
    # BCAdat=list()
    # BCAdat$beta=BCA_beta$beta
    # names(BCAdat$beta)=selsnpinfo$snp
    # BCAdat$varbeta=BCA_beta$varbeta
    # names(BCAdat$varbeta)=selsnpinfo$snp
    # BCAdat$N=N_BCA
    # BCAdat$type="cc"
    # BCAdat$snp=selsnpinfo$snp
    # BCAdat$postion=selsnpinfo$hg19_position
    
    # #compute LD
    # chr=unlist(strsplit(selsnps[1],":"))[1]
    # #read data from 1000 genome v3 EUR
    # thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/" #to keep 1000 genome V3 data
    # 
    # 
    # if (oldchr1!=chr)
    # {
    #   refbim=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.bim"))
    #   #refsnps=intersect(rsid[!is.na(rsid)],refbim$V2)
    #   refraw=fread(paste0(thousanddir,"ALL.chr",chr,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.traw"))
    #   oldchr1=chr
    # }else
    # 
    # #refbimstr=paste0(refbim$V1,":",refbim$V4)
    # refbimstr1=paste0(refbim$V1,":",refbim$V4,"_",refbim$V5,"_",refbim$V6)
    # refbimstr2=paste0(refbim$V1,":",refbim$V4,"_",refbim$V6,"_",refbim$V5)
    # 
    # Z1=Z2=NULL
    # tmp1=intersect(selsnpinfo$snphg19,refbimstr1)
    # tmp2=intersect(selsnpinfo$snphg19,refbimstr2)
    # idx1=match(tmp1,refbimstr1)
    # if (length(tmp1)>0)
    # {
    #   Z1=as.data.frame(refraw[idx1,7:ncol(refraw),drop=F])
    #   rownames(Z1)=tmp1
    # }
    # idx2=match(tmp2,refbimstr2)
    # if (length(tmp2)>0)
    # {
    #   Z2=as.data.frame(refraw[idx2,7:ncol(refraw),drop=F])
    #   #change id
    #   rownames(Z2)=tmp2
    #   Z2=2-Z2
    # }
    # Z3=as.data.frame(t(rbind(Z1,Z2)))
    # tmp=selsnpinfo$snphg19[selsnpinfo$snphg19 %in% colnames(Z3)]
    # if (length(tmp)<length(selsnps)) print(paste0(gene,":",length(selsnps)-length(tmp)," snps are missing in 1000genome control"))
    # Z3=Z3[,match(tmp,colnames(Z3))]
    # tmp=rowMeans(Z3)
    # Z3=Z3[!is.na(tmp),]
    # LD=cor(Z3)
    # idx=match(colnames(Z3),selsnpinfo$snphg19)
    # rownames(LD)=colnames(LD)=selsnpinfo$snp[idx]
    # 
    # tmp=selsnpinfo$snp[idx]
    # idx=match(tmp,BCAdat$snp)
    # BCAdat$beta=BCAdat$beta[idx]
    # BCAdat$varbeta=BCAdat$varbeta[idx]
    # BCAdat$snp=BCAdat$snp[idx]
    # BCAdat$postion=BCAdat$postion[idx]
    # BCAdat$LD=LD
    # 
    # idx=match(tmp,GTExdat$snp)
    # GTExdat$beta=GTExdat$beta[idx]
    # GTExdat$varbeta=GTExdat$varbeta[idx]
    # GTExdat$snp=GTExdat$snp[idx]
    # GTExdat$postion=GTExdat$postion[idx]
    # GTExdat$LD=LD
    # 
    # check_dataset(GTExdat)
    # check_dataset(BCAdat)
    # my.res <- coloc.abf(dataset1=GTExdat, dataset2=BCAdat)
    # print(my.res)
    # pdf(paste0("../result/coloc1_",gene,"_abf.pdf"),width=12)
    # sensitivity(my.res,"H4 > 0.8")
    # dev.off()
    # 
    # S3=runsusie(GTExdat)
    # S4=runsusie(BCAdat)
    # susie.res=NULL
    # if(requireNamespace("susieR",quietly=TRUE)) {
    #   susie.res=coloc.susie(S3,S4)
    #   print(susie.res$summary)
    # }
    # if (!is.null(susie.res$summary))
    # {
    #   if(requireNamespace("susieR",quietly=TRUE) & nrow(susie.res$summary)>1) {
    #     pdf(paste0("../result/coloc1_",gene,"_susie.pdf"),width=12)
    #     sensitivity(susie.res,"H4 > 0.8",row=1,dataset1=GTExdat,dataset2=Bonndat)
    #     sensitivity(susie.res,"H4 > 0.8",row=2,dataset1=GTExdat,dataset2=Bonndat)
    #     dev.off()
    #   }
    # }
    # 
    
    #validation Bonn data
    cat("Bonn..")
    type=traits[i]
    Bonn_N=NA
    if (type=="BE")
    {
      sumdat=BE_Bonnsummarydat
      Bonn_N=1037+3537
    }
    if (type=="EA")
    {
      sumdat=EA_Bonnsummarydat
      Bonn_N=1609+3537
    }
    if (type=="BEEA")
    {
      sumdat=BEEA_Bonnsummarydat
      Bonn_N=1037+1609+3537
    }
    idx1=which(selsnpinfo$snphg19 %in% sumdat$snpname)
    idx2=which(selsnpinfo$snphg19 %in% sumdat$snpname1)
    selsnpinfo$Bonnflip=F
    selsnpinfo$Bonnflip[selsnpinfo$snphg19 %in% sumdat$snpname1]=T
    Bonn_beta=data.frame(snp=selsnps,beta=NA,varbeta=NA,pvalue=NA)
    Bonn_beta$freq_non_effect_allele_cases=Bonn_beta$freq_non_effect_allele_controls=NA
    for (j in 1:length(selsnps))
    {
      idx1=which(sumdat$snpname==selsnpinfo$snphg19[j])
      if (length(idx1)>0)
      {
        Bonn_beta$beta[j]=sumdat$BETA[idx1]
        Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
        Bonn_beta$pvalue[j]=sumdat$P[idx1]
        Bonn_beta$freq_non_effect_allele_cases[j]=sumdat$freq_non_effect_allele_cases[idx1]
        Bonn_beta$freq_non_effect_allele_controls[j]=sumdat$freq_non_effect_allele_controls[idx1]
      }
      idx1=which(sumdat$snpname1==selsnpinfo$snphg19[j])
      if (length(idx1)>0)
      {
        Bonn_beta$beta[j]=-sumdat$BETA[idx1]
        Bonn_beta$varbeta[j]=sumdat$SE[idx1]^2
        Bonn_beta$pvalue[j]=sumdat$P[idx1]
        Bonn_beta$freq_non_effect_allele_cases[j]=1-sumdat$freq_non_effect_allele_cases[idx1]
        Bonn_beta$freq_non_effect_allele_controls[j]=1-sumdat$freq_non_effect_allele_controls[idx1]
      }
    }
    
    # Bonndat=list()
    # Bonndat$beta=Bonn_beta$beta
    # names(Bonndat$beta)=selsnpinfo$snp
    # Bonndat$varbeta=Bonn_beta$varbeta
    # names(Bonndat$varbeta)=selsnpinfo$snp
    # Bonndat$N=Bonn_N
    # Bonndat$type="cc"
    # Bonndat$snp=selsnpinfo$snp
    # Bonndat$postion=selsnpinfo$hg19_position
    # 
    
    allbeta=data.frame(snp=selsnpinfo$snp,chr=chr,poistion=selsnpinfo$hg19_position,effect_allele=NA,non_effect_allel=NA,
                       Beta_BCA=BCA_beta$beta,Pvalue_BCA=BCA_beta$pvalue,Beta_Bonn=Bonn_beta$beta,Pvalue_Bonn=Bonn_beta$pvalue,selected=F,
                       freq_non_effect_allele_cases_Bonn=Bonn_beta$freq_non_effect_allele_cases,freq_non_effect_allele_controls_Bonn=Bonn_beta$freq_non_effect_allele_controls,
                       freq_non_effect_allele_cases_BCA=BCA_beta$freq_non_effect_allele_cases,freq_non_effect_allele_controls_BCA=BCA_beta$freq_non_effect_allele_controls,
                       BCAflip=selsnpinfo$BCAflip,Bonnflip=selsnpinfo$Bonnflip)
    for (j in 1:nrow(allbeta))
    {
      tmp=unlist(strsplit(selsnpinfo$selsnps[j],"_"))
      allbeta$effect_allele[j]=tmp[2]
      allbeta$non_effect_allele[j]=tmp[3]
    }
    idx=which(rownames(res_min)==gene)
    selsnps1=unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T))
    idx=match(selsnps1,selsnps)
    allbeta$selected[idx]=T
    allbeta$gene=gene
    allbeta$trait=trait
    allbeta0=rbind(allbeta0,allbeta)
  }

  return(allbeta0)
}
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
idx=order(gtexv8_ge_anno$Chromosome,gtexv8_ge_anno$start)
gtexv8_ge_anno=gtexv8_ge_anno[idx,]
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

idx=which(proteingenes=="ZNF641")
proteingenes[(idx-2):(idx+2)]
#"ASB8"     "CCDC184"  "OR10AD1"  "H1FNT"    "ZNF641"   "ANP32D"   "C12orf54" "OR8S1"    "LALBA"  
genes=c("ASB8","CCDC184","OR10AD1","H1FNT","ZNF641","ANP32D","C12orf54","OR8S1","LALBA")
genes=c("ZNF641")
sum(genes %in% rownames(res_min))
genes[genes %in% rownames(res_min)]
# "ASB8"     "CCDC184"  "OR10AD1"  "H1FNT"    "C12orf54" "OR8S1"

genes=genes[genes %in% rownames(res_min)]
traits=rep("BEEA",length(genes))
allbeta_chr12_2=NULL
tmp=check_beta(genes=genes,traits=traits)
allbeta_chr12_2=rbind(allbeta_chr12_2,tmp)
save(allbeta_chr12_2,file="../result/allbeta_chr12_2.RData")
traits=rep("BE",length(genes))
tmp=check_beta(genes=genes,traits=traits)
allbeta_chr12_2=rbind(allbeta_chr12_2,tmp)
save(allbeta_chr12_2,file="../result/allbeta_chr12_2.RData")
traits=rep("EA",length(genes))
tmp=check_beta(genes=genes,traits=traits)
allbeta_chr12_2=rbind(allbeta_chr12_2,tmp)
save(allbeta_chr12_2,file="../result/allbeta_chr12_2.RData")

# allbeta_chr12_3=NULL
# traits=rep("EA",length(genes))
# tmp=check_beta(genes=genes,traits=traits)
# allbeta_chr12_3=rbind(allbeta_chr12_3,tmp)
# save(allbeta_chr12_3,file="../result/allbeta_chr12_3.RData")
# idx=which(allbeta_chr12_3$gene=="ZNF641")
# allbeta=allbeta_chr12_3[idx,]
# plot(allbeta$Beta_BCA,allbeta$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
# abline(0,1)
# plot(allbeta$freq_non_effect_allele_cases_BCA-allbeta$freq_non_effect_allele_controls_BCA,allbeta$freq_non_effect_allele_cases_Bonn-allbeta$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
# abline(0,1)

load("../result/allbeta_chr12_2.RData")

#only keep snps without flipping allels, seems the pattern is different in BCA and Bonn
idx=which(allbeta_chr12_2$trait=="BEEA" & !allbeta_chr12_2$BCAflip & !allbeta_chr12_2$Bonnflip)
allbeta=allbeta_chr12_2[idx,]
plot(allbeta$Beta_BCA,allbeta$Beta_Bonn,xlab="Beta (BCA)",ylab="Beta (Bonn)",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
plot(allbeta$freq_non_effect_allele_cases_BCA-allbeta$freq_non_effect_allele_controls_BCA,allbeta$freq_non_effect_allele_cases_Bonn-allbeta$freq_non_effect_allele_controls_Bonn,xlab="Frequency non_effect_allel diff (case-control) BCA",ylab="Frequency non_effect_allel diff (case-control) Bonn",cex.lab=1.3,cex.axis=1.3)
abline(0,1)
