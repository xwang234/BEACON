#!usr/bin/env Rscript
#check the results based on 1se rule

library(data.table)
compute_fwer_fdr=function(dat=skat_min_code,cutoff=0.05)
{
  compute_fdr=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method="fdr")
    if (sum(tmp<cutoff1,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff)]
      res=res[order(res)]
    }
    return(res)
  }
  compute_fdr_addpvalue=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method="fdr")
    if (sum(tmp<cutoff1,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff)]
      res=res[order(res)]
      tmp1=names(res)
      idx=match(tmp1,names(pvalues))
      res=pvalues[idx]
    }
    return(res)
  }
  compute_fdr_addfwer=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method="fdr")
    if (sum(tmp<cutoff1,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff)]
      res=res[order(res)]
      tmp1=names(res)
      tmp=p.adjust(pvalues,method="bonferroni")
      idx=match(tmp1,names(tmp))
      res=tmp[idx]
    }
    return(res)
  }
  compute_fwer=function(dat,namecol="BE_p",cutoff1=cutoff)
  {
    res=NULL
    pvalues=dat[,which(colnames(dat)==namecol)]
    names(pvalues)=rownames(dat)
    tmp=p.adjust(pvalues,method = "bonferroni")
    if (sum(tmp<cutoff,na.rm = T)>0)
    {
      res=tmp[which(tmp<cutoff1)]
      res=res[order(res)]
    }
    return(res)
  }
  BE_fdr=compute_fdr(dat,namecol="BE_p")
  BE_fdr_p=compute_fdr_addpvalue(dat,namecol="BE_p")
  BE_fdr_fwer=compute_fdr_addfwer(dat,namecol="BE_p")
  BE_fwer=compute_fwer(dat,namecol="BE_p")
  EA_fdr=compute_fdr(dat,namecol="EA_p")
  EA_fdr_p=compute_fdr_addpvalue(dat,namecol="EA_p")
  EA_fdr_fwer=compute_fdr_addfwer(dat,namecol="EA_p")
  EA_fwer=compute_fwer(dat,namecol="EA_p")
  BEA_fdr=compute_fdr(dat,namecol="BEA_p")
  BEA_fwer=compute_fwer(dat,namecol="BEA_p")
  BEEA_fdr=compute_fdr(dat,namecol="BEEA_p")
  BEEA_fdr_p=compute_fdr_addpvalue(dat,namecol="BEEA_p")
  BEEA_fdr_fwer=compute_fdr_addfwer(dat,namecol="BEEA_p")
  BEEA_fwer=compute_fwer(dat,namecol="BEEA_p")
  return(list(BE_fdr=BE_fdr,BE_fdr_p=BE_fdr_p,BE_fdr_fwer=BE_fdr_fwer,BE_fwer=BE_fwer,EA_fdr=EA_fdr,EA_fdr_p=EA_fdr_p,EA_fdr_fwer=EA_fdr_fwer,EA_fwer=EA_fwer,
              BEA_fdr=BEA_fdr,BEEA_fdr_p=BEEA_fdr_p,BEA_fwer=BEA_fwer,BEEA_fdr=BEEA_fdr,BEEA_fdr_fwer=BEEA_fdr_fwer,BEEA_fwer=BEEA_fwer))
}

count_validated_gene=function(dat=GTEx_BE_oxford_res1)
{
  res=data.frame(gene=unique(dat$gene),nsnp=NA,nvalidatedsnp=NA,propvalidatedsnp=NA,stringsAsFactors = F)
  for (i in 1:nrow(res))
  {
    res$nsnp[i]=sum(dat$gene==res$gene[i])
    res$nvalidatedsnp[i]=sum(dat$gene==res$gene[i]& dat$validated==T)
    res$propvalidatedsnp[i]=res$nvalidatedsnp[i]/res$nsnp[i]
  }
  return(res)
}

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

# allsnps=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt"))
# allsnps=allsnps[allsnps$`#chrom` %in% paste0("chr",c(1:22,"X","Y")),]
# allsnps$`#chrom`=gsub("chr","",allsnps$`#chrom`)

library("biomaRt")
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")

validate_snp=function(model=res_min,fdrres=fdr_fwer_res$BE_fdr,pres=fdr_fwer_res$BE_fdr_p,fwerres=fdr_fwer_res$BE_fdr_fwer,
                      summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt"),
                      beaconfile=paste0(beaconfolder,"imp_TCGA_BE_CO.assoc.logistic"))
{
  allres=NULL
  
  if(length(fdrres)>0)
  {
    #summary stat from validation
    summarydat=as.data.frame(fread(summaryfile,header=T))
    tmp=intersect("SNP",colnames(summarydat))
    if (length(tmp)==0)
    {
      colnames(summarydat)[which(colnames(summarydat)=="rsid")]="SNP"
      colnames(summarydat)[which(colnames(summarydat)=="pvalue")]="P"
      colnames(summarydat)[which(colnames(summarydat)=="beta")]="BETA"
      colnames(summarydat)[which(colnames(summarydat)=="effect-allele")]="effect_allele"
      colnames(summarydat)[which(colnames(summarydat)=="non-effect-allele")]="non_effect_allele"
    }
    #summary from beacon
    beacondat=as.data.frame(fread(beaconfile,header=T))
    #for each gene
    for (i in 1:length(fdrres))
    {
      idx=which(rownames(model) %in% names(fdrres[i]))
      selsnps=unique(unlist(strsplit(model$selectedsnps[idx],"|",fixed=T)))
      tmp=unlist(strsplit(selsnps[1],":"))
      chr=tmp[1]
      pos=rsid=allele1=allele2=rep(NA,length(selsnps))
      #selsnps="19:19259262_C_A"
      for (j in 1:length(selsnps))
      {
        tmp=unlist(strsplit(selsnps[j],":"))
        tmp=unlist(strsplit(tmp[2],"_"))
        allalleles=tmp[2:3]
        allele1[j]=tmp[2] #minor allele
        allele2[j]=tmp[3]
        pos[j]=as.numeric(tmp[1])
        #find rsid based on position and allles
        #tmp=listAttributes(snpmart)
        tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart)
        if (nrow(tmp)>0)
        {
          for (k in 1:nrow(tmp))
          {
            myalleles=unlist(strsplit(tmp$allele[k],"/",fixed=T))
            if (all(allalleles %in% myalleles))
            {
              rsid[j]=tmp$refsnp_id[k]
              break
            }
          }
        }
      }
      n=length(selsnps)
      res=data.frame(snp=selsnps,rsid=rsid,chr=rep(chr,n),pos=pos,gene=rep(names(fdrres[i]),n),twas_allele1=allele1,twas_allele2=allele2,
                     beacon_beta=NA,twas_p=rep(pres[i],n),twas_fdr=rep(fdrres[i],n),twas_fwer=rep(fwerres[i],n),
                     validation_p=NA,validation_beta=NA,validation_effective=NA,validation_noneffective=NA,validation_info=NA)
      tmp=intersect(rsid[!is.na(rsid)],summarydat$SNP)
      if (length(tmp)>0)
      {
        idx1=match(tmp,res$rsid)
        idx2=match(tmp,summarydat$SNP)
        snpid=as.character(res$snp[match(tmp,res$rsid)])
        idx3=match(snpid,beacondat$SNP)
        idx4=which(is.na(idx3)) #snp not found in discovery data,may be due to allele def
        if (length(idx4)>0)
        {
          for (i in idx4)
          {
            mysnpid=snpid[i]
            tmp1=unlist(strsplit(mysnpid,":"))
            tmp2=unlist(strsplit(tmp1[2],"_"))
            mysnpid1=paste0(tmp1[1],":",tmp2[1],"_",tmp2[3],"_",tmp2[2])
            idx5=which(beacondat$SNP==mysnpid1)
            if (length(idx5)>0)
            {
              beacondat$SNP[idx5]=mysnpid
              beacondat$OR[idx5]=exp(-log(beacondat$OR[idx5])) #change sign of beta
            }
          }
        }
        idx3=match(snpid,beacondat$SNP)
        res$beacon_beta[idx1]=log(beacondat$OR[idx3]) #beta=log(OR)
        res$validation_p[idx1]=summarydat$P[idx2]
        res$validation_beta[idx1]=summarydat$BETA[idx2]
        #res$validation_or[idx1]=summarydat$OR[idx2]
        res$validation_effective[idx1]=summarydat$effect_allele[idx2]
        res$validation_noneffective[idx1]=summarydat$non_effect_allele[idx2]
        res$validation_info[idx1]=summarydat$info[idx2]
      }
      res$validated=F
      idx=which(res$validation_p<0.05)
      if (length(idx)>0)
      {
        for (i in idx)
        {
          if (!is.na(res$beacon_beta[i]))
          {
            if (res$twas_allele1[i]==res$validation_effective[i] & res$validation_beta[i]*res$beacon_beta[i]>0)
            {
              res$validated[i]=T
            }else
            {
              if (res$twas_allele2[i]==res$validation_effective[i] & res$validation_beta[i]*res$beacon_beta[i]<0)
              {
                res$validated[i]=T
              }
            }
          }
        }
      }
      allres=rbind(allres,res)
    }
  }
  print(paste0("total num of snp: ",nrow(allres)))
  idx=which(is.na(allres$rsid))
  if (length(idx)>0) print(paste0(length(idx)," snp can't be validated, not found rsid"))
  idx=which(!is.na(allres$rsid) & is.na(allres$validation_p))
  if (length(idx)>0) print(paste0(length(idx)," snp can't be validated, not found in validation summary"))
  idx=which(allres$validation_p<0.05 & is.na(allres$beacon_beta))
  if (length(idx)>0) print(paste0(length(idx)," snp can't be validated, not found in beacon summary"))
  idx=which(allres$validation_p<0.05)
  if (length(idx)>0) print(paste0(length(idx)," snp been validated only based on p-value"))

  idx=which(allres$validated==T)
  if (length(idx)>0) print(paste0(length(idx)," snp been validated"))
  
  return(allres)
}

#GTEx derived TWAS
prefix="dist500K_GTEx_PC4"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/skat_res.RData"))
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
skat_1se_code=skat_1se[rownames(skat_1se) %in% proteingenes,]
fdr_fwer_res=compute_fwer_fdr(dat=skat_1se_code)
fdr_fwer_res=compute_fwer_fdr(dat=skat_1se)
beaconfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_GTExmodel/" #beacon association results for GTEx selected snps

GTEx_BE_Bonn_res1=validate_snp(model=res_1se,beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic"))
# [1] "total num of snp: 196"
# [1] "28 snp can't be validated, not found in validation summary"
# [1] "46 snp been validated only based on p-value"
# [1] "43 snp been validated"
length(intersect(GTEx_BE_Bonn_res$snp[GTEx_BE_Bonn_res$validated==T],
                 GTEx_BE_Bonn_res1$snp[GTEx_BE_Bonn_res1$validated==T])) #19 were found in min model
length(unique(GTEx_BE_Bonn_res1$gene)) #3
length(unique(GTEx_BE_Bonn_res$gene)) #10
sum(unique(GTEx_BE_Bonn_res1$gene) %in% unique(GTEx_BE_Bonn_res$gene)) #1 coded gene found in min model
GTEx_BE_oxford_res1=validate_snp(model=res_1se,summaryfile=paste0(summaryfolder,"BE_oxford_autosomes.txt"),
                                beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic"))
# [1] "total num of snp: 196"
# [1] "19 snp can't be validated, not found in validation summary"
# [1] "45 snp been validated only based on p-value"
# [1] "40 snp been validated"
sum(GTEx_BE_Bonn_res1$snp[GTEx_BE_Bonn_res1$validated==T] %in% GTEx_BE_oxford_res1$snp[GTEx_BE_oxford_res1$validated==T]) #18

GTEx_EA_Bonn_res1=validate_snp(model=res_1se,fdrres=fdr_fwer_res$EA_fdr,pres=fdr_fwer_res$EA_fdr_p,fwerres=fdr_fwer_res$EA_fdr_fwer,summaryfile=paste0(summaryfolder,"EA_Bonn_autosomes.txt"),
                              beaconfile=paste0(beaconfolder,"imp_GTEx_EA_CO.assoc.logistic"))

GTEx_BEEA_Bonn_res1=validate_snp(model=res_1se,fdrres=fdr_fwer_res$BEEA_fdr,pres=fdr_fwer_res$BEEA_fdr_p,fwerres=fdr_fwer_res$BEEA_fdr_fwer,summaryfile=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"),
                                beaconfile=paste0(beaconfolder,"imp_GTEx_BEEA_CO.assoc.logistic"))
# [1] "total num of snp: 145"
# [1] "13 snp can't be validated, not found in validation summary"
# [1] "21 snp been validated only based on p-value"
# [1] "18 snp been validated"


save(GTEx_BE_Bonn_res1,GTEx_BE_oxford_res1,GTEx_EA_Bonn_res1,GTEx_BEEA_Bonn_res1,
     file="../result/validate_twas2.RData")

#problems:
#didn't find rsid for some snps based on position. there were not included in the biomart database
#tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(6,29912030,29912030),mart = snpmart)
