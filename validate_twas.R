#!usr/bin/env Rscript
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

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

# allsnps=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt"))
# allsnps=allsnps[allsnps$`#chrom` %in% paste0("chr",c(1:22,"X","Y")),]
# allsnps$`#chrom`=gsub("chr","",allsnps$`#chrom`)

library("biomaRt")
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")

validate_snp=function(fdrres=fdr_fwer_res$BE_fdr,pres=fdr_fwer_res$BE_fdr_p,fwerres=fdr_fwer_res$BE_fdr_fwer,summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt"))
{
  allres=NULL
  
  if(length(fdrres)>0)
  {
    summarydat=as.data.frame(fread(summaryfile,header=T))
    tmp=intersect("SNP",colnames(summarydat))
    if (length(tmp)==0)
    {
      colnames(summarydat)[which(colnames(summarydat)=="rsid")]="SNP"
      colnames(summarydat)[which(colnames(summarydat)=="pvalue")]="P"
    }
    for (i in 1:length(fdrres))
    {
      idx=which(rownames(res_min) %in% names(fdrres[i]))
      selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
      tmp=unlist(strsplit(selsnps[1],":"))
      chr=tmp[1]
      pos=rsid=rep(NA,length(selsnps))
      for (j in 1:length(selsnps))
      {
        tmp=unlist(strsplit(selsnps[j],":"))
        tmp=unlist(strsplit(tmp[2],"_"))
        allallels=tmp[2:3]
        pos[j]=as.numeric(tmp[1])
        #find rsid based on position and allles
        tmp=getBM(attributes=c('refsnp_id','allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart)
        if (nrow(tmp)>0)
        {
          if (nrow(tmp)==1)
          {
            rsid[j]=tmp$refsnp_id
          }else
          {
            tmp=tmp[tmp$chrom_end-tmp$chrom_start==0,]
            if (nrow(tmp)>0)
            {
              if (nrow(tmp)==1)
              {
                rsid[j]=tmp$refsnp_id
              }else
              {
                #remove del
                tmp=tmp[!grepl("-",tmp$allele),]
                if (nrow(tmp)>0)
                {
                  if (nrow(tmp)==1)
                  {
                    rsid[j]=tmp$refsnp_id
                  }else
                  {
                     #for multivariates SNPS
                      for (k in 1:nrow(tmp))
                      {
                        alllels=unlist(strsplit(tmp$allele[k],"/",fixed=T))
                        if (all(allallels %in% alllels))
                        {
                          rsid[j]=tmp$refsnp_id[k]
                          break
                        }
                      }
                  }
                }
              }
            }
          }
        }
      }
      
      res=data.frame(snp=selsnps,rsid=rsid,chr=chr,pos=pos,gene=names(fdrres[i]),twas_p=pres[i],twas_fdr=fdrres[i],twas_fwer=fwerres[i],gwas_p=NA,gwas_info=NA)
      tmp=intersect(rsid[!is.na(rsid)],summarydat$SNP)
      if (length(tmp)>0)
      {
        idx1=match(tmp,res$rsid)
        idx2=match(tmp,summarydat$SNP)
        res$gwas_p[idx1]=summarydat$P[idx2]
        #res$gwas_or[idx1]=summarydat$OR[idx2]
        res$gwas_info[idx1]=summarydat$info[idx2]
      }
      
      allres=rbind(allres,res)
    }
    
  }
  return(allres)
}

#TCGA derived TWAS
#load TWAS result (tryscat1.R):
prefix="dist500K_TCGA_PC4"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/skat_res.RData"))
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
fdr_fwer_res=compute_fwer_fdr()

BE_Bonn_res=validate_snp()
BE_oxford_res=validate_snp(summaryfile=paste0(summaryfolder,"BE_oxford_autosomes.txt"))
BE_oxford_res[which(BE_oxford_res$gwas_p<0.05),]
# snp        rsid chr      pos  gene       twas_p   twas_fdr  twas_fwer     gwas_p gwas_info
# DDX49 19:18817903_A_C  rs10423674  19 18817903 DDX49 1.781413e-06 0.01750238 0.01750238 0.00203924  1.000000
# 10    19:18962186_A_G rs114786762  19 18962186 MEF2B 3.595918e-06 0.01766495 0.03532989 0.01621030  0.904966
# 18    19:18829770_G_A  rs11668587  19 18829770 MEF2B 3.595918e-06 0.01766495 0.03532989 0.00226376  0.993737
# 22    19:18826975_T_C  rs10221489  19 18826975 MEF2B 3.595918e-06 0.01766495 0.03532989 0.00260866  0.994527
# 37    19:18944025_G_A   rs2239371  19 18944025 MEF2B 3.595918e-06 0.01766495 0.03532989 0.03746450  1.000000

EA_Bonn_res=validate_snp(fdrres=fdr_fwer_res$EA_fdr,pres=fdr_fwer_res$EA_fdr_p,fwerres=fdr_fwer_res$EA_fdr_fwer,summaryfile=paste0(summaryfolder,"EA_Bonn_autosomes.txt"))
nrow(EA_Bonn_res) #40
EA_Bonn_res[which(EA_Bonn_res$gwas_p<0.05),]
# snp      rsid chr      pos  gene       twas_p   twas_fdr  twas_fwer    gwas_p gwas_info
# 14 19:18801813_A_G rs2385178  19 18801813 MEF2B 4.171275e-06 0.02533991 0.04099112 0.0455728  0.990243
# 35 19:18835115_A_G rs2051815  19 18835115 MEF2B 4.171275e-06 0.02533991 0.04099112 0.0290132  0.947020

BEEA_Bonn_res=validate_snp(fdrres=fdr_fwer_res$BEEA_fdr,pres=fdr_fwer_res$BEEA_fdr_p,fwerres=fdr_fwer_res$BEEA_fdr_fwer,summaryfile=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"))
BEEA_Bonn_res[which(BEEA_Bonn_res$gwas_p<0.05),]
# snp       rsid chr      pos  gene      twas_p     twas_fdr  twas_fwer    gwas_p gwas_info
# 12 19:19240470_T_A rs75031662  19 19240470 MEF2B 1.70529e-07 0.0008382352 0.00167647 0.0341526  0.569197
# 14 19:18801813_A_G  rs2385178  19 18801813 MEF2B 1.70529e-07 0.0008382352 0.00167647 0.0429329  0.990419

#GTEx derived TWAS
prefix="dist500K_GTEx_PC4"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
load(paste0(outfolder,"/skat_res.RData"))
skat_min_code=skat_min[rownames(skat_min) %in% proteingenes,]
fdr_fwer_res=compute_fwer_fdr()

GTEx_BE_Bonn_res=validate_snp()
nrow(GTEx_BE_Bonn_res) #1060
GTEx_BE_Bonn_res[which(GTEx_BE_Bonn_res$gwas_p<0.05),]
nrow(GTEx_BE_Bonn_res[which(GTEx_BE_Bonn_res$gwas_p<0.05),]) #98

GTEx_BE_oxford_res=validate_snp(summaryfile=paste0(summaryfolder,"BE_oxford_autosomes.txt"))
nrow(GTEx_BE_oxford_res) #1060
GTEx_BE_oxford_res[which(GTEx_BE_oxford_res$gwas_p<0.05),]
nrow(GTEx_BE_oxford_res[which(GTEx_BE_oxford_res$gwas_p<0.05),]) #127

GTEx_EA_Bonn_res=validate_snp(fdrres=fdr_fwer_res$EA_fdr,pres=fdr_fwer_res$EA_fdr_p,fwerres=fdr_fwer_res$EA_fdr_fwer,summaryfile=paste0(summaryfolder,"EA_Bonn_autosomes.txt"))
nrow(GTEx_EA_Bonn_res) #25
GTEx_EA_Bonn_res[which(GTEx_EA_Bonn_res$gwas_p<0.05),]
# snp       rsid chr      pos gene       twas_p     twas_fdr    twas_fwer     gwas_p gwas_info
# 2  19:18688985_C_T  rs4808821  19 18688985 KXD1 4.290802e-08 0.0004479597 0.0004479597 0.04620170  0.734253
# 5  19:18669987_A_G  rs3803916  19 18669987 KXD1 4.290802e-08 0.0004479597 0.0004479597 0.01289870  0.784282
# 7  19:18796177_A_T  rs8102768  19 18796177 KXD1 4.290802e-08 0.0004479597 0.0004479597 0.02076700  0.973371
# 9  19:18669987_A_G  rs3803916  19 18669987 COMP 6.851423e-06 0.0357644268 0.0715288536 0.01289870  0.784282
# 23 19:18886415_A_G rs56033554  19 18886415 COMP 6.851423e-06 0.0357644268 0.0715288536 0.00190007  0.730930

GTEx_BEEA_Bonn_res=validate_snp(fdrres=fdr_fwer_res$BEEA_fdr,pres=fdr_fwer_res$BEEA_fdr_p,fwerres=fdr_fwer_res$BEEA_fdr_fwer,summaryfile=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt"))
nrow(GTEx_BEEA_Bonn_res) #1333
GTEx_BEEA_Bonn_res[which(GTEx_BEEA_Bonn_res$gwas_p<0.05),]
nrow(GTEx_BEEA_Bonn_res[which(GTEx_BEEA_Bonn_res$gwas_p<0.05),]) #171

save(BE_Bonn_res,BE_oxford_res,EA_Bonn_res,BEEA_Bonn_res,GTEx_BE_Bonn_res,GTEx_BE_oxford_res,GTEx_EA_Bonn_res,GTEx_BEEA_Bonn_res,
     file="../result/validate_twas.RData")
