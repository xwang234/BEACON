#!/usr/bin/env Rscript
library(data.table)
qqplot=function(pvalue=NULL,fwer=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (-log base 10)",
       ylab="Observed p-value (-log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
  pvalue_order=pvalue[order(pvalue)]
  if (is.null(fwer))
  {
    fwer=p.adjust(pvalue_order,method="bonferroni")
    idx=which(fwer<0.05)
    if (length(idx)>0)
    {
      points(-log((1:n)/n,base=10)[idx],-log(pvalue[order(pvalue)],base=10)[idx],pch=16,col="red")
      print(paste0("#FWER<0.05:",length(idx)))
      legend("topleft",legend = "FWER<0.05",pch=16,col="red")
    }
  }else
  {
    idx1=sum(fwer<0.05)
    if (length(idx1)>0)
    {
      idx=1:idx1
      points(-log((1:n)/n,base=10)[idx],-log(pvalue[order(pvalue)],base=10)[idx],pch=16,col="red")
      print(paste0("#FWER<0.05:",length(idx)))
      legend("topleft",legend = "FWER<0.05",pch=16,col="red")
    }
  }
}
#V7 version
lookuptable=fread("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")
eqtldir="/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/GTEx_Analysis_v7_eQTL/"

#Junction results---------------
jointeqtl=as.data.frame(fread(paste0(eqtldir,"Esophagus_Gastroesophageal_Junction.v7.signif_variant_gene_pairs.txt")))
save(jointeqtl,file="../result/GTEx_junction_V7_eQTL_gene_pairs.RData")
all(jointeqtl$variant_id %in% lookuptable$variant_id)
idx=match(jointeqtl$variant_id,lookuptable$variant_id)
jointeqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

jointeqtl$variant_id1=gsub("_b37","",jointeqtl$variant_id)
tmp=unlist(strsplit(jointeqtl$variant_id1,"_"))
jointeqtl$Chr=tmp[seq(1,length(tmp),4)]
jointeqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
jointeqtl$Alt=tmp[seq(3,length(tmp),4)]
jointeqtl$Ref=tmp[seq(4,length(tmp),4)]
jointeqtl$name=paste0(jointeqtl$Chr,"_",jointeqtl$Position,"_",jointeqtl$Alt,"_",jointeqtl$Ref)
jointeqtl$namepos=paste0(jointeqtl$Chr,"_",jointeqtl$Position)

library(qvalue)

jointalleqtl=as.data.frame(fread(paste0(eqtldir,"Esophagus_Gastroesophageal_Junction.allpairs.txt")))
jointalleqtl$variant_id1=gsub("_b37","",jointalleqtl$variant_id)
tmp=unlist(strsplit(jointalleqtl$variant_id1,"_"))
jointalleqtl$Chr=tmp[seq(1,length(tmp),4)]
jointalleqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
jointalleqtl$Alt=tmp[seq(3,length(tmp),4)]
jointalleqtl$Ref=tmp[seq(4,length(tmp),4)]
jointalleqtl$name=paste0(jointalleqtl$Chr,"_",jointalleqtl$Position,"_",jointalleqtl$Alt,"_",jointalleqtl$Ref)
jointalleqtl$namepos=paste0(jointalleqtl$Chr,"_",jointalleqtl$Position)
all(jointalleqtl$variant_id %in% lookuptable$variant_id)

idx=match(jointalleqtl$variant_id,lookuptable$variant_id)
jointalleqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

tmp=qvalue(jointalleqtl$pval_nominal)
jointalleqtl$qvalue=tmp$qvalues
idx=which(jointalleqtl$qvalue<0.05)
length(unique(jointalleqtl$variant_id[idx])) #733707
jointalleqtl$fwer=p.adjust(jointalleqtl$pval_nominal,method="bonferroni")
idx=which(jointalleqtl$fwer<0.05)
length(unique(jointalleqtl$variant_id[idx])) #155606
#qqplot(jointalleqtl$pval_nominal)

#check SNP from literatures
dong23snp=read.table("../data/Dong23snp.txt",header=T,stringsAsFactors = F,sep="\t")
colnames(dong23snp)[colnames(dong23snp)=="Effect.allele"]="Alt"
dong23snp$Ref=" "
puya14snp=read.table("../data/Puya14snp.txt",header=T,stringsAsFactors = F,sep=" ",fill=T)
colnames(puya14snp)[which(colnames(puya14snp)=="Tested_allele")]="Alt"
colnames(puya14snp)[which(colnames(puya14snp)=="Other_allele")]="Ref"
levine12snp=read.table("../data/Levine12snp.txt",header=T,stringsAsFactors = F)
check_snp_eqtl=function(dat=puya14snp,opt=3,gtextable=jointeqtl)
{
  if (opt==1) #use alleles, can't match SNPs on different strands (e.g. T/G to A/C)
  {
    if ("Ref" %in% colnames(dat))
    {
      str11=paste0(dat$Chr,"_",dat$Position,"_",dat$Alt,"_",dat$Ref)
      str12=paste0(dat$Chr,"_",dat$Position,"_",dat$Ref,"_",dat$Alt)
      #sum(gtextable$name %in% str11)
      #sum(gtextable$name %in% str12)
      tmp=c(str11[str11 %in% gtextable$name],
            str12[str12 %in% gtextable$name])
      print(paste0("#SNPS found:",length(tmp)))
      print(paste0(tmp))
    }
  }
  if (opt==2) #use position
  {
    dat$name=paste0(dat$Chr,"_",dat$Position,"_",dat$Alt,"_",dat$Ref)
    allpairs=NULL
    dat$GTEx_minp=dat$GTEx_maf=dat$GTEx_name=dat$GTEx_SNP=NA
    str1=paste0(dat$Chr,"_",dat$Position)
    tmp=str1[str1 %in% gtextable$namepos]
    print(paste0("#SNPS found:",length(tmp)))
    print(paste0(tmp))
    if (length(tmp)>0)
    {
      idx=gtextable$namepos %in% tmp
      allpairs=gtextable[idx,]
      dat$GTEx_minp=dat$GTEx_maf=dat$GTEx_name=dat$GTEx_SNP=NA
      idx=match(tmp,str1)
      for (i in idx)
      {
        idx1=which(gtextable$namepos==str1[i])
        dat$GTEx_minp[i]=min(gtextable$pval_nominal[idx1])
        dat$GTEx_maf[i]=gtextable$maf[idx1[1]]
        dat$GTEx_name[i]=gtextable$name[idx1[1]]
        dat$GTEx_SNP[i]=gtextable$SNP[idx1[1]]
      }
    }
    return(list(allpairs=allpairs,snptable=dat))
  }
  
  if (opt==3) #use SNPid
  {
    tmp=intersect(gtextable$SNP,dat$SNP)
    print(paste0("#SNPS found:",length(tmp)))
    print(paste0(tmp))
    allpairs=NULL
    dat$GTEx_minp=dat$GTEx_maf=NA
    if (length(tmp)>0)
    {
      idx=gtextable$SNP %in% tmp
      allpairs=gtextable[idx,]
      dat$GTEx_minp=dat$GTEx_maf=NA
      idx=match(tmp,dat$SNP)
      for (i in idx)
      {
        idx1=which(gtextable$SNP==dat$SNP[i])
        dat$GTEx_minp[i]=min(gtextable$pval_nominal[idx1])
        dat$GTEx_maf[i]=gtextable$maf[idx1[1]]
      }
    }
    return(list(allpairs=allpairs,snptable=dat))
  }
}

all(puya14snp$SNP %in% dong23snp$SNP) #T
sum(levine12snp$SNP %in% dong23snp$SNP) #3
levine12snp[levine12snp$SNP %in% dong23snp$SNP,]

puyares=check_snp_eqtl(gtextable = jointeqtl)
# [1] "#SNPS found:3"
# [1] "rs7255"     "rs9257809"  "rs10108511"
puyares1=check_snp_eqtl(gtextable = jointeqtl,opt=2)
#View(puyares$snptable)
#View(puyares$eqtl)
dongres=check_snp_eqtl(gtextable = jointeqtl,dat=dong23snp)
# [1] "#SNPS found:5"
# [1] "rs3072"     "rs7255"     "rs9257809"  "rs10108511" "rs10423674"
dongres1=check_snp_eqtl(gtextable =jointeqtl,dat=dong23snp,opt=2)
levineres=check_snp_eqtl(gtextable =jointeqtl,dat=levine12snp)
# [1] "#SNPS found:1"
# [1] "rs10423674"
levineres1=check_snp_eqtl(gtextable =jointeqtl,dat=levine12snp,opt=2)

puyaresall=check_snp_eqtl(gtextable = jointalleqtl)
# [1] "#SNPS found:10"
# [1] "rs7255"     "rs2687202"  "rs9257809"  "rs17451754" "rs17749155" "rs10108511" "rs7852462"  "rs1247942" 
# [9] "rs2464469"  "rs1979654" 
puyaresall1=check_snp_eqtl(gtextable = jointalleqtl,opt=2)
# [1] "#SNPS found:11"
# [1] "2_20878820"   "15_58362025"  "7_117256712"  "8_10068073"   "8_11435516"   "3_70929983"   "12_114673723"
# [8] "6_29356331"   "9_100310501"  "2_200045039"  "16_86396835"
which(puyaresall1$snptable$SNP!=puyaresall1$snptable$GTEx_SNP)
#12 this is an example one SNP has 2 rsid (in different dbsnp version)
puyaresall1$eqtl=puyares1$snptable$SNP[!is.na(puyares1$snptable$GTEx_minp)]
puyaresall1$eqtlpairs=puyares1$allpairs

dongresall=check_snp_eqtl(dat=dong23snp,gtextable = jointalleqtl)
# [1] "#SNPS found:16"
# [1] "rs3072"     "rs7255"     "rs2687202"  "rs9823696"  "rs9257809"  "rs2188554"  "rs17451754" "rs17749155"
# [9] "rs10108511" "rs11789015" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"  "rs10419226" "rs10423674"
dongresall1=check_snp_eqtl(dat=dong23snp,gtextable = jointalleqtl,opt=2)
# [1] "#SNPS found:19"
# [1] "2_20878820"   "2_200045039"  "2_20878406"   "3_70929983"   "3_183783353"  "6_29356331"   "6_62391538"  
# [8] "7_117256712"  "7_117040117"  "8_11435516"   "8_10068073"   "9_100310501"  "9_96716028"   "12_114673723"
# [15] "15_58362025"  "15_58267416"  "16_86396835"  "19_18803172"  "19_18817903"
dongresall1$eqtl=dongres1$snptable$SNP[!is.na(dongres1$snptable$GTEx_minp)]
dongresall1$eqtlpairs=dongres1$allpairs

levineresall=check_snp_eqtl(dat=levine12snp,gtextable = jointalleqtl)
# [1] "#SNPS found:12"
# [1] "rs2687201"  "rs11789015" "rs6479527"  "rs1490865"  "rs3111601"  "rs9936833"  "rs1728400"  "rs3950627" 
# [9] "rs2178146"  "rs13332095" "rs10419226" "rs10423674"
levineresall1=check_snp_eqtl(dat=levine12snp,gtextable = jointalleqtl,opt=2)
levineresall1$eqtl=levineres1$snptable$SNP[!is.na(levineres1$snptable$GTEx_minp)]
levineresall1$eqtlpairs=levineres1$allpairs

save(dongresall1,puyaresall1,levineresall1,file="../result/Dong_Puya_Levine_GTEx_junction.RData")

#Mucosa results--------------------
mucosaeqtl=as.data.frame(fread(paste0(eqtldir,"Esophagus_Mucosa.v7.signif_variant_gene_pairs.txt")))
save(mucosaeqtl,file="../result/GTEx_mucosa_V7_eQTL_gene_pairs.RData")
all(mucosaeqtl$variant_id %in% lookuptable$variant_id)
idx=match(mucosaeqtl$variant_id,lookuptable$variant_id)
mucosaeqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

mucosaeqtl$variant_id1=gsub("_b37","",mucosaeqtl$variant_id)
tmp=unlist(strsplit(mucosaeqtl$variant_id1,"_"))
mucosaeqtl$Chr=tmp[seq(1,length(tmp),4)]
mucosaeqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
mucosaeqtl$Alt=tmp[seq(3,length(tmp),4)]
mucosaeqtl$Ref=tmp[seq(4,length(tmp),4)]
mucosaeqtl$name=paste0(mucosaeqtl$Chr,"_",mucosaeqtl$Position,"_",mucosaeqtl$Alt,"_",mucosaeqtl$Ref)
mucosaeqtl$namepos=paste0(mucosaeqtl$Chr,"_",mucosaeqtl$Position)

mucosaalleqtl=as.data.frame(fread(paste0(eqtldir,"Esophagus_Mucosa.allpairs.txt")))
mucosaalleqtl$variant_id1=gsub("_b37","",mucosaalleqtl$variant_id)
tmp=unlist(strsplit(mucosaalleqtl$variant_id1,"_"))
mucosaalleqtl$Chr=tmp[seq(1,length(tmp),4)]
mucosaalleqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
mucosaalleqtl$Alt=tmp[seq(3,length(tmp),4)]
mucosaalleqtl$Ref=tmp[seq(4,length(tmp),4)]
mucosaalleqtl$name=paste0(mucosaalleqtl$Chr,"_",mucosaalleqtl$Position,"_",mucosaalleqtl$Alt,"_",mucosaalleqtl$Ref)
mucosaalleqtl$namepos=paste0(mucosaalleqtl$Chr,"_",mucosaalleqtl$Position)
all(mucosaalleqtl$variant_id %in% lookuptable$variant_id)
idx=match(mucosaalleqtl$variant_id,lookuptable$variant_id)
mucosaalleqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

puyares=check_snp_eqtl(gtextable = mucosaeqtl)
# [1] "#SNPS found:4"
# [1] "rs7255"     "rs9257809"  "rs17749155" "rs10108511"
puyares1=check_snp_eqtl(gtextable = mucosaeqtl,opt=2)
#View(puyares$snptable)
#View(puyares$eqtl)
dongres=check_snp_eqtl(gtextable = mucosaeqtl,dat=dong23snp)
# [1] "#SNPS found:5"
# [1] "rs3072"     "rs7255"     "rs9257809"  "rs17749155" "rs10108511"
dongres1=check_snp_eqtl(gtextable = mucosaeqtl,dat=dong23snp,opt=2)

levineres=check_snp_eqtl(gtextable = mucosaeqtl,dat=levine12snp)
# [1] "#SNPS found:0"
levineres1=check_snp_eqtl(gtextable = mucosaeqtl,dat=levine12snp,opt=2)

puyaresall=check_snp_eqtl(gtextable = mucosaalleqtl)
# [1] "#SNPS found:10"
# [1] "rs7255"     "rs2687202"  "rs9257809"  "rs17451754" "rs17749155" "rs10108511" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"
puyaresall1=check_snp_eqtl(gtextable = mucosaalleqtl,opt=2)
# [1] "#SNPS found:11"
# [1] "2_20878820"   "15_58362025"  "7_117256712"  "8_10068073"   "8_11435516"   "3_70929983"   "12_114673723"
# [8] "6_29356331"   "9_100310501"  "2_200045039"  "16_86396835"
which(puyaresall1$snptable$SNP!=puyaresall1$snptable$GTEx_SNP)
#12 this is an example one SNP has 2 rsid (in different dbsnp version)
puyaresall1$eqtl=puyares1$snptable$SNP[!is.na(puyares1$snptable$GTEx_minp)]
puyaresall1$eqtlpairs=puyares1$allpairs

dongresall=check_snp_eqtl(dat=dong23snp,gtextable = mucosaalleqtl)
# [1] "#SNPS found:16"
# [1] "rs3072"     "rs7255"     "rs2687202"  "rs9823696"  "rs9257809"  "rs2188554"  "rs17451754" "rs17749155"
# [9] "rs10108511" "rs11789015" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"  "rs10419226" "rs10423674"
dongresall1=check_snp_eqtl(dat=dong23snp,gtextable = mucosaalleqtl,opt=2)
# "#SNPS found:18"
# [1] "2_20878820"   "2_200045039"  "2_20878406"   "3_70929983"   "3_183783353"  "6_29356331"   "7_117256712"  "7_117040117" 
# [9] "8_11435516"   "8_10068073"   "9_100310501"  "9_96716028"   "12_114673723" "15_58362025"  "15_58267416"  "16_86396835" 
# [17] "19_18803172"  "19_18817903" 
dongresall1$eqtl=dongres1$snptable$SNP[!is.na(dongres1$snptable$GTEx_minp)]
dongresall1$eqtlpairs=dongres1$allpairs

levineresall=check_snp_eqtl(dat=levine12snp,gtextable = mucosaalleqtl)
# [1] "#SNPS found:12"
# [1] "rs2687201"  "rs11789015" "rs6479527"  "rs1490865"  "rs3111601"  "rs9936833"  "rs1728400"  "rs3950627" 
# [9] "rs2178146"  "rs13332095" "rs10419226" "rs10423674"
levineresall1=check_snp_eqtl(dat=levine12snp,gtextable = mucosaalleqtl,opt=2)
levineresall1$eqtl=levineres1$snptable$SNP[!is.na(levineres1$snptable$GTEx_minp)]
levineresall1$eqtlpairs=levineres1$allpairs

save(dongresall1,puyaresall1,levineresall1,file="../result/Dong_Puya_Levine_GTEx_mucosa.RData")

#Stomach--------
stomacheqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Stomach.v7.signif_variant_gene_pairs.txt.gz")))
save(stomacheqtl,file="../result/GTEx_stomach_V7_eQTL_gene_pairs.RData")
all(stomacheqtl$variant_id %in% lookuptable$variant_id)
idx=match(stomacheqtl$variant_id,lookuptable$variant_id)
stomacheqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

stomacheqtl$variant_id1=gsub("_b37","",stomacheqtl$variant_id)
tmp=unlist(strsplit(stomacheqtl$variant_id1,"_"))
stomacheqtl$Chr=tmp[seq(1,length(tmp),4)]
stomacheqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
stomacheqtl$Alt=tmp[seq(3,length(tmp),4)]
stomacheqtl$Ref=tmp[seq(4,length(tmp),4)]
stomacheqtl$name=paste0(stomacheqtl$Chr,"_",stomacheqtl$Position,"_",stomacheqtl$Alt,"_",stomacheqtl$Ref)
stomacheqtl$namepos=paste0(stomacheqtl$Chr,"_",stomacheqtl$Position)

stomachalleqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Stomach.allpairs.txt.gz")))
stomachalleqtl$variant_id1=gsub("_b37","",stomachalleqtl$variant_id)
tmp=unlist(strsplit(stomachalleqtl$variant_id1,"_"))
stomachalleqtl$Chr=tmp[seq(1,length(tmp),4)]
stomachalleqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
stomachalleqtl$Alt=tmp[seq(3,length(tmp),4)]
stomachalleqtl$Ref=tmp[seq(4,length(tmp),4)]
stomachalleqtl$name=paste0(stomachalleqtl$Chr,"_",stomachalleqtl$Position,"_",stomachalleqtl$Alt,"_",stomachalleqtl$Ref)
stomachalleqtl$namepos=paste0(stomachalleqtl$Chr,"_",stomachalleqtl$Position)
all(stomachalleqtl$variant_id %in% lookuptable$variant_id)
idx=match(stomachalleqtl$variant_id,lookuptable$variant_id)
stomachalleqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

puyares=check_snp_eqtl(gtextable = stomacheqtl)
# [1] "#SNPS found:3"
# [1] "rs7255"     "rs9257809"  "rs10108511"
puyares1=check_snp_eqtl(gtextable = stomacheqtl,opt=2)
#View(puyares$snptable)
#View(puyares$eqtl)
dongres=check_snp_eqtl(gtextable = stomacheqtl,dat=dong23snp)
# [1] "#SNPS found:4"
# [1] "rs3072"     "rs7255"     "rs9257809"  "rs10108511"
dongres1=check_snp_eqtl(gtextable = stomacheqtl,dat=dong23snp,opt=2)

levineres=check_snp_eqtl(gtextable = stomacheqtl,dat=levine12snp)
# [1] "#SNPS found:0"
levineres1=check_snp_eqtl(gtextable = stomacheqtl,dat=levine12snp,opt=2)

puyaresall=check_snp_eqtl(gtextable = stomachalleqtl)
# [1] "#SNPS found:10"
# [1] [1] "rs7255"     "rs2687202"  "rs9257809"  "rs17451754" "rs17749155" "rs10108511" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654" 
puyaresall1=check_snp_eqtl(gtextable = stomachalleqtl,opt=2)
# [1] "#SNPS found:11"
# [1] "2_20878820"   "15_58362025"  "7_117256712"  "8_10068073"   "8_11435516"   "3_70929983"   "12_114673723"
# [8] "6_29356331"   "9_100310501"  "2_200045039"  "16_86396835"
which(puyaresall1$snptable$SNP!=puyaresall1$snptable$GTEx_SNP)
#12 this is an example one SNP has 2 rsid (in different dbsnp version)
puyaresall1$eqtl=puyares1$snptable$SNP[!is.na(puyares1$snptable$GTEx_minp)]
puyaresall1$eqtlpairs=puyares1$allpairs

dongresall=check_snp_eqtl(dat=dong23snp,gtextable = stomachalleqtl)
# [1] "#SNPS found:16"
# [1] "rs3072"     "rs7255"     "rs2687202"  "rs9823696"  "rs9257809"  "rs2188554"  "rs17451754" "rs17749155"
# [9] "rs10108511" "rs11789015" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"  "rs10419226" "rs10423674"
dongresall1=check_snp_eqtl(dat=dong23snp,gtextable = stomachalleqtl,opt=2)
# [1] "#SNPS found:19"
# [1] "2_20878820"   "2_200045039"  "2_20878406"   "3_70929983"   "3_183783353"  "6_29356331"   "6_62391538"   "7_117256712" 
# [9] "7_117040117"  "8_11435516"   "8_10068073"   "9_100310501"  "9_96716028"   "12_114673723" "15_58362025"  "15_58267416" 
# [17] "16_86396835"  "19_18803172"  "19_18817903"
dongresall1$eqtl=dongres1$snptable$SNP[!is.na(dongres1$snptable$GTEx_minp)]
dongresall1$eqtlpairs=dongres1$allpairs

levineresall=check_snp_eqtl(dat=levine12snp,gtextable = stomachalleqtl)
# [1] "#SNPS found:12"
# [1] "rs2687201"  "rs11789015" "rs6479527"  "rs1490865"  "rs3111601"  "rs9936833"  "rs1728400"  "rs3950627" 
# [9] "rs2178146"  "rs13332095" "rs10419226" "rs10423674"
levineresall1=check_snp_eqtl(dat=levine12snp,gtextable = stomachalleqtl,opt=2)
levineresall1$eqtl=levineres1$snptable$SNP[!is.na(levineres1$snptable$GTEx_minp)]
levineresall1$eqtlpairs=levineres1$allpairs

save(dongresall1,puyaresall1,levineresall1,file="../result/Dong_Puya_Levine_GTEx_stomach.RData")

#Esophagus_Muscularis------------------------------------------
musculariseqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Esophagus_Muscularis.v7.signif_variant_gene_pairs.txt.gz")))
save(musculariseqtl,file="../result/GTEx_muscularis_V7_eQTL_gene_pairs.RData")
all(musculariseqtl$variant_id %in% lookuptable$variant_id)
idx=match(musculariseqtl$variant_id,lookuptable$variant_id)
musculariseqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

musculariseqtl$variant_id1=gsub("_b37","",musculariseqtl$variant_id)
tmp=unlist(strsplit(musculariseqtl$variant_id1,"_"))
musculariseqtl$Chr=tmp[seq(1,length(tmp),4)]
musculariseqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
musculariseqtl$Alt=tmp[seq(3,length(tmp),4)]
musculariseqtl$Ref=tmp[seq(4,length(tmp),4)]
musculariseqtl$name=paste0(musculariseqtl$Chr,"_",musculariseqtl$Position,"_",musculariseqtl$Alt,"_",musculariseqtl$Ref)
musculariseqtl$namepos=paste0(musculariseqtl$Chr,"_",musculariseqtl$Position)

muscularisalleqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Esophagus_Muscularis.allpairs.txt.gz")))
muscularisalleqtl$variant_id1=gsub("_b37","",muscularisalleqtl$variant_id)
tmp=unlist(strsplit(muscularisalleqtl$variant_id1,"_"))
muscularisalleqtl$Chr=tmp[seq(1,length(tmp),4)]
muscularisalleqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
muscularisalleqtl$Alt=tmp[seq(3,length(tmp),4)]
muscularisalleqtl$Ref=tmp[seq(4,length(tmp),4)]
muscularisalleqtl$name=paste0(muscularisalleqtl$Chr,"_",muscularisalleqtl$Position,"_",muscularisalleqtl$Alt,"_",muscularisalleqtl$Ref)
muscularisalleqtl$namepos=paste0(muscularisalleqtl$Chr,"_",muscularisalleqtl$Position)
all(muscularisalleqtl$variant_id %in% lookuptable$variant_id)
idx=match(muscularisalleqtl$variant_id,lookuptable$variant_id)
muscularisalleqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

puyares=check_snp_eqtl(gtextable = musculariseqtl)
# [1] "#SNPS found:3"
# [1] "rs7255"     "rs9257809" "rs10108511"
puyares1=check_snp_eqtl(gtextable = musculariseqtl,opt=2)
#View(puyares$snptable)
#View(puyares$eqtl)
dongres=check_snp_eqtl(gtextable = musculariseqtl,dat=dong23snp)
# [1] "#SNPS found:4"
# [1] "rs3072"     "rs7255"     "rs9257809"  "rs10108511"
dongres1=check_snp_eqtl(gtextable = musculariseqtl,dat=dong23snp,opt=2)

levineres=check_snp_eqtl(gtextable = musculariseqtl,dat=levine12snp)
# [1] "#SNPS found:0"
levineres1=check_snp_eqtl(gtextable = musculariseqtl,dat=levine12snp,opt=2)

puyaresall=check_snp_eqtl(gtextable = muscularisalleqtl)
# [1] "#SNPS found:10"
# [1] "rs7255"     "rs2687202"  "rs9257809"  "rs17451754" "rs17749155" "rs10108511" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"
puyaresall1=check_snp_eqtl(gtextable = muscularisalleqtl,opt=2)
# [1] "#SNPS found:11"
# [1] "2_20878820"   "15_58362025"  "7_117256712"  "8_10068073"   "8_11435516"   "3_70929983"   "12_114673723"
# [8] "6_29356331"   "9_100310501"  "2_200045039"  "16_86396835"
which(puyaresall1$snptable$SNP!=puyaresall1$snptable$GTEx_SNP)
#12 this is an example one SNP has 2 rsid (in different dbsnp version)
puyaresall1$eqtl=puyares1$snptable$SNP[!is.na(puyares1$snptable$GTEx_minp)]
puyaresall1$eqtlpairs=puyares1$allpairs

dongresall=check_snp_eqtl(dat=dong23snp,gtextable = muscularisalleqtl)
# [1] "#SNPS found:16"
# [1] "rs3072"     "rs7255"     "rs2687202"  "rs9823696"  "rs9257809"  "rs2188554"  "rs17451754" "rs17749155"
# [9] "rs10108511" "rs11789015" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"  "rs10419226" "rs10423674"
dongresall1=check_snp_eqtl(dat=dong23snp,gtextable = muscularisalleqtl,opt=2)
# [1] "#SNPS found:19"
# [1] "2_20878820"   "2_200045039"  "2_20878406"   "3_70929983"   "3_183783353"  "6_29356331"   "6_62391538"   "7_117256712"  "7_117040117" 
# [10] "8_11435516"   "8_10068073"   "9_100310501"  "9_96716028"   "12_114673723" "15_58362025"  "15_58267416"  "16_86396835"  "19_18803172" 
# [19] "19_18817903" 
dongresall1$eqtl=dongres1$snptable$SNP[!is.na(dongres1$snptable$GTEx_minp)]
dongresall1$eqtlpairs=dongres1$allpairs

levineresall=check_snp_eqtl(dat=levine12snp,gtextable = muscularisalleqtl)
# [1] "#SNPS found:12"
# [1] "rs2687201"  "rs11789015" "rs6479527"  "rs1490865"  "rs3111601"  "rs9936833"  "rs1728400"  "rs3950627" 
# [9] "rs2178146"  "rs13332095" "rs10419226" "rs10423674"
levineresall1=check_snp_eqtl(dat=levine12snp,gtextable = muscularisalleqtl,opt=2)
levineresall1$eqtl=levineres1$snptable$SNP[!is.na(levineres1$snptable$GTEx_minp)]
levineresall1$eqtlpairs=levineres1$allpairs

save(dongresall1,puyaresall1,levineresall1,file="../result/Dong_Puya_Levine_GTEx_muscularis.RData")

#Whole blood------------------------------------------------
bloodeqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Whole_Blood.v7.signif_variant_gene_pairs.txt.gz")))
save(bloodeqtl,file="../result/GTEx_blood_V7_eQTL_gene_pairs.RData")
all(bloodeqtl$variant_id %in% lookuptable$variant_id)
idx=match(bloodeqtl$variant_id,lookuptable$variant_id)
bloodeqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

bloodeqtl$variant_id1=gsub("_b37","",bloodeqtl$variant_id)
tmp=unlist(strsplit(bloodeqtl$variant_id1,"_"))
bloodeqtl$Chr=tmp[seq(1,length(tmp),4)]
bloodeqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
bloodeqtl$Alt=tmp[seq(3,length(tmp),4)]
bloodeqtl$Ref=tmp[seq(4,length(tmp),4)]
bloodeqtl$name=paste0(bloodeqtl$Chr,"_",bloodeqtl$Position,"_",bloodeqtl$Alt,"_",bloodeqtl$Ref)
bloodeqtl$namepos=paste0(bloodeqtl$Chr,"_",bloodeqtl$Position)

bloodalleqtl=as.data.frame(fread(paste0("gunzip -cq ",eqtldir,"Whole_Blood.allpairs.txt.gz")))
bloodalleqtl$variant_id1=gsub("_b37","",bloodalleqtl$variant_id)
tmp=unlist(strsplit(bloodalleqtl$variant_id1,"_"))
bloodalleqtl$Chr=tmp[seq(1,length(tmp),4)]
bloodalleqtl$Position=as.numeric(tmp[seq(2,length(tmp),4)])
bloodalleqtl$Alt=tmp[seq(3,length(tmp),4)]
bloodalleqtl$Ref=tmp[seq(4,length(tmp),4)]
bloodalleqtl$name=paste0(bloodalleqtl$Chr,"_",bloodalleqtl$Position,"_",bloodalleqtl$Alt,"_",bloodalleqtl$Ref)
bloodalleqtl$namepos=paste0(bloodalleqtl$Chr,"_",bloodalleqtl$Position)
all(bloodalleqtl$variant_id %in% lookuptable$variant_id)
idx=match(bloodalleqtl$variant_id,lookuptable$variant_id)
bloodalleqtl$SNP=lookuptable$rs_id_dbSNP147_GRCh37p13[idx]

puyares=check_snp_eqtl(gtextable = bloodeqtl)
# [1] "#SNPS found:3"
# [1] "rs7255"     "rs9257809" "rs10108511"
puyares1=check_snp_eqtl(gtextable = bloodeqtl,opt=2)
#View(puyares$snptable)
#View(puyares$eqtl)
dongres=check_snp_eqtl(gtextable = bloodeqtl,dat=dong23snp)
# [1] "#SNPS found:5"
# [1] "rs3072"     "rs7255"     "rs9257809"  "rs10108511" "rs10419226"
dongres1=check_snp_eqtl(gtextable = bloodeqtl,dat=dong23snp,opt=2)

levineres=check_snp_eqtl(gtextable = bloodeqtl,dat=levine12snp)
# [1] "#SNPS found:1"  "rs10419226"
levineres1=check_snp_eqtl(gtextable = bloodeqtl,dat=levine12snp,opt=2)

puyaresall=check_snp_eqtl(gtextable = bloodalleqtl)
# [1] "#SNPS found:10"
# [1] "rs7255"     "rs2687202"  "rs9257809"  "rs17451754" "rs17749155" "rs10108511" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"
puyaresall1=check_snp_eqtl(gtextable = bloodalleqtl,opt=2)
# [1] "#SNPS found:11"
# [1] "2_20878820"   "15_58362025"  "7_117256712"  "8_10068073"   "8_11435516"   "3_70929983"   "12_114673723"
# [8] "6_29356331"   "9_100310501"  "2_200045039"  "16_86396835"
which(puyaresall1$snptable$SNP!=puyaresall1$snptable$GTEx_SNP)
#12 this is an example one SNP has 2 rsid (in different dbsnp version)
puyaresall1$eqtl=puyares1$snptable$SNP[!is.na(puyares1$snptable$GTEx_minp)]
puyaresall1$eqtlpairs=puyares1$allpairs

dongresall=check_snp_eqtl(dat=dong23snp,gtextable = bloodalleqtl)
# [1] "#SNPS found:16"
# [1] "rs3072"     "rs7255"     "rs2687202"  "rs9823696"  "rs9257809"  "rs2188554"  "rs17451754" "rs17749155"
# [9] "rs10108511" "rs11789015" "rs7852462"  "rs1247942"  "rs2464469"  "rs1979654"  "rs10419226" "rs10423674"
dongresall1=check_snp_eqtl(dat=dong23snp,gtextable = bloodalleqtl,opt=2)
# [1] "#SNPS found:19"
# [1] "2_20878820"   "2_200045039"  "2_20878406"   "3_70929983"   "3_183783353"  "6_29356331"   "6_62391538"   "7_117256712"  "7_117040117" 
# [10] "8_11435516"   "8_10068073"   "9_100310501"  "9_96716028"   "12_114673723" "15_58362025"  "15_58267416"  "16_86396835"  "19_18803172" 
# [19] "19_18817903" 
dongresall1$eqtl=dongres1$snptable$SNP[!is.na(dongres1$snptable$GTEx_minp)]
dongresall1$eqtlpairs=dongres1$allpairs

levineresall=check_snp_eqtl(dat=levine12snp,gtextable = bloodalleqtl)
# [1] "#SNPS found:12"
# [1] "rs2687201"  "rs11789015" "rs6479527"  "rs1490865"  "rs3111601"  "rs9936833"  "rs1728400"  "rs3950627" 
# [9] "rs2178146"  "rs13332095" "rs10419226" "rs10423674"
levineresall1=check_snp_eqtl(dat=levine12snp,gtextable = bloodalleqtl,opt=2)
levineresall1$eqtl=levineres1$snptable$SNP[!is.na(levineres1$snptable$GTEx_minp)]
levineresall1$eqtlpairs=levineres1$allpairs

save(dongresall1,puyaresall1,levineresall1,file="../result/Dong_Puya_Levine_GTEx_blood.RData")

dongeqtl=dong23snp
dongeqtl$junction=dongeqtl$mucosa=dongeqtl$musculari=dongeqtl$stomach=dongeqtl$blood=F
load("../result/Dong_Puya_Levine_GTEx_junction.RData")
dongeqtl$junction[dongeqtl$SNP %in% dongresall1$eqtl]=T
load("../result/Dong_Puya_Levine_GTEx_mucosa.RData")
dongeqtl$mucosa[dongeqtl$SNP %in% dongresall1$eqtl]=T
load("../result/Dong_Puya_Levine_GTEx_muscularis.RData")
dongeqtl$musculari[dongeqtl$SNP %in% dongresall1$eqtl]=T
load("../result/Dong_Puya_Levine_GTEx_stomach.RData")
dongeqtl$stomach[dongeqtl$SNP %in% dongresall1$eqtl]=T
load("../result/Dong_Puya_Levine_GTEx_blood.RData")
dongeqtl$blood[dongeqtl$SNP %in% dongresall1$eqtl]=T 
dongeqtl$sum=rowSums(dongeqtl[,9:13])
idx=order(dongeqtl$sum,decreasing = T)
dongeqtl=dongeqtl[idx,]
for (i in 9:13)
{
  dongeqtl[,i]=as.numeric(dongeqtl[,i])
}
tmp=dongeqtl[,9:13]
rownames(tmp)=dongeqtl$SNP
par(mar=c(6,2,2,5))
library(RColorBrewer)
# redgreen <- c("red", "green")
# pal <- colorRampPalette(redgreen)(100)
heatmap(as.matrix(tmp),Rowv=NA,Colv = NA,scale="none",margins = c(7,7))
