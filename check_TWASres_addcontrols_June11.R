#!/usr/bin/env Rscript
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
  BE_fwer=compute_fwer(dat,namecol="BE_p")
  EA_fdr=compute_fdr(dat,namecol="EA_p")
  EA_fwer=compute_fwer(dat,namecol="EA_p")
  BEA_fdr=compute_fdr(dat,namecol="BEA_p")
  BEA_fwer=compute_fwer(dat,namecol="BEA_p")
  BEEA_fdr=compute_fdr(dat,namecol="BEEA_p")
  BEEA_fwer=compute_fwer(dat,namecol="BEEA_p")
  return(list(BE_fdr=BE_fdr,BE_fwer=BE_fwer,EA_fdr=EA_fdr,EA_fwer=EA_fwer,
              BEA_fdr=BEA_fdr,BEA_fwer=BEA_fwer,BEEA_fdr=BEEA_fdr,BEEA_fwer=BEEA_fwer))
}


#load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11" #models stored for GTEx, came from prediction_michigan_models6_GTExV8.R
check_skatres=function(dat=skat_min1,opt="plot")
{
  colnames(dat)=c("BE_p","EA_p","BEA_p","BEEA_p")
  dat_code=dat[rownames(dat) %in% proteingenes,]
  tmp=compute_fwer_fdr(dat=dat_code)
  if (opt=="plot")
  {
    for(i in 1:length(tmp))
    {
      if (!is.null(tmp[[i]]))
      {
        print(paste0(names(tmp)[i],":"))
        tmp1=tmp[[i]]
        if (length(tmp1)<=7)
        {
          print(tmp1)
        }else
        {
          j=seq(1,length(tmp1),7)
          j[2:length(j)]=j[2:length(j)]-1
          j=c(j,length(tmp1))
          for (k in 1:(length(j)-1))
          {
            tmp2=tmp1[j[k]:j[k+1]]
            print(tmp2)
          }
        }
        #print(tmp[[i]])
      }
    }
    par(mfrow=c(2,2))
    par(mar=c(6,6,3,1))
    qqplot(dat_code$BE_p,main="BE")
    qqplot(dat_code$EA_p,main="EA")
    qqplot(dat_code$BEA_p,main="BE vs EA")
    qqplot(dat_code$BEEA_p,main="BEEA")
  }
  return(tmp)
}

check_model=function(mod=res_min,r2cut=0)
{
  mod=mod[which(mod$r2>=r2cut),]
  print(nrow(mod))
  print(quantile(mod$r2))
  print(quantile(mod$numselectedsnp))
}

scinumber=function(num)
{
  formatC(num,format="e",digits = 2)
}

check_res=function(r2cutoff=0)
{
  print(paste0("r2cutoff=",r2cutoff))
  genesr2=rownames(res_min_code)[res_min_code$r2>r2cutoff]
  idx=which(rownames(assoc_min_code) %in% genesr2)
  tmp=compute_fwer_fdr(dat=assoc_min_code[idx,])
  print(paste0("regular TWAS:"))
  print(tmp)
  par(mfrow=c(2,2),mar=c(5,5.5,2,1))
  qqplot(assoc_min_code$BE_p[idx],main="BE")
  qqplot(assoc_min_code$EA_p[idx],main="EA")
  qqplot(assoc_min_code$BEA_p[idx],main="BEA") #EA vs BE
  qqplot(assoc_min_code$BEEA_p[idx],main="BEEA")
  idx=which(rownames(skat_min2_code) %in% genesr2)
  print("skat-TWAS result:")
  tmp=p.adjust(skat_min2_code$BEA_p[idx],method="fdr")
  print(paste0("BE/EA,fdr:",min(tmp)))
  tmp=check_skatres(dat=skat_min2_code[idx,])
  genes=NULL
  for (i in 1:length(tmp))
  {
    genes=unique(c(genes,names(tmp[[i]])))
  }
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
    idx2=which(rownames(res_min)==gene)
    print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                 ",R2=",scinumber(res_min$r2[idx2]),",#SNP=",res_min$numselectedsnp[idx2]))
    
  }
  return(skat_min2_code[idx,])
}

print_genes=function(genes=NULL)
{
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
    idx2=which(rownames(res_min)==gene)
    print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                 ",R2=",scinumber(res_min$r2[idx2]),",#SNP=",res_min$numselectedsnp[idx2]))
    
  }
}

#GTEx junction------
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)

# #check v7 genes
# genes=c("KXD1","UBAC1","ISYNA1")
# for (gene in genes)
# {
#   idx=which(rownames(skat_min2_code)==gene)
#   print(paste0(gene,": BE-p ",scinumber(skat_min2_code$BE_p[idx])))
#   print(paste0(gene,": EA-p ",scinumber(skat_min2_code$EA_p[idx])))
#   print(paste0(gene," BEEA-p ",scinumber(skat_min2_code$BEEA_p[idx])))
# }

#GTEx stomach------
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)

#GTEx blood------
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8blooddata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
print_genes(genes="TAPT1")
#"TAPT1 chr4:16160505,R2=6.58e-02,#SNP=63"
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)

# load(paste0(outfolder,"/skat_3gwas_res.RData"))
# colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
# skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
# tmp1=check_res()
# idx=match(genes,rownames(tmp))
# tmp1[idx,]
# # BE_p      EA_p     BEA_p    BEEA_p
# # KLHL26 0.4915183 0.9735564 0.2844699 0.8912872
# # PGPEP1 0.4032439 0.8299318 0.3935246 0.5317461
# # LSM4   0.3838819 0.9321208 0.7118257 0.4740115
# # GDF15  0.5593960 0.7814989 0.3670895 0.5893278
# # HOMER3 0.7362117 0.5671278 0.2337186 0.7806328

#GTEx mucosa------
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)


# load(paste0(outfolder,"/skat_3gwas_res.RData"))
# colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
# skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
# tmp1=check_res()
# idx=match(genes,rownames(tmp))
# tmp1[idx,]
# # BE_p       EA_p     BEA_p     BEEA_p
# # SLC25A42 0.2745835 0.03689584 0.2225576 0.10134638
# # TMEM161A 0.1393658 0.04392426 0.1934317 0.06722125



#meta analysis
load(paste0(outfolder,"/skat_meta_res.RData"))
skat_min2=meta_p$resp
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p")
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)

#compare p-values-
outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
outfolder2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
snppos1=snppos
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_addcontrols_for_prediction.RData")
snppos2=snppos
load(paste0(outfolder1,"/bca_assoc.RData"))
assoc_min_code1=assoc_min[rownames(assoc_min) %in% proteingenes,]
load(paste0(outfolder2,"/bca_assoc.RData"))
assoc_min_code2=assoc_min[rownames(assoc_min) %in% proteingenes,]
comgenes=intersect(rownames(assoc_min_code1),rownames(assoc_min_code2))
idx1=match(comgenes,rownames(assoc_min_code1))
idx2=match(comgenes,rownames(assoc_min_code2))
plot(-log10(assoc_min_code1$BEEA_p[idx1]),-log10(assoc_min_code2$BEEA_p[idx2]),xlab="Old p-value (-log10)",ylab="New p-value (-log10)")

load(paste0(outfolder1,"/skat_res.RData"))
skat_min2_code1=skat_min2[rownames(skat_min2) %in% proteingenes,]
load(paste0(outfolder2,"/skat_res.RData"))
skat_min2_code2=skat_min2_pc4[rownames(skat_min2_pc4) %in% proteingenes,]
colnames(skat_min2_code1)=colnames(skat_min2_code2)=c("BE","EA","BEA","BEEA")
comgenes=intersect(rownames(skat_min2_code1),rownames(skat_min2_code2))
idx1=match(comgenes,rownames(skat_min2_code1))
idx2=match(comgenes,rownames(skat_min2_code2))
plot(-log10(skat_min2_code1$V4[idx1]),-log10(skat_min2_code2$V4[idx2]),xlab="Old p-value (-log10)",ylab="New p-value (-log10)")
abline(0,1,col="red")
cor1=round(cor(-log10(skat_min2_code1$V4[idx1]),-log10(skat_min2_code2$V4[idx2])),2)

plot(-log10(skat_min2_code1$V2[idx1]),-log10(skat_min2_code2$V2[idx2]),xlab="Old p-value (-log10)",ylab="New p-value (-log10)")
abline(0,1,col="red")
cor1=round(cor(-log10(skat_min2_code1$V2[idx1]),-log10(skat_min2_code2$V2[idx2])),2)

plot(-log10(skat_min2_code1$V1[idx1]),-log10(skat_min2_code2$V1[idx2]),xlab="Old p-value (-log10)",ylab="New p-value (-log10)")
abline(0,1,col="red")
cor1=round(cor(-log10(skat_min2_code1$V1[idx1]),-log10(skat_min2_code2$V1[idx2])),2)
which(-log10(skat_min2_code1$V1[idx1])>6 & -log10(skat_min2_code2$V1[idx2])<2) #10787

load(paste0(outfolder1,"/preidiction_michigan_model.RData"))
res_min1=res_min
load(paste0(outfolder2,"/preidiction_michigan_model.RData"))
res_min2=res_min
load(paste0(outfolder1,"/bca_predict_geneexp.RData"))
predict_min10=predict_min[,1:2]
predict_min1=predict_min[,3:ncol(predict_min)]
geneexpsamplenames=strsplit(colnames(predict_min1),"_") #use localid
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_min1)=geneexpsamplenames
load(paste0(outfolder2,"/bca_predict_geneexp.RData"))
predict_min20=predict_min[,1:2]
predict_min2=predict_min[,3:ncol(predict_min)]

load(load(paste0(outfolder1,"/bca_extractgenotype.RData")))
bcagenotype1=bcagenotype
geneexpsamplenames=strsplit(colnames(bcagenotype1),"_") #use localid
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(bcagenotype1)=geneexpsamplenames
load(load(paste0(outfolder2,"/bca_extractgenotype.RData")))
bcagenotype2=bcagenotype
get_genmodel=function(gene="PPP2R5C",genemodel=res_min1)
{
  idx=which(rownames(genemodel)==gene)
  res=data.frame(snp=rep(NA,genemodel$numselectedsnp[idx]),coeff=NA)
  res$snp=unlist(strsplit(genemodel$selectedsnps[idx],"|",fixed=T))
  res$coeff=as.numeric(unlist(strsplit(genemodel$selectedsnps_coeff[idx],"|",fixed=T)))
  print(paste0("r2=",genemodel$r2[idx]))
  return(res)
}

get_bcagenotype=function(gene="PPP2R5C",genemodel=res_min1,genotype=bcagenotype1)
{
  tmp=get_genmodel(gene=gene,genemodel = genemodel)
  tmp=intersect(tmp$snp,rownames(genotype))
  idx=match(tmp,rownames(genotype))
  res=genotype[idx,]
  return(res)
}
get_geneexp=function(gene="PPP2R5C",geneexp=predict_min1)
{
  idx=which(rownames(geneexp)==gene)
  res=data.frame(geneexp=unlist(geneexp[idx,]))
  return(res)
}
check_gene=function(gene="PPP2R5C")
{
  genemodel1=get_genmodel(gene=gene,genemodel = res_min1)
  genemodel2=get_genmodel(gene=gene,genemodel = res_min2)
  print(paste0(sum(genemodel1$snp %in% genemodel2$snp)," common snps")) #25
  sum(genemodel1$snp %in% rownames(snppos2))
  sum(genemodel2$snp %in% rownames(snppos1))
  geneexp1=get_geneexp(gene=gene,geneexp = predict_min1)
  geneexp2=get_geneexp(gene=gene,geneexp = predict_min2)
  comsamples=intersect(rownames(geneexp1),rownames(geneexp2))
  idx1=match(comsamples,rownames(geneexp1))
  idx2=match(comsamples,rownames(geneexp2))
  plot(geneexp1$geneexp[idx1],geneexp2$geneexp[idx2],xlab="Predicted expr (old)",ylab="Predicted expr (new)")
  abline(0,1,col="red")
  print(paste0("correlation between predicted expr: ",cor(geneexp1$geneexp[idx1],geneexp2$geneexp[idx2])))
  idx1=which(rownames(assoc_min_code1)==gene)
  idx2=which(rownames(assoc_min_code2)==gene)
  print(assoc_min_code1[idx1,])
  print(assoc_min_code2[idx2,])
  genotype1=get_bcagenotype(gene=gene,genemodel = res_min1,genotype=bcagenotype1)
  genotype2=get_bcagenotype(gene=gene,genemodel = res_min2,genotype=bcagenotype2)
  comsnps=intersect(rownames(genotype1),rownames(genotype2))
  comsamples=intersect(colnames(genotype1),colnames(genotype2))
  genotype1_=genotype1[match(comsnps,rownames(genotype1)),match(comsamples,colnames(genotype1))]
  genotype2_=genotype2[match(comsnps,rownames(genotype2)),match(comsamples,colnames(genotype2))]
  sum(genotype1_!=genotype2_)/nrow(genotype1_)/ncol(genotype1_)
  skat_min2_code1[which(rownames(skat_min2_code1)==gene),]
  skat_min2_code2[which(rownames(skat_min2_code2)==gene),]
}
gene="TMEM161A"
gene="SLC25A42"
gene="LDAH"
gene="JUND"
skat_min2_code1[idx1[10787],]
#GTEx muscularis------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8muscularisdata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
print_genes(genes="HOXB8")

idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)

#GTEx adipose------
#outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8adiposedata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_Jan19"
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code=res_min[rownames(res_min) %in% proteingenes,]
par(mfrow=c(1,1))
hist(res_min_code$r2[res_min_code$glmflag==1],main="",xlab="R-squared",col="blue")
check_model(mod=res_min_code)
check_model(mod=res_min_code,r2cut=0.05)
#Regular TWAS result
load(paste0(outfolder,"/bca_assoc.RData"))
assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
assoc_min_pc6_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]

#SKAT TWAS result
load(paste0(outfolder,"/skat_res.RData"))
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()

idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)
tmp1=p.adjust(tmp$EA_p,method="fdr")
min(tmp1)

#for pc6
assoc_min_code=assoc_min_pc6_code
skat_min2=skat_min2_pc6
colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") 
skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
tmp=check_res()
print_genes(genes="C10orf67")
idx=order(tmp$EA_p)
genes=rownames(tmp)[idx[1:2]]
tmp=check_res(r2cutoff=0.05)
tmp=check_res(r2cutoff=0.1)


get_genelist=function(tissues=c("junction","stomach","mucosa","muscularis","blood","adipose"),
                      r2cutoff=0)
{
  genetable=NULL
  for (i in 1:length(tissues))
  {
    tissue=tissues[i]
    outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_June11")
    if (tissue=="junction")
    {
      outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11")
    }
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/skat_res.RData"))
    colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    genesr2=rownames(res_min)[which(res_min$r2>r2cutoff)]
    skat_min2_code=skat_min2_code[rownames(skat_min2_code) %in% genesr2,]
    res=check_skatres(dat=skat_min2_code,opt="noplot")
    for (j in 1:(length(res)/2))
    {
      if (j==1)
      {
        comparison="BEvsCO"
      }
      if (j==2)
      {
        comparison="EAvsCO"
      }
      if (j==3)
      {
        comparison="EAvsBE"
      }
      if (j==4)
      {
        comparison="BEEAvsCO"
      }
      k=(j-1)*2+1 #fdr
      if (length(res[[k]])>0)
      {
        genesfdr=names(res[[k]])
        genesfwer=names(res[[k+1]])
        genesfdr=genesfdr[!genesfdr %in% genesfwer]
        if (length(genesfwer)>0)
        {
          for (l in 1:length(genesfwer))
          {
            gene=genesfwer[l]
            idx1=which(rownames(res_min)==gene)
            idx2=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
            generes=data.frame(gene=gene,type=comparison,cutoff="FWER",value=res[[k+1]][l],
                               tissue=tissue,R2=res_min$r2[idx1],numsnp=res_min$numselectedsnp[idx1],
                               chr=gtexv8_ge_anno$Chromosome[idx2],position=gtexv8_ge_anno$start[idx2],stringsAsFactors = F)
            genetable=rbind(genetable,generes)
          }
        }
        if (length(genesfdr)>0)
        {
          for (l in 1:length(genesfdr))
          {
            gene=genesfdr[l]
            idx1=which(rownames(res_min)==gene)
            idx2=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
            generes=data.frame(gene=gene,type=comparison,cutoff="FDR",value=res[[k]][l],
                               tissue=tissue,R2=res_min$r2[idx1],numsnp=res_min$numselectedsnp[idx1],
                               chr=gtexv8_ge_anno$Chromosome[idx2],position=gtexv8_ge_anno$start[idx2],stringsAsFactors = F)
            genetable=rbind(genetable,generes)
          }
        }
      }
    }
  }
  return(genetable)
}
genelist=get_genelist()
genelist05=get_genelist(r2cutoff = 0.05)
genelist1=get_genelist(r2cutoff = 0.1)

#list result for each gene
get_genetable=function(genelist)
{
  genes=unique(genelist$gene)
  genetable=data.frame(gene=genes,tissue=NA,type=NA,
                       chr=NA,position=NA,numtissue=NA,numtype=NA,stringsAsFactors = F)
  for (i in 1:nrow(genetable))
  {
    idx=which(genelist$gene==genes[i])
    genetable$chr[i]=genelist$chr[idx[1]]
    genetable$position[i]=genelist$position[idx[1]]
    genetable$tissue[i]=paste0(genelist$tissue[idx],collapse = "|")
    genetable$type[i]=paste0(genelist$type[idx],collapse = "|")
    genetable$numtissue[i]=length(unique(genelist$tissue[idx]))
    genetable$numtype[i]=length(unique(genelist$type[idx]))
  }
  return(genetable)
}
genetable=get_genetable(genelist)
genetable05=get_genetable(genelist05)
genetable1=get_genetable(genelist1)
write.csv(genelist,file="../result/SKAT_TWAS_genes_r2_0_.csv",row.names = F)
write.csv(genelist05,file="../result/SKAT_TWAS_genes_r2_0.05_.csv",row.names = F)
write.csv(genelist1,file="../result/SKAT_TWAS_genes_r2_0.1_.csv",row.names = F)

#check the overlap of eQTLs 
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_junction=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_blood=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
sum(snp_junction %in% snp_blood) #194059
