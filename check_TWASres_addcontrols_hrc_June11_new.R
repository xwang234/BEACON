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

#mark FDR genes with blue
qqplot=function(pvalue=NULL,fwer=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (-log base 10)",
       ylab="Observed p-value (-log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
  pvalue_order=pvalue[order(pvalue)]
  fwer=p.adjust(pvalue_order,method="bonferroni")
  fdr=p.adjust(pvalue_order,method="fdr")
  mycolor=rep("black",length(pvalue))
  mycolor[fdr<0.05]="blue"
  mycolor[fwer<0.05]="red"
  points(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),pch=16,col=mycolor)
  print(paste0("#FWER<0.05:",sum(mycolor=="red")))
  print(paste0("#FDR<0.05:",sum(mycolor=="blue")))
  if (sum(mycolor=="red")>0 & sum(mycolor=="blue")>0)
  legend("topleft",legend = c("FWER<0.05","FDR<0.05"),pch=c(16,16),col=c("red","blue"))
  if (sum(mycolor=="red")==0 & sum(mycolor=="blue")>0)
    legend("topleft",legend = c("FDR<0.05"),pch=c(16),col=c("blue"))
}

compute_fdr=function(dat,namecol="BE_p",cutoff=0.05)
{
  res=NULL
  pvalues=dat[,which(colnames(dat)==namecol)]
  names(pvalues)=rownames(dat)
  tmp=p.adjust(pvalues,method="fdr")
  if (sum(tmp<cutoff,na.rm = T)>0)
  {
    res=tmp[which(tmp<cutoff)]
    res=res[order(res)]
  }
  return(res)
}
compute_fwer=function(dat,namecol="BE_p",cutoff=0.05)
{
  res=NULL
  pvalues=dat[,which(colnames(dat)==namecol)]
  names(pvalues)=rownames(dat)
  tmp=p.adjust(pvalues,method = "bonferroni")
  if (sum(tmp<cutoff,na.rm = T)>0)
  {
    res=tmp[which(tmp<cutoff)]
    res=res[order(res)]
  }
  return(res)
}

compute_fwer_fdr=function(dat=skat_min_code,cutoff=0.05)
{
  
  BE_fdr=compute_fdr(dat,namecol="BE_p",cutoff = cutoff)
  BE_fwer=compute_fwer(dat,namecol="BE_p",cutoff = cutoff)
  EA_fdr=compute_fdr(dat,namecol="EA_p",cutoff = cutoff)
  EA_fwer=compute_fwer(dat,namecol="EA_p",cutoff = cutoff)
  BEA_fdr=compute_fdr(dat,namecol="BEA_p",cutoff = cutoff)
  BEA_fwer=compute_fwer(dat,namecol="BEA_p",cutoff = cutoff)
  BEEA_fdr=compute_fdr(dat,namecol="BEEA_p",cutoff = cutoff)
  BEEA_fwer=compute_fwer(dat,namecol="BEEA_p",cutoff = cutoff)
  return(list(BE_fdr=BE_fdr,BE_fwer=BE_fwer,EA_fdr=EA_fdr,EA_fwer=EA_fwer,
              BEA_fdr=BEA_fdr,BEA_fwer=BEA_fwer,BEEA_fdr=BEEA_fdr,BEEA_fwer=BEEA_fwer))
}

print_compute_fwer_fdr=function(dat=compute_fwer_fdr_res)
{
  for(i in 1:length(dat))
  {
    if (!is.null(dat[[i]]))
    {
      print(paste0(names(dat)[i],":"))
      dat1=dat[[i]]
      if (length(dat1)<=7)
      {
        print(dat1)
      }else
      {
        j=seq(1,length(dat1),7)
        j[2:length(j)]=j[2:length(j)]-1
        j=c(j,length(dat1))
        for (k in 1:(length(j)-1))
        {
          dat2=dat1[j[k]:j[k+1]]
          print(dat2)
        }
      }
    }
  }
}
#load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtexv8_ge_anno.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])

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

check_res=function(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,r2cutoff=0,plotopt="T")
{
  print(paste0("r2cutoff=",r2cutoff))
  genesr2=rownames(genemodel)[genemodel$r2>r2cutoff]
  idx=which(rownames(assocres) %in% genesr2)
  twasreg=compute_fwer_fdr(dat=assocres[idx,])
  print(paste0("regular TWAS:"))
  print_compute_fwer_fdr(twasreg)
  if (plotopt)
  {
    par(mfrow=c(2,2),mar=c(5,5.5,2,1))
    qqplot(assocres$BE_p[idx],main="BE")
    qqplot(assocres$EA_p[idx],main="EA")
    qqplot(assocres$BEA_p[idx],main="BEA") #EA vs BE
    qqplot(assocres$BEEA_p[idx],main="BEEA")
  }
  
  #print gene info
  genes=NULL
  for (i in 1:length(twasreg))
  {
    genes=unique(c(genes,names(twasreg[[i]])))
  }
  if (length(genes)>0)
  {
    for (i in 1:length(genes))
    {
      gene=genes[i]
      idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
      idx2=which(rownames(genemodel)==gene)
      print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                   ",R2=",scinumber(genemodel$r2[idx2]),",#SNP=",genemodel$numselectedsnp[idx2]))
    }
  }
  
  idx=which(rownames(skatres) %in% genesr2)
  print("skat-TWAS result:")
  # tmp=p.adjust(skatres$BEA_p[idx],method="fdr")
  # print(paste0("BE/EA,fdr:",min(tmp)))
  twasskat=compute_fwer_fdr(dat=skatres[idx,])
  print_compute_fwer_fdr(twasskat)
  if (plotopt)
  {
    qqplot(skatres$BE_p[idx],main="BE")
    qqplot(skatres$EA_p[idx],main="EA")
    qqplot(skatres$BEA_p[idx],main="BEA") #EA vs BE
    qqplot(skatres$BEEA_p[idx],main="BEEA")
  }
  
  genes=NULL
  for (i in 1:length(twasskat))
  {
    genes=unique(c(genes,names(twasskat[[i]])))
  }
  if (length(genes)>0)
  {
    for (i in 1:length(genes))
    {
      gene=genes[i]
      idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
      idx2=which(rownames(genemodel)==gene)
      print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                   ",R2=",scinumber(genemodel$r2[idx2]),",#SNP=",genemodel$numselectedsnp[idx2]))
    }
  }
  
  return(list(twasreg=twasreg,twasskat=twasskat))
}

check_res_heritability=function(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,heritability=NULL,plotopt="T")
{
  print(paste0("use heritability p<0.05"))
  #genesr2=rownames(genemodel)[genemodel$r2>r2cutoff]
  genesh=rownames(heritability)[heritability$V3<0.05]
  idx=which(rownames(assocres) %in% genesh)
  twasreg=compute_fwer_fdr(dat=assocres[idx,])
  print(paste0("regular TWAS:"))
  print_compute_fwer_fdr(twasreg)
  if (plotopt)
  {
    par(mfrow=c(2,2),mar=c(5,5.5,2,1))
    qqplot(assocres$BE_p[idx],main="BE")
    qqplot(assocres$EA_p[idx],main="EA")
    qqplot(assocres$BEA_p[idx],main="BEA") #EA vs BE
    qqplot(assocres$BEEA_p[idx],main="BEEA")
  }
  
  #print gene info
  genes=NULL
  for (i in 1:length(twasreg))
  {
    genes=unique(c(genes,names(twasreg[[i]])))
  }
  if (length(genes)>0)
  {
    for (i in 1:length(genes))
    {
      gene=genes[i]
      idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
      idx2=which(rownames(genemodel)==gene)
      print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                   ",R2=",scinumber(genemodel$r2[idx2]),",#SNP=",genemodel$numselectedsnp[idx2]))
    }
  }
  
  idx=which(rownames(skatres) %in% genesh)
  print("skat-TWAS result:")
  # tmp=p.adjust(skatres$BEA_p[idx],method="fdr")
  # print(paste0("BE/EA,fdr:",min(tmp)))
  twasskat=compute_fwer_fdr(dat=skatres[idx,])
  print_compute_fwer_fdr(twasskat)
  if (plotopt)
  {
    qqplot(skatres$BE_p[idx],main="BE")
    qqplot(skatres$EA_p[idx],main="EA")
    qqplot(skatres$BEA_p[idx],main="BEA") #EA vs BE
    qqplot(skatres$BEEA_p[idx],main="BEEA")
  }
  
  genes=NULL
  for (i in 1:length(twasskat))
  {
    genes=unique(c(genes,names(twasskat[[i]])))
  }
  if (length(genes)>0)
  {
    for (i in 1:length(genes))
    {
      gene=genes[i]
      idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
      idx2=which(rownames(genemodel)==gene)
      print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                   ",R2=",scinumber(genemodel$r2[idx2]),",#SNP=",genemodel$numselectedsnp[idx2]))
    }
  }
  
  return(list(twasreg=twasreg,twasskat=twasskat))
}

print_genes=function(genes=NULL,genemodel=res_min_code)
{
  for (i in 1:length(genes))
  {
    gene=genes[i]
    idx1=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
    idx2=which(rownames(genemodel)==gene)
    print(paste0(gene," ",gtexv8_ge_anno$Chromosome[idx1],":",gtexv8_ge_anno$start[idx1],
                 ",R2=",scinumber(genemodel$r2[idx2]),",#SNP=",genemodel$numselectedsnp[idx2]))
    
  }
}

check_organ=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC",
                     outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)
{
  print(outfolder)
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  #check gene models
  check_model(mod=res_min_code)
  check_model(mod=res_min_code,r2cut=0.01)
  #heritability
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  print("heritability:")
  check_model(mod=res_min_code[rownames(res_min_code) %in% rownames(heritability)[heritability$V3<0.05],])
  print(quantile(heritability$V1[rownames(heritability) %in% rownames(res_min_code) & heritability$V3<0.05]))
  load(paste0(outfolder,"/bca_assoc.RData"))
  load(paste0(outfolder,"/skat_res.RData"))
  #if (!exists("skat_min2"))
  #{
    skat_min2=skat_min2_pc6
  #}
  colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p")
  
  res1=res2=res3=NULL
  if (opt=="pc4")
  {
    assoc_min_code=assoc_min[rownames(assoc_min) %in% proteingenes,]
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    res1=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt)
    res2=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt,r2cutoff=0.05)
    res3=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt,r2cutoff=0.1)
  }
  if (opt=="pc6")
  {
    assoc_min_code=assoc_min_pc6[rownames(assoc_min_pc6) %in% proteingenes,]
    colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p")
    skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]
    #res1=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt)
    res1=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt,r2cutoff=0.01)
    # res2=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt,r2cutoff=0.05)
    # res3=check_res(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,plotopt = plotopt,r2cutoff=0.1)
    #res4=check_res_heritability(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,heritability = heritability,plotopt = plotopt)
  }
  if (opt=="heritability")
  {
    res1=check_res_heritability(genemodel=res_min_code,assocres=assoc_min_code,skatres=skat_min2_code,heritability = heritability,plotopt = plotopt)
  }
  #return(list(res_min_code=res_min_code,assoc_min_code=assoc_min_code,skat_min2_code=skat_min2_code,res1=res1,res2=res2,res3=res3))
  #return(list(res_min_code=res_min_code,assoc_min_code=assoc_min_code,skat_min2_code=skat_min2_code,res1=res1))
}
#GTEx junction------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"
# junction_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F) #this is the result before using modified removehighcorr
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar"
# junction_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
junction_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

#GTEx stomach------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_June11"
# stomach_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor_allcovar"
# stomach_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
stomach_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

#GTEx blood------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
# blood_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8blooddata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor_allcovar"
# blood_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
blood_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

#GTEx mucosa------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_June11"
# mucosa_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor_allcovar"
# mucosa_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
mucosa_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

#GTEx muscularis------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_June11"
# muscularis_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8muscularisdata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor_allcovar"
# muscularis_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
muscularis_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

#GTEx adipose------
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_June11"
# adipose_old_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8adiposedata_ambiguous_TPM_for_prediction.RData")
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor"
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar"
# adipose_new_pc4=check_organ(outfolder=outfolder,opt="pc4",plotopt=F)
adipose_new_pc6=check_organ(outfolder=outfolder,outfolder1="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor",opt="pc6",plotopt=T)

plot_2pvalues=function(pvalue1=junction_new_pc4$assoc_min_code,pvalue2=junction_new_pc6$assoc_min_code,coln=1,opt="log10",xlab="",ylab="",main="")
{
  p1=pvalue1[,coln]
  names(p1)=rownames(pvalue1)
  p2=pvalue2[,coln]
  names(p2)=rownames(pvalue2)
  comgenes=intersect(names(p1),names(p2))
  idx1=match(comgenes,names(p1))
  idx2=match(comgenes,names(p2))
  p1=p1[idx1]
  p2=p2[idx2]
  par(mar=c(6,6,2,1))
  #par(mfrow=c(1,1))
  if (opt=="log10")
  {
    plot(-log10(p1),-log10(p2),xlab=paste0(xlab,": -log10(p-value)"),ylab=paste0(ylab,": -log10(p-value)"),main=main,cex.axis=1.3,cex.lab=1.3)
  }else
  {
    plot(p1,p2,xlab=xlab,ylab=ylab,main=main,cex.axis=1.3,cex.lab=1.3)
  }
  abline(0,1,col="red")
}

#regular twas
#rown=2,coln=1;rown=4,coln=3,rown=6,coln=5,rown=8,coln=7
#skat twas
#rown=2,coln=1;rown=4,coln=2,rown=6,coln=3,rown=8,coln=4
compres=function(dat,rown=2,coln=1,pre1,pre2,source1,source2)
{
  tmp1=unlist(strsplit(dat[rown,1],";"))
  tmp2=unlist(strsplit(dat[rown,2],";"))
  tmp3=unique(c(tmp1,tmp2))
  res=NULL
  if (length(tmp3)>0)
  {
    res=data.frame(matrix(nrow=length(tmp3),ncol=6))
    rownames(res)=tmp3
    colnames(res)=c(pre1,pre2,"pvalue1","pvalue2","chr","position")
    res[,1]=res[,2]=F
    for (i in 1:length(tmp3))
    {
      gene=tmp3[i]
      if (gene %in% tmp1)
      {
        res[i,1]=T
      }
      if (gene %in% tmp2)
      {
        res[i,2]=T
      }
      idx=which(rownames(source1)==gene)
      res$pvalue1[i]=source1[idx,coln]
      idx=which(rownames(source2)==gene)
      res$pvalue2[i]=source2[idx,coln]
      idx=which(gtexv8_ge_anno$Symbol==gene & gtexv8_ge_anno$V3=="gene")[1]
      res$chr[i]=gtexv8_ge_anno$Chromosome[idx]
      res$position[i]=gtexv8_ge_anno$start[idx]
    }
  }
  if (!is.null(res))
  {
    tmp=rownames(res)[res[,1] & res[,2]]
    if (length(tmp)>0)
    {
      print(paste0("genes in ",pre1," and ",pre2," (n=",length(tmp),") : ",paste0(tmp,collapse = " ")))
    }
    tmp=rownames(res)[res[,1] & !res[,2]]
    if (length(tmp)>0)
    {
      print(paste0("genes only in ",pre1," (n=",length(tmp),") : ",paste0(tmp,collapse = " ")))
    }
    tmp=rownames(res)[!res[,1] & res[,2]]
    if (length(tmp)>0)
    {
      print(paste0("genes only in ",pre2," (n=",length(tmp),") : ",paste0(tmp,collapse = " ")))
    }
  }
  return(res)
  
}
#regular twas:colns=c(1,3,5,7)
# compare_fwer_fdr=function(dat1=junction_new_pc4$res1$twasskat,dat2=junction_new_pc6$res2$twasskat,
#                  source1=junction_new_pc4$skat_min2_code,source2=junction_new_pc6$skat_min2_code,
#                  colns=c(1,2,3,4))
compare_fwer_fdr=function(dat1=NULL,dat2=NULL,
                          source1=NULL,source2=NULL,
                          colns=c(1,2,3,4),pre1=NULL,pre2=NULL)
{
  dat=data.frame(matrix(nrow=length(dat1),ncol=2))
  rownames(dat)=names(dat1)
  for(i in 1:nrow(dat))
  {
    dat[i,1]=paste0(names(dat1[[i]]),collapse=";")
    dat[i,2]=paste0(names(dat2[[i]]),collapse=";")
  }
  #only compare fwer
  print("BE:")
  BE_fwer=compres(dat=dat,rown=2,coln=colns[1],pre1=pre1,pre2=pre2,source1=source1,source2=source2)
  print("EA:")
  EA_fwer=compres(dat=dat,rown=4,coln=colns[2],pre1=pre1,pre2=pre2,source1=source1,source2=source2)
  print("BEA:")
  BEA_fwer=compres(dat=dat,rown=6,coln=colns[3],pre1=pre1,pre2=pre2,source1=source1,source2=source2)
  print("BEEA:")
  BEEA_fwer=compres(dat=dat,rown=8,coln=colns[4],pre1=pre1,pre2=pre2,source1=source1,source2=source2)
  
  return(list(BE_fwer=BE_fwer,EA_fwer=EA_fwer,BEA_fwer=BEA_fwer,BEEA_fwer=BEEA_fwer))
  
}

compare_2res=function(res1=junction_new_pc4,res2=junction_new_pc6,pre1="PC4",pre2="PC6")
{
  #compare regular twas
  par(mfrow=c(2,2))
  plot_2pvalues(pvalue1=res1$assoc_min_code,pvalue2=res2$assoc_min_code,coln = 1,main="BE",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$assoc_min_code,pvalue2=res2$assoc_min_code,coln = 3,main="EA",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$assoc_min_code,pvalue2=res2$assoc_min_code,coln = 5,main="BEA",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$assoc_min_code,pvalue2=res2$assoc_min_code,coln = 7,main="BEEA",xlab=pre1,ylab=pre2)
  print("compare twasreg:")
  compare_twasreg=compare_fwer_fdr(dat1=res1$res1$twasreg,dat2=res2$res1$twasreg,
                                   source1=res1$assoc_min_code,source2=res2$assoc_min_code,
                                   colns=c(1,3,5,7),pre1=pre1,pre2=pre2)
  #compare skat twas
  par(mfrow=c(2,2))
  plot_2pvalues(pvalue1=res1$skat_min2_code,pvalue2=res2$skat_min2_code,coln = 1,main="BE",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$skat_min2_code,pvalue2=res2$skat_min2_code,coln = 2,main="EA",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$skat_min2_code,pvalue2=res2$skat_min2_code,coln = 3,main="BEA",xlab=pre1,ylab=pre2)
  plot_2pvalues(pvalue1=res1$skat_min2_code,pvalue2=res2$skat_min2_code,coln = 4,main="BEEA",xlab=pre1,ylab=pre2)
  print("compare twasskat:")
  compare_twasskat=compare_fwer_fdr(dat1=res1$res1$twasskat,dat2=res2$res1$twasskat,
                                    source1=res1$skat_min2_code,source2=res2$skat_min2_code,
                                    colns=c(1,2,3,4),pre1=pre1,pre2=pre2)
  return(list(compare_twasreg=compare_twasreg,compare_twasskat=compare_twasskat))
  
}

compare_2res()
compare_2res(res1=mucosa_new_pc4,res2=mucosa_new_pc6)
compare_2res(res1=stomach_new_pc4,res2=stomach_new_pc6)
compare_2res(res1=blood_new_pc4,res2=blood_new_pc6)
compare_2res(res1=muscularis_new_pc4,res2=muscularis_new_pc6)
compare_2res(res1=adipose_new_pc4,res2=adipose_new_pc6)

compare_2res(res1=junction_new_pc6,res2=junction_old_pc4,pre1="NewPC6",pre2="Old")
compare_2res(res1=mucosa_new_pc6,res2=mucosa_old_pc4,pre1="NewPC6",pre2="Old")
compare_2res(res1=stomach_new_pc6,res2=stomach_old_pc4,pre1="NewPC6",pre2="Old")
compare_2res(res1=blood_new_pc6,res2=blood_old_pc4,pre1="NewPC6",pre2="Old")
compare_2res(res1=muscularis_new_pc6,res2=muscularis_old_pc4,pre1="NewPC6",pre2="Old")
compare_2res(res1=adipose_new_pc6,res2=adipose_old_pc4,pre1="NewPC6",pre2="Old")


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
outfolder2="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC"
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


check_skatres=function(dat=skat_min1,opt="noplot")
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
dong23snp=read.table("../data/Dong23snp.txt",header=T,sep="\t",stringsAsFactors = F)
extra3snp=read.table("../result/Extra3SNPs.txt",header = T)
dong26snp=rbind(dong23snp[,1:3],extra3snp)
#write.table(dong26snp,file="../data/Dong26snp.txt",row.names = F,col.names = T,quote=F,sep="\t")
tmp=read.table("../result/Dong26SNPs_postion.txt")
tmp$pos=NA
for (i in 1:nrow(tmp))
{
  tmp$pos[i]=unlist(strsplit(tmp$V1[i],"-"))[2]
}
all(tmp$pos==dong26snp$Position) #T
library(rtracklayer)
library(GenomicRanges)
chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg19ToHg38.over.chain")

gr_dong26snp=GRanges(seqnames = paste0("chr",dong26snp$Chr),ranges=IRanges(start=dong26snp$Position,width=1))
tmp=liftOver(gr_dong26snp,chain)
dong26snp$pos38=NA
for (i in 1:length(tmp))
{
  tmp1=tmp[i]
  if (length(tmp1)>0)
    dong26snp$pos38[i]=as.numeric(start(tmp1))
}
#use hg38
gr_dong26snp=GRanges(seqnames = paste0("chr",dong26snp$Chr),ranges=IRanges(start=dong26snp$pos38,width=1))
get_genelist=function(tissues=c("junction","stomach","mucosa","muscularis","blood","adipose"),
                      r2cutoff=0,opt="old")
{
  genetable=NULL
  for (i in 1:length(tissues))
  {
    tissue=tissues[i]
    if (opt=="old")
    {
      outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_June11")
      if (tissue=="junction")
      {
        outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11")
      }
    }
    if (opt %in% c("newpc4","newpc6"))
    {
      #outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_HRC")
      outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_HRC_MAF005_rmhighcor")
    }
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/skat_res.RData"))
    if (opt %in% c("old","newpc4"))
    {
      colnames(skat_min2)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
      skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    }
    if (opt=="newpc6")
    {
      colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
      skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]
    }
  
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
            generes=data.frame(gene=gene,type=comparison,cutoff="FWER",r2cutoff=r2cutoff,
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
            generes=data.frame(gene=gene,type=comparison,cutoff="FDR",r2cutoff=r2cutoff,
                               tissue=tissue,R2=res_min$r2[idx1],numsnp=res_min$numselectedsnp[idx1],
                               chr=gtexv8_ge_anno$Chromosome[idx2],position=gtexv8_ge_anno$start[idx2],stringsAsFactors = F)
            genetable=rbind(genetable,generes)
          }
        }
      }
    }
  }
  gr_genetable=GRanges(seqnames = genetable$chr,ranges=IRanges(start=genetable$position,width=1))
  genetable$dist_gwassnp=genetable$gwassnp=NA
  for (i in 1:nrow(genetable))
  {
    tmp=distance(gr_genetable[i],gr_dong26snp)
    if (any(!is.na(tmp)))
    {
      idx=which.min(tmp)
      genetable$gwassnp[i]=dong26snp$SNP[idx]
      genetable$dist_gwassnp[i]=tmp[idx]
    }
  }
  idx=order(genetable$type,genetable$tissue,genetable$dist_gwassnp)
  genetable=genetable[idx,]
  idx=which(is.na(genetable$dist_gwassnp) | genetable$dist_gwassnp>10000000)
  genetable$new=F
  genetable$new[idx]=T
  print(paste0("total genes:",nrow(genetable),"; ","new genes:",sum(genetable$new)))
  return(genetable)
}
genelist_old=get_genelist()
genelist05_old=get_genelist(r2cutoff = 0.05)
genelist1_old=get_genelist(r2cutoff = 0.1)

genelist_newpc4=get_genelist(opt="newpc4")
genelist05_newpc4=get_genelist(r2cutoff = 0.05,opt="newpc4")
genelist1_newpc4=get_genelist(r2cutoff = 0.1,opt="newpc4")

genelist_newpc6=get_genelist(opt="newpc6")
genelist01_newpc6=get_genelist(r2cutoff = 0.01, opt="newpc6")
genelist05_newpc6=get_genelist(r2cutoff = 0.05,opt="newpc6")
genelist1_newpc6=get_genelist(r2cutoff = 0.1,opt="newpc6")

# #list result for each gene
# get_genetable=function(genelist)
# {
#   genes=unique(genelist$gene)
#   genetable=data.frame(gene=genes,tissue=NA,type=NA,
#                        chr=NA,position=NA,numtissue=NA,numtype=NA,stringsAsFactors = F)
#   for (i in 1:nrow(genetable))
#   {
#     idx=which(genelist$gene==genes[i])
#     genetable$chr[i]=genelist$chr[idx[1]]
#     genetable$position[i]=genelist$position[idx[1]]
#     genetable$tissue[i]=paste0(genelist$tissue[idx],collapse = "|")
#     genetable$type[i]=paste0(genelist$type[idx],collapse = "|")
#     genetable$numtissue[i]=length(unique(genelist$tissue[idx]))
#     genetable$numtype[i]=length(unique(genelist$type[idx]))
#   }
#   return(genetable)
# }

get_genelist_heritability=function(tissues=c("junction","stomach","mucosa","muscularis","blood","adipose"))
{
  genetable=NULL
  
  for (i in 1:length(tissues))
  {
    tissue=tissues[i]
    outfolder_heritability=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_HRC_MAF005_rmhighcor")
    heritfile=paste0(outfolder_heritability,"/heritability_","1grm",".txt")
    heritability=read.table(heritfile,header = T)
    
    outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",tissue,"_HRC_MAF005_rmhighcor")
    
    load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
    load(paste0(outfolder,"/skat_res.RData")) #skat_min2,skat_min2_pc6

    colnames(skat_min2_pc6)=c("BE_p","EA_p","BEA_p","BEEA_p") #weights.beta=c(1,1), the one picked before
    skat_min2_code=skat_min2_pc6[rownames(skat_min2_pc6) %in% proteingenes,]
  
    idx=match(rownames(heritability),rownames(res_min))
    if (sum(is.na(idx))>0) warning("Some genes in heritability are missing")
    res_min$heritability=res_min$heritability_p=NA
    res_min$heritability[idx]=heritability$V1
    res_min$heritability_p[idx]=heritability$V3
    genesr2=rownames(res_min)[which(res_min$heritability_p<0.05)]
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
            generes=data.frame(gene=gene,type=comparison,cutoff="FWER",
                               tissue=tissue,R2=res_min$r2[idx1],numsnp=res_min$numselectedsnp[idx1],heritability=res_min$heritability[idx1],heritability_pvalue=res_min$heritability_p[idx1],
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
            generes=data.frame(gene=gene,type=comparison,cutoff="FDR",
                               tissue=tissue,R2=res_min$r2[idx1],numsnp=res_min$numselectedsnp[idx1],heritability=res_min$heritability[idx1],heritability_pvalue=res_min$heritability_p[idx1],
                               chr=gtexv8_ge_anno$Chromosome[idx2],position=gtexv8_ge_anno$start[idx2],stringsAsFactors = F)
            genetable=rbind(genetable,generes)
          }
        }
      }
    }
  }
  gr_genetable=GRanges(seqnames = genetable$chr,ranges=IRanges(start=genetable$position,width=1))
  genetable$dist_gwassnp=genetable$gwassnp=NA
  for (i in 1:nrow(genetable))
  {
    tmp=distance(gr_genetable[i],gr_dong26snp)
    if (any(!is.na(tmp)))
    {
      idx=which.min(tmp)
      genetable$gwassnp[i]=dong26snp$SNP[idx]
      genetable$dist_gwassnp[i]=tmp[idx]
    }
  }
  idx=order(genetable$type,genetable$tissue,genetable$dist_gwassnp)
  genetable=genetable[idx,]
  idx=which(is.na(genetable$dist_gwassnp) | genetable$dist_gwassnp>10000000)
  genetable$new=F
  genetable$new[idx]=T
  print(paste0("total genes:",nrow(genetable),"; ","new genes:",sum(genetable$new)))
  return(genetable)
}

genelist_heritability=get_genelist_heritability()
write.csv(genelist_heritability,file="../result/SKAT_TWAS_newpc6_genes_heritability.csv",row.names = F)
get_genetable1=function(genelist,cutoffopt="fwer")
{
  if (cutoffopt=="fwer")
  {
    genelist=genelist[genelist$cutoff=="FWER",]
  }
  
  alltypes=unique(genelist$type)
  allgenetable=NULL
  for (i in 1:length(alltypes))
  {
    idx=which(genelist$type==alltypes[i])
    genes=unique(genelist$gene[idx])
    genetable=data.frame(gene=genes,tissue=NA,type=NA,
                         chr=NA,position=NA,numtissue=NA,stringsAsFactors = F)
    genelist1=genelist[idx,]
    for (j in 1:nrow(genetable))
    {
      idx=which(genelist1$gene==genes[j])
      genetable$chr[j]=genelist1$chr[idx[1]]
      genetable$position[j]=genelist1$position[idx[1]]
      genetable$tissue[j]=paste0(genelist1$tissue[idx],collapse = "|")
      genetable$type[j]=alltypes[i]
      genetable$numtissue[j]=length(unique(genelist1$tissue[idx]))
    }
    allgenetable=rbind(allgenetable,genetable)
  }
  allgenetable$new=allgenetable$dist_gwassnp=allgenetable$gwassnp=NA
  for (i in 1:nrow(allgenetable))
  {
    idx=which(genelist$gene==allgenetable$gene[i])
    allgenetable$gwassnp[i]=genelist$gwassnp[idx[1]]
    allgenetable$dist_gwassnp[i]=genelist$dist_gwassnp[idx[1]]
    allgenetable$new[i]=genelist$new[idx[1]]
  }
  print(paste0("total genes:",nrow(allgenetable),"; ","new genes:",sum(allgenetable$new)))
  print(table(allgenetable$type))
  return(allgenetable)
}

genetable_old=get_genetable1(genelist_old)
# [1] "total genes:40; new genes:5"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 25       13        2 
genetable05_old=get_genetable1(genelist05_old)
# [1] "total genes:18; new genes:2"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 10        5        3 
genetable1_old=get_genetable1(genelist1_old)
# [1] "total genes:8; new genes:1"
# 
# BEEAvsCO   BEvsCO   EAvsBE 
# 4        3        1 
write.csv(genelist_old,file="../result/SKAT_TWAS_genes_r2_0_.csv",row.names = F)
write.csv(genelist05_old,file="../result/SKAT_TWAS_genes_r2_0.05_.csv",row.names = F)
write.csv(genelist1_old,file="../result/SKAT_TWAS_genes_r2_0.1_.csv",row.names = F)
write.csv(genetable_old,file="../result/SKAT_TWAS_genetable_r2_0_.csv",row.names = F)
write.csv(genetable05_old,file="../result/SKAT_TWAS_genetable_r2_0.05_.csv",row.names = F)
write.csv(genetable1_old,file="../result/SKAT_TWAS_genetable_r2_0.1_.csv",row.names = F)

genetable_newpc4=get_genetable1(genelist_newpc4)
# [1] "total genes:27; new genes:5"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 20        3        4 
genetable05_newpc4=get_genetable1(genelist05_newpc4)
# [1] "total genes:13; new genes:3"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 9        2        2 
genetable1_newpc4=get_genetable1(genelist1_newpc4)
# [1] "total genes:7; new genes:1"
# 
# BEEAvsCO   EAvsCO 
# 4        3 
write.csv(genelist_newpc4,file="../result/SKAT_TWAS_newpc4_genes_r2_0_.csv",row.names = F)
write.csv(genelist05_newpc4,file="../result/SKAT_newpc4_TWAS_genes_r2_0.05_.csv",row.names = F)
write.csv(genelist1_newpc4,file="../result/SKAT_newpc4_TWAS_genes_r2_0.1_.csv",row.names = F)
write.csv(genetable_newpc4,file="../result/SKAT_TWAS_newpc4_genetable_r2_0_.csv",row.names = F)
write.csv(genetable05_newpc4,file="../result/SKAT_TWAS_newpc4_genetable_r2_0.05_.csv",row.names = F)
write.csv(genetable1_newpc4,file="../result/SKAT_TWAS_newpc4_genetable_r2_0.1_.csv",row.names = F)

genetable_newpc6=get_genetable1(genelist_newpc6)
# [1] "total genes:30; new genes:4"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 21        6        3 
tmp1=table(genetable_newpc6$chr)
# chr15 chr16 chr19  chr2 chr22  chr6 
# 2     2    21     1     1     3 
genetable05_newpc6=get_genetable1(genelist05_newpc6)
# [1] "total genes:13; new genes:2"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 10        2        1 
tmp2=table(genetable05_newpc6$chr)
# chr14 chr15 chr16 chr19  chr2 
# 2     3     2     5     1
genetable1_newpc6=get_genetable1(genelist1_newpc6)
# [1] "total genes:8; new genes:1"
# 
# BEEAvsCO   BEvsCO   EAvsCO 
# 4        2        2 
tmp3=table(genetable1_newpc6$chr)
# chr11 chr15 chr16 chr19  chr2 
# 1     2     2     2     1
write.csv(genelist_newpc6,file="../result/SKAT_TWAS_newpc6_genes_r2_0_.csv",row.names = F)
write.csv(genelist01_newpc6,file="../result/SKAT_TWAS_newpc6_genes_r2_01_.csv",row.names = F)
write.csv(genelist05_newpc6,file="../result/SKAT_newpc6_TWAS_genes_r2_0.05_.csv",row.names = F)
write.csv(genelist1_newpc6,file="../result/SKAT_newpc6_TWAS_genes_r2_0.1_.csv",row.names = F)
write.csv(genetable_newpc6,file="../result/SKAT_TWAS_newpc6_genetable_r2_0_.csv",row.names = F)
write.csv(genetable05_newpc6,file="../result/SKAT_TWAS_newpc6_genetable_r2_0.05_.csv",row.names = F)
write.csv(genetable1_newpc6,file="../result/SKAT_TWAS_newpc6_genetable_r2_0.1_.csv",row.names = F)

genetable_heritability=get_genetable1(genelist_heritability)
# [1] "total genes:18; new genes:4"
# BEEAvsCO   BEvsCO 
# 13        5 
write.csv(genetable_heritability,file="../result/SKAT_TWAS_newpc6_genetable_heritability.csv",row.names = F)
chrs=unique(c(names(tmp1),names(tmp2),names(tmp3)))
#[1] "chr15" "chr16" "chr19" "chr2"  "chr22" "chr6"  "chr14" "chr11"
chrs=gsub("chr","",chrs)
usedsnps=dong26snp$SNP[dong26snp$Chr %in% chrs]
read_gwassnp=function(genotypefolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/")
{
  gwasdat=read.table(paste0(genotypefolder,"Dong26SNPs_genotype.txt"))
  idx=match(gwasdat$V2,dong26snp$Position)
  rownames(gwasdat)=dong26snp$SNP[idx]
  gwasdat=gwasdat[rownames(gwasdat)%in%usedsnps,]
  tmp=gwasdat
  tmp1=data.frame(matrix(NA,nrow=ncol(gwasdat)-9+4,ncol=nrow(gwasdat)))
  colnames(tmp1)=rownames(gwasdat)
  tmp1[1,]=tmp[,1]
  tmp1[2,]=tmp[,2]
  tmp1[3,]=tmp[,4]
  tmp1[4,]=tmp[,5]
  
  for (i in 1:nrow(gwasdat))
  {
    idx=which(grepl("0|0",tmp[i,10:ncol(tmp)],fixed = T))
    tmp1[idx+4,i]=0
    idx=which(grepl("0|1",tmp[i,10:ncol(tmp)],fixed = T))
    tmp1[idx+4,i]=1
    idx=which(grepl("1|0",tmp[i,10:ncol(tmp)],fixed = T))
    tmp1[idx+4,i]=1
    
    idx=which(grepl("1|1",tmp[i,10:ncol(tmp)],fixed = T))
    tmp1[idx+4,i]=2
  }
  #check sample name
  tmp=data.table::fread(paste0(genotypefolder,"chr6.dose.vcf.gz.head"),skip = 18)
  tmp=colnames(tmp)[10:ncol(tmp)]
  tmp=gsub("SEP","",tmp)
  samples=strsplit(tmp,"_") #use localid
  samplenames=sapply(1:length(samples),function(x){
    tmp=samples[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
  all(samplenames %in% colnames(bcagenotype))
  rownames(tmp1)=c("chr","position","ref","alt",samplenames)
  gwasdat=tmp1[5:nrow(tmp1),]
  gwassnp=tmp1[1:4,]
  return(list(gwasdat=gwasdat,gwassnp=gwassnp))
}

gwasdat_beacon=read_gwassnp()
gwasdat_cambridge=read_gwassnp(genotypefolder ="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/")
all(colnames(gwasdat_beacon$gwasdat)==colnames(gwasdat_cambridge$gwasdat)) #T
allgwasdat=rbind(gwasdat_beacon$gwasdat,gwasdat_cambridge$gwasdat)
all(colnames(bcagenotype) %in% rownames(allgwasdat)) #T
idx=match(colnames(bcagenotype),rownames(allgwasdat))
allgwasdat=allgwasdat[idx,]
write.table(allgwasdat,file="../result/dong26snp_addcontrols_genotypedat.txt",sep="\t",quote=F)
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC"
load(paste0(outfolder,"/bca_extractgenotype.RData"))
sum(tmp %in% colnames(bcagenotype))

#add pvaules/fwer and check if SKAT gene is significant after adjusting for gwas snps
update_genelist=function(genelist=genelist_newpc6,opt="fwer",adjustopt="PC6",r2cutoff=0)
{
  if (opt=="fwer")
  {
    genelist=genelist[genelist$cutoff=="FWER",]
  }
  genelist$fwer_addgwas=genelist$pvalue_addgwas=genelist$BEEA_fwer=genelist$BEEA_pvalue=genelist$BEA_fwer=genelist$BEA_pvalue=genelist$EA_fwer=genelist$EA_pvalue=genelist$BE_fwer=genelist$BE_pvalue=NA
  for (i in 1:nrow(genelist))
  {
    if (genelist$tissue[i]=="adipose")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC"
    }
    if (genelist$tissue[i]=="blood")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC"
    }
    if (genelist$tissue[i]=="junction")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC"
    }
    if (genelist$tissue[i]=="mucosa")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC"
    }
    if (genelist$tissue[i]=="muscularis")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC"
    }
    if (genelist$tissue[i]=="stomach")
    {
      outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC"
    }
    #work on pvalues
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/skat_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    idx=which(rownames(skat_min2_code)==genelist$gene[i])
    genelist$BEEA_pvalue[i]=skat_min2_code$BEEA_pvalue[idx]
    genelist$BEEA_fwer[i]=skat_min2_code$BEEA_fwer[idx]
    genelist$BEA_pvalue[i]=skat_min2_code$BEA_pvalue[idx]
    genelist$BEA_fwer[i]=skat_min2_code$BEA_fwer[idx]
    genelist$EA_pvalue[i]=skat_min2_code$EA_pvalue[idx]
    genelist$EA_fwer[i]=skat_min2_code$EA_fwer[idx]
    genelist$BE_pvalue[i]=skat_min2_code$BE_pvalue[idx]
    genelist$BE_fwer[i]=skat_min2_code$BE_fwer[idx]
    
    #work on gwas snp adjusted pvalues
    load(paste0(outfolder,"/skat_gwas_res.RData"))
    if (adjustopt=="PC6")
    {
      skat_min2=skat_min2_pc6
    }
    colnames(skat_min2)=c("BE_pvalue","EA_pvalue","BEA_pvalue","BEEA_pvalue")
    skat_min2_code=skat_min2[rownames(skat_min2) %in% proteingenes,]
    genesr2=rownames(res_min)[which(res_min$r2>r2cutoff)]
    skat_min2_code=skat_min2_code[rownames(skat_min2_code) %in% genesr2,]
    skat_min2_code$BEEA_fwer=skat_min2_code$BEA_fwer=skat_min2_code$EA_fwer=skat_min2_code$BE_fwer=NA
    skat_min2_code$BEEA_fwer=p.adjust(skat_min2_code$BEEA_pvalue,method="bonferroni")
    skat_min2_code$BEA_fwer=p.adjust(skat_min2_code$BEA_pvalue,method="bonferroni")
    skat_min2_code$EA_fwer=p.adjust(skat_min2_code$EA_pvalue,method="bonferroni")
    skat_min2_code$BE_fwer=p.adjust(skat_min2_code$BE_pvalue,method="bonferroni")
    idx=which(rownames(skat_min2_code)==genelist$gene[i])
    if (genelist$type[i]=="BEEAvsCO")
    {
      genelist$pvalue_addgwas[i]=skat_min2_code$BEEA_pvalue[idx]
      genelist$fwer_addgwas[i]=skat_min2_code$BEEA_fwer[idx]
    }
    if (genelist$type[i]=="EAvsBE")
    {
      genelist$pvalue_addgwas[i]=skat_min2_code$BEA_pvalue[idx]
      genelist$fwer_addgwas[i]=skat_min2_code$BEA_fwer[idx]
    }
    if (genelist$type[i]=="EAvsCO")
    {
      genelist$pvalue_addgwas[i]=skat_min2_code$EA_pvalue[idx]
      genelist$fwer_addgwas[i]=skat_min2_code$EA_fwer[idx]
    }
    if (genelist$type[i]=="BEvsCO")
    {
      genelist$pvalue_addgwas[i]=skat_min2_code$BE_pvalue[idx]
      genelist$fwer_addgwas[i]=skat_min2_code$BE_fwer[idx]
    }
  }
  return(genelist)
}

genelist_newpc6_full=update_genelist()
dim(genelist_newpc6_full) #36 21
length(unique(genelist_newpc6_full$gene))
genelist05_newpc6_full=update_genelist(genelist=genelist05_newpc6,r2cutoff = 0.05)
dim(genelist05_newpc6_full) #17 21
genelist1_newpc6_full=update_genelist(genelist=genelist1_newpc6,r2cutoff = 0.1)
dim(genelist1_newpc6_full) #11 21
new_genelist=function(oldgenelist=genelist_newpc6_full,genelist=genelist05_newpc6_full)
{
  tmp1=paste0(oldgenelist$gene,"_",oldgenelist$type,"_",oldgenelist$cutoff,"_",oldgenelist$tissue)
  tmp2=paste0(genelist$gene,"_",genelist$type,"_",genelist$cutoff,"_",genelist$tissue)
  newgenelist=NULL
  idx=which(!tmp2 %in% tmp1)
  if (length(idx)>0)
  {
    newgenelist=genelist[idx,]
  }
  return(newgenelist)
}
genelist05_newpc6_full_new=new_genelist()
dim(genelist05_newpc6_full_new) #5 21
genelist1_newpc6_full_new=new_genelist(oldgenelist = rbind(genelist_newpc6_full,genelist05_newpc6_full_new),genelist=genelist1_newpc6_full)
dim(genelist1_newpc6_full_new) #4 21
#how many gene_type
tmp1=paste0(genelist_newpc6_full$gene,"_",genelist_newpc6_full$type,"_",genelist_newpc6_full$new)
tmp2=paste0(genelist05_newpc6_full_new$gene,"_",genelist05_newpc6_full_new$type,"_",genelist05_newpc6_full_new$new)
tmp3=paste0(genelist1_newpc6_full_new$gene,"_",genelist1_newpc6_full_new$type,"_",genelist1_newpc6_full_new$new)
tmp11=paste0(genelist_newpc6_full$gene,"_",genelist_newpc6_full$type,"_",genelist_newpc6_full$r2cutoff,"_",genelist_newpc6_full$new)
tmp22=paste0(genelist05_newpc6_full_new$gene,"_",genelist05_newpc6_full_new$type,"_",genelist05_newpc6_full_new$r2cutoff,"_",genelist05_newpc6_full_new$new)
tmp33=paste0(genelist1_newpc6_full_new$gene,"_",genelist1_newpc6_full_new$type,"_",genelist1_newpc6_full_new$r2cutoff,"_",genelist1_newpc6_full_new$new)
tmp=unique(c(tmp1,tmp2,tmp3))
idx1=match(intersect(tmp,tmp1),tmp1)
idx2=match(intersect(tmp,tmp2),tmp2)
idx3=match(intersect(tmp,tmp3),tmp3)
tmp=c(tmp11[idx1],tmp22[idx2],tmp33[idx3])
uniqe_genelist=data.frame(gene=rep(NA,length(tmp)),type=NA,r2cutoff=NA,new=NA)
for (i in 1:length(tmp))
{
  tmp1=unlist(strsplit(tmp[i],"_"))
  uniqe_genelist$gene[i]=tmp1[1]
  uniqe_genelist$type[i]=tmp1[2]
  uniqe_genelist$r2cutoff[i]=tmp1[3]
  uniqe_genelist$new[i]=tmp1[4]
}
tmp=paste0(uniqe_genelist$gene,"_",uniqe_genelist$type)
which(duplicated(tmp))
uniqe_genelist=uniqe_genelist[!duplicated(tmp),]
length(unique(uniqe_genelist$gene)) #25
idx=order(uniqe_genelist$type,uniqe_genelist$new,decreasing = T)
uniqe_genelist=uniqe_genelist[idx,]
write.csv(genelist_newpc6_full,file="../result/SKAT_TWAS_newpc6_genes_r2_0_.csv",row.names = F)
write.csv(genelist05_newpc6_full_new,file="../result/SKAT_TWAS_newpc6_genes_r2_0.05_.csv",row.names = F)
write.csv(genelist1_newpc6_full_new,file="../result/SKAT_TWAS_newpc6_genes_r2_0.1_.csv",row.names = F)

#check the overlap of eQTLs 
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_junction=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_June11"
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
snp_blood=unique(unlist(strsplit(res_min$selectedsnps,"|",fixed=T)))
sum(snp_junction %in% snp_blood) #194059
