#!/usr/bin/evn Rscript
#The code is used to predict gene expression based on GTEx model we have.


#based on mucosa models
prefix="dist500K_GTEx_mucosa_June11"
#snp,snppos
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8mucosadata_ambiguous_TPM_for_prediction.RData")
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))

#save results here
prefix="dist500K_GTEx_mucosa_beacon_Jan15"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
#extract snps from the prediction models

extract_snp2=function(dat=rbind(res_min),prefix="dist500K_GTEx_mucosa_beacon_Jan15")
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
  dat=dat[dat$glmflag==1,]
  snps=unique(unlist(strsplit(dat$selectedsnps,"|",fixed=T)))
  res=data.frame(snp=snps,chr=NA,pos=NA,stringsAsFactors = F)

  probenames=rownames(snppos)
  idx=match(res$snp,probenames)
  res$chr=snppos$chr[idx]
  res$pos=snppos$pos[idx]
  idx=order(res$chr,res$pos)
  res=res[idx,]
  res=res[order(res[,2],res[,3]),]
  tmp=paste0(res[,2],"_",res[,3])
  idx=duplicated(tmp) #multi-allel
  res=res[!idx,]
  #create regions file used to extact genotype data
  for (i in 1:length(unique(res$chr)))
  {
    chr=unique(res$chr)[i]
    idx=which(res$chr==chr)
    tmp=ceiling(length(idx)/2)
    chrres=paste0(chr,":",res$pos[idx[1:tmp]],-res$pos[idx[1:tmp]])
    chrres=paste0(chrres,collapse = "\n")
    filename=paste0(outfolder,"/prediction_snps_tabix_chr",chr,"_1.txt")
    fileCon=file(filename)
    writeLines(chrres,fileCon)
    close(fileCon)

    chrres=paste0(chr,":",res$pos[idx[(tmp+1):length(idx)]],-res$pos[idx[(tmp+1):length(idx)]])
    chrres=paste0(chrres,collapse = "\n")
    filename=paste0(outfolder,"/prediction_snps_tabix_chr",chr,"_2.txt")
    fileCon=file(filename)
    writeLines(chrres,fileCon)
    close(fileCon)
  }
  return(res)
}

#BEACON
tmp=extract_snp2()
#CAMBRIDGE
tmp=extract_snp2(prefix = "dist500K_GTEx_mucosa_cambridge_Jan15")
tmp=extract_snp2(prefix="dist500K_GTEx_mucosa_mergeqc_Jan15")

prefix="dist500K_GTEx_mucosa_beacon_Jan15"
dataset="beacondbgapcontrol_qc_1000g"

prefix="dist500K_GTEx_mucosa_cambridge_Jan15"
dataset="cambridgewtccc_qc_1000g"

prefix="dist500K_GTEx_mucosa_mergeqc_Jan15"
dataset="merge_beacon_cambridge_qc_1000g"
#run bca_extract_genotype2_v8_ambiguousSNP_addnewcontrols.sh
cmd=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/code/bca_extract_genotype2_V8_ambiguousSNP_addnewcontrols.sh ",prefix,dataset)
system(cmd,wait = T)

print("save BCA genotypte---")
library(data.table)
extractgenotype=function(prefix="dist500K_GTEx_mucosa_beacon_Jan15")
{
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
  for (i in 1:22)
  {
    cat(i,'..')
    tmp1=fread(paste0(outfolder,"/chr",i,"_select.bim"))
    tmp1=as.data.frame(tmp1)
    tmp=fread(paste0(outfolder,"/chr",i,"_select.raw"),header=T)
    tmp=as.data.frame(tmp)
    idx=duplicated(tmp$IID) #check this
    tmp=tmp[!idx,]
    rownames(tmp)=tmp$IID
    tmp=tmp[,7:ncol(tmp)]
    colnames(tmp)=paste0(tmp1[,1],":",tmp1[,4],"_",tmp1[,5],"_",tmp1[,6])
    if (i==1)
    {
      bcagenotype=tmp
      bcabim=tmp1
    }else
    {
      bcagenotype=cbind.data.frame(bcagenotype,tmp)
      bcabim=rbind(bcabim,tmp1)
    }
  }
  mycolnames=colnames(bcagenotype)
  idx=duplicated(mycolnames)
  bcagenotype=bcagenotype[,!idx]
  bcagenotype=as.data.frame(t(bcagenotype))
  bcabim=bcabim[!idx,]
  save(bcagenotype,bcabim,file=paste0(outfolder,"/bca_extractgenotype.RData"))
  Sys.time()
}

extractgenotype(prefix="dist500K_GTEx_mucosa_mergeqc_Jan15")

print("predict gene expression---")
predict_geneexp=function(i=1,modeltable=res_min)
{
  idx=modeltable$glmflag==1
  #pick the lasso selected genes
  modeltable=modeltable[idx,]
  predicted_geneexp=data.frame(matrix(NA,nrow=1,ncol=2+ncol(bcagenotype_chunk)))
  rownames(predicted_geneexp)=rownames(modeltable)[i]
  colnames(predicted_geneexp)=c("n_totalsnp","n_avaisnp",colnames(bcagenotype_chunk))
  predicted_geneexp[,1]=modeltable$numselectedsnp[i]
  selectedsnps=unlist(strsplit(modeltable$selectedsnps[i],"|",fixed=T))
  selectedcoeff=as.numeric(unlist(strsplit(modeltable$selectedsnps_coeff[i],"|",fixed=T)))
  
  #if some imputed snps need to to flipped
  
  correctedsnps=NULL
  if (length(intersect(selectedsnps,rownames(bcagenotype_chunk)))<length(selectedsnps))
  {
    missingsnps=selectedsnps[!selectedsnps %in% rownames(bcagenotype_chunk)]
    for (j in 1:length(missingsnps))
    {
      tmp=unlist(strsplit(missingsnps[j],"_"))
      tmp1=paste0(tmp[c(1,3,2)],collapse = "_")  #change the order of allele
      idx=which(rownames(bcagenotype_chunk)==tmp1)
      if (length(idx)>0)
      {
        correctedsnps=c(correctedsnps,tmp1)
        idx1=which(selectedsnps==missingsnps[j])
        selectedsnps[idx1]=tmp1 #change the snp name to make it consistent with bca
      }
      
    }
  }
  idx=match(selectedsnps,rownames(bcagenotype_chunk))
  navaisnp=sum(!is.na(idx))
  if (navaisnp>0)
  {
    predicted_geneexp[,2]=navaisnp
    idx=selectedsnps %in% rownames(bcagenotype_chunk)
    selectedsnps=selectedsnps[idx]
    selectedcoeff=selectedcoeff[idx]
    idx1=match(selectedsnps,rownames(bcagenotype_chunk))
    availmat=bcagenotype_chunk[idx1,] #create a small matrix to avoid too large memory usage
    if (length(correctedsnps)>0)
    {
      idx2=match(correctedsnps,rownames(availmat))
      availmat[idx2,]=2-availmat[idx2,]
    }
    geneexp=as.matrix(t(availmat)) %*% selectedcoeff
    predicted_geneexp[,3:ncol(predicted_geneexp)]=geneexp
  }
  return(predicted_geneexp)
}

#ml fhR
#export OMP_NUM_THREADS=1
#salloc --constraint=gizmok -n 21 -t 3-1 mpirun -n 1 R --interactive
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(res_min)
#run for each prefix
prefix="dist500K_GTEx_mucosa_cambridge_Jan15"

prefix="dist500K_GTEx_mucosa_beacon_Jan15"

prefix="dist500K_GTEx_mucosa_mergeqc_Jan15"

outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_extractgenotype.RData"))

mpi.bcast.Robj2slave(outfolder)
#mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
mpi.bcast.Robj2slave(bcabim)
mpi.bcast.Robj2slave(predict_geneexp)
#mpi.bcast.Robj2slave(snppos)
mpi_predict_geneexp=function(modeltable=res_min)
{
  modeltable=modeltable[modeltable$glmflag==1,]
  rows=1:nrow(modeltable)
  snps_bcagenotype=unlist(strsplit(rownames(bcagenotype),"_"))
  snps_bcagenotype=snps_bcagenotype[seq(1,length(snps_bcagenotype),3)]
  modeltable$cumsumsnp=cumsum(modeltable$numselectedsnp)
  nchunks=ceiling(sum(modeltable$numselectedsnp)/1500) #each chunk picks ~1500 snps
  joblables=cut(modeltable$cumsumsnp,nchunks)
  res=NULL
  # n=njobs
  # nchunks=ceiling(length(rows)/n)
  print(paste0("number of total:",nchunks))
  for (i in 1:nchunks)
  {
    if (i %% 5==0) cat(i,"..")
    #if (i %% 10==0) mpi.remote.exec(gc())
    seq=which(joblables==levels(joblables)[i])
    selectedsnps=unique(unlist(strsplit(modeltable$selectedsnps[seq],"|",fixed=T)))
    snps_selectedsnps=unlist(strsplit(selectedsnps,"_"))
    snps_selectedsnps=snps_selectedsnps[seq(1,length(snps_selectedsnps),3)]
    idx=snps_bcagenotype %in% snps_selectedsnps #extract genotypes for each chunk (use chr:pos)
    bcagenotype_chunk=bcagenotype[idx,]
    mpi.bcast.Robj2slave(bcagenotype_chunk)
    tmp=mpi.parSapply(X=seq,FUN=predict_geneexp,modeltable=modeltable,job.num=min(c(njobs,length(seq))))
    res1=as.data.frame(matrix(unlist(tmp),ncol=2+ncol(bcagenotype),byrow = T))
    if (length(unlist(tmp)) %% (2+ncol(bcagenotype)) !=0) stop(i)
    res=rbind(res,res1)
  }
  res=as.data.frame(res)
  rownames(res)=rownames(modeltable)
  colnames(res)=c("n_totalsnp","n_avaisnp",colnames(bcagenotype))
  return(res)
}
predict_min=mpi_predict_geneexp()
save(predict_min,file=paste0(outfolder,"/bca_predict_geneexp.RData"))
Sys.time()

#after got the predict expression, get the outcomes
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

prefix="dist500K_GTEx_mucosa_cambridge_Jan15"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
predict_cambridge=predict_min[!is.na(predict_min$n_avaisnp),]
colnames(predict_cambridge)=gsub("SEP","",colnames(predict_cambridge))

prefix="dist500K_GTEx_mucosa_beacon_Jan15"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
predict_beacon=predict_min[!is.na(predict_min$n_avaisnp),]
colnames(predict_beacon)=gsub("SEP","",colnames(predict_beacon))

comgenes=intersect(rownames(predict_cambridge),rownames(predict_beacon))
idx=match(comgenes,rownames(predict_cambridge))
predict_cambridge=predict_cambridge[idx,]
idx=match(comgenes,rownames(predict_beacon))
predict_beacon=predict_beacon[idx,]

get_output=function(dat=predict_cambridge)
{
  res=data.frame(phenoBE=rep(-9,ncol(dat)),phenoEA=-9,phenoEABE=-9)
  rownames(res)=colnames(dat)
  idx=colnames(dat) %in% sampletable$localid
  res$phenoBE[!idx]=res$phenoEA[!idx]=res$phenoEABE[!idx]=1
  avaisamples=colnames(dat)[colnames(dat) %in% sampletable$localid]
  idx1=match(avaisamples,colnames(dat))
  idx2=match(avaisamples,sampletable$localid)
  res$phenoBE[idx1]=sampletable$phenoBE_bca[idx2]
  res$phenoEA[idx1]=sampletable$phenoEA_bca[idx2]
  res$phenoEABE[idx1]=sampletable$phenoEABE_bca[idx2]
  return(res)
}

cambridge_output=get_output()
beacon_output=get_output(dat=predict_beacon)

imputedqcBEACON=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_1000g/chr1_filter_hg19tohg38_flip.fam")
imputedqcBEACONsamples=gsub("SEP","",imputedqcBEACON$V2)

imputedqcCambridge=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/cambridgewtccc_qc_1000g/chr17_filter_hg19tohg38_flip.fam")
imputedqcCambridgesamples=gsub("SEP","",imputedqcCambridge$V2)

prefix="dist500K_GTEx_mucosa_mergeqc_Jan15"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
sum(!is.na(predict_min$n_avaisnp))
predict_min=predict_min[!is.na(predict_min$n_avaisnp),]

# quantile(predict_min$n_avaisnp/predict_min$n_totalsnp,c(0,0.01,0.05,0.1,0.15,0.25,0.5,1))
# 0%        1%        5%       10%       15%       25%       50%      100% 
# 0.2000000 0.7500000 0.8823529 0.9152542 0.9310345 0.9500000 0.9750000 1.0000000 
#tmp=predict_min$n_avaisnp/predict_min$n_totalsnp
colnames(predict_min)=gsub("SEP","",colnames(predict_min))
predict_beacon=predict_min[,colnames(predict_min) %in% imputedqcBEACONsamples]
predict_cambridge=predict_min[,colnames(predict_min) %in% imputedqcCambridgesamples]
cambridge_output=get_output(dat=predict_cambridge1)
beacon_output=get_output(dat=predict_beacon1)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
sum(rownames(predict_cambridge) %in% proteingenes)
idx=match(rownames(predict_min),rownames(res_min))
info=data.frame(r2=res_min$r2[idx])
rownames(info)=rownames(res_min)[idx]
idx=match(rownames(info),rownames(predict_min))
info$totalsnp=predict_min$n_totalsnp[idx]
info$avaisnp=predict_min$n_avaisnp[idx]
quantile(predict_min$n_avaisnp[idx]/predict_min$n_totalsnp[idx],c(0,0.01,0.05,0.1,0.15,0.25,0.5,1))
# 0%        1%        5%       10%       15%       25%       50%      100% 
#   0.2000000 0.7500000 0.8823529 0.9152542 0.9310345 0.9500000 0.9750000 1.0000000 
predict_beacon_protein=predict_beacon[rownames(predict_beacon) %in% proteingenes,]
save(predict_cambridge,predict_beacon,cambridge_output,beacon_output,proteingenes,info,file="../result/Beacon_Cambridge_predictgenexp.RData")
for ( i in 1:3)
{
  idx=which(beacon_output[,i]==-9)
  beacon_output[idx,i]=NA
  idx=which(cambridge_output[,i]==-9)
  cambridge_output[idx,i]=NA
}

gene="SSBP4"
gene="KXD1"
gene="ARMC6"
idx=which(rownames(predict_beacon)==gene)
info[idx,]
fit=glm(I(beacon_output[,1]==2)~as.numeric(predict_beacon[idx,]),family = "binomial")
fit=glm(I(cambridge_output[,1]==2)~as.numeric(predict_cambridge[idx,]),family = "binomial")
summary(fit)
