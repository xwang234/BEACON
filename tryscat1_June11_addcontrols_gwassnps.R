#!/usr/bin/env Rscript

#salloc -t 1-1 --constraint=gizmok -n 155 mpirun -n 1 R --interactive
library(SKAT)
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
# skat_p=function(Z,Covariate,idx1,idx2)
# {
#   y=c(rep(1,length(idx1)),rep(0,length(idx2)))
#   
#   Z1=Z[c(idx1,idx2),,drop=F]
#   Covariate1=as.matrix(Covariate[c(idx1,idx2),])
#   tmp=rowSums(Z1,na.rm=T)
#   if (sum(tmp!=0) >1)
#   {
#     obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
#     out2<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #0.0004209424
#     res=out2$p.value
#     return(res)
#   }else
#   {
#     return(NA)
#   }
# }
skat_p=function(Z,Covariate,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=Covariate[c(idx1,idx2),]
  tmp=rowSums(Z1,na.rm=T)
  if (sum(tmp!=0) >1)
  {
    obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
    out2<-SKAT(Z1, obj.s,weights.beta=c(1,1),r.corr=0,is_dosage=T) #0.0004209424
    res=out2$p.value
    return(res)
  }else
  {
    return(NA)
  }
}

#genemodel,bcagenotype,Covariate were loaded,used for mpi
compute_p_arow=function(i)
{
  genes=rownames(genemodel)[genemodel$glmflag==1]
  gene=genes[i]
  res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA,stringsAsFactors = F)
  rownames(res)=gene
  selectedsnps=unlist(strsplit(genemodel$selectedsnps[which(rownames(genemodel)==gene)],"|",fixed=T))
  idx=match(selectedsnps,rownames(bcagenotype))
  idxtocorrect=which(is.na(idx))
  
  if (length(idxtocorrect)>0)
  {
    for (j in idxtocorrect)
    {
      tmp=unlist(strsplit(selectedsnps[j],"_"))
      selectedsnps[j]=paste0(tmp[c(1,3,2)],collapse = "_")
    }
    idx=match(selectedsnps,rownames(bcagenotype))
  }
  
  Z=t(bcagenotype[idx,,drop=F])
  if (length(idxtocorrect)>0)
  {
    Z[idxtocorrect,]=2-Z[idxtocorrect,]
  }
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  res$BE_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoEA_bca==1)
  res$EA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoEA_bca==2)
  idx2=which(covariatetable$phenoBE_bca==2)
  res$BEA_p=skat_p(Z,Covariate,idx1,idx2)
  idx1=which(covariatetable$phenoBE_bca==2 | covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  res$BEEA_p=skat_p(Z,Covariate,idx1,idx2)
  return(res)
}
mpi.bcast.Robj2slave(skat_p)
mpi.bcast.Robj2slave(compute_p_arow)
mpi.remote.exec(library(SKAT))

outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_Jan19")
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_Jan19"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_Jan19"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_Jan19"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_Jan19"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_Jan19"
# outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_Jan19"
outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor")

outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005_rmhighcor_allcovar",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005_rmhighcor_allcovar",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005_rmhighcor_allcovar",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005_rmhighcor_allcovar",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005_rmhighcor_allcovar",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005_rmhighcor_allcovar")

#first run PC4, then run PC6
opt="PC4"
opt="PC6"
for (ii in 1:length(outfolders))
{
  set.seed(1000)
  outfolder=outfolders[ii]
  assocfile=paste0(outfolder,"/bca_assoc.RData")
  #if (file.exists(assocfile) & !file.exists(paste0(outfolder,"/skat_gwas_res.RData")))
  {
    print(outfolder)
    load(paste0(outfolder,"/preidiction_michigan_model.RData"))
    load(paste0(outfolder,"/bca_extractgenotype.RData"))
    library(readxl)
    sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
    sampletable=as.data.frame(sampletable)
    for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
    
    #add covariate table (pc1-pc4)
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
      return(tmp)
    }
    
    covariatetable=readeigenstrat()
    rownames(covariatetable)=gsub("SEP","",rownames(covariatetable))
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
    
    allsamples=intersect(colnames(bcagenotype),rownames(covariatetable))
    idx=match(allsamples,colnames(bcagenotype))
    bcagenotype=bcagenotype[,idx]
    idx=match(allsamples,rownames(covariatetable))
    covariatetable=covariatetable[idx,]
    
    #try different number of PCs
    if (opt=="PC4")
      Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
    if (opt=="PC6")
      Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]
    allgwassnps=read.table("../result/dong26snp_addcontrols_genotypedat.txt")
    #allgwassnps=read.table("../result/dong26snp_addcontrols_genotypedat_24snp.txt") #skat_gwas_res_old.RData
    
    idx=match(rownames(Covariate),rownames(allgwassnps))
    tmp=allgwassnps[idx,]
    
    Covariate=cbind.data.frame(Covariate,tmp)
    
    mpi.bcast.Robj2slave(outfolder)
    #mpi.remote.exec(load(paste0(outfolder,"/preidiction_michigan_model.RData")))
    #mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
    mpi.bcast.Robj2slave(sampletable)
    mpi.bcast.Robj2slave(allsamples)
    mpi.bcast.Robj2slave(Covariate)
    mpi.bcast.Robj2slave(covariatetable)
    
    #add chr
    get_chr_model=function(genemodel=res_min)
    {
      genemodel=genemodel[genemodel$glmflag==1,]
      genemodel$chr=genemodel$pos=NA
      for (i in 1:nrow(genemodel))
      {
        tmp=unlist(strsplit(genemodel$selectedsnps[i],"|",fixed=T))
        tmp=unlist(strsplit(tmp,":"))
        genemodel$chr[i]=tmp[1]
        tmp=unlist(strsplit(tmp[2],"_"))
        genemodel$pos[i]=as.numeric(tmp[1])
      }
      return(genemodel)
    }
    genemodel=get_chr_model()
    
    library(GenomicRanges)
    bcabim$V1=as.character(bcabim$V1)
    gr_allbcagenotype=GRanges(seqnames = bcabim$V1,ranges = IRanges(start=bcabim$V4,width=1))
    
    mpi_compute_p=function(genemodel=genemodel)
    {
      genemodel=genemodel[genemodel$glmflag==1,]
      bcabim$V1=as.character(bcabim$V1)
      allgenemodel=genemodel
      allbcagenotype=bcagenotype
      resall=NULL
      chrs=unique(allgenemodel$chr)
      for (chr in chrs)
      {
        res=NULL
        print(paste0(chr,"---"))
        genemodel=allgenemodel[allgenemodel$chr==chr,]
        idx=order(genemodel$chr,genemodel$pos)
        genemodel=genemodel[idx,]
        gr_genemodel=GRanges(seqnames = genemodel$chr,ranges=IRanges(start=genemodel$pos,width = 1))
        mpi.bcast.Robj2slave(genemodel)
        genes=rownames(genemodel)
        rows=1:nrow(genemodel)
        #
        n=200 #work on a subset of genes each time, to save memory
        nchunks=ceiling(length(rows)/n)
        print(paste0("number of total:",nchunks))
        for (i in 1:nchunks)
        {
          #if (i %% 5==0) cat(i,"..")
          if (i<nchunks)
          {
            seq=rows[((i-1)*n+1):(i*n)]
          }else
          {
            seq=rows[((i-1)*n+1):length(rows)]
          }
          gr_chunk=GRanges(seqnames = chr,ranges=IRanges(start=max(c(genemodel$pos[seq[1]]-1e6,1)),end=genemodel$pos[seq[length(seq)]]+1e6)) #
          tmp=distance(gr_allbcagenotype,gr_chunk)
          idx=which(tmp==0)
          bcagenotype=allbcagenotype[idx,]
          mpi.bcast.Robj2slave(bcagenotype)
          tmp=mpi.parSapply(X=seq,FUN=compute_p_arow,job.num=njobs)
          if (length(tmp) %% 4 !=0) print(i)
          res1=as.data.frame(matrix(unlist(tmp),ncol=4,byrow = T))
          #rownames(res1)=rownames(genemodel)[seq]
          res=rbind(res,res1)
        }
        rownames(res)=rownames(genemodel)
        resall=rbind(resall,res)
      }
      
      return(resall)
    }
    #load(paste0(outfolder,"/skat_gwas_res.RData"))
    if (opt=="PC4")
    {
      skat_min2=mpi_compute_p(genemodel=genemodel)
      #save skat result
      save(skat_min2,file=paste0(outfolder,"/skat_gwas_res.RData"))
    }
    if (opt=="PC6")
    {
      load(paste0(outfolder,"/skat_gwas_res.RData"))
      skat_min2_pc6=mpi_compute_p(genemodel=genemodel)
      #save skat result
      save(skat_min2,skat_min2_pc6,file=paste0(outfolder,"/skat_gwas_res.RData"))
    }
  }
}


mpi.close.Rslaves()
mpi.quit()
quit()