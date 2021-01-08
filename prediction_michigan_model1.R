#!/usr/bin/env Rscript
# #get all the data
# # #load copynumber
# load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp.RData")
# load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGA_EAC_michigan_genotype.RData")
# #geneexp
# phenotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_norm.txt"
# phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
# covariatefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_COVA_GE_PEER_clinical_used.txt"
# snp=genotype
# snppos=genotypepos
# colnames(snppos)=c("chr","pos","mag","min")
# snppos$chr=gsub(23,"X",snppos$chr)
# idx=is.na(unlist(snp[1,]))
# phenotype=fread(phenotypefile,header=T,sep="\t")
# phenotype=as.data.frame(phenotype)
# rownames(phenotype)=phenotype[,1]
# phenotype=phenotype[,-1]
# phenotypepos=fread(phenotypeposfile,header=T,sep="\t")
# phenotypepos=as.data.frame(phenotypepos)
# rownames(phenotypepos)=phenotypepos[,1]
# #all(rownames(phenotypepos)==rownames(phenotype)) #T
# covariate=as.data.frame(fread(covariatefile))
# rownames(covariate)=covariate[,1]
# covariate=covariate[,-1]
# covariate=t(covariate)
# covariate=as.data.frame(covariate)
# covariate$stage[is.na(covariate$stage)]=3 #glmnet not allow NA
# #1 sample "TCGA-2H-A9GG" doesn't have copynumber data
# cosample=intersect(colnames(snp),colnames(copynumber_new))
# cosample=intersect(colnames(phenotype),cosample)
# idx=match(cosample,colnames(snp))
# snp=snp[,idx]
# idx=match(cosample,colnames(phenotype))
# phenotype=phenotype[,idx]
# idx=match(cosample,colnames(copynumber_new))
# copynumber=copynumber_new[,idx]
# idx=match(cosample,rownames(covariate))
# covariate=covariate[idx,]
# #sum(rownames(covariate)==colnames(phenotype))
# #sum(colnames(phenotype)==colnames(snp))
# #sum(colnames(copynumber)==colnames(snp))
# #sum(rownames(copynumber)==rownames(phenotype))
# idx=match(cosample,colnames(mutation))
# mutation=mutation[,idx]
# sum(rownames(mutation) %in% rownames(phenotype))
# rownames(mutation)[!rownames(mutation) %in% rownames(phenotype)][1:3]
# mutation1=mutation[rownames(mutation) %in% rownames(phenotype),]
# mutation_new=data.frame(matrix(0,nrow=nrow(phenotype),ncol=ncol(phenotype)))
# rownames(mutation_new)=rownames(phenotype)
# colnames(mutation_new)=colnames(phenotype)
# idx=match(rownames(mutation1),rownames(mutation_new))
# mutation_new[idx,]=mutation1
# mutation=mutation_new
# all(colnames(mutation_new)==colnames(phenotype)) #T
#save(snp,snppos,phenotype,phenotypepos,copynumber,mutation,covariate,file="../result/TCGAdatafor_prediction_michigan.RData")


#create model and compute r-squared for a gene
compute_cor_arow=function(i,opt="1se",ncv=5,distcutoff=1e6)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff) 
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  if (length(idx)>1)
  {
    X=t(snp[idx,])
    mutationV=unlist(mutation[i,])
    if (all(mutationV==0)) #no mutation
    {
      Xall=as.matrix(cbind(X,cn=unlist(copynumber[i,]),covariate))
      covariateall=as.matrix(cbind(cn=unlist(copynumber[i,]),covariate))
    }else
    {
      Xall=as.matrix(cbind(X,cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
      covariateall=as.matrix(cbind(cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
    }
    
    set.seed(i+10000)
    if (opt=="min")
    {
      penalty=rep(1,ncol(Xall))
      penalty[(length(idx)+1):length(penalty)]=0 #force the covariates to be included in the model
      #cvfit <- cv.glmnet(as.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty)
      #fit=glmnet(as.matrix(Xall),Y,nlambda = 100,penalty.factor=penalty)
      cvfit <- cv.glmnet(as.matrix(Xall),Y,nlambda=100,nfolds=10)
      fit=glmnet(as.matrix(Xall),Y,nlambda = 100)
      lamba_sel=cvfit$lambda.min
    }else
    {
      #if use 1se and include all the covariates, very few snps will be selected, so beter to not force the covariates to be included 
      cvfit <- cv.glmnet(as.matrix(Xall),Y,nlambda=100,nfolds=10)
      fit=glmnet(as.matrix(Xall),Y,nlambda = 100)
      lamba_sel=cvfit$lambda.1se
    }
    lambda_1varsel=fit$lambda[which(fit$df>0)[1]] #lambda value when 1 snp variable is selected should be greater than this
    
    # allcor=sapply(1:ncol(X),function(x){
    #   cor(Y,Xall[,x])
    # })
    if (lamba_sel<lambda_1varsel) #if glmnet selected some variables
    {
      glmcoeff=as.matrix(coef(fit,s=lamba_sel))
      glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
      idx1=match(glmleftsnp,rownames(glmcoeff))
      glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snp by effect size
      numvar=length(glmleftsnp)
      if (numvar>0)
      {
        idx1=match(glmleftsnp,colnames(X))
        Xsel=X[,idx1]
        if (numvar>1) #check if number of covariate is greater than sample size
        {
          nummaxvar=min(nrow(X)-ncol(covariateall)-1,numvar)
          numvar=nummaxvar
          Xsel=Xsel[,1:nummaxvar]
        }
        if (numvar==1) #keep Xsel as in matrix form
        {
          Xsel=matrix(Xsel,ncol=1)
          colnames(Xsel)=glmleftsnp
        }
        Xall1=as.matrix(cbind(Xsel,covariateall))
        #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
        fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
        #summary(fit1)$r.squared
        lmcoeff=summary(fit1)$coefficients
        #align up coeff with Xsel
        rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
        lmleftsnp=rownames(lmcoeff)[rownames(lmcoeff) %in% colnames(Xsel)] 
        numvar=length(lmleftsnp)
        if (numvar>0)
        {
          glmflag=1
          idx1=match(lmleftsnp,rownames(lmcoeff))
          selectedsnps=paste0(rownames(lmcoeff)[idx1],collapse = "|")
          selectedsnps_coeff=paste0(lmcoeff[idx1,1],collapse = "|")
          idx1=match(lmleftsnp,colnames(Xsel))
          Xsel=Xsel[,idx1]
          if (numvar==1)
          {
            Xsel=matrix(Xsel,ncol=1)
            colnames(Xsel)=lmleftsnp
          }
          fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
          r2=cor(fitted,Y)^2
        }
        
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff))
}

#predicted geneexp using crossvalidation
fitted_cv=function(Xsel,covariateall,Y,ncv=5)
{
  Xall=as.matrix(cbind(Xsel,covariateall))
  maxnumvar=floor(length(Y)*(1-1/ncv))
  if (ncol(Xall)>maxnumvar-1) #number of covariates is greater than sample size, select subset of covariates
  {
    lmfit1=lm(Y~Xall)
    lmcoeff1=summary(lmfit1)$coefficients
    rownames(lmcoeff1)=gsub("Xall","",rownames(lmcoeff1))
    lmleftsnp1=rownames(lmcoeff1)[rownames(lmcoeff1) %in% colnames(Xsel)] 
    idx1=match(lmleftsnp1,rownames(lmcoeff1))
    lmleftsnp1=lmleftsnp1[order(abs(lmcoeff1[idx1,1]),decreasing = T)]
    idx1=match(lmleftsnp1,colnames(Xsel))
    Xsel=Xsel[,idx1]
    Xsel=Xsel[,1:(maxnumvar-ncol(covariateall)-1)]
    
    Xall=as.matrix(cbind(Xsel,covariateall))
  }
  
  fitted1=rep(0,length(Y))
  set.seed(10000)
  permutidx=sample(1:length(Y))
  idxs=as.integer(seq(1,length(Y)+1,length.out = ncv+1)) #boundary points of cv segments
  
  for (ii in 1:ncv)
  {
    idx_predict=rep(F,length(Y))
    idx_predict[idxs[ii]:(idxs[ii+1]-1)]=T
    trainfm=lm(Y[permutidx[!idx_predict]]~Xall[permutidx[!idx_predict],])
    traincoeff=summary(trainfm)$coefficients
    rownames(traincoeff)=gsub("Xall[permutidx[!idx_predict], ]","",rownames(traincoeff),fixed = T)
    trainleftsnps=rownames(traincoeff)[rownames(traincoeff) %in% colnames(Xsel)]
    numvar=length(trainleftsnps)
    idx1=match(trainleftsnps,colnames(Xsel))
    Xsel1=Xsel[,idx1]
    if (numvar==1)
    {
      Xsel1=matrix(Xsel1,ncol=1)
      colnames(Xsel1)=trainleftsnps
    }
    fitted1[permutidx[idx_predict]]=rep(traincoeff[1,1],sum(idx_predict)) #intercept term
    idx1=match(trainleftsnps,rownames(traincoeff))
    if (numvar>0) #to add each selected snp term
    {
      for (j in 1:numvar)
      {
        fitted1[permutidx[idx_predict]]=fitted1[permutidx[idx_predict]]+Xsel1[permutidx[idx_predict],j]*traincoeff[idx1[j],1]
      }
    }
  }
  return(fitted1)
}

#another version, only pick snps when cvfit is concave curve,  i.e. applying cvfit$lambda.1se selected some snps
compute_cor_arow=function(i,opt="1se",ncv=10,distcutoff=1e6)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff) 
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  if (length(idx)>1)
  {
    X=t(snp[idx,])
    mutationV=unlist(mutation[i,])
    if (all(mutationV==0)) #no mutation
    {
      Xall=as.matrix(cbind(X,cn=unlist(copynumber[i,]),covariate))
      covariateall=as.matrix(cbind(cn=unlist(copynumber[i,]),covariate))
    }else
    {
      Xall=as.matrix(cbind(X,cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
      covariateall=as.matrix(cbind(cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
    }
    
    set.seed(i+10000)
    cvfit <- cv.glmnet(as.matrix(Xall),Y,nlambda=100,nfolds=10)
    fit=glmnet(as.matrix(Xall),Y,nlambda = 100)
    lambda_1varsel=fit$lambda[which(fit$df>0)[1]] #lambda value when 1 snp variable is selected should be greater than this
    
    if (cvfit$lambda.1se<lambda_1varsel) #if glmnet 1se selected some variables
    {
      glmcoeff=as.matrix(coef(fit,s=cvfit$lambda.1se))
      if (sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))>0) #if SNPs were selected
      {
        if (opt=="min")
        {
          lamba_sel=cvfit$lambda.min
        }else
        {
          lamba_sel=cvfit$lambda.1se
        }
        glmcoeff=as.matrix(coef(fit,s=lamba_sel))
        glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
        idx1=match(glmleftsnp,rownames(glmcoeff))
        glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snp by effect size
        numvar=length(glmleftsnp)
        if (numvar>0)
        {
          idx1=match(glmleftsnp,colnames(X))
          Xsel=X[,idx1]
          if (numvar>1) #check if number of covariate is greater than sample size
          {
            nummaxvar=min(nrow(X)-ncol(covariateall)-1,numvar)
            numvar=nummaxvar
            Xsel=Xsel[,1:nummaxvar]
          }
          if (numvar==1) #keep Xsel as in matrix form
          {
            Xsel=matrix(Xsel,ncol=1)
            colnames(Xsel)=glmleftsnp
          }
          Xall1=as.matrix(cbind(Xsel,covariateall))
          #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
          fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
          #summary(fit1)$r.squared
          lmcoeff=summary(fit1)$coefficients
          #align up coeff with Xsel
          rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
          lmleftsnp=rownames(lmcoeff)[rownames(lmcoeff) %in% colnames(Xsel)] 
          numvar=length(lmleftsnp)
          if (numvar>0)
          {
            glmflag=1
            idx1=match(lmleftsnp,rownames(lmcoeff))
            selectedsnps=paste0(rownames(lmcoeff)[idx1],collapse = "|")
            selectedsnps_coeff=paste0(lmcoeff[idx1,1],collapse = "|")
            idx1=match(lmleftsnp,colnames(Xsel))
            Xsel=Xsel[,idx1]
            if (numvar==1)
            {
              Xsel=matrix(Xsel,ncol=1)
              colnames(Xsel)=lmleftsnp
            }
            fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
            r2=cor(fitted,Y)^2
          }
          
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff))
}

#load data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp


#compute
test=compute_cor_arow(i=15859,opt="1se",ncv=10,distcutoff = 1e6)

# #parallel
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData")
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.remote.exec(load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData"))
#mpi.remote.exec(sum(is.na(snp)))
# mpi.bcast.Robj2slave(snppos)
# mpi.bcast.Robj2slave(phenotype)
# mpi.bcast.Robj2slave(phenotypepos)
# mpi.bcast.Robj2slave(copynumber)
# mpi.bcast.Robj2slave(covariate)
mpi.bcast.Robj2slave(gr_snp)
mpi.bcast.Robj2slave(gr_pos)
mpi.bcast.cmd(library(GenomicRanges))
mpi.bcast.cmd(library(glmnet))
mpi.bcast.Robj2slave(compute_cor_arow)
mpi.bcast.Robj2slave(fitted_cv)


mpi_compute_cor=function(rows=1:nrow(phenotypepos),opt="min",ncv=10,distcutoff=1e6)
{
  res=NULL
  n=1000
  nchunks=ceiling(length(rows)/n)
  print(paste0("number of total:",nchunks))
  for (i in 1:nchunks)
  {
    cat(i,"..")
    if (i<nchunks)
    {
      seq=rows[((i-1)*n+1):(i*n)]
    }else
    {
      seq=rows[((i-1)*n+1):length(rows)]
    }
    tmp=mpi.parSapply(X=seq,FUN=compute_cor_arow,opt=opt,ncv=ncv,distcutoff=distcutoff,job.num=njobs)
    res1=matrix(unlist(tmp),ncol=6,byrow = T)
    res=rbind(res,res1)
  }
  colnames(res)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff")
  res=as.data.frame(res)
  rownames(res)=rownames(phenotypepos)
  res[,1]=as.numeric(as.character(res[,1]))
  res[,2]=as.numeric(as.character(res[,2]))
  res[,3]=as.numeric(as.character(res[,3]))
  res[,4]=as.numeric(as.character(res[,4]))
  res[,5]=as.character(res[,5])
  res[,6]=as.character(res[,6])
  return(res)
}

prefix="dist1M"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
if (!dir.exists(outfolder)) dir.create(outfolder)
#rum mpi function and save the result
res_min=mpi_compute_cor()
res_1se=mpi_compute_cor(opt="1se")
save(res_min,res_1se,file=paste0(outfolder,"/preidiction_michigan_model.RData"))
#extract snps from the prediction models
extract_snp=function(dat=rbind(res_min,res_1se))
{
  dat=dat[dat$glmflag==1,]
  snps=unique(unlist(strsplit(dat$selectedsnps,"|",fixed=T)))
  res=data.frame(snp=snps,chr=NA,pos=NA,stringsAsFactors = F)
  probenames=paste0(rownames(snppos),"_",snppos$mag)
  idx=match(res$snp,probenames)
  res$chr=snppos$chr[idx]
  res$pos=snppos$pos[idx]
  res=res[order(res[,2],res[,3]),]
  tmp=paste0(res[,2],"_",res[,3])
  idx=duplicated(tmp)
  res=res[!idx,]
  #create regions file used to extact genotype data
  for (i in 1:length(unique(res$chr)))
  {
    chr=unique(res$chr)[i]
    idx=which(res$chr==chr)
    chrres=paste0(chr,":",res$pos[idx],-res$pos[idx])
    chrres=paste0(chrres,collapse = "\n")
    filename=paste0(outfolder,"/prediction_snps_tabix_chr",chr,".txt")
    #fileCon=file(filename)
    #writeLines(chrres,fileCon)
    #close(fileCon)
  }
  return(res)
}
tmp=extract_snp()

cmd=paste0("./bca_extract_genotype.sh ",prefix)
system(cmd,wait = T)
#after run extract_genotype.sh, read the extracted genotype
library(data.table)
for (i in 1:22)
{
  cat(i,'..')
  tmp=fread(paste0(outfolder,"/chr",i,"_select.raw"),header=T)
  tmp=as.data.frame(tmp)
  rownames(tmp)=tmp$IID
  tmp=tmp[,7:ncol(tmp)]
  tmp1=fread(paste0(outfolder,"/chr",i,"_select.bim"))
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
bcagenotype=as.data.frame(t(bcagenotype))
save(bcagenotype,bcabim,file=paste0(outfolder,"/bca_extractgenotype.RData"))
# read_genotype=function(genotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/imputation_michigan/chr3_filter.raw",
#                        snpfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/imputation_michigan/chr3_filter.bim",
#                        snpname="3:196519878_C")
# {
#   genotype1=fread(genotypefile)
#   genotype1=as.data.frame(genotype1)
#   idx=which(colnames(genotype1)==snpname)
#   res=genotype1[,idx]
# }
# tcga_res=read_genotype()
# table(tcga_res)
# # 0  1  2 
# # 57 75 53
# bca_res=read_genotype(genotypefile = "/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr3_select.raw",snpname ="3:196519878_T")
# table(bca_res)
# # 0    1    2 
# # 3173 4430 1600 
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
      idx1=which(rownames(snppos)==tmp[1])
      idx=which(bcabim$V2==tmp[1] & bcabim$V6==tmp[2]& bcabim$V5==snppos[idx1,4])[1]
      if (length(idx)>0)
      {
        idx1=which(selectedsnps==missingsnps[j])
        selectedsnps[idx1]=paste0(bcabim$V2[idx],"_",bcabim$V5[idx])
        correctedsnps=c(correctedsnps,paste0(bcabim$V2[idx],"_",bcabim$V5[idx]))
        #idx2=which(rownames(bcagenotype_chunk)==selectedsnps[idx1])
        #bcagenotype_chunk[idx2,]=2-bcagenotype_chunk[idx2,]
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



##salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
if (!exists("phenotypepos")) load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData") #phenotypepos
load(paste0(outfolder,"/preidiction_michigan_model.RData")) #models
load(paste0(outfolder,"/bca_extractgenotype.RData")) #extracted genotype
# library(Rmpi)
# njobs=mpi.universe.size() - 1
print(njobs)
# mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(res_min)
mpi.bcast.Robj2slave(res_1se)
mpi.bcast.Robj2slave(outfolder)
#mpi.remote.exec(load(paste0(outfolder,"/bca_extractgenotype.RData")))
mpi.bcast.Robj2slave(bcabim)
mpi.bcast.Robj2slave(predict_geneexp)
mpi.bcast.Robj2slave(snppos)
mpi_predict_geneexp=function(modeltable=res_min)
{
  modeltable=modeltable[modeltable$glmflag==1,]
  rows=1:nrow(modeltable)
  snps_bcagenotype=unlist(strsplit(rownames(bcagenotype),"_"))
  snps_bcagenotype=snps_bcagenotype[seq(1,length(snps_bcagenotype),2)]
  
  res=NULL
  n=njobs
  nchunks=ceiling(length(rows)/n)
  print(paste0("number of total:",nchunks))
  for (i in 1:nchunks)
  {
    if (i %% 5==0) cat(i,"..")
    if (i<nchunks)
    {
      seq=rows[((i-1)*n+1):(i*n)]
    }else
    {
      seq=rows[((i-1)*n+1):length(rows)]
    }
    selectedsnps=unlist(strsplit(modeltable$selectedsnps[seq],"|",fixed=T))
    snps_selectedsnps=unlist(strsplit(selectedsnps,"_"))
    snps_selectedsnps=snps_selectedsnps[seq(1,length(snps_selectedsnps),2)]
    idx=snps_bcagenotype %in% snps_selectedsnps #extract genotypes for each chunk
    bcagenotype_chunk=bcagenotype[idx,]
    mpi.bcast.Robj2slave(bcagenotype_chunk)
    tmp=mpi.parSapply(X=seq,FUN=predict_geneexp,modeltable=modeltable,job.num=njobs)
    res1=matrix(unlist(tmp),ncol=2+ncol(bcagenotype),byrow = T)
    res=rbind(res,res1)
  }
  res=as.data.frame(res)
  rownames(res)=rownames(modeltable)
  colnames(res)=c("n_totalsnp","n_avaisnp",colnames(bcagenotype))
  return(res)
}
predict_min=mpi_predict_geneexp()
predict_1se=mpi_predict_geneexp(modeltable=res_1se)
save(predict_min,predict_1se,file=paste0(outfolder,"/bca_predict_geneexp.RData"))
#after run mpi_predict_geneexp


#check predicted geneexp
tmp=rowMeans(predict_min[,3:ncol(predict_min)],na.rm=T)
sum(tmp==0,na.rm=T) #21
sum(is.na(tmp)) #6
which(tmp==0)
# LOC728788 ANKRD20A4      BCL8    BMS1P5  C1orf152     CBWD6   CXADRP2    FAM35B   FOXD4L6 
# 18       547      1156      1227      1741      2422      3516      4764      5175 
# GPRIN2 HIST1H2BL  HIST1H4E  HIST1H4K  KIAA0368 LOC645166 LOC646214 LOC648691 LOC654342 
# 5745      6084      6099      6102      6905      7736      7741      7751      7760 
# PGM5P2    SPDYE1     SYT15 
# 10062     12874     13234 

#association
if (!exists("sampletable"))
{
  #library(xlsx)
  #sampletable=read.xlsx("/fh/fast/dai_j/BEACON/BeagessCambridgeAmos/bc_AT_JD1.xlsx",sheetIndex = 1,stringsAsFactors=F)
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data//bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  sampletable$site[sampletable$site=="NA"]=NA
  geneexpsamplenames=strsplit(colnames(predict_min)[3:ncol(predict_min)],"_")
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}

for (i in unique(sampletable$site))
{
  cat(i,"..")
  idx=which(sampletable$site==i)
  print(table(sampletable$phenoBE_bc[idx],sampletable$phenoEA_bc[idx]))
}
readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  if (opt==1)
  {
    eigsamples=read.table(eigsampfile,stringsAsFactors = F)
    eigsamples=eigsamples$V2
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=thesamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}

eigenstratmatrix=readeigenstrat()

pvalue_arow=function(i,idx1,idx2)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  x=as.numeric(predict_bcageneexp[i,c(idx1,idx2)])
  BE_p=NA
  if (sum(is.na(x))<0.6*length(x))
  {
    fm=glm(y~x+age+sex+pc1+pc2+pc3,data=covariates,family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      BE_p=summary(fm)$coefficients[2,4]
    }
  }
  return(BE_p)
}


qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

load(paste0(outfolder,"/bca_predict_geneexp.RData"))
recodesite=function(covariates1=covariates,idx1,idx2)
{
  caseidx=1:length(idx1)
  coidx=(length(idx1)+1):nrow(covariates1)
  covariates1$site1=NA
  for (mysite in unique(covariates1$site))
  {
    idx=which(covariates1$site==mysite)
    if (sum(idx %in% caseidx)>0 && sum(idx %in% coidx)>0)
    {
      covariates1$site1[idx]=mysite
    }else
    {
      covariates1$site1[idx]="other"
    }
  }
  return(covariates1)
}
mpi_association=function(predict_bcageneexp=predict_min)
{
  predict_bcageneexp=predict_bcageneexp[,3:ncol(predict_bcageneexp)]
  geneexpsamplenames=colnames(predict_bcageneexp)
  geneexpsamplenames=strsplit(geneexpsamplenames,"_")
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
  colnames(predict_bcageneexp)=geneexpsamplenames
  geneexpsamplenames=intersect(sampletable$localid,geneexpsamplenames)
  idx=match(geneexpsamplenames,colnames(predict_bcageneexp))
  predict_bcageneexp=predict_bcageneexp[,idx]
  idx=match(geneexpsamplenames,sampletable$localid)
  sampletable1=sampletable[idx,]
  for (i in 1:ncol(sampletable1)) sampletable1[,i][sampletable1[,i]==-9]=NA
  idx=match(geneexpsamplenames,colnames(eigenstratmatrix))
  eigenstratmatrix=eigenstratmatrix[,idx]
  sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:3,]))
  #BE
  mpi.bcast.Robj2slave(pvalue_arow)
  mpi.bcast.Robj2slave(recodesite)
  mpi.bcast.Robj2slave(predict_bcageneexp)
  mpi.bcast.Robj2slave(sampletable1)
  idx1=which(sampletable1$phenoBE_bc==2) #case
  length(idx1) #3288
  idx2=which(sampletable1$phenoBE_bc==1)
  length(idx2) #2176
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  seq=1:nrow(predict_bcageneexp)
  BE_p=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
  #EA
  idx1=which(sampletable1$phenoEA_bc==2) #case
  length(idx1) #2514
  idx2=which(sampletable1$phenoEA_bc==1)
  length(idx2) #2179
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  EA_p=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
  #EA vs # BE
  idx1=which(sampletable1$phenoEA_bc==2) # EAcase
  length(idx1) #2514
  idx2=which(sampletable1$phenoBE_bc==2) #BEcase
  length(idx2) #3288
  length(intersect(idx1,idx2))
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  BEA_p=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
  res=data.frame(BE_p=BE_p,BE_fdr=p.adjust(BE_p),EA_p=EA_p,EA_fdr=p.adjust(EA_p),BEA_p=BEA_p,BEA_fdr=p.adjust(BEA_p))
  rownames(res)=rownames(predict_bcageneexp)
  return(res)
}
assoc_min=mpi_association()
assoc_1se=mpi_association(predict_bcageneexp = predict_1se)
save(assoc_min,assoc_1se,file=paste0(outfolder,"/bca_assoc.RData"))

#check RBBP5
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist1M/bca_predict_geneexp.RData")
predict_bcageneexp=predict_1se
predict_bcageneexp=predict_bcageneexp[,3:ncol(predict_bcageneexp)]
geneexpsamplenames=colnames(predict_bcageneexp)
geneexpsamplenames=strsplit(geneexpsamplenames,"_")
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_bcageneexp)=geneexpsamplenames
geneexpsamplenames=intersect(sampletable$localid,geneexpsamplenames)
idx=match(geneexpsamplenames,colnames(predict_bcageneexp))
predict_bcageneexp=predict_bcageneexp[,idx]
idx=match(geneexpsamplenames,sampletable$localid)
sampletable1=sampletable[idx,]
for (i in 1:ncol(sampletable1)) sampletable1[,i][sampletable1[,i]==-9]=NA
idx=match(geneexpsamplenames,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:3,]))
gene="RBBP5"
idx1=which(sampletable1$phenoEA_bc==2) # EAcase
length(idx1) #2514
idx2=which(sampletable1$phenoBE_bc==2) #BEcase
length(idx2) #3288
covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                      bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                      pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                      stringsAsFactors = F)
covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
idxgene=which(rownames(predict_bcageneexp)==gene) #4642
pvalue_arow_sites=function(i,idx1,idx2)
{
  caseidx=1:length(idx1)
  coidx=(length(idx1)+1):nrow(covariates1)
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  x=as.numeric(predict_bcageneexp[i,c(idx1,idx2)])
  n=length(unique(covariates$site1))
  res=data.frame(matrix(NA,nrow=n,ncol=4))
  colnames(res)=c("p","coeff","ncase","nco")
  rownames(res)=unique(covariates$site1)
  for (i in 1:n)
  {
    idx=which(covariates$site1==unique(covariates$site1)[i])
    res$ncase[i]=sum(idx %in% caseidx)
    res$nco[i]=sum(idx %in% coidx)
    fm=glm(y[idx]~x[idx]+age+sex+pc1+pc2+pc3,data=covariates[idx,],family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x[idx]" %in% rownames(tmp))
    {
      res$p[i]=summary(fm)$coefficients[2,4]
      res$coeff[i]=summary(fm)$coefficients[2,1]
    }
  }
  return(res)
}

#check FLJ10661
gene="FLJ10661"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K/bca_predict_geneexp.RData")
predict_bcageneexp=predict_min
predict_bcageneexp=predict_bcageneexp[,3:ncol(predict_bcageneexp)]
geneexpsamplenames=colnames(predict_bcageneexp)
geneexpsamplenames=strsplit(geneexpsamplenames,"_")
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_bcageneexp)=geneexpsamplenames
geneexpsamplenames=intersect(sampletable$localid,geneexpsamplenames)
idx=match(geneexpsamplenames,colnames(predict_bcageneexp))
predict_bcageneexp=predict_bcageneexp[,idx]
idx=match(geneexpsamplenames,sampletable$localid)
sampletable1=sampletable[idx,]
for (i in 1:ncol(sampletable1)) sampletable1[,i][sampletable1[,i]==-9]=NA
idx=match(geneexpsamplenames,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:3,]))
idx1=which(sampletable1$phenoBE_bc==2) #case
length(idx1) #3288
idx2=which(sampletable1$phenoBE_bc==1)
length(idx2) #2176
covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                      bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                      pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                      stringsAsFactors = F)
covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
idxgene=which(rownames(predict_bcageneexp)==gene) #2182
FLJ10661res=pvalue_arow_sites(i=idxgene,idx1,idx2)

idx=which(res_1se_500k$r2>0.05)
idx=which(rownames(assoc_1se_500K) %in% rownames(res_1se_500k)[idx])
qqplot(assoc_1se_500K$BE_p[idx])

idx=which(res_1se_1M$r2>0.05)
idx=which(rownames(assoc_1se_1M) %in% rownames(res_1se_1M)[idx])
qqplot(assoc_1se_1M$BE_p[idx])

