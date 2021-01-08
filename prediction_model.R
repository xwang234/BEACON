#!/usr/bin/env Rscript
# #get all the data
# # #load copynumber
# load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp.RData")
# snpfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_score9_SNP_GE.txt" #withi info>0.9
# snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_score9_SNP_POS.txt"
# #geneexp
# phenotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_norm.txt"
# phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
# covariatefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_COVA_GE_PEER_clinical_used.txt"
# snp=fread(snpfile)
# snp=as.data.frame(snp)
# idxdup=which(duplicated(snp$id)) #remove duplicated snp (a snp appear twice with the same coding or different coding)
# idxrm=snp$id %in% snp$id[idxdup]
# sum(idxrm) #1744
# #View(snp[idxrm,])
# snp=snp[!idxrm,]
# rownames(snp)=snp$id
# snp=snp[,-1]
# snppos=fread(snpposfile)
# snppos=as.data.frame(snppos)
# snppos=snppos[!idxrm,]
# rownames(snppos)=snppos$snp
# all(rownames(snp)==rownames(snppos)) #T
# snppos=snppos[,-1]
# snppos$chr=gsub(23,"X",snppos$chr)
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
# cosample=colnames(snp)[colnames(snp) %in% colnames(copynumber_new)]
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
# save(snp,snppos,phenotype,phenotypepos,copynumber,mutation,covariate,file="../result/TCGAdatafor_prediction.RData")

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

#using glmnet on snps and other covariates, compute fitted using 5fold cross-validation
compute_cor_arow=function(i,opt="1se",ncv=5)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<5e5) #select snp within 500K of geneexp
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
#load data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp

#compute
test=compute_cor_arow(i=1)
res=data.frame(matrix(NA,nrow=nrow(phenotype),ncol=6))
rownames(res)=rownames(phenotype)
colnames(res)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff")
for (i in 1:nrow(res))
{
  if (i %%100==0) cat(i,'..')
  tmp=compute_cor_arow(i=i,opt="min")
  res[i,]=unlist(tmp)
}
res_1se=data.frame(matrix(NA,nrow=nrow(phenotype),ncol=6))
rownames(res_1se)=rownames(phenotype)
colnames(res_1se)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff")
for (i in 1:nrow(res_1se))
{
  if (i %%100==0) cat(i,'..')
  tmp=compute_cor_arow(i=i,opt="1se")
  res_1se[i,]=unlist(tmp)
}

res_10cv=data.frame(matrix(NA,nrow=nrow(phenotype),ncol=6))
rownames(res_10cv)=rownames(phenotype)
colnames(res_10cv)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff")
for (i in 1:nrow(res_10cv))
{
  if (i %%100==0) cat(i,'..')
  tmp=compute_cor_arow(i=i,opt="min",ncv=10)
  res_10cv[i,]=unlist(tmp)
}

res_1se_10cv=data.frame(matrix(NA,nrow=nrow(phenotype),ncol=6))
rownames(res_1se_10cv)=rownames(phenotype)
colnames(res_1se_10cv)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff")
for (i in 1:nrow(res_1se_10cv))
{
  if (i %%100==0) cat(i,'..')
  tmp=compute_cor_arow(i=i,opt="1se",ncv=10)
  res_1se_10cv[i,]=unlist(tmp)
}

#save(res_old,res_new,file="../result/prediction_model.RData") #res_new is the one use 1se
#save(res,res_10cv,res_1se,res_1se_10cv,file="../result/prediction_model_Jan4.RData") #Jan4 force all covariates to be included
save(res,res_10cv,res_1se,res_1se_10cv,file="../result/prediction_model_Jan7.RData")
hist(as.numeric(res$r2),col="blue")
quantile(as.numeric(res$numtotalsnp),probs=c(0,0.01,0.1,0.25,0.5,0.75,1))
# 0%   25%   50%   75%  100% 
# 0  2121  2676  3223 21490 
quantile(as.numeric(res$numselectedsnp))
# 0%  25%  50%  75% 100% 
# 0    4   16   29   91
sum(is.na(res_old$r2))/nrow(res_old) # [1] 0.1808857 glmnet failed to select snps
sum(is.na(res_new$r2))/nrow(res_new) #0.6578065 glmnet failed
quantile(as.numeric(res_old$r2),na.rm=T)
# 0%          25%          50%          75%         100% 
# 8.921500e-09 3.262870e-02 7.387062e-02 1.284006e-01 9.569630e-01 
quantile(as.numeric(res_new$r2),na.rm=T)
# 0%          25%          50%          75%         100% 
# 3.668102e-10 9.560576e-03 2.743567e-02 5.653237e-02 9.501451e-01 

# #parallel
#salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 /app/easybuild/software/R/3.3.3-foss-2016b-fh1/bin/R --interactive
# library(Rmpi)
# njobs=mpi.universe.size() - 1
# print(njobs)
# mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
# rownames(snp)=1:nrow(snp)
# #mpi.bcast.Robj2slave(snp) #snp is too big to bcast
# #mpi.bcast.cmd(snp<-rbind(snp1,snp2,snp3,snp4,snp5)) #one way to get around the problem
# #mpi.bcast.cmd(load("result/tmp_snp_for_r2.RData")) #another way
# mpi.remote.exec(load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction.RData"))
# #mpi.remote.exec(sum(is.na(snp)))
# mpi.bcast.Robj2slave(snppos)
# mpi.bcast.Robj2slave(phenotype)
# mpi.bcast.Robj2slave(phenotypepos)
# mpi.bcast.Robj2slave(copynumber)
# mpi.bcast.Robj2slave(covariate)
# mpi.bcast.Robj2slave(gr_snp)
# mpi.bcast.Robj2slave(gr_pos)
# mpi.bcast.cmd(library(GenomicRanges))
# mpi.bcast.cmd(library(glmnet))
# mpi.bcast.Robj2slave(compute_cor_arow)
# mpi.bcast.Robj2slave(fitted_cv)
# 
# 
# mpi_compute_cor=function(rows,opt="min",ncv=5)
# {
#   res=NULL
#   n=1000
#   nchunks=ceiling(length(rows)/n)
#   print(paste0("number of total:",nchunks))
#   for (i in 1:nchunks)
#   {
#     cat(i,"..")
#     if (i<nchunks)
#     {
#       seq=rows[((i-1)*n+1):(i*n)]
#     }else
#     {
#       seq=rows[((i-1)*n+1):length(rows)]
#     }
#     tmp=mpi.parSapply(X=seq,FUN=compute_cor_arow,opt=opt,ncv=ncv,job.num=njobs)
#     res1=matrix(unlist(tmp),ncol=6,byrow = T)
#     res=rbind(res,res1)
#   }
#   return(res)
# }

#extract snps from the prediction models
extract_snp=function(dat=res_new)
{
  dat=dat[dat$glmflag==1,]
  snps=unique(unlist(strsplit(dat$selectedsnps,"|",fixed=T)))
  res=data.frame(snp=snps,chr=NA,pos=NA,stringsAsFactors = F)  
  idx=match(res$snp,rownames(snppos))
  res$chr=snppos$chr[idx]
  res$pos=snppos$pos[idx]
  res=res[order(res[,2],res[,3]),]
  return(res)
}

