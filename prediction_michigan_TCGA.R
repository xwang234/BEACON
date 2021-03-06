#!/usr/bin/env Rscript
#modified from model7.R
#USE 87 TCGA EAC samples and imputed genotypes
# #get all the data
# # #load copynumber
load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp.RData")
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGA_EAC_imp_michigan_genotype.RData")
#geneexp
phenotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_norm.txt"
phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
covariatefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_COVA_GE_PEER_clinical_pc4_used.txt"
phenotype=fread(phenotypefile,header=T,sep="\t")
phenotype=as.data.frame(phenotype)
rownames(phenotype)=phenotype[,1]
phenotype=phenotype[,-1]
phenotypepos=fread(phenotypeposfile,header=T,sep="\t")
phenotypepos=as.data.frame(phenotypepos)
rownames(phenotypepos)=phenotypepos[,1]
#all(rownames(phenotypepos)==rownames(phenotype)) #T
covariate=as.data.frame(fread(covariatefile))
rownames(covariate)=covariate[,1]
covariate=covariate[,-1]
covariate=t(covariate)
covariate=as.data.frame(covariate)
covariate$stage[is.na(covariate$stage)]=3 #glmnet not allow NA
idx=match(rownames(covariate),tcgaclinical$bcr_patient_barcode)
covariate$gender=tcgaclinical$gender[idx]
covariate$gender[covariate$gender=="MALE"]=1
covariate$gender[covariate$gender=="FEMALE"]=0
covariate$disease=NA #"Adenomas"
tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Adenomas_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
idx=match(tmp$submitter_id,rownames(covariate))
covariate$disease[idx]=0 #"Adenomas"
tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Squmous_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
idx=match(tmp$submitter_id,rownames(covariate))
covariate$disease[idx]=1 #"Squamous"
#1 sample "TCGA-2H-A9GG" doesn't have copynumber data
cosample=intersect(colnames(snp),colnames(copynumber_new))
cosample=intersect(colnames(phenotype),cosample)
cosample=intersect(rownames(covariate)[covariate$disease==0],cosample)
idx=match(cosample,colnames(snp))
snp=snp[,idx]
idx=match(cosample,colnames(phenotype))
phenotype=phenotype[,idx]
idx=match(cosample,colnames(copynumber_new))
copynumber=copynumber_new[,idx]
idx=match(cosample,rownames(covariate))
covariate=covariate[idx,]
covariate=covariate[,-ncol(covariate)]
#sum(rownames(covariate)==colnames(phenotype))
#sum(colnames(phenotype)==colnames(snp))
#sum(colnames(copynumber)==colnames(snp))
#sum(rownames(copynumber)==rownames(phenotype))
idx=match(cosample,colnames(mutation))
mutation=mutation[,idx]
sum(rownames(mutation) %in% rownames(phenotype))
rownames(mutation)[!rownames(mutation) %in% rownames(phenotype)][1:3]
mutation1=mutation[rownames(mutation) %in% rownames(phenotype),]
mutation_new=data.frame(matrix(0,nrow=nrow(phenotype),ncol=ncol(phenotype)))
rownames(mutation_new)=rownames(phenotype)
colnames(mutation_new)=colnames(phenotype)
idx=match(rownames(mutation1),rownames(mutation_new))
mutation_new[idx,]=mutation1
mutation=mutation_new
all(colnames(mutation_new)==colnames(phenotype)) #T
#save(snp,snppos,phenotype,phenotypepos,copynumber,mutation,covariate,file="../result/TCGAdatafor_prediction_michigan.RData")



#predicted geneexp using crossvalidation
#predicted geneexp using crossvalidation
fitted_cv=function(Xsel,covariateall,Y,ncv=5)
{
  Xall=data.matrix(cbind(Xsel,covariateall))
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
    Xsel=Xsel[,idx1,drop=F]
    Xsel=Xsel[,1:min(maxnumvar-ncol(covariateall)-1,ncol(Xsel)),drop=F]
    
    Xall=data.matrix(cbind(Xsel,covariateall))
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
#include p-values of disease types
compute_cor_arow=function(i,opt="1se",ncv=10,distcutoff=1e6,alpha=1)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  p_cn=NA
  p_mutation=NA
  #p_disease=NA
  p_gender=NA
  p_stage=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  if (length(idx)>1)
  {
    X=t(snp[idx,])
    
    mutationV=unlist(mutation[i,])
    if (all(mutationV==0)) #no mutation
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),covariate))
      covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),covariate))
    }else
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
      covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
    }
    
    penalty=rep(1,ncol(Xall))
    penalty[(length(idx)+1):length(penalty)]=0 #force the covariates to be included in the model
    set.seed(i+10000)
    cvfit=tryCatch(
      {
        cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=alpha)
      },
      error=function(e)
      {
        return(F)
      }
    )
    if (is.list(cvfit))
    {
      fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=alpha)
      lambda_1varsel=fit$lambda[which(fit$df>ncol(covariateall))[1]] #lambda value when 1 snp variable is selected should be greater than this
      if (cvfit$lambda.min<=lambda_1varsel) #if glmnet min selected some variables
      {
        glmcoeff=as.matrix(coef(fit,s=cvfit$lambda.min))
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
          idx2=match(glmleftsnp,rownames(glmcoeff))
          glmleftsnpcoeff=glmcoeff[idx2]
          numvar=length(glmleftsnp)
          if (numvar>0)
          {
            idx1=match(glmleftsnp,colnames(X))
            Xsel=X[,idx1,drop=F]
            Xall1=data.matrix(cbind(Xsel,covariateall))
            #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
            fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
            #summary(fit1)$r.squared
            lmcoeff=summary(fit1)$coefficients
            rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
            if (sum(rownames(lmcoeff)=="cn")>0) p_cn=lmcoeff[which(rownames(lmcoeff)=="cn"),4]
            if (sum(rownames(lmcoeff)=="mutation")>0) p_mutation=lmcoeff[which(rownames(lmcoeff)=="mutation"),4]
            if (sum(rownames(lmcoeff)=="stage")>0) p_stage=lmcoeff[which(rownames(lmcoeff)=="stage"),4]
            if (sum(rownames(lmcoeff)=="age")>0) p_age=lmcoeff[which(rownames(lmcoeff)=="age"),4]
            #p_disease=lmcoeff[which(rownames(lmcoeff)=="disease"),4]
            if (sum(rownames(lmcoeff)=="gender")>0) p_gender=lmcoeff[which(rownames(lmcoeff)=="gender"),4]
            glmflag=1
            selectedsnps=paste0(glmleftsnp,collapse = "|")
            selectedsnps_coeff=paste0(glmleftsnpcoeff,collapse = "|")
            # if (numvar==1)
            # {
            #   Xsel=matrix(Xsel,ncol=1)
            #   colnames(Xsel)=glmleftsnp
            # }
            fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
            r2=cor(fitted,Y)^2
          }
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              p_cn=p_cn,p_mutation=p_mutation,p_stage=p_stage,p_age=p_age,p_gender=p_gender))
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
test=compute_cor_arow(i=4777,opt="min",ncv=10,distcutoff = 5e5)

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


mpi_compute_cor=function(rows=1:nrow(phenotypepos),opt="min",ncv=10,distcutoff=1e6,alpha=1)
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
    tmp=mpi.parSapply(X=seq,FUN=compute_cor_arow,opt=opt,ncv=ncv,distcutoff=distcutoff,alpha=alpha,job.num=njobs)
    if (length(unlist(tmp)) %% 11 !=0) print (i)
    res1=matrix(unlist(tmp),ncol=11,byrow = T)
    res=rbind(res,res1)
  }
  colnames(res)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff","p_cn","p_mutation","p_stage","p_age","p_gender")
  res=as.data.frame(res)
  rownames(res)=rownames(phenotypepos)
  res[,5]=as.character(res[,5])
  res[,6]=as.character(res[,6])
  for (i in c(1:4,7:ncol(res))) res[,i]=as.numeric(as.character(res[,i]))
  return(res)
}

#optionA
prefix="glmnet_dist500K_EAC2"
#optionB
prefix="elnet_dist500K_EAC2"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
if (!dir.exists(outfolder)) dir.create(outfolder)
#rum mpi function and save the result
#optionA
res_min=mpi_compute_cor(distcutoff = 5e5)
res_1se=mpi_compute_cor(distcutoff = 5e5,opt="1se")
#optionB
res_min=mpi_compute_cor(distcutoff = 5e5,alpha = 0.5)
res_1se=mpi_compute_cor(distcutoff = 5e5,opt="1se",alpha = 0.5)
# res_min=mpi_compute_cor()
# res_1se=mpi_compute_cor(opt="1se")
save(res_min,res_1se,file=paste0(outfolder,"/preidiction_michigan_model.RData"))
#extract snps from the prediction models
extract_snp=function(dat=rbind(res_min,res_1se))
{
  dat=dat[dat$glmflag==1,]
  snps=unique(unlist(strsplit(dat$selectedsnps,"|",fixed=T)))
  res=data.frame(snp=snps,chr=NA,pos=NA,stringsAsFactors = F)
  
  probenames=rownames(snppos)
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
    fileCon=file(filename)
    writeLines(chrres,fileCon)
    close(fileCon)
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
  tmp1=fread(paste0(outfolder,"/chr",i,"_select.bim"))
  tmp1=as.data.frame(tmp1)
  tmp=fread(paste0(outfolder,"/chr",i,"_select.raw"),header=T)
  tmp=as.data.frame(tmp)
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
  snps_bcagenotype=snps_bcagenotype[seq(1,length(snps_bcagenotype),3)]
  
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
    snps_selectedsnps=snps_selectedsnps[seq(1,length(snps_selectedsnps),3)]
    idx=snps_bcagenotype %in% snps_selectedsnps #extract genotypes for each chunk (use chr:pos)
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
  geneexpsamplenames=strsplit(colnames(predict_1se)[3:ncol(predict_1se)],"_")
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
    fm=glm(y~x+age+sex+pc1+pc2+pc3+pc4,data=covariates,family=binomial)
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
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log10)",
       ylab="Observed p-value (log10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
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
  sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:4,]))
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
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],pc4=sampletable1$pc4[c(idx1,idx2)],
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
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],pc4=sampletable1$pc4[c(idx1,idx2)],
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
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],pc4=sampletable1$pc4[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  BEA_p=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
  
  #EA and BE vs control
  idx1=which(sampletable1$phenoBE_bc==2 | sampletable1$phenoEA_bc==2) #BE or EAcase
  length(idx1) #5802
  idx2=which(sampletable1$phenoBE_bc==1 |sampletable1$phenoEA_bc==1) #Control
  length(idx2) #2176
  length(intersect(idx1,idx2))
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],pc4=sampletable1$pc4[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  BEEA_p=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
  res=data.frame(BE_p=BE_p,BE_fdr=p.adjust(BE_p),EA_p=EA_p,EA_fdr=p.adjust(EA_p),BEA_p=BEA_p,BEA_fdr=p.adjust(BEA_p),BEEA_p=BEEA_p,BEEA_fdr=p.adjust(BEEA_p))
  rownames(res)=rownames(predict_bcageneexp)
  return(res)
}
assoc_min=mpi_association()
assoc_1se=mpi_association(predict_bcageneexp = predict_1se)
save(assoc_min,assoc_1se,file=paste0(outfolder,"/bca_assoc.RData"))

#check RBBP5
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC/bca_predict_geneexp.RData")
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
#gene="RBBP5"
idx1=which(sampletable1$phenoEA_bc==2) # EAcase
length(idx1) #2514
idx2=which(sampletable1$phenoBE_bc==2) #BEcase
length(idx2) #3288knowndat2=knowndat2[idx,
covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                      bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                      pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                      stringsAsFactors = F)
covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
idxgene=which(rownames(predict_bcageneexp)==gene) #4642
pvalue_arow_sites=function(i,idx1,idx2)
{
  caseidx=1:length(idx1)
  coidx=(length(idx1)+1):nrow(covariates)
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

#check DDX49
gene="DDX49"
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC/bca_predict_geneexp.RData")
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
idx1=which(sampletable1$phenoBE_bc==2) #case
length(idx1) #3288
idx2=which(sampletable1$phenoBE_bc==1)
length(idx2) #2176
idx1=which(sampletable1$phenoBE_bc==2 | sampletable1$phenoEA_bc==2) #BE or EAcase
length(idx1) #5802
covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                      bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                      pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                      stringsAsFactors = F)
covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
idxgene=which(rownames(predict_bcageneexp)==gene) #2182
DDX49res=pvalue_arow_sites(i=idxgene,idx1,idx2)

idx=which(res_1se_500k$r2>0.05)
idx=which(rownames(assoc_1se_500K) %in% rownames(res_1se_500k)[idx])
qqplot(assoc_1se_500K$BE_p[idx])
min(p.adjust(assoc_1se_500K$BE_p[idx]),na.rm=T)
rownames(assoc_1se_500K)[idx][which.min(p.adjust(assoc_1se_500K$BE_p[idx]))]
qqplot(assoc_1se_500K$EA_p[idx])
min(p.adjust(assoc_1se_500K$EA_p[idx]),na.rm=T)
rownames(assoc_1se_500K)[idx][which.min(p.adjust(assoc_1se_500K$EA_p[idx]))]
qqplot(assoc_1se_500K$BEA_p[idx])
min(p.adjust(assoc_1se_500K$BEA_p[idx]),na.rm=T)
rownames(assoc_1se_500K)[idx][which.min(p.adjust(assoc_1se_500K$BEA_p[idx]))]

idx=which(res_1se_1M$r2>0.05)
idx=which(rownames(assoc_1se_1M) %in% rownames(res_1se_1M)[idx])
qqplot(assoc_1se_1M$BE_p[idx])
min(p.adjust(assoc_1se_1M$BE_p[idx]),na.rm=T)
rownames(assoc_1se_1M)[idx][which.min(p.adjust(assoc_1se_1M$BE_p[idx]))]
qqplot(assoc_1se_1M$EA_p[idx])
min(p.adjust(assoc_1se_1M$EA_p[idx]),na.rm=T)
rownames(assoc_1se_1M)[idx][which.min(p.adjust(assoc_1se_1M$EA_p[idx]))]
qqplot(assoc_1se_1M$BEA_p[idx])
min(p.adjust(assoc_1se_1M$BEA_p[idx]),na.rm=T)
rownames(assoc_1se_1M)[idx][which.min(p.adjust(assoc_1se_1M$BEA_p[idx]))]




idx1=which(res_1se_500k$glmflag==1)
idx2=which(res_1se_1M$glmflag==1)
length(intersect(idx1,idx2))
idx=intersect(idx1,idx2)
sum(res_1se_1M$r2[idx]>res_1se_500k$r2[idx])
sum(res_1se_1M$selectedsnps[idx]==res_1se_500k$selectedsnps[idx])

load(paste0(outfolder,"/bca_extractgenotype.RData"))
pvalue_assoc_gene=function(gene="DDX49",idx1,idx2)
{
  bcagenotype1=bcagenotype
  # idx=grepl("\\w*_\\w*_\\w*",colnames(bcagenotype1)) #remove those NA samples with more than 1 _
  # bcagenotype1=bcagenotype1[,!idx]
  # tmp=unlist(strsplit(colnames(bcagenotype1),"_"))
  tmp=rep(NA,ncol(bcagenotype1))
  for (i in 1:length(tmp))
  {
    tmp1=unlist(strsplit(colnames(bcagenotype1)[i],"_"))
    tmp[i]=paste0(tmp1[2:length(tmp1)],collapse = "_")
  }
  colnames(bcagenotype1)=tmp
  idx=match(colnames(predict_bcageneexp),colnames(bcagenotype1))
  bcagenotype1=bcagenotype1[,idx]
  tmp=compute_cor_arow(i=which(rownames(res_1se_500k)==gene),opt="1se",ncv=10,distcutoff = 5e5)
  selectedsnps=unlist(strsplit(tmp$selectedsnps,"|",fixed=T))
  pvalues=rep(NA,length(selectedsnps))
  names(pvalues)=selectedsnps
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  for (i in 1:length(selectedsnps))
  {
    x=as.numeric(bcagenotype1[which(rownames(bcagenotype1)==selectedsnps[i]),c(idx1,idx2)])
    fm=glm(y~x+age+sex+pc1+pc2+pc3,data=covariates,family=binomial)
    #fm=glm(y~x+age+sex+site1+pc1+pc2+pc3,data=covariates,family=binomial) #including site
    tmp=summary(fm)$coefficients
    if ("x" %in% rownames(tmp))
    {
      pvalues[i]=summary(fm)$coefficients[2,4]
    }
  }
  return(pvalues)
}

computeeQTL=function(selectedsnps,gene="DDX49")
{
  res=rep(NA,length(selectedsnps))
  j=which(rownames(phenotypepos)==gene)
  Y=unlist(phenotype[j,])
  for (i in 1:length(selectedsnps))
  {
    X=t(snp[which(rownames(snp)==selectedsnps[i]),])
    mutationV=unlist(mutation[j,])
    if (all(mutationV==0)) #no mutation
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[j,]),covariate))
      covariateall=data.matrix(cbind(cn=unlist(copynumber[j,]),covariate))
    }else
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[j,]),mutation=mutationV,covariate))
      covariateall=data.matrix(cbind(cn=unlist(copynumber[j,]),mutation=mutationV,covariate))
    }
    fit1=lm(Y~Xall)
    res[i]=summary(fit1)$coefficients[2,4]
  }
}
idx=res_1se_500k$glmflag==1
quantile(res_1se_500k$r2[idx])
quantile(res_1se_500k$numselectedsnp[idx])
idx=res_1se_1M$glmflag==1
quantile(res_1se_1M$r2[idx])
quantile(res_1se_1M$numselectedsnp[idx])

load("../result/Dong23SNPs.RData")
checksnps=c("rs199620551","rs10419226","rs10423674")
knowndat=allgenotypedat[5:nrow(allgenotypedat),which(colnames(allgenotypedat) %in% checksnps)]
# tmp=rep(NA,nrow(knowndat))
# for (i in 1:length(tmp))
# {
#   tmp1=unlist(strsplit(rownames(knowndat)[i],"_"))
#   tmp[i]=paste0(tmp1[2:length(tmp1)],collapse = "_")
# }
# rownames(knowndat)=tmp
seldat=bcagenotype[which(rownames(bcagenotype) %in% selectedsnps),]
seldat=as.data.frame(t(seldat))
idx=match(rownames(seldat),rownames(knowndat))
knowndat=knowndat[idx,]
test=cor(data.matrix(cbind(seldat,knowndat)))
knowndat1=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/bca_1000g/chr19.3knownsnps.raw",header=T,stringsAsFactors = F)
rownames(knowndat1)=knowndat1$IID
knowndat1=knowndat1[,7:9]
idx=match(rownames(seldat),rownames(knowndat1))
test1=cor(data.matrix(cbind(seldat,knowndat1)))
write.csv(test1,file="../result/Cormatrix_4selectedSNPS_3knownSNPS_DDX49.csv")
alldat=data.matrix(cbind(seldat,knowndat1))
for (i in 1:ncol(alldat))
{
  print((sum(alldat[,i]==1)+sum(alldat[,i]==2)*2)/2/nrow(alldat))
}
knowndat2=knowndat1
tmp=rep(NA,nrow(knowndat2))
for (i in 1:length(tmp))
{
  tmp1=unlist(strsplit(rownames(knowndat2)[i],"_"))
  tmp[i]=paste0(tmp1[2:length(tmp1)],collapse = "_")
}
rownames(knowndat2)=tmp
idx=match(geneexpsamplenames,rownames(knowndat2))
knowndat2=knowndat2[idx,]
y=c(rep(1,length(idx1)),rep(0,length(idx2)))
x=as.numeric(predict_bcageneexp[idxgene,c(idx1,idx2)])
colnames(knowndat2)=c("snp1","snp2","snp3")
fm=glm(y~x+age+sex+pc1+pc2+pc3+snp1+snp2+snp3,data=cbind.data.frame(covariates,knowndat2[c(idx1,idx2),]),family=binomial)
fm=glm(y~x+age+sex+pc1+pc2+pc3,data=cbind.data.frame(covariates,knowndat2[c(idx1,idx2),]),family=binomial)
fm=glm(y~snp1+age+sex+pc1+pc2+pc3,data=cbind.data.frame(covariates,knowndat2[c(idx1,idx2),]),family=binomial)
fm=glm(y~snp2+age+sex+pc1+pc2+pc3,data=cbind.data.frame(covariates,knowndat2[c(idx1,idx2),]),family=binomial)
fm=glm(y~snp3+age+sex+pc1+pc2+pc3,data=cbind.data.frame(covariates,knowndat2[c(idx1,idx2),]),family=binomial)

gene="CRTC1"
knownsnps=colnames(knowndat1)
knownsnps=gsub("X","",knownsnps)
knownsnps=gsub(".",":",knownsnps,fixed = T)

TCGAknownsnps=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/chr19_3knownsnps.txt",stringsAsFactors = F)
TCGAknownsnps1=data.frame(matrix(NA,nrow=185,ncol=3))
colnames(TCGAknownsnps1)=c("rs199620551","rs10419226","rs10423674")
tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/imputationJune2014/plink/TCGAtumors_chr19_flip.fam",stringsAsFactors = F)
rownames(TCGAknownsnps1)=tmp$V1
for (i in 1:185)
{
  tmp=TCGAknownsnps[,5+(i-1)*3+2]+2*TCGAknownsnps[,5+(i-1)*3+3]
  TCGAknownsnps1[i,]=tmp
}
idx=match(colnames(phenotype),rownames(TCGAknownsnps1))
TCGAknownsnps1=TCGAknownsnps1[idx,]
Y=unlist(phenotype[which(rownames(phenotype)==gene),]) 
mutationV=unlist(mutation[which(rownames(phenotype)==gene),])
covariateall=data.matrix(cbind(cn=unlist(copynumber[which(rownames(phenotype)==gene),]),mutation=mutationV,covariate))
for (i in 1:3)
{
  X1=cbind.data.frame(x=TCGAknownsnps1[,i],covariateall)
  fm=glm(as.formula("Y~."),data=X1)
  print(summary(fm)$coefficients[2,4])
}
which(rownames(snp) %in% c("19:18803172_T","19:18804294_TG","19:18817903_A"))
rownames(snp)[which(rownames(snp) %in% c("19:18803172_T","19:18804294_TG","19:18817903_A"))] #"19:18803172_T"
X1=cbind.data.frame(x=unlist(snp[6035350,]),covariateall)
which(rownames(res_1se_500k)==gene)
res_1se_500k[4266,]
res_1se_1M[4266,]
