#!/usr/bin/env Rscript
# #USE 87 TCGA EAC samples and imputed genotypes, work on one gene
# # #get all the data
# # # #load copynumber
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
# idx=match(rownames(covariate),tcgaclinical$bcr_patient_barcode)
# covariate$gender=tcgaclinical$gender[idx]
# covariate$gender[covariate$gender=="MALE"]=1
# covariate$gender[covariate$gender=="FEMALE"]=0
# covariate$disease=NA #"Adenomas"
# tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Adenomas_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
# idx=match(tmp$submitter_id,rownames(covariate))
# covariate$disease[idx]=0 #"Adenomas"
# tmp=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_ESCA_Squmous_clinical.tsv",header=T,stringsAsFactors = F,sep="\t")
# idx=match(tmp$submitter_id,rownames(covariate))
# covariate$disease[idx]=1 #"Squamous"
# #1 sample "TCGA-2H-A9GG" doesn't have copynumber data
# cosample=intersect(colnames(snp),colnames(copynumber_new))
# cosample=intersect(colnames(phenotype),cosample)
# cosample=intersect(rownames(covariate)[covariate$disease==0],cosample)
# idx=match(cosample,colnames(snp))
# snp=snp[,idx]
# idx=match(cosample,colnames(phenotype))
# phenotype=phenotype[,idx]
# idx=match(cosample,colnames(copynumber_new))
# copynumber=copynumber_new[,idx]
# idx=match(cosample,rownames(covariate))
# covariate=covariate[idx,]
# covariate=covariate[,-ncol(covariate)]
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
# #save(snp,snppos,phenotype,phenotypepos,copynumber,mutation,covariate,file="../result/TCGAdatafor_prediction_michigan.RData")


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

##build model and compute r-squared for a gene, opt:1se or min rule
compute_cor_arow=function(gene="DDX49",opt="min",ncv=10,distcutoff=5e5)
{
  i=which(rownames(phenotypepos)==gene)
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected SNPs for the gene
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA 
  selectedsnps_coeff=NA #Beta
  p_cn=NA
  p_mutation=NA
  #p_disease=NA
  p_gender=NA
  p_stage=NA
  p_age=NA
  if (length(idx)>1)
  {
    X=t(snp[idx,]) #cis-SNPs
    mutationV=unlist(mutation[i,]) #mutation
    if (all(mutationV==0)) #no mutation
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),covariate)) #genotype+covariates
      covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),covariate)) #covariates
    }else
    {
      Xall=data.matrix(cbind(X,cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
      covariateall=data.matrix(cbind(cn=unlist(copynumber[i,]),mutation=mutationV,covariate))
    }
    
    penalty=rep(1,ncol(Xall))
    penalty[(length(idx)+1):length(penalty)]=0 #force the covariates to be included in the model
    set.seed(i+10000)
    cvfit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=1)
    fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=1)
    lambda_1varsel=fit$lambda[which(fit$df>ncol(covariateall))[1]] #lambda value when 1 snp is selected
    
    if (cvfit$lambda.min<lambda_1varsel) #if glmnet min selected some snps
    {
      glmcoeff=as.matrix(coef(fit,s=cvfit$lambda.min)) #glmnet coefficients
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
        # idx1=which(glmcoeff[,1]!=0)
        # View(glmcoeff[idx1,])
        glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
        idx1=match(glmleftsnp,rownames(glmcoeff))
        glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snps by their effect size
        numvar=length(glmleftsnp)
        if (numvar>0)
        {
          idx1=match(glmleftsnp,colnames(X))
          Xsel=X[,idx1]
          if (numvar>1) #number of selected SNPs is bounded by sample size, can't be greater than nummaxvar (64) 
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
          fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to multyple selected SNPS having identical genotype, those SNPs are usually close to each other.
          #summary(fit1)$r.squared
          lmcoeff=summary(fit1)$coefficients
          #align up coeff with Xsel
          rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
          if (sum(rownames(lmcoeff)=="cn")>0) p_cn=lmcoeff[which(rownames(lmcoeff)=="cn"),4]
          if (sum(rownames(lmcoeff)=="mutation")>0) p_mutation=lmcoeff[which(rownames(lmcoeff)=="mutation"),4]
          if (sum(rownames(lmcoeff)=="stage")>0) p_stage=lmcoeff[which(rownames(lmcoeff)=="stage"),4]
          if (sum(rownames(lmcoeff)=="age")>0) p_age=lmcoeff[which(rownames(lmcoeff)=="age"),4]
          #p_disease=lmcoeff[which(rownames(lmcoeff)=="disease"),4]
          if (sum(rownames(lmcoeff)=="gender")>0) p_gender=lmcoeff[which(rownames(lmcoeff)=="gender"),4]
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
            fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv) #get the predicted geneexp using CV
            r2=cor(fitted,Y)^2
          }
          
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              p_cn=p_cn,p_mutation=p_mutation,p_stage=p_stage,p_age=p_age,p_gender=p_gender))
}

#load TCGA data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/TCGAdatafor_prediction_michigan.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp

#get the prediction model for a gene
gene="DDX49" #all gene names are stored in rownames(phenotypepos)
genemodel=compute_cor_arow(gene=gene,opt="min",ncv=10,distcutoff = 5e5)
genemodel
# $r2
# [1] 0.004223446
# $glmflag
# [1] 1
# $numvar
# [1] 3
# $numsnp
# [1] 1259
# $selectedsnps
# [1] "19:18782197_C|19:18799388_G|19:18774810_G"
# $selectedsnps_coeff
# [1] "0.123771669097826|0.103060701491264|0.00812414396175863"
# $p_cn
# [1] 1.536283e-11
# $p_mutation
# [1] NA
# $p_stage
# [1] 0.6186303
# $p_age
# [1] 0.5355178
# $p_gender
# [1] 0.4494572
#extract SNPs within a gene in BCA data
prefix=paste0("BCA_",gene)
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
if (!dir.exists(outfolder)) dir.create(outfolder)
#extract BCA SNPs close to the gene (1MB) 
if (! file.exists(paste0(outfolder,"/dose.raw")))
{
  idx=which(rownames(phenotypepos)==gene)
  chr=phenotypepos[idx,2]
  start_loc=max(1,phenotypepos[idx,3]-1e6)
  end_loc=phenotypepos[idx,4]+1e6
  cmd=paste("./bca_extract_genotype_foragene.sh",prefix,chr,start_loc,end_loc,collapse = " ")
  system(cmd,wait = T)
}
#after run extract_genotyp_foragenee.sh, read the extracted genotype
library(data.table)
bcagenotype=fread(paste0(outfolder,"/dose.raw"),header=T)
bcagenotype=as.data.frame(bcagenotype)
rownames(bcagenotype)=bcagenotype$IID
bcagenotype=bcagenotype[,7:ncol(bcagenotype)]
bcabim=fread(paste0(outfolder,"/dose.bim"))
bcabim=as.data.frame(bcabim)
save(bcagenotype,bcabim,file=paste0(outfolder,"/bca_extractgenotype.RData"))

predict_geneexp=function()
{
  predicted_geneexp=data.frame(matrix(NA,nrow=1,ncol=2+nrow(bcagenotype)))
  rownames(predicted_geneexp)=gene
  colnames(predicted_geneexp)=c("n_totalsnp","n_avaisnp",rownames(bcagenotype))
  predicted_geneexp[,1]=genemodel$numvar
  selectedsnps=unlist(strsplit(genemodel$selectedsnps,"|",fixed=T))
  selectedcoeff=as.numeric(unlist(strsplit(genemodel$selectedsnps_coeff,"|",fixed=T)))
  
  #if some imputed snps need to to be flipped
  correctedsnps=NULL
  if (length(intersect(selectedsnps,colnames(bcagenotype)))<length(selectedsnps))
  {
    missingsnps=selectedsnps[!selectedsnps %in% colnames(bcagenotype)]
    for (j in 1:length(missingsnps))
    {
      tmp=unlist(strsplit(missingsnps[j],"_"))
      idx1=which(rownames(snppos)==tmp[1])
      if (length(idx1)>0)
      {
        idx=which(bcabim$V2==tmp[1] & bcabim$V6==tmp[2]& bcabim$V5==snppos[idx1,4])
        if (length(idx)>0)
        {
          idx=idx[1]
          idx1=which(selectedsnps==missingsnps[j])
          selectedsnps[idx1]=paste0(bcabim$V2[idx],"_",bcabim$V5[idx])
          correctedsnps=c(correctedsnps,paste0(bcabim$V2[idx],"_",bcabim$V5[idx]))
        }
      }
    }
  }
  idx=match(selectedsnps,colnames(bcagenotype))
  navaisnp=sum(!is.na(idx))
  if (navaisnp>0)
  {
    predicted_geneexp[,2]=navaisnp
    idx=selectedsnps %in% colnames(bcagenotype)
    selectedsnps=selectedsnps[idx]
    selectedcoeff=selectedcoeff[idx]
    idx1=match(selectedsnps,colnames(bcagenotype))
    availmat=bcagenotype[,idx1] #create a small matrix to avoid too large memory usage
    if (length(correctedsnps)>0)
    {
      idx2=match(correctedsnps,colnames(availmat))
      availmat[,idx2]=2-availmat[,idx2]
    }
    geneexp=as.matrix(availmat) %*% selectedcoeff #apply the linear model
    predicted_geneexp[,3:ncol(predicted_geneexp)]=geneexp
  }
  return(predicted_geneexp)
}

predicted_geneexp=predict_geneexp()

#association analysis
if (!exists("sampletable"))
{
  sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data//bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
  sampletable$site[sampletable$site=="NA"]=NA
  for (i in 1:ncol(sampletable)) sampletable[,i][sampletable[,i]==-9]=NA
  geneexpsamplenames=strsplit(colnames(predicted_geneexp)[3:ncol(predicted_geneexp)],"_")
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}

#PCs
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

pvalue_arow=function(idx1,idx2,covariates,predict_bcageneexp)
{
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  x=as.numeric(predict_bcageneexp[1,c(idx1,idx2)])
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

#load(paste0(outfolder,"/bca_predict_geneexp.RData"))
#combine sites
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
association_gene=function(predict_bcageneexp=predicted_geneexp)
{
  predict_bcageneexp=predict_bcageneexp[,3:ncol(predict_bcageneexp),drop=F]
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
  
  idx=match(geneexpsamplenames,colnames(eigenstratmatrix))
  eigenstratmatrix=eigenstratmatrix[,idx]
  sampletable1=cbind.data.frame(sampletable1,t(eigenstratmatrix[1:3,]))
  #Compare BE with CO
  idx1=which(sampletable1$phenoBE_bc==2) #case
  length(idx1) #3288
  idx2=which(sampletable1$phenoBE_bc==1)
  length(idx2) #2176
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  BE_p=pvalue_arow(idx1=idx1,idx2=idx2,covariates = covariates,predict_bcageneexp=predict_bcageneexp)
  #EA with CO
  idx1=which(sampletable1$phenoEA_bc==2) #case
  length(idx1) #2514
  idx2=which(sampletable1$phenoEA_bc==1)
  length(idx2) #2179
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  EA_p=pvalue_arow(idx1=idx1,idx2=idx2,covariates = covariates,predict_bcageneexp=predict_bcageneexp)
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
  BEA_p=pvalue_arow(idx1=idx1,idx2=idx2,covariates = covariates,predict_bcageneexp=predict_bcageneexp)
  
  #EA and BE vs control
  idx1=which(sampletable1$phenoBE_bc==2 | sampletable1$phenoEA_bc==2) #BE or EAcase
  length(idx1) #5802
  idx2=which(sampletable1$phenoBE_bc==1 |sampletable1$phenoEA_bc==1) #Control
  length(idx2) #2176
  length(intersect(idx1,idx2))
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$pc1[c(idx1,idx2)],pc2=sampletable1$pc2[c(idx1,idx2)],pc3=sampletable1$pc3[c(idx1,idx2)],
                        stringsAsFactors = F)
  covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  BEEA_p=pvalue_arow(idx1=idx1,idx2=idx2,covariates = covariates,predict_bcageneexp=predict_bcageneexp)
  res=data.frame(BE_p=BE_p,EA_p=EA_p,BEA_p=BEA_p,BEEA_p=BEEA_p)
  rownames(res)=rownames(predict_bcageneexp)
  return(res)
}
assoc_gene_res=association_gene()
assoc_gene_res
# BE_p        EA_p     BEA_p       BEEA_p
# DDX49 4.173502e-07 4.10393e-06 0.7449309 4.607661e-08

