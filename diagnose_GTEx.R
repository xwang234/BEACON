
##############################################################################
### diagnosing the differences of three GTEx datasets in prediction models ###
### June 11, 2020
##############################################################################

rm(list=ls())

removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}

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
    Xsel=Xsel[,idx1]
    Xsel=Xsel[,1:(maxnumvar-ncol(covariateall)-1)]
    
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


library(glmnet)
library(GenomicRanges)

genename <- "HSP90AA1"
genename <- "ISYNA1"
genename <- "CERS1"
genename <- "KXD1"
genename <- "UBAC1"
genename <- "FOXF1"
genename <- "ZPLD1"

genename <- "RAB34"

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")

snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp


#compute
#test=compute_cor_arow(i=31038,opt="min",ncv=10,distcutoff = 5e5)
#which(row.names(phenotype)=="THAP6")
#[1] 5657
#which(row.names(phenotype)=="CERS1")
#[1] 31036
#which(row.names(phenotype)=="KXD1")
#i <- 31038
#which(row.names(phenotype)=="HSP90AA1")
#i <- 25105
#which(row.names(phenotype)=="UBAC1")
#i <- 17803
#which(row.names(phenotype)=="FBP2")
#i <- 17243

#which(row.names(phenotype)=="HSP90AA1")
#genename <- "HSP90AA1"
#genename <- "ISYNA1"

#genename <- "ZPLD1"

  i <- which(row.names(phenotype)==genename)
  ncv=10
  distcutoff = 5e5
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
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  
  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  }  
  hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  dim(X1)
  X <- removehighcorr(X1,0.9)
  hist(ucor[colnames(X1)%in% colnames(X),1])
  hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      #I used standardize
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  
  #cverr <- matrix(NA,length(cvfit$lambda),100)
  #cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  #rownames(cverr)=cvfit$lambda
  #for (l in 1:100) {
  #  set.seed(l+100)
  #  fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #  alllambda=intersect(cvfit$lambda,fit$lambda)
  #  idx1=match(alllambda,cvfit$lambda)
  #  idx2=match(alllambda,fit$lambda)
  #  cverr[idx1,l] <- fit$cvm[idx2]
  #  cvse[idx1,l] <- fit$cvsd[idx2]
  #}
  
  #merr <- apply(cverr,1,mean,na.rm=T)
  #mcvse <- sqrt(apply(cvse^2,1,mean))/10
  
  #lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  #lambda.min <- cvfit$lambda[which.min(merr)]
  
  #plot(log(cvfit$lambda),merr)
  #abline(v=log(lambda.1se),lty=2)
  #abline(v=log(lambda.min),lty=3)
  
  #plot(cvfit$lambda,merr)
  
  glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.min))
  sum(glmcoeff[,1]!=0)
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) #else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  Xsel <- X[,colnames(X)%in%selsnps]
  pred.Y <- fitted_cv(Xsel,covariate,Y,ncv=5)
  cor(Y,pred.Y)
  
  gene1 <- Y
  genotype1 <- X1
  cov1 <- covariate
  
  
  
  
  load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondatafor_prediction.RData")
  
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  i <- which(row.names(phenotype)==genename)
  ncv=10
  distcutoff = 5e5
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
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  
  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  }  
  hist(ucor[,1])

  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  dim(X1)
  X <- removehighcorr(X1,0.9)
  hist(ucor[colnames(X1)%in% colnames(X),1])
  hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      #I used standardize
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  #cverr <- matrix(NA,length(cvfit$lambda),100)
  #cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  #rownames(cverr)=cvfit$lambda
  #for (l in 1:100) {
  #  set.seed(l+100)
  #  fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
  #  alllambda=intersect(cvfit$lambda,fit$lambda)
  #  idx1=match(alllambda,cvfit$lambda)
  #  idx2=match(alllambda,fit$lambda)
  #  cverr[idx1,l] <- fit$cvm[idx2]
  #  cvse[idx1,l] <- fit$cvsd[idx2]
  #}
  
  #merr <- apply(cverr,1,mean,na.rm=T)
  #mcvse <- sqrt(apply(cvse^2,1,mean))/10
  
  #lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  #lambda.min <- cvfit$lambda[which.min(merr)]
  
  #plot(log(cvfit$lambda),merr)
  #abline(v=log(lambda.1se),lty=2)
  #abline(v=log(lambda.min),lty=3)
  
  #plot(cvfit$lambda,merr)
  
  glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.min))
  sum(glmcoeff[,1]!=0)
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) #else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  Xsel <- X[,colnames(X)%in%selsnps]
  pred.Y <- fitted_cv(Xsel,covariate,Y,ncv=5)
  cor(Y,pred.Y)
  
  
  gene2 <- Y
  genotype2 <- X1
  cov2 <- covariate
  
  
  
  load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")
  
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  i <- which(row.names(phenotype)==genename)
  ncv=10
  distcutoff = 5e5
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
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  
  X1=t(snp[idx,])
  ucor <- matrix(0,ncol(X1),2)
  for (l in 1:ncol(X1)){
    ucov<- data.matrix(cbind(X1[,l],covariate))
    ufit <- lm(Y~ucov)
    ucor[l,1] <- summary(ufit)$coef[2,4]
    ucor[l,2] <- cor(Y,X1[,l])
  }  
  hist(ucor[,1])
  
  pcor <- ucor[,1]
  X1 <- X1[,order(pcor,decreasing=T)]
  dim(X1)
  X <- removehighcorr(X1,0.9)
  hist(ucor[colnames(X1)%in% colnames(X),1])
  hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  covariate <- covariate[,-(24:27)]
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      #I used standardize
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  
  cverr <- matrix(NA,length(cvfit$lambda),100)
  cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  rownames(cverr)=cvfit$lambda
  for (l in 1:100) {
    set.seed(l+100)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    alllambda=intersect(cvfit$lambda,fit$lambda)
    idx1=match(alllambda,cvfit$lambda)
    idx2=match(alllambda,fit$lambda)
    cverr[idx1,l] <- fit$cvm[idx2]
    cvse[idx1,l] <- fit$cvsd[idx2]
  }
  
  merr <- apply(cverr,1,mean,na.rm=T)
  mcvse <- sqrt(apply(cvse^2,1,mean))/10
  
  lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  lambda.min <- cvfit$lambda[which.min(merr)]
  
  plot(log(cvfit$lambda),merr)
  abline(v=log(lambda.1se),lty=2)
  abline(v=log(lambda.min),lty=3)
  
  glmcoeff=as.matrix(coef(cvfit,s=lambda.min))
  sum(glmcoeff[,1]!=0)
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) #else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  Xsel <- X[,colnames(X)%in%selsnps]
  Xsel <- as.matrix(Xsel)
  pred.Y <- fitted_cv(Xsel,covariate,Y,ncv=5)
  cor(Y,pred.Y,method="spearman")
  
  gene3 <- Y
  genotype3 <- X1
  cov3 <- covariate
  
  
  ### first investigate the differences of gene expression
  mean(names(gene3)==names(gene2))
  mean(names(gene1)%in%names(gene2))  ## about ~20 sample names not shown up in later sampes ##
  mean(names(gene2)%in%names(gene1))

  hist(gene1)
  hist(gene2)
  hist(gene3)
  
  plot(gene2,gene3)  
  idx <- match(names(gene1),names(gene2))
  mean(names(gene2[idx[!is.na(idx)]])==names(gene1[!is.na(idx)]))
  
  plot(gene2[idx[!is.na(idx)]],gene1[!is.na(idx)])      ## observed outliers in gene1
  plot(gene2[idx[!is.na(idx)]],gene1[!is.na(idx)],ylim=c(-1,2))
  abline(0,1)
  plot(gene3[idx[!is.na(idx)]],gene1[!is.na(idx)])
  plot(gene3[idx[!is.na(idx)]],gene1[!is.na(idx)],ylim=c(-2,2),xlim=c(-2,2))
  abline(0,1)
  
  
  cor(gene2[idx[!is.na(idx)]],gene1[!is.na(idx)])
  cor(gene2[idx[!is.na(idx)]],gene1[!is.na(idx)],method="spearman")
  
  
  cor(gene3[idx[!is.na(idx)]],gene1[!is.na(idx)])
  cor(gene3[idx[!is.na(idx)]],gene1[!is.na(idx)],method="spearman")
  
  ## now investigate the differences of genotypes ##
  
  dim(genotype1)
  dim(genotype2)
  dim(genotype3)
  mean(colnames(genotype1) %in% colnames(genotype2))  ## genotype1 names did not match, likely due to build
  #[1] 0
  mean(colnames(genotype2) %in% colnames(genotype3))
  #[1] 0.6015038
  ### now investigate the covariate ###

  dim(cov1)
  dim(cov2)
  dim(cov3)
  
  colnames(cov1)
  colnames(cov2)
  colnames(cov3)  ## cov3 messed up, 4 repeat columns
  
  mean(rownames(cov1)%in%rownames(cov2))
  mean(rownames(cov2)==rownames(cov3))

  idx <- match(rownames(cov2),rownames(cov1))
  cor(cov2[!is.na(idx),]$pc1,cov1[idx[!is.na(idx)],]$pc1)
  
  cor(cov2[!is.na(idx),]$factor1,cov1[idx[!is.na(idx)],]$factor1)
  cor(cov2[!is.na(idx),]$factor5,cov1[idx[!is.na(idx)],]$factor5)
  cor(cov2[!is.na(idx),]$factor2,cov1[idx[!is.na(idx)],]$factor2)
  cor(cov2[!is.na(idx),]$factor3,cov1[idx[!is.na(idx)],]$factor3)
  cor(cov2[!is.na(idx),]$factor4,cov1[idx[!is.na(idx)],]$factor4)
  
  cor(cov2[!is.na(idx),]$factor15,cov1[idx[!is.na(idx)],]$factor15)
  mean(row.names(cov2[!is.na(idx),])==row.names(cov1[idx[!is.na(idx)],]))
  
  mean(cov2[!is.na(idx),]$gender==cov1[idx[!is.na(idx)],]$gender)
  mean(cov2[!is.na(idx),]$age==cov1[idx[!is.na(idx)],]$age)
  
  
  
  
  
  
  