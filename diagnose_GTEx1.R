
##############################################################################
### diagnosing the differences of three GTEx datasets in prediction models ###
### June 11, 2020
##############################################################################

rm(list=ls())

#check overlap of V7 SNPs and V8 SNPs
#hg38 --->hg19---
hg38tohg19=function(snpnames=rownames(genotypeV8))
{
  library(rtracklayer)
  library(GenomicRanges)
  chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
  tmp=unlist(strsplit(snpnames,":"))
  chr0=tmp[1]
  chr=paste0("chr",tmp[1])
  tmp=tmp[seq(2,length(tmp),2)]
  tmp=unlist(strsplit(tmp,"_"))
  pos=as.integer(tmp[seq(1,length(tmp),3)])
  alt=tmp[seq(2,length(tmp),3)]
  ref=tmp[seq(3,length(tmp),3)]
  gr_dat=GRanges(seqnames = chr,ranges=IRanges(start=pos,width=1))
  tmp=liftOver(gr_dat,chain)
  newsnpnames=newpos=rep(NA,length(tmp))
  for (i in 1:length(tmp))
  {
    tmp1=unlist(tmp[i])
    if (length(tmp1)==0)
    {
      warning(paste0(snpnames[i]," not transformed"))
    }else
    {
      if (length(tmp1)==1)
      {
        newpos[i]=start(tmp1)
      }else
      {
        warning(paste0(snpnames[i]," transformed to ",length(tmp1)," snps"))
        newpos[i]=start(tmp1)[1]
      }
    }
  }
  newsnpnames=paste0(chr0,":",newpos,"_",alt,"_",ref)
  res=data.frame(snphg38=snpnames,snphg19=newsnpnames,stringsAsFactors = F)
  return(res)
}

#check overlap of V8 snps and V7 snps, put V8 first---
overlap_snp1_snp2=function(snpnames1=rownames(genotypeV8),snpnames2=rownames(genotypeV7))
{
  newsnpnames1=hg38tohg19(snpnames=snpnames1)
  newsnpnames1=newsnpnames1[!is.na(newsnpnames1$snphg19),]
  newsnpnames1_hg19=newsnpnames1$snphg19
  sum(newsnpnames1_hg19 %in% snpnames2)
  #some thimes the above not working, they may sit on different strands (A/G and T/C), to use positon to match
  posnames1=rep(NA,length(newsnpnames1_hg19))
  posnames2=rep(NA,length(snpnames2))
  get_posname=function(snpnames=newsnpnames1_hg19)
  {
    tmp=unlist(strsplit(snpnames,"_"))
    posnames=tmp[seq(1,length(tmp),3)]
    return(posnames)
  }
  posnames1=get_posname()
  posnames2=get_posname(snpnames=snpnames2)
  newsnpnames1$posnames=posnames1
  newsnpnames2=data.frame(snphg19=snpnames2,posnames=posnames2,stringsAsFactors = F)
  comsnps=intersect(posnames1,posnames2)
  print(paste0("snp1: ",length(snpnames1),", snp2: ",length(snpnames2),", common snp: ",length(comsnps)))
  newsnpnames1$overlap=F
  idx=match(comsnps,posnames1)
  newsnpnames1$overlap[idx]=T
  idx=match(comsnps,posnames2)
  newsnpnames2$overlap=F
  newsnpnames2$overlap[idx]=T
  return(list(snpnames1=newsnpnames1,snpnames2=newsnpnames2))
}
#extract overlap snp data, put V8 first,used to compare snp data---
extract_overlapsnpdat=function(snpnames1=rownames(genotypeV8),snpnames2=rownames(genotypeV7))
{
  olapsnpnames=overlap_snp1_snp2(snpnames1,snpnames2)
  idx=which(olapsnpnames$snpnames1$overlap==T)
  idx1=match(olapsnpnames$snpnames1$snphg38[idx],rownames(snpV8))
  snpdat1=snpV8[idx1,]
  rownames(snpdat1)=olapsnpnames$snpnames1$posnames[idx]
  idx=which(olapsnpnames$snpnames2$overlap==T)
  idx1=match(olapsnpnames$snpnames2$snphg19[idx],rownames(snpV7))
  snpdat2=snpV7[idx1,]
  rownames(snpdat2)=olapsnpnames$snpnames2$posnames[idx]
  idx=match(rownames(snpdat1),rownames(snpdat2))
  snpdat2=snpdat2[idx,]
  tmp=intersect(colnames(snpdat1),colnames(snpdat2))
  idx1=match(tmp,colnames(snpdat1))
  idx2=match(tmp,colnames(snpdat2))
  snpdat1=snpdat1[,idx1]
  snpdat2=snpdat2[,idx2]
  return(list(snpdat1=snpdat1,snpdat2=snpdat2))
}
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


compute_cor_arow=function(i,ncv=10,distcutoff=5e5,corcutoff=0.9,phenotype=NA,snp=NA,covariate=NA,gr_snp=NA,gr_pos=NA)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numlowcorsnp=0 #number of SNP after removing high correlation SNPs
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  selectcov=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  if (length(idx)>1)
  {
    #too many highly correlated SNPs are included for glmnet
    X1=t(snp[idx,])
    ucor <- matrix(0,ncol(X1),2)
    for (l in 1:ncol(X1)){
      ucov<- data.matrix(cbind(X1[,l],covariate))
      ufit <- lm(Y~ucov)
      ucor[l,1] <- summary(ufit)$coef[2,4]
      ucor[l,2] <- cor(Y,X1[,l])
    } 
    
    #hist(ucor[,1])
    pcor <- ucor[,1]
    ## I order the snps based on their correlation with gene expression, first delete those low correlation-with-gene SNPs
    X1 <- X1[,order(pcor,decreasing=T)]
    X <- removehighcorr(X1,corcutoff)
    if (class(X)=="numeric") #only 1 snp left
    {
      X=matrix(X,nrow=length(X),ncol=1)
    }
    numlowcorsnp=ncol(X)
    #hist(ucor[colnames(X1)%in% colnames(X),1])
    #hist(ucor[colnames(X1)%in% colnames(X),2])
    #dim(X)
    
    Xall=data.matrix(cbind(X,covariate))
    #not include covariate
    #Xall=X
    covariateall=covariate
    
    penalty=rep(1,ncol(Xall))
    #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
    set.seed(i+10000)
    cvfit=tryCatch(
      {
        ### I change alpha to 0.5 for better variable selection when highly correlated features
        cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
      },
      error=function(e)
      {
        return(F)
      }
    )
    plot(cvfit)
    if (is.list(cvfit))
    {
      ## do 100 times cv.glmnet and take average for cverr
      ## the number of cvfit$lambda may be less than 100 sometimes even you specified 100
      cverr <- matrix(NA,length(cvfit$lambda),100)
      rownames(cverr)=cvfit$lambda
      for (l in 1:100) {
        set.seed(l+100)
        fit=tryCatch(
          {
            ### I change alpha to 0.5 for better variable selection when highly correlated features
            cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
          },
          error=function(e)
          {
            return(F)
          }
        )
        if (is.list(fit)) #Error in predmat[which, seq(nlami)] = preds : replacement has length zero
        {
          #fit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
          alllambda=intersect(cvfit$lambda,fit$lambda)
          idx1=match(alllambda,cvfit$lambda)
          idx2=match(alllambda,fit$lambda)
          cverr[idx1,l] <- fit$cvm[idx2]
        }
      }
      merr <- apply(cverr,1,mean,na.rm=T)
      #plot(log(cvfit$lambda),merr)
      lambda.best <- cvfit$lambda[which.min(merr)]
      #fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=0.5,standardize=T)
      glmcoeff=as.matrix(coef(cvfit,s=lambda.best))
      #sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
      #glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
      #the selected covariate
      selectcov=selectcovariate=rownames(glmcoeff)[!rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0]
      selectcovariate=selectcovariate[selectcovariate!="(Intercept)"]
      selectcovariate=covariateall[,colnames(covariateall) %in% selectcovariate,drop=F]
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
          nummaxvar=min(nrow(X)-ncol(selectcovariate)-1,numvar)
          numvar=nummaxvar
          Xsel=Xsel[,1:nummaxvar]
        }
        if (numvar==1) #keep Xsel as in matrix form
        {
          Xsel=matrix(Xsel,ncol=1)
          colnames(Xsel)=glmleftsnp
        }
        Xall1=data.matrix(cbind(Xsel,selectcovariate))
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
          fitted=fitted_cv(Xsel,selectcovariate,Y,ncv=ncv)
          r2=cor(fitted,Y)^2
          #tmp=abs(cor(Y,Xsel))
          #hist(tmp,main=paste0(numvar," snps"),xlab="correlation between snps and expr")
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              numlowcorsnp=numlowcorsnp,allsnps=paste0(rownames(snp)[idx],collapse = "|"),
              selectcov=paste0(selectcov,collapse="|")))
}

genename <- "KXD1"
genename <- "HSP90AA1"
genename <- "ISYNA1"
genename <- "CERS1"
genename <- "UBAC1"
genename <- "FOXF1"


load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData") #V7
snpV7=snp
phenotypeV7=phenotype
snppos$chr[snppos$chr==23]="X"
snpposV7=snppos
phenotypeposV7=phenotypepos
covariateV7=covariate
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snpV7=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_posV7=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
V7=compute_cor_arow(i=which(rownames(phenotype)==genename),ncv=10,distcutoff = 5e5,corcutoff = 0.9)
V7[c("r2","numvar","numsnp","numlowcorsnp")]
# $r2
# [1] 0.07076597
# 
# $numvar
# [1] 48
# 
# $numsnp
# [1] 2006
# 
# $numlowcorsnp
# [1] 420

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8junctiondata_ambiguous_TPM_for_prediction.RData")
snpV8=snp
phenotypeV8=phenotype
snpposV8=snppos
phenotypeposV8=phenotypepos
covariateV8=covariate
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snpV8=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_posV8=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/tmp_junction_SNP.RData")
snpV8_clean=snp1
snpposV8_clean=snppos1
gr_snpV8_clean=gr_snp=GRanges(seqnames = snpposV8_clean$chr,ranges=IRanges(start=snpposV8_clean$pos,width = 1)) #SNP

genes=c("KXD1","HSP90AA1","ISYNA1","CERS1","UBAC1","FOXF1")


V8res=data.frame(tpm_snp=rep(NA,6),tpm_r2=NA,rpkm_snp=NA,rpkm_r2=NA,com_snp=NA,stringsAsFactors = F)
rownames(V8res)[1:6]=genes
for (i in 1:length(genes))
{
  gene=genes[i]
  print(gene)
  #use tpm+tmm gene expr
  phenotype=phenotypeV8
  covariate=covariateV8
  # V8=compute_cor_arow(i=which(rownames(phenotype)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,)
  # V8res$tpm_snp[i]=V8$numvar
  # V8res$tpm_r2[i]=V8$r2
  #use quantile normalization gene expr
  V8_q=compute_cor_arow(i=which(rownames(phenotype)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,
                        phenotype=phenotype_quantile,snp=snpV8,covariate=covariate_quantile,
                        gr_snp=gr_snpV8,gr_pos=gr_posV8)
  V8res$rpkm_snp[i]=V8_q$numvar
  V8res$rpkm_r2[i]=V8_q$r2
  # V8res$com_snp[i]=sum(unlist(strsplit(V8$selectedsnps,"|",fixed = T)) %in% unlist(strsplit(V8_q$selectedsnps,"|",fixed = T)))
}

V7V8res=data.frame(V7_snp=rep(NA,6),V7_r2=NA,V7_totalsnp=NA,V8_snp=NA,V8_r2=NA,V8_totalsnp=NA,com_snp=NA,com_totalsnp=NA,cor_expr=NA,stringsAsFactors = F)
rownames(V7V8res)=genes
for (i in 1:length(genes))
{
  gene=genes[i]
  print(gene)
  V8=compute_cor_arow(i=which(rownames(phenotypeV8)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,
                        phenotype=phenotypeV8,snp=snpV8,covariate=covariateV8,
                        gr_snp=gr_snpV8,gr_pos=gr_posV8)
  V7=compute_cor_arow(i=which(rownames(phenotypeV7)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,
                      phenotype=phenotypeV7,snp=snpV7,covariate=covariateV7,
                      gr_snp=gr_snpV7,gr_pos=gr_posV7)
  
  V7V8res$V7_snp[i]=V7$numvar
  V7V8res$V7_r2[i]=V7$r2
  V7V8res$V7_totalsnp[i]=V7$numsnp
  V7V8res$V8_snp[i]=V8$numvar
  V7V8res$V8_r2[i]=V8$r2
  V7V8res$V8_totalsnp[i]=V8$numsnp
  comsamples=intersect(colnames(phenotypeV7),colnames(phenotypeV8))
  idx1=which(rownames(phenotypeV7)==gene)
  idx2=which(rownames(phenotypeV8)==gene)
  V7V8res$cor_expr[i]=cor(unlist(phenotypeV7[idx1,match(comsamples,colnames(phenotypeV7))]),
                          unlist(phenotypeV8[idx2,match(comsamples,colnames(phenotypeV8))]))
  if (V7$numvar*V8$numvar>0)
  {
    olapsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$selectedsnps,"|",fixed = T)),
                                 snpnames2=unlist(strsplit(V7$selectedsnps,"|",fixed = T)))
    V7V8res$com_snp[i]=sum(olapsnpnames$snpnames1$overlap==T)
    olapallsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$allsnps,"|",fixed = T)),
                                      snpnames2=unlist(strsplit(V7$allsnps,"|",fixed = T)))
    V7V8res$com_totalsnp[i]=sum(olapallsnpnames$snpnames1$overlap==T)
  }
  
}
V7V8res_withcovriate=V7V8res
V7V8res_without=V7V8res

V7V8res=data.frame(V7_snp=rep(NA,6),V7_r2=NA,V7_totalsnp=NA,V8_snp=NA,V8_r2=NA,V8_totalsnp=NA,com_snp=NA,com_totalsnp=NA,cor_expr=NA,stringsAsFactors = F)
rownames(V7V8res)=genes
for (i in 1:length(genes))
{
  gene=genes[i]
  print(gene)
  V8=compute_cor_arow(i=which(rownames(phenotypeV8)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,
                      phenotype=phenotypeV8,snp=snpV8_clean,covariate=covariateV8,
                      gr_snp=gr_snpV8_clean,gr_pos=gr_posV8)
  V7=compute_cor_arow(i=which(rownames(phenotypeV7)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9,
                      phenotype=phenotypeV7,snp=snpV7,covariate=covariateV7,
                      gr_snp=gr_snpV7,gr_pos=gr_posV7)
  
  V7V8res$V7_snp[i]=V7$numvar
  V7V8res$V7_r2[i]=V7$r2
  V7V8res$V7_totalsnp[i]=V7$numsnp
  V7V8res$V8_snp[i]=V8$numvar
  V7V8res$V8_r2[i]=V8$r2
  V7V8res$V8_totalsnp[i]=V8$numsnp
  comsamples=intersect(colnames(phenotypeV7),colnames(phenotypeV8))
  idx1=which(rownames(phenotypeV7)==gene)
  idx2=which(rownames(phenotypeV8)==gene)
  V7V8res$cor_expr[i]=cor(unlist(phenotypeV7[idx1,match(comsamples,colnames(phenotypeV7))]),
                          unlist(phenotypeV8[idx2,match(comsamples,colnames(phenotypeV8))]))
  if (V7$numvar*V8$numvar>0)
  {
    olapsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$selectedsnps,"|",fixed = T)),
                                   snpnames2=unlist(strsplit(V7$selectedsnps,"|",fixed = T)))
    V7V8res$com_snp[i]=sum(olapsnpnames$snpnames1$overlap==T)
    olapallsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$allsnps,"|",fixed = T)),
                                      snpnames2=unlist(strsplit(V7$allsnps,"|",fixed = T)))
    V7V8res$com_totalsnp[i]=sum(olapallsnpnames$snpnames1$overlap==T)
  }
}

#check overlap of SNPs, for selected snps
olapsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$selectedsnps,"|",fixed = T)),
                               snpnames2=unlist(strsplit(V7$selectedsnps,"|",fixed = T)))
#[1] "snp1: 33, snp2: 48, common snp: 9"
#check overlap of SNPs, for all snps
olapallsnpnames=overlap_snp1_snp2(snpnames1=unlist(strsplit(V8$allsnps,"|",fixed = T)),
                               snpnames2=unlist(strsplit(V7$allsnps,"|",fixed = T)))
#[1] "snp1: 1721, snp2: 2006, common snp: 1700"

allV7V8snpdat=extract_overlapsnpdat(snpnames1=unlist(strsplit(V8$allsnps,"|",fixed = T)),
                                snpnames2=unlist(strsplit(V7$allsnps,"|",fixed = T)))
agreement=rep(0,nrow(allV7V8snpdat$snpdat1))
for (i in 1:length(agreement))
{
  agreement[i]=sum(as.numeric(allV7V8snpdat$snpdat1[i,])==as.numeric(allV7V8snpdat$snpdat2[i,]))/ncol(allV7V8snpdat$snpdat1)
}
quantile(agreement,c(0,0.05,0.1,0.25,1))
# 0%        5%       10%       25%      100% 
# 0.4325843 0.9775281 0.9887640 1.0000000 1.0000000

comsamples=intersect(colnames(snpV7),colnames(snpV8))
idx=match(comsamples,colnames(snpV7))
com_snpV7=snpV7[,idx]
com_phenotypeV7=phenotypeV7[,idx]
com_covariateV7=covariateV7[idx,]
idx=match(comsamples,colnames(snpV8))
com_snpV8=snpV8[,idx]
com_phenotypeV8=phenotypeV8[,idx]
com_covariateV8=covariateV8[idx,]

snp=com_snpV7
phenotype=com_phenotypeV7
snppos=snpposV7
phenotypepos=phenotypeposV7
covariate=com_covariateV7
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
com_V7=compute_cor_arow(i=which(rownames(phenotype)==genename),ncv=10,distcutoff = 5e5,corcutoff = 0.9)

snp=com_snpV8
phenotype=com_phenotypeV8
snppos=snpposV8
phenotypepos=phenotypeposV8
covariate=com_covariateV8
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
com_V8=compute_cor_arow(i=which(rownames(phenotype)==genename),ncv=10,distcutoff = 5e5,corcutoff = 0.9)

#check covariate
all(com_covariateV7$age==com_covariateV8$age)
all(com_covariateV7$gender==com_covariateV8$gender)
#top pc
tmp=cor(com_covariateV7[,1:4],com_covariateV8[,1:4])
heatmap(abs(tmp),Colv = NA,Rowv = NA)
diag(tmp)
# pc1         pc2         pc3         pc4 
# 0.94947572  0.83352734 -0.04392324 -0.02670200 
tmp[3,4]
tmp=cor(com_covariateV7[,5:19],com_covariateV8[,5:19])
heatmap(abs(tmp),Colv = NA,Rowv = NA)
diag(tmp)
# factor1     factor2     factor3     factor4     factor5     factor6     factor7     factor8     factor9 
# -0.97150447 -0.62768052 -0.45287291 -0.11345059  0.04277765  0.28300305  0.03158802 -0.34103435 -0.39802461 
# factor10    factor11    factor12    factor13    factor14    factor15 
# -0.46266729  0.12912521 -0.23609532  0.15313051  0.18516262  0.01869393 


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
  
  
  
  #plot(cvfit$lambda,merr)
  
  glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.min))
  sum(glmcoeff[,1]!=0)
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) #else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  Xsel <- X[,colnames(X)%in%selsnps]
  Xsel <- as.matrix(Xsel)
  pred.Y <- fitted_cv(Xsel,covariate,Y,ncv=5)
  cor(Y,pred.Y)
  
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
  
  
  
  
  
  
  