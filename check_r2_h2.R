#to check difference of R-squared and H^2

setwd("/fh/fast/dai_j/BEACON/BEACON_GRANT/code")
#organ="stomach"
#gene="TMEM30A"


organ="junction"
gene="COX7A2"

organ="adipose"
gene="LDAH"

organ="muscularis"
gene="GGA1"

organ="junction"
gene="BORCS8"

organ="mucosa"
gene="SLC25A42"

#compute r-squared-------------

rdata=paste0("../result/GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
load(rdata)

removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  ## first remove redundant SNPs corr=1 ##
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= 1))]
  ## next remove high correlation SNPs but not exactly same ##
  tmp <- cor(datnew)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  dlist <- (1:ncol(datnew))[apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  rlist <- (1:ncol(datnew))[!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  clist <- NULL
  for (l in dlist) {
    clist <- c(clist, which(abs(tmp[,l])>=corcutoff))
  } 
  
  ## this is to make sure every deleted SNP has one tagSNP retained
  
  l<-1
  while (l<=length(dlist)) {
    ll <- dlist[l]
    tt <- which(abs(tmp[,ll])>=corcutoff)
    if (sum(tt %in% rlist)==0) {
      tcount <- rep(0,length(tt))
      for (k in 1:length(tt)) {
        tcount[k] <- sum(clist==tt[k])
      }
      rlist <- c(rlist,tt[which.max(tcount)])
      dlist <- dlist[-which(dlist==tt[which.max(tcount)])]
    } 
    l <- l+1
  } 
  rlist <- rlist[order(rlist)]
  datnew <- datnew[,rlist]
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

#new minimum rule,run 100 cvfit,distcutoff:dist to pick cis SNPs, corcutoff:to remove high correlated SNPs
# So I went through diagnosis: here are my choices
# 
# Use the last version of datasets with ambiguous data and : TPM+TMM+standardize
# Include covariates into model selection,
# Use standardize=T in glmnet
# Still do 100 CV to stabilize the selection
# Use correlation filtering threshold 0.9
compute_cor_arow=function(i,ncv=10,distcutoff=5e5,corcutoff=0.9)
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
  numsnp=NA
  selectedsnps=NA
  selectedsnps_coeff=NA
  numlowcorsnp=NA
  
  if (length(idx)>1)
  {
    tmp=quantile(Y,probs=c(0.15,0.85))
    if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
    
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
    if (class(X)[1]=="numeric") #only 1 snp left
    {
      X=matrix(X,nrow=length(X),ncol=1)
    }
    numlowcorsnp=ncol(X)
    #hist(ucor[colnames(X1)%in% colnames(X),1])
    #hist(ucor[colnames(X1)%in% colnames(X),2])
    #dim(X)
    
    Xall=data.matrix(cbind(X,covariate))
    
    ## compute the univariate eQTL association ##
    upval <- matrix(0,1,ncol(X))
    for (i in 1:ncol(X)){
      genotype <- X[,i]
      uXall <- data.frame(cbind(Y,genotype,covariate))
      ufit <- glm(Y~.,data=uXall) 
      upval[1,i] <- summary(ufit)$coef[2,4]
    }  
    
   
    upval <- data.frame(upval)
    colnames(upval)=colnames(X)
    
    
    upval <-  matrix(0,1,ncol(X))
    for (i in 1:ncol(X)){
      genotype <- X[,i]
      uXall <- data.frame(cbind(Y,genotype))
      ufit <- glm(Y~.,data=uXall) 
      upval[i] <- summary(ufit)$coef[2,4]
    }  
    
    upval <- data.frame(upval)
    colnames(upval)=colnames(X)
    
    
    covariateall=covariate
    
    penalty=rep(1,ncol(Xall))
    penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
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
      for (l in 1:ncol(cverr)) {
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
      selectcovariate=rownames(glmcoeff)[!rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0]
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
          #fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
          r2=cor(fitted,Y)^2
          #tmp=abs(cor(Y,Xsel))
          #hist(tmp,main=paste0(numvar," snps"),xlab="correlation between snps and expr")
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              numlowcorsnp=numlowcorsnp,fitted=fitted))
}

library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp

res1=compute_cor_arow(i=which(rownames(phenotype)==gene))
r2=res1$r2


#compute H2---------------------------------------------

plink="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink" ##
gctafolder="/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.93.2beta/" ##


genotypefolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/" ##
gtexprefix="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_flip_chr"  ## this is to include all the GTEx genotype
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor/") ##

#generate covariates,categorical data
generate_cov=function(outfile=paste0(outfolder,"GCTA.covar"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=5))
  res[,1]=res[,2]=rownames(covariate)
  res[,3]=as.character(covariate$gender)
  res[res[,3]==1,3]="M"
  res[res[,3]==2,3]="F"
  res[,4]=as.character(covariate$pcr)
  res[res[,4]==0,4]="PCR0"
  res[res[,4]==1,4]="PCR1"
  res[,5]=as.character(covariate$platform)
  res[res[,5]==0,5]="Platform0"
  res[res[,5]==1,5]="Platform1"
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}

#qov,numeric data,pcs+peer factors+age 
generate_qcov=function(outfile=paste0(outfolder,"GCTA.qcovar"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=3))
  res[,1]=res[,2]=rownames(covariate)
  res[,3]=covariate$age
  for (i in 1:4)
  {
    res=cbind(res,covariate[,paste0("pc",i)])
  }
  #no peer
  for (i in 1:15)
  {
    res=cbind(res,covariate[,paste0("factor",i)])
  }
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}

#generate fam select files to keep samples from the organ
generate_fam=function(outfile=paste0(outfolder,"GCTA.select.fam"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=2))
  res[,1]=res[,2]=rownames(covariate)
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}
generate_cov()
generate_qcov()
generate_fam()


#read result
readhsq=function(hsqfile=paste0(outfolder,"LDAH.hsq"),opt="mgrm")
{
  tmp1=NULL
  if (file.exists(hsqfile))
  {
    tmp=read.table(hsqfile,sep="\t",fill=T,stringsAsFactors = F,header = T)
    if (opt=="1grm")
    {
      idx1=which(tmp$Source=="V(G)/Vp")
      idx2=which(tmp$Source=="Pval")
      if (length(idx1)>0 & length(idx2)>0)
        tmp1=data.frame(H=tmp$Variance[c(idx1)],SE=tmp$SE[c(idx1)],pvalue=tmp$Variance[idx2],stringsAsFactors = F)
    }else #multiple grms
    {
      idx=which(grepl("Sum",tmp[,1]))
      if(length(idx)>0)
        tmp1=data.frame(H=tmp$Variance[idx],SE=tmp$SE[idx],pvalue=tmp$Variance[nrow(tmp)-1],stringsAsFactors = F)
    }
  }
  return(tmp1)
}

#estimate heritarbility
estimate_heri=function(genename="TMEM30A",distcutoff=5e5,opt="1grm")
{
  allres=NULL
  #phenotype file
  phenofile=paste0(outfolder,genename,".pheno")
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=3))
  res[,1]=res[,2]=rownames(covariate)
  idx=which(rownames(phenotype)==genename)
  res[,3]=unlist(phenotype[idx,])
  write.table(res,file=phenofile,sep=" ",row.names = F,col.names = F,quote=F)
  
  idx=which(rownames(phenotypepos)==genename)
  tmp=distance(gr_snp,gr_pos[idx])
  idx1=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx1,]))
  idx1=idx1[tmp!=0] #remove all 0 genotypes
  if (length(idx1)>10)
  {
    bimfile=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx],".bim")
    bim=as.data.frame(data.table::fread(bimfile))
    idx2=match(snppos$pos[idx1],bim$V4)
    snpfilename=paste0(outfolder,genename,".snpname")
    tmp=data.frame(snps=bim$V2[idx2])
    write.table(tmp,file=snpfilename,row.names = F,col.names = F,quote=F)
    famfile=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx],".fam")
    fam=read.table(famfile,stringsAsFactors = F)
    select.fam=paste0(outfolder,"GCTA.select.fam")
    
    prefix1=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx])
    genofilename=paste0(outfolder,genename)
    cmd=paste0(plink," --bfile ",prefix1," --keep ",select.fam," --extract ",snpfilename," --make-bed --out ", genofilename)
    system(cmd)
    tmp=read.table(paste0(genofilename,".fam"))
    if(any(tmp$V1!=rownames(covariate))) warning("the fam file has different sample order!")
    
    #run gcta, use mgrm
    if (opt=="mgrm")
    {
      cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --ld-score-region 200 --thread-num 12 --out ",genofilename)
      system(cmd)
      cmd=paste0("./stratify_GCTAsnps.R ",genofilename)
      system(cmd)
      if (file.exists(paste0(genofilename,".score.ld")))
      {
        if (file.size(paste0(genofilename,"_group1.txt"))>0 & file.size(paste0(genofilename,"_group2.txt"))>0 &file.size(paste0(genofilename,"_group3.txt"))>0 & file.size(paste0(genofilename,"_group4.txt"))>0)
        {
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group1.txt --make-grm --thread-num 12 --out ",genofilename,"_group1")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group2.txt --make-grm --thread-num 12 --out ",genofilename,"_group2")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group3.txt --make-grm --thread-num 12 --out ",genofilename,"_group3")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group4.txt --make-grm --thread-num 12 --out ",genofilename,"_group4")
          system(cmd)
          fileConn<-file(paste0(genofilename,"_mult_GRMs.txt"))
          writeLines(paste0(genofilename,"_group1\n",genofilename,"_group2\n",genofilename,"_group3\n",genofilename,"_group4"), fileConn)
          close(fileConn)
          cmd=paste0(gctafolder,"gcta64 --reml --mgrm ",genofilename,"_mult_GRMs.txt --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --reml-maxit 1000 --out ",genofilename)
          system(cmd)
          
          allres=readhsq(hsqfile = paste0(outfolder,genename,".hsq"))
          if (!is.null(allres)) rownames(allres)=genename
        }
      }
      #remove intermediate files
      cmd=paste0("rm ",outfolder,genename,".*")
      system(cmd)
      cmd=paste0("rm ",outfolder,genename,"_*.*")
      system(cmd)
    }
    
    #run gcta, use REML
    if (opt=="1grm")
    {
      cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --make-grm --thread-num 12 --out ",genofilename)
      system(cmd)
      #use AI by setting reml-alg 0, it is the default value
      cmd=paste0(gctafolder,"gcta64 --reml --reml-alg 0 --grm ",genofilename," --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --out ",genofilename)
      system(cmd)
      allres=readhsq(hsqfile = paste0(outfolder,genename,".hsq"),opt="1grm")
      if (!is.null(allres)) rownames(allres)=genename
    }
    cmd=paste0("rm ",outfolder,genename,".*")
    system(cmd)
  }
  
  return(allres)
}

h2=estimate_heri(genename = gene,opt="1grm")
#h2m=estimate_heri(genename = gene,opt="mgrm")

#get H2/R-squared table for all the results
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
organs=c("adipose","blood","junction","mucosa","muscularis","stomach")
tmp=data.frame(matrix(NA,nrow=20000,ncol=12))
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  tmp1=rownames(res_min)[rownames(res_min) %in% proteingenes & !is.na(res_min$r2)]
  tmp[1:length(tmp1),(i-1)*2+1]=tmp1
  outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  tmp1=rownames(heritability)[rownames(heritability) %in% proteingenes]
  tmp[1:length(tmp1),(i-1)*2+2]=tmp1
}
allgenes=tmp[,1]
for (i in 2:12)
{
  allgenes=unique(c(allgenes,tmp[,i]))
}
allgenes=allgenes[!is.na(allgenes)]
tmp=data.frame(matrix(NA,nrow=length(allgenes),ncol=12))
colnames(tmp)[seq(1,ncol(tmp),2)]=paste0(organs,"_R2")
colnames(tmp)[seq(2,ncol(tmp),2)]=paste0(organs,"_H2")
rownames(tmp)=allgenes
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  tmp1=intersect(allgenes,rownames(res_min))
  idx1=match(tmp1,allgenes)
  idx2=match(tmp1,rownames(res_min))
  tmp[idx1,(i-1)*2+1]=res_min$r2[idx2]
  outfolder1=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor")
  heritfile=paste0(outfolder1,"/heritability_","1grm",".txt")
  heritability=read.table(heritfile,header = T)
  tmp1=intersect(allgenes,rownames(heritability))
  idx1=match(tmp1,allgenes)
  idx2=match(tmp1,rownames(heritability))
  tmp[idx1,(i-1)*2+2]=heritability$V1[idx2]
}
tmp1=sapply(1:nrow(tmp), function(x) sum(!is.na(tmp[x,])))
idx=order(tmp1,decreasing = T)
write.csv(tmp[idx,],file="../result/All6organs_R2_H2.csv")
