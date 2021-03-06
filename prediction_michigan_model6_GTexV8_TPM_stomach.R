#!/usr/bin/env Rscript
#USE GTEx stomach V8 data, including ambiguous SNPs, use TPM+TMM gene expression

# #get all the data
#generated by read_GTEx_EC.R
#snp,snppos,rownames:1:1234_C_T
#to get the results based on 1se rule
#use the code from James 4/18/2020

#salloc -t 6-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 Rscript ./prediction_michigan_model6_GTexV8_TPM_stomach.R &
#regular TWAS-step 1----------------
removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
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
  selectedsnps=NA
  selectedsnps_coeff=NA
  
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
    #plot(cvfit)
    if (is.list(cvfit))
    {
      ## do 10 times cv.glmnet and take average for cverr
      ## the number of cvfit$lambda may be less than 100 sometimes even you specified 100
      cverr <- matrix(NA,length(cvfit$lambda),10)
      rownames(cverr)=cvfit$lambda
      for (l in 1:10) {
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
          r2=cor(fitted,Y)^2
          #tmp=abs(cor(Y,Xsel))
          #hist(tmp,main=paste0(numvar," snps"),xlab="correlation between snps and expr")
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              numlowcorsnp=numlowcorsnp))
}

# compute_cor_arow_fusion=function(i,ncv=10,distcutoff=5e5)
# {
#   
#   Y=unlist(phenotype[i,]) #geneexp
#   r2=NA
#   glmflag=0 #if glm selected variables
#   tmp=distance(gr_snp,gr_pos[i])
#   idx=which(tmp<distcutoff)
#   tmp=rowSums(data.matrix(snp[idx,]))
#   idx=idx[tmp!=0] #remove all 0 genotypes
#   numvar=0 #number of snp selected by glmnet
#   selectedsnps=NA
#   selectedsnps_coeff=NA
#   p_gender=NA
#   p_age=NA
#   tmp=quantile(Y,probs=c(0.15,0.85))
#   if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
#   if (length(idx)>1)
#   {
#     #too many highly correlated SNPs are included for glmnet
#     X1=t(snp[idx,])
#     ucor <- matrix(0,ncol(X1),2)
#     for (l in 1:ncol(X1)){
#       ucov<- data.matrix(cbind(X1[,l],covariate))
#       ufit <- lm(Y~ucov)
#       ucor[l,1] <- summary(ufit)$coef[2,4]
#       ucor[l,2] <- cor(Y,X1[,l])
#     } 
#     
#     #hist(ucor[,1])
#     pcor <- ucor[,1]
#     ## I order the snps based on their correlation with gene expression, first delete those low correlation-with-gene SNPs
#     X1 <- X1[,order(pcor,decreasing=T)]
#     X <- removehighcorr(X1,0.8)
#     if (class(X)=="numeric") #only 1 snp left
#     {
#       X=matrix(X,nrow=length(X),ncol=1)
#     }
#     
#     #fusion add this:----
#     X2=scale(X)
#     #regress covariates out of the genotypes as well (this is more accurate but slower)
#     for (i in 1:ncol(X))
#     {
#       dat=data.frame(y=X[,i],covariate)
#       X2[,i]= summary(lm(y~.,data= dat ))$resid
#     }
#     X=scale(X2)
#     #remove variance in phenotype explained by covariates
#     dat=data.frame(y=Y,covariate)
#     Y2 = summary(lm(y~., data=dat ))$resid
#     Y=scale(Y2)
#     #end adding
#     
#     #Xall=data.matrix(cbind(X,covariate))
#     
#     #covariateall=covariate
#     
#     #penalty=rep(1,ncol(Xall))
#     #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
#     set.seed(i+10000)
#     cvfit=tryCatch(
#       {
#         ### I change alpha to 0.5 for better variable selection when highly correlated features
#         cv.glmnet(data.matrix(X),Y,nlambda=100,nfolds=10,alpha=0.5,standardize=T)
#       },
#       error=function(e)
#       {
#         return(F)
#       }
#     )
#     
#     if (is.list(cvfit))
#     {
#       ## do 100 times cv.glmnet and take average for cverr
#       ## the number of cvfit$lambda may be less than 100 sometimes even you specified 100
#       cverr <- matrix(NA,length(cvfit$lambda),100)
#       rownames(cverr)=cvfit$lambda
#       for (l in 1:100) {
#         set.seed(l+100)
#         fit=tryCatch(
#           {
#             ### I change alpha to 0.5 for better variable selection when highly correlated features
#             cv.glmnet(data.matrix(X),Y,nlambda=100,nfolds=10,alpha=0.5,standardize=T)
#           },
#           error=function(e)
#           {
#             return(F)
#           }
#         )
#         if (is.list(fit)) #Error in predmat[which, seq(nlami)] = preds : replacement has length zero
#         {
#           #fit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
#           alllambda=intersect(cvfit$lambda,fit$lambda)
#           idx1=match(alllambda,cvfit$lambda)
#           idx2=match(alllambda,fit$lambda)
#           cverr[idx1,l] <- fit$cvm[idx2]
#         }
#       }
#       merr <- apply(cverr,1,mean,na.rm=T)
#       #plot(log(cvfit$lambda),merr)
#       lambda.best <- cvfit$lambda[which.min(merr)]
#       #fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=0.5)
#       glmcoeff=as.matrix(coef(fit,s=lambda.best))
#       #sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
#       #glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
#       glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
#       idx1=match(glmleftsnp,rownames(glmcoeff))
#       glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snp by effect size
#       numvar=length(glmleftsnp)
#       if (numvar>0)
#       {
#         idx1=match(glmleftsnp,colnames(X))
#         Xsel=X[,idx1]
#         if (numvar>1) #check if number of covariate is greater than sample size
#         {
#           nummaxvar=min(nrow(X)-ncol(covariateall)-1,numvar)
#           numvar=nummaxvar
#           Xsel=Xsel[,1:nummaxvar]
#         }
#         if (numvar==1) #keep Xsel as in matrix form
#         {
#           Xsel=matrix(Xsel,ncol=1)
#           colnames(Xsel)=glmleftsnp
#         }
#         Xall1=data.matrix(cbind(Xsel,covariateall))
#         #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
#         fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
#         #summary(fit1)$r.squared
#         lmcoeff=summary(fit1)$coefficients
#         #align up coeff with Xsel
#         rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
#         if (sum(rownames(lmcoeff)=="age")>0) p_age=lmcoeff[which(rownames(lmcoeff)=="age"),4]
#         #p_disease=lmcoeff[which(rownames(lmcoeff)=="disease"),4]
#         if (sum(rownames(lmcoeff)=="gender")>0) p_gender=lmcoeff[which(rownames(lmcoeff)=="gender"),4]
#         lmleftsnp=rownames(lmcoeff)[rownames(lmcoeff) %in% colnames(Xsel)] 
#         numvar=length(lmleftsnp)
#         if (numvar>0)
#         {
#           glmflag=1
#           idx1=match(lmleftsnp,rownames(lmcoeff))
#           selectedsnps=paste0(rownames(lmcoeff)[idx1],collapse = "|")
#           selectedsnps_coeff=paste0(lmcoeff[idx1,1],collapse = "|")
#           idx1=match(lmleftsnp,colnames(Xsel))
#           Xsel=Xsel[,idx1]
#           if (numvar==1)
#           {
#             Xsel=matrix(Xsel,ncol=1)
#             colnames(Xsel)=lmleftsnp
#           }
#           fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
#           r2=cor(fitted,Y)^2
#           tmp=abs(cor(Y,Xsel))
#           hist(tmp,main=paste0(numvar," snps"),xlab="correlation between snps and expr")
#         }
#       }
#     }
#   }
#   return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
#               p_age=p_age,p_gender=p_gender))
# }

#load data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdata_ambiguous_TPM_for_prediction.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_allsnp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_allpos=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
allsnp=snp
allsnppos=snppos
allphenotype=phenotype
allphenotypepos=phenotypepos


#compute

Sys.time()
prefix="dist500K_GTEx_stomach_June11"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
if (!dir.exists(outfolder)) dir.create(outfolder)
sink(paste0(outfolder,"/TWAS.log"))

#KXD1
gene="UBAC1"
gene="S100A11"
gene="KXD1"
test=compute_cor_arow(i=which(rownames(phenotype)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.99)
test1=compute_cor_arow(i=which(rownames(phenotype)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.9)
length(intersect(unlist(strsplit(test$selectedsnps,"|",fixed = T)),unlist(strsplit(test1$selectedsnps,"|",fixed = T))))
#test2=compute_cor_arow(i=which(rownames(phenotype)==gene),ncv=10,distcutoff = 5e5,corcutoff = 0.5)

Sys.time()
# #parallel
#salloc -t 6-1  -n 99 mpirun -n 1 R --interactive
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdatafor_prediction.RData")
library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
#mpi.remote.exec(load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExstomachdatafor_prediction.RData"))
#mpi.remote.exec(sum(is.na(snp)))
# mpi.bcast.Robj2slave(snppos)
# mpi.bcast.Robj2slave(phenotype)
# mpi.bcast.Robj2slave(phenotypepos)
# mpi.bcast.Robj2slave(copynumber)
# mpi.bcast.Robj2slave(covariate)
#mpi.bcast.Robj2slave(gr_snp)
#mpi.bcast.Robj2slave(gr_pos)
mpi.bcast.cmd(library(GenomicRanges))
mpi.bcast.cmd(library(glmnet))
mpi.bcast.Robj2slave(compute_cor_arow)
mpi.bcast.Robj2slave(fitted_cv)
mpi.bcast.Robj2slave(removehighcorr)
mpi.bcast.Robj2slave(covariate) 

#run on each chr
mpi_compute_cor=function(ncv=10,distcutoff=5e5,corcutoff=0.9)
{
  res=NULL
  n=njobs
  for (chr in 1:22)
  {
    print(paste0(chr,"---"))
    idx=which(allphenotypepos$chr==chr)
    phenotype=allphenotype[idx,]
    phenotypepos=allphenotypepos[idx,]
    gr_pos=gr_allpos[idx]
    idx=which(allsnppos$chr==chr)
    snp=allsnp[idx,]
    gr_snp=gr_allsnp[idx]
    mpi.bcast.Robj2slave(gr_snp)
    mpi.bcast.Robj2slave(gr_pos)
    mpi.bcast.Robj2slave(phenotype)
    mpi.bcast.Robj2slave(snp)
    rows=1:nrow(phenotype)
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
      tmp=mpi.parSapply(X=seq,FUN=compute_cor_arow,ncv=ncv,distcutoff=distcutoff,corcutoff=corcutoff,job.num=njobs)
      if (length(unlist(tmp)) %% 7 !=0) print (i)
      res1=matrix(unlist(tmp),ncol=7,byrow = T)
      rownames(res1)=rownames(phenotypepos)[seq]
      res=rbind(res,res1)
    }
    save(chr,res,file=paste0(outfolder,"/TWAStmpresult.RData"))
    Sys.time()
  }
  colnames(res)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff","numlowcorsnp")
  res=as.data.frame(res)
  
  res[,5]=as.character(res[,5])
  res[,6]=as.character(res[,6])
  for (i in c(1:4,7:ncol(res))) res[,i]=as.numeric(as.character(res[,i]))
  return(res)
}



print("start to build models---")
Sys.time()
#rum mpi function and save the result
res_min=mpi_compute_cor(distcutoff = 5e5)
save(res_min,file=paste0(outfolder,"/preidiction_michigan_model.RData"))
Sys.time()

#extract snps from the prediction models


extract_snp2=function(dat=rbind(res_min))
{
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

tmp=extract_snp2()
print("extract BCA genotypte---")
Sys.time()
cmd=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/code/bca_extract_genotype2_V8_ambiguousSNP.sh ",prefix)
system(cmd,wait = T)
#after run extract_genotype.sh, read the extracted genotype
print("save BCA genotypte---")
library(data.table)
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

print("predict BCA gene expression---")
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

##NOTE: this step requires large nodes. It is not reliable to run on regular nodes (hang frequently)
##salloc -t 1-1 --mem-per-cpu 32G -n 21 --partition=largenode mpirun -n 1 R --interactive
if (!exists("phenotypepos")) load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8stomachdatafor_prediction.RData") #phenotypepos
if (!exists("res_min")) load(paste0(outfolder,"/preidiction_michigan_model.RData")) #models
if (!exists("bcagenotype")) load(paste0(outfolder,"/bca_extractgenotype.RData")) #extracted genotype
# library(Rmpi)
# njobs=mpi.universe.size() - 1
# print(njobs)
# mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
mpi.bcast.Robj2slave(res_min)
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
#predict_1se=mpi_predict_geneexp(modeltable=res_1se)
save(predict_min,file=paste0(outfolder,"/bca_predict_geneexp.RData"))
Sys.time()



#step2 association---
print("start association analysis---")
library(readxl)

sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
geneexpsamplenames=strsplit(colnames(predict_min)[3:ncol(predict_min)],"_") #use localid
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})

for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  sampletable[idx,i]=NA
}


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

if (!exists("predict_min")) load(paste0(outfolder,"/bca_predict_geneexp.RData"))

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
  all_predict_bcageneexp=predict_bcageneexp=predict_bcageneexp[,idx]
  idx=match(geneexpsamplenames,sampletable$localid)
  sampletable1=sampletable[idx,]
  for (i in 1:ncol(sampletable1)) sampletable1[,i][sampletable1[,i]==-9]=NA
  #BE
  mpi.bcast.Robj2slave(pvalue_arow)
  #mpi.bcast.Robj2slave(recodesite)
  #mpi.bcast.Robj2slave(predict_bcageneexp)
  mpi.bcast.Robj2slave(sampletable1)
  idx1=which(sampletable1$phenoBE_bca==2) #case
  length(idx1) #3288
  idx2=which(sampletable1$phenoBE_bca==1)
  length(idx2) #3195
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$ev1_bca[c(idx1,idx2)],pc2=sampletable1$ev2_bca[c(idx1,idx2)],pc3=sampletable1$ev3_bca[c(idx1,idx2)],pc4=sampletable1$ev4_bca[c(idx1,idx2)],
                        stringsAsFactors = F)
  #covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  n=1000 #1000 snps one time
  rows=1:nrow(all_predict_bcageneexp)
  nchunks=ceiling(nrow(all_predict_bcageneexp)/n)
  print(paste0("BE:number of total:",nchunks))
  BE_p=NULL
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
    predict_bcageneexp=all_predict_bcageneexp[seq,]
    mpi.bcast.Robj2slave(predict_bcageneexp)
    seq=1:nrow(predict_bcageneexp)
    tmp=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
    BE_p=c(BE_p,tmp)
  }
  
  
  #EA
  idx1=which(sampletable1$phenoEA_bca==2) #case
  length(idx1) #2514
  idx2=which(sampletable1$phenoEA_bca==1)
  length(idx2) #3198
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$ev1_bca[c(idx1,idx2)],pc2=sampletable1$ev2_bca[c(idx1,idx2)],pc3=sampletable1$ev3_bca[c(idx1,idx2)],pc4=sampletable1$ev4_bca[c(idx1,idx2)],
                        stringsAsFactors = F)
  #covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  print(paste0("EA:number of total:",nchunks))
  EA_p=NULL
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
    predict_bcageneexp=all_predict_bcageneexp[seq,]
    mpi.bcast.Robj2slave(predict_bcageneexp)
    seq=1:nrow(predict_bcageneexp)
    tmp=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
    EA_p=c(EA_p,tmp)
  }
  
  #EA vs # BE
  idx1=which(sampletable1$phenoEA_bca==2) # EAcase
  length(idx1) #2514
  idx2=which(sampletable1$phenoBE_bca==2) #BEcase
  length(idx2) #3288
  length(intersect(idx1,idx2))
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$ev1_bca[c(idx1,idx2)],pc2=sampletable1$ev2_bca[c(idx1,idx2)],pc3=sampletable1$ev3_bca[c(idx1,idx2)],pc4=sampletable1$ev4_bca[c(idx1,idx2)],
                        stringsAsFactors = F)
  #covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  print(paste0("BEA:number of total:",nchunks))
  BEA_p=NULL
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
    predict_bcageneexp=all_predict_bcageneexp[seq,]
    mpi.bcast.Robj2slave(predict_bcageneexp)
    seq=1:nrow(predict_bcageneexp)
    tmp=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
    BEA_p=c(BEA_p,tmp)
  }
  
  #EA and BE vs control
  idx1=which(sampletable1$phenoBE_bca==2 | sampletable1$phenoEA_bca==2) #BE or EAcase
  length(idx1) #5802
  idx2=which(sampletable1$phenoBE_bca==1 |sampletable1$phenoEA_bca==1) #Control
  length(idx2) #3198
  length(intersect(idx1,idx2))
  covariates=data.frame(age=sampletable1$age[c(idx1,idx2)],sex=sampletable1$sex[c(idx1,idx2)],site=sampletable1$site[c(idx1,idx2)],
                        bmi=sampletable1$bmi_recent_healthy[c(idx1,idx2)],cig_smk=sampletable1$cig_smk_ever[c(idx1,idx2)],cig_pack_yrs=sampletable1$cig_pack_yrs[c(idx1,idx2)],
                        pc1=sampletable1$ev1_bca[c(idx1,idx2)],pc2=sampletable1$ev2_bca[c(idx1,idx2)],pc3=sampletable1$ev3_bca[c(idx1,idx2)],pc4=sampletable1$ev4_bca[c(idx1,idx2)],
                        stringsAsFactors = F)
  #covariates=recodesite(covariates1=covariates,idx1=idx1,idx2=idx2)
  mpi.bcast.Robj2slave(idx1)
  mpi.bcast.Robj2slave(idx2)
  mpi.bcast.Robj2slave(covariates)
  print(paste0("BEEA:number of total:",nchunks))
  BEEA_p=NULL
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
    predict_bcageneexp=all_predict_bcageneexp[seq,]
    mpi.bcast.Robj2slave(predict_bcageneexp)
    seq=1:nrow(predict_bcageneexp)
    tmp=mpi.parSapply(X=seq,FUN=pvalue_arow,idx1=idx1,idx2=idx2,job.num=njobs)
    BEEA_p=c(BEEA_p,tmp)
  }
  
  res=data.frame(BE_p=BE_p,BE_fdr=p.adjust(BE_p),EA_p=EA_p,EA_fdr=p.adjust(EA_p),BEA_p=BEA_p,BEA_fdr=p.adjust(BEA_p),BEEA_p=BEEA_p,BEEA_fdr=p.adjust(BEEA_p))
  rownames(res)=rownames(all_predict_bcageneexp)
  return(res)
}
assoc_min=mpi_association()
#assoc_1se=mpi_association(predict_bcageneexp = predict_1se)
save(assoc_min,file=paste0(outfolder,"/bca_assoc.RData"))
load(paste0(outfolder,"/bca_assoc.RData"))
quantile(assoc_min[,1],na.rm=T)
quantile(assoc_min[,3],na.rm=T)
quantile(assoc_min[,5],na.rm=T)
quantile(assoc_min[,7],na.rm=T)
Sys.time()
print("Done")
sink()
## quit program
mpi.close.Rslaves()
mpi.quit()
quit()
