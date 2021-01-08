
############################################################################
### This code modify Kevin's code and check eQTL prediction model and    ###
### perform TWAS analysis in BEACON data                                 ###
############################################################################

rm(list=ls())


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


compute_cor_arow=function(i,opt="1se",ncv=10,distcutoff=1e6)
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
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  if (length(idx)>1)
  {
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
    X <- removehighcorr(X1,0.8)
    hist(ucor[colnames(X1)%in% colnames(X),1])
    hist(ucor[colnames(X1)%in% colnames(X),2])
    dim(X)
    Xall=data.matrix(cbind(X,covariate))
    covariateall=covariate
    
    penalty=rep(1,ncol(Xall))
    penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
    #set.seed(i+10000)
    cvfit=tryCatch(
      {
        cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
      },
      error=function(e)
      {
        return(F)
      }
    )
    plot(cvfit)
    cverr <- matrix(0,length(cvfit$lambda),100) 
    for (l in 1:100) {
      fit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
      alllambda=intersect(cvfit$lambda,fit$lambda)
      idx1=match(alllambda,cvfit$lambda)
      idx2=match(alllambda,fit$lambda)
      cverr[idx1,l] <- fit$cvm[idx2]
    }  
    merr <- apply(cverr,1,mean)
    plot(log(fit$lambda),merr)
    
    lambda.best <- fit$lambda[which.min(merr)]
    glmcoeff=as.matrix(coef(fit,s=lambda.best))
    sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
    glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
    
    
    if (is.list(cvfit))
    {
      fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=1)
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
            Xall1=data.matrix(cbind(Xsel,covariateall))
            #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
            fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
            #summary(fit1)$r.squared
            lmcoeff=summary(fit1)$coefficients
            #align up coeff with Xsel
            rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
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
              fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
              r2=cor(fitted,Y)^2
            }
            
          }
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              p_age=p_age,p_gender=p_gender))
}


removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}

library(SKAT)


prefix="dist500K_GTEx_April18"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(file=paste0(outfolder,"/bca_extractgenotype.RData"))
tmp=colnames(bcagenotype)
if (grepl("_",tmp[1]))
{
  tmp=strsplit(tmp,"_")
  tmp1=sapply(1:length(tmp),function(x){
    tmp1=tmp[[x]]
    paste0(tmp1[2:length(tmp1)],collapse = "_")
  })
  colnames(bcagenotype)=tmp1 #use localid
}



library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)

sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
#eigenstratmatrix=readeigenstrat()
allsamples=intersect(colnames(bcagenotype),sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid
idx=match(allsamples,colnames(bcagenotype))
bcagenotype=bcagenotype[,idx]
# idx=match(allsamples,colnames(eigenstratmatrix))
# eigenstratmatrix=eigenstratmatrix[,idx]
# sampletable=cbind(sampletable,t(eigenstratmatrix[1:4,]))
Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","pc1","pc2","pc3","pc4")]
Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)


beaconfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GWAS/imputed_GTExmodel/"

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])

# allsnps=as.data.frame(fread("/fh/fast/stanford_j/Xiaoyu/Tools/annotation/All_SNPs.txt"))
# allsnps=allsnps[allsnps$`#chrom` %in% paste0("chr",c(1:22,"X","Y")),]
# allsnps$`#chrom`=gsub("chr","",allsnps$`#chrom`)

library("biomaRt")
mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")

summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_autosomes.txt")

beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic")

library(data.table)
#summary stat from validation
summarydat=as.data.frame(fread(summaryfile,header=T))
summarydat1=as.data.frame(fread(summaryfile1,header=T))

tmp=intersect("SNP",colnames(summarydat))
if (length(tmp)==0)
{
  colnames(summarydat)[which(colnames(summarydat)=="rsid")]="SNP"
  colnames(summarydat)[which(colnames(summarydat)=="pvalue")]="P"
  colnames(summarydat)[which(colnames(summarydat)=="beta")]="BETA"
  colnames(summarydat)[which(colnames(summarydat)=="effect-allele")]="effect_allele"
  colnames(summarydat)[which(colnames(summarydat)=="non-effect-allele")]="non_effect_allele"
}

tmp=intersect("SNP",colnames(summarydat1))
if (length(tmp)==0)
{
  colnames(summarydat1)[which(colnames(summarydat1)=="rsid")]="SNP"
  colnames(summarydat1)[which(colnames(summarydat1)=="pvalue")]="P"
  colnames(summarydat1)[which(colnames(summarydat1)=="beta")]="BETA"
  colnames(summarydat1)[which(colnames(summarydat1)=="effect-allele")]="effect_allele"
  colnames(summarydat1)[which(colnames(summarydat1)=="non-effect-allele")]="non_effect_allele"
}

#summary from beacon
beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))




#load data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp


#compute
#test=compute_cor_arow(i=31038,opt="min",ncv=10,distcutoff = 5e5)
which(row.names(phenotype)=="DDX49")
#[1] 31038
which(row.names(phenotype)=="CERS1")
#[1] 31036
which(row.names(phenotype)=="KXD1")
#i <- 31038
which(row.names(phenotype)=="HSP90AA1")
#i <- 25105
which(row.names(phenotype)=="UBAC1")
#i <- 17803
which(row.names(phenotype)=="FBP2")
#i <- 17243


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
  X <- removehighcorr(X1,0.8)
  hist(ucor[colnames(X1)%in% colnames(X),1])
  hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  cvfit$lambda.1se
  
  cverr <- matrix(NA,length(cvfit$lambda),100)
  cvse  <- matrix(NA,length(cvfit$lambda),100)
  
  rownames(cverr)=cvfit$lambda
  for (l in 1:100) {
    set.seed(l+100)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5)
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
  
  #plot(cvfit$lambda,merr)
  
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  glmcoeff=as.matrix(coef(cvfit,s=lambda.min))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]


  mean(selsnps %in% row.names(bcagenotype))
  
  #if some imputed snps need to to flipped
  
  correctedsnps=NULL
  if (length(intersect(selsnps,rownames(bcagenotype)))<length(selsnps))
  {
    missingsnps=selsnps[!selsnps %in% rownames(bcagenotype)]
    for (j in 1:length(missingsnps))
    {
      tmp=unlist(strsplit(missingsnps[j],"_"))
      tmp1=paste0(tmp[c(1,3,2)],collapse = "_")  #change the order of allele
      idx=which(rownames(bcagenotype)==tmp1)
      if (length(idx)>0)
      {
        correctedsnps=c(correctedsnps,tmp1)
        idx1=which(selsnps==missingsnps[j])
        selsnps[idx1]=tmp1 #change the snp name to make it consistent with bca
      }
      
    }
  } 
  
  ### still missing snps from bca genotype data ###
  mean(selsnps %in% row.names(bcagenotype))
  
  selectedsnps <- selsnps
  
  
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
  ## there a few SNPs not found in bcagenotype
  sum(is.na(idx))
  
  selsnps <- selsnps[!is.na(idx)]
  selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
  #if (length(idxtocorrect)>0)
  #{
  #  Z[idxtocorrect,]=2-Z[idxtocorrect,]
  #}
  rownames(Z)=colnames(bcagenotype)
  idx=match(allsamples,rownames(Z))
  Z=Z[idx,,drop=F]
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoBE_bc==2) #case
  idx2=which(sampletable$phenoBE_bc==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
 
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  out$p.value
   
    
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
 
      
   
      
      selsnps <- colnames(Z1)
      tmp=unlist(strsplit(selsnps[1],":"))
      chr=tmp[1]
      pos=rsid=allele1=allele2=p_uni=p_mul=cor_uni=rep(NA,length(selsnps))
        #selsnps="19:19259262_C_A"
      for (j in 1:length(selsnps))
        {
          cat(j,"\n")
          tmp=unlist(strsplit(selsnps[j],":"))
          tmp=unlist(strsplit(tmp[2],"_"))
          allalleles=tmp[2:3]
          allele1[j]=tmp[2] #minor allele
          allele2[j]=tmp[3]
          pos[j]=as.numeric(tmp[1])
          #find rsid based on position and allles
          #tmp=listAttributes(snpmart)
          tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
          if (nrow(tmp)>0)
          {
            for (k in 1:nrow(tmp))
            {
              myalleles=unlist(strsplit(tmp$allele[k],"/",fixed=T))
              if (all(allalleles %in% myalleles))
              {
                rsid[j]=tmp$refsnp_id[k]
                break
              }
            }
          }
      }
      
        ### validation in Bonn ###
      
        tmp=intersect(rsid[!is.na(rsid)],summarydat$SNP)
        valsnps1 <- selsnps[which(rsid %in% tmp)]
        rsid1 <- rsid[rsid %in% tmp]
        val1 <- summarydat[summarydat$SNP %in% tmp,]
        val1 <- val1[match(rsid1,val1$SNP),]
        uu <- val1$BETA
        vv <- as.numeric(val1$SE)
        ZZ <- Z1[,colnames(Z1) %in% valsnps1]
        rr <- cor(ZZ,use="pairwise.complete.obs")
        VV <- diag(vv) %*% rr %*% diag(vv)
        lamb <- eigen(VV)$values
        Q <- drop(t(uu) %*% uu)
        mean(val1$P<0.05)
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
        
        
        
        
        ### validation in Oxford ###
        
          
        tmp=intersect(rsid[!is.na(rsid)],summarydat1$SNP)
        valsnps2 <- selsnps[which(rsid %in% tmp)]
        rsid1 <- rsid[rsid %in% tmp]
        val2 <- summarydat1[summarydat1$SNP %in% tmp,]
        val2 <- val2[match(rsid1,val2$SNP),]
        uu <- val2$BETA
        vv <- as.numeric(val2$se)
        ZZ <- Z1[,colnames(Z1) %in% valsnps2]
        rr <- cor(ZZ,use="pairwise.complete.obs")
        VV <- diag(vv) %*% rr %*% diag(vv)
        lamb <- eigen(VV)$values
        Q <- drop(t(uu) %*% uu)
        mean(val1$P<0.05)
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
        #sum(lamb*(1-pchisq(Q,df=1)))
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
        
        
        #### combine the two validation datasets ####
        
        val1 <- val1[val1$SNP %in% val2$SNP,]
        valsnps1 <- valsnps1[valsnps1%in%valsnps2]
        
        val2 <- val2[val2$SNP %in% val1$SNP,]
        valsnps2 <- valsnps2[valsnps2%in%valsnps1]
        
        uu1 <- val1$BETA
        vv1 <- as.numeric(val1$SE)
        uu2 <- val2$BETA
        vv2 <- as.numeric(val2$se)
        
        uu <- uu1+uu2
        ZZ <- Z1[,colnames(Z1) %in% valsnps1]
        rr <- cor(ZZ,use="pairwise.complete.obs")
        VV1 <- diag(vv1) %*% rr %*% diag(vv1)
        VV2 <- diag(vv2) %*% rr %*% diag(vv2)
        VV <- VV1 + VV2
        lamb <- eigen(VV)$values
        Q <- drop(t(uu) %*% uu)
        mean(val1$P<0.05)
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
        #sum(lamb*(1-pchisq(Q,df=1)))
        pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")

  
        
        
        
        
        
        
        
        
        
        
        
        
  
which(row.names(phenotype)=="CERS1")
#[1] 31036
which(row.names(phenotype)=="PLEKHF2")
#[1] 15931
which(row.names(phenotype)=="GATA4")

which(row.names(phenotype)=="CTSB")

which(row.names(phenotype)=="SGK223")

which(row.names(phenotype)=="KXD1")

which(row.names(phenotype)=="PPP1R3B")

which(row.names(phenotype)=="COMP")

which(row.names(phenotype)=="OR11A1")

which(row.names(phenotype)=="SPAG11A")

which(row.names(phenotype)=="MEF2B")
