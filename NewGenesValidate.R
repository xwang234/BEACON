
###############################################################################
##### Discover and Validate several NEW genes in summary stat data      #######
##### May 2021                                                          #######
###############################################################################



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



#load BCA covariate tables
library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable=as.data.frame(sampletable)
for (i in 1:ncol(sampletable))
{
  idx=which(sampletable[,i]==-9)
  if (length(idx)>0)
    sampletable[idx,i]=NA
}

#BCA covariate table, principle components
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
#Covariate is used for SKAT
opt="PC6"
#try different number of PCs
if (opt=="PC4")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","sex")] #pc4
# Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","sex")]
if (opt=="PC6")
  Covariate=covariatetable[,colnames(covariatetable) %in% c("pc1","pc2","pc3","pc4","pc5","pc6","sex")]




TWAS_SKAT_gene <- function(genename,cvmin=1,phenotype,phenotypepos,snp,snppos,bcagenotype,covariate,verbose=0) {
  
  #idx2=which(codinggenetable$Symbol==genename)
  #if (verbose>0) print(paste0("work on ",genename,". Its closest gwas snp: ",codinggenetable$gwas_snp[idx2], ", distance:",round(codinggenetable$dist_to_gwas_snp[idx2]/1000000,2),"MB"))
  library(SKAT)
  library(glmnet)
  library(GenomicRanges)
  snppos$chr[snppos$chr==23]="X"
  phenotypepos$chr[phenotypepos$chr==23]="X"
  gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
  gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
  
  result <- data.frame(matrix(0,8,4))  
  if (verbose>0) print("get the gene model---")
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
  dim(X1)
  #ucor <- matrix(0,ncol(X1),2)
  #for (l in 1:ncol(X1)){
  #  ucov<- data.matrix(cbind(X1[,l],covariate))
  #  ufit <- lm(Y~ucov)
  #  ucor[l,1] <- summary(ufit)$coef[2,4]
  #  ucor[l,2] <- cor(Y,X1[,l])
  #}  
  #hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  #pcor <- ucor[,1]
  #X1 <- X1[,order(pcor,decreasing=T)]
  X <- removehighcorr(X1,0.9)
  #hist(ucor[colnames(X1)%in% colnames(X),1])
  #hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  max(cvfit$lambda[cvfit$cvm < min(cvfit$cvm) + cvfit$cvsd[which.min(cvfit$cvm)]])
  cvfit$lambda.1se
  cvfit$lambda.min
  
  
  cverr <- matrix(NA,length(cvfit$lambda),100)
  cvse  <- matrix(NA,length(cvfit$lambda),100)
  rownames(cverr)=cvfit$lambda
  for (l in 1:ncol(cverr)) {
    set.seed(l+100)
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    alllambda=intersect(cvfit$lambda,fit$lambda)
    idx1=match(alllambda,cvfit$lambda)
    idx2=match(alllambda,fit$lambda)
    cverr[idx1,l] <- fit$cvm[idx2]
    cvse[idx1,l] <- fit$cvsd[idx2]
  }
  
  merr <- apply(cverr,1,mean,na.rm=T)
  mcvse <- sqrt(apply(cvse^2,1,mean))
  
  lambda.1se <- cvfit$lambda[min(which(merr< (min(merr) + mcvse[which.min(merr)])))]
  lambda.min <- cvfit$lambda[which.min(merr)]
  
  #plot(log(cvfit$lambda),merr)
  #abline(v=log(lambda.1se),lty=2)
  #abline(v=log(lambda.min),lty=3)
  #plot(cvfit$lambda,merr)
  #glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  lambda.min <- cvfit$lambda.min
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  
  ## check the marginal eQTL association ##
  outp <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat <- data.matrix(cbind(X[,colnames(X)==selsnps[i]],covariate))
    fit <- lm(Y~Xmat)
    outp[i] <- summary(fit)$coef[2,4]
  }  
  outp
  
  mean(selsnps %in% row.names(bcagenotype))
  ## here we should have a function to get the genotypes for SNPs being selected in selsnps ##
  bcagenotype=extractBCAgenotype(selsnps = selsnps) #takes 2 minutes for 38 snps
  
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
  
  if (verbose>0) print("compute skat result----")
  idx1=match(colnames(bcagenotype),rownames(covariatetable))
  covariatetable=covariatetable[idx1,]
  idx1=match(colnames(bcagenotype),rownames(Covariate))
  Covariate=Covariate[idx1,]
  if (any(rownames(covariatetable)!=colnames(bcagenotype))) #covariatetable keeps BCA covariates
    warning("bcagenotype doesn't match covariatetable")
  
  idx=match(selectedsnps,rownames(bcagenotype))
  ## there a few SNPs not found in bcagenotype
  if (sum(is.na(idx))>0)
    print(paste0(sum(is.na(idx))," out of ", length(idx)," selected snps not been found in genotype data"))
  
  selsnps <- selsnps[!is.na(idx)]
  #print(length(selsnps))
  selcoeff <- selcoeff[!is.na(idx)]
  idx <- idx[!is.na(idx)]
  Z=t(bcagenotype[idx,,drop=F])
  if (length(correctedsnps)>0) #flip snps
  {
    idxtocorrect=match(correctedsnps,colnames(Z))
    Z[,idxtocorrect]=2-Z[,idxtocorrect]
  }
  #colnames(Z) <- row.names(bcagenotype)[idx]
  mean(colnames(Z)==selsnps)
  
  #rownames(Z)=colnames(bcagenotype)
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoBE_bca==2) #case
  idx2=which(covariatetable$phenoBE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[6,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[6,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[6,4] <- out$p.value
  rownames(result)[6]="BE_SKAT"
  result[6,]
  # for (i in 1:ncol(Z1)) {
  #   Covariate2 <- as.data.frame(cbind(Z1[,i],Covariate1))
  #   outfit <- glm(y~.,data=Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEA_bca==2) #case
  idx2=which(covariatetable$phenoEA_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[7,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[7,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[7,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[7,4] <- out$p.value
  rownames(result)[7]="EA_SKAT"
  result[7,]
  #for (i in 1:ncol(Z1)) {
  #   Covariate2 <- cbind(Z1[,i],Covariate1)
  #   outfit <- glm(y~Covariate2,family="binomial")
  #   print(summary(outfit)$coef[2,4])
  # }  
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  outp1
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(covariatetable$phenoEABE_bca==2) #case
  idx2=which(covariatetable$phenoEABE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.data.frame(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[8,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(as.formula("y ~."), data=Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T) 
  result[8,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),method="SKATO",is_dosage=T)
  result[8,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),is_dosage=T)
  result[8,4] <- out$p.value
  rownames(result)[8]="BEEA_SKAT"
  result[8,]
  
  outp1 <- rep(0,length(selsnps))
  for (i in 1:length(selsnps)) {
    Xmat<- data.matrix(cbind(Z1[,colnames(Z1)==selsnps[i]],Covariate1))
    fit <- glm(y~Xmat,family=binomial) 
    outp1[i] <- summary(fit)$coef[2,4]
  }  
  outp1
  
  
  if (verbose>0) print("start validation---")
  selsnps <- colnames(Z1)
  tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
  idx=which(is.na(tmp$snp))
  if (length(idx)>0) #use biomart, slower method
  {
    pos=allele1=allele2=rep(NA,length(selsnps))
    for (j in idx)
    {
      #cat(j,"\n")
      tmp1=unlist(strsplit(tmp$selsnps[j],":"))
      chr=tmp1[1]
      tmp1=unlist(strsplit(tmp1[2],"_"))
      allalleles=tmp1[2:3]
      allele1[j]=tmp1[2] #minor allele
      allele2[j]=tmp1[3]
      pos[j]=as.numeric(tmp1[1])
      #find rsid based on position and allles
      #tmp=listAttributes(snpmart)
      tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
      if (nrow(tmp2)>0)
      {
        for (k in 1:nrow(tmp2))
        {
          myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
          if (all(allalleles %in% myalleles))
          {
            tmp$snp[j]=tmp2$refsnp_id[k]
            break
          }
        }
      }
    }
  }
  rsid=tmp$snp
  if (sum(is.na(rsid))) print(paste0(sum(is.na(rsid)," snps can't find snp rsid")))
  library(survey)
  ### validation in Bonn EA ###
  
  tmp=intersect(rsid[!is.na(rsid)],EA_Bonnsummarydat$SNP)
  
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- EA_Bonnsummarydat[EA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[1,1] <- length(tmp)
  result[1,2] <- mean(val1$P<0.05)
  result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[1]="EA_Bonn"  
  result[1,]
  
  ### validation in Bonn EA/BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BEEA_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BEEA_Bonnsummarydat[BEEA_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[2,1] <- length(tmp)
  result[2,2] <- mean(val1$P<0.05)
  result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[2]="BEEA_Bonn"
  result[2,]
  ### validation in Bonn BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Bonnsummarydat$SNP)
  valsnps1 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val1 <- BE_Bonnsummarydat[BE_Bonnsummarydat$SNP %in% tmp,]
  val1 <- val1[match(rsid1,val1$SNP),]
  uu <- val1$BETA
  vv <- as.numeric(val1$SE)
  ZZ <- Z1[,match(valsnps1,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  result[3,1] <- length(tmp)
  result[3,2] <- mean(val1$P<0.05)
  result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[3]="BE_Bonn"
  result[3,]
  ### validation in Oxford BE ###
  
  tmp=intersect(rsid[!is.na(rsid)],BE_Oxfordsummarydat$SNP)
  valsnps2 <- selsnps[which(rsid %in% tmp)]
  rsid1 <- rsid[which(rsid %in% tmp)]
  val2 <- BE_Oxfordsummarydat[BE_Oxfordsummarydat$SNP %in% tmp,]
  val2 <- val2[match(rsid1,val2$SNP),]
  uu <- val2$BETA
  vv <- as.numeric(val2$se)
  ZZ <- Z1[,match(valsnps2,colnames(Z1))]
  rr <- cor(ZZ,use="pairwise.complete.obs")
  VV <- diag(vv) %*% rr %*% diag(vv)
  lamb <- eigen(VV)$values
  Q <- drop(t(uu) %*% uu)
  
  result[4,1] <- length(tmp)
  result[4,2] <- mean(val2$P<0.05)
  result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
  #sum(lamb*(1-pchisq(Q,df=1)))
  result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
  rownames(result)[4]="BE_Oxford"
  result[4,]
  #### combine the two validation datasets ####
  
  val1 <- val1[val1$SNP %in% val2$SNP,]
  valsnps1 <- valsnps1[valsnps1%in%valsnps2]
  
  val2 <- val2[val2$SNP %in% val1$SNP,]
  valsnps2 <- valsnps2[valsnps2%in%valsnps1]
  if (nrow(val1)>1)  {
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
    result[5,1] <- nrow(val1)
    result[5,2] <- mean(val1$P<0.05 & val2$P <0.05)
    result[5,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
    result[5,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
    rownames(result)[5]="BE_Combined"
  } 
  
  return(result)  
}

extractBCAgenotype=function(selsnps=NULL)
{
  allgenotypedat=NULL
  dataset="merge_beacon_cambridge_hrc_maf001_snp" #for HRC maf001
  chr=unlist(strsplit(selsnps[1],":"))[1]
  genotypefile=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/",dataset,"/","chr",chr,"_filter_hg19tohg38_flip.traw")
  bimfile=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/",dataset,"/","chr",chr,"_filter_hg19tohg38_flip.bim")
  bim=as.data.frame(data.table::fread(bimfile))
  columnname=as.data.frame(data.table::fread(genotypefile,nrows=1))
  #sample name
  samplename=colnames(columnname)[7:ncol(columnname)]
  samplename=unlist(strsplit(samplename,"_"))
  samplename=samplename[seq(2,length(samplename),2)] 
  samplename=gsub("SEP","",samplename)
  for (i in 1:length(selsnps))
  {
    mypos=unlist(strsplit(selsnps[i],"_"))[1]
    mypos=unlist(strsplit(mypos,":"))[2]
    #where is the snp
    nline=which(bim$V4==mypos)
    
    if (length(nline)>0) #read data
    {
      genotypedat=as.data.frame(data.table::fread(genotypefile,nrows=length(nline),skip=nline))
      rownames(genotypedat)=paste0(genotypedat[,1],":",genotypedat[,4],"_",genotypedat[,5],"_",genotypedat[,6])
      genotypedat=genotypedat[,7:ncol(genotypedat)]
      colnames(genotypedat)=samplename
      allgenotypedat=rbind(allgenotypedat,genotypedat)
    }
  }
  
  return(allgenotypedat)
}



organidx <- 2
opt <- "HRC"

organs=c("mucosa","junction","stomach","muscularis","adipose","blood")
outfolders=c("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_mucosa_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_junction_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_stomach_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_muscularis_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_adipose_HRC_MAF005",
             "/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_blood_HRC_MAF005")

outfolder=outfolders[organidx]

#load(paste0(outfolder,"/preidiction_michigan_model.RData")) #the saved model res_min
#no need to load bcagenotype data
#load(paste0(outfolder,"/bca_extractgenotype.RData")) #the saved genotype data, bca_genotype
#load(paste0(outfolder,"/bca_predict_geneexp.RData")) #predicted geneexp, predict_min
#load(paste0(outfolder,"/skat_res.RData")) #saved skat res, skat_min2_pc6,skat_min2
#load GTEx gene expression data, snp,phenotype
if (opt=="HRC")
{
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
}


genename <- "IL2RB"
genename <- "COX7A2"
genename <- "FILIP1"
genename <- "HSP90AA1"
genename <- "KRTAP5-8"
genename <- "FOXF1"
genename <- "LDAH"
genename <- "ISYNA1"
genename <- "UBAC1"

cvmin=1
verbose=T
eqtl <- vector("list", 6)
for (organidx in 1:6){
  
  outfolder=outfolders[organidx]
  load(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExV8",organs[organidx],"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData"))
  
  
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
  #ucor <- matrix(0,ncol(X1),2)
  #for (l in 1:ncol(X1)){
  #  ucov<- data.matrix(cbind(X1[,l],covariate))
  #  ufit <- lm(Y~ucov)
  #  ucor[l,1] <- summary(ufit)$coef[2,4]
  #  ucor[l,2] <- cor(Y,X1[,l])
  #}  
  #hist(ucor[,1])
  #X<- X[,ucor[,1]<0.2]
  #pcor<- ucor[ucor[,1]<0.2,1]
  #pcor <- ucor[,1]
  #X1 <- X1[,order(pcor,decreasing=T)]
  X <- removehighcorr(X1,0.9)
  #hist(ucor[colnames(X1)%in% colnames(X),1])
  #hist(ucor[colnames(X1)%in% colnames(X),2])
  dim(X)
  Xall=data.matrix(cbind(X,covariate))
  covariateall=covariate
  
  penalty=rep(1,ncol(Xall))
  #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
  set.seed(i+10000)
  cvfit=tryCatch(
    {
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=T)
    },
    error=function(e)
    {
      return(F)
    }
  )
  plot(cvfit)
  
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=cvfit$lambda.1se))
  
  
  sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
  selcoeff <- glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
  selsnps <- rownames(glmcoeff)[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X)]
  eqtl[[organidx]] <- selsnps
  print(organs[organidx])
  print(length(Y))
  print(eqtl[[organidx]])
}  
  

