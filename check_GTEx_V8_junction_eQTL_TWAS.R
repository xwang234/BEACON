
############################################################################
### This code modify Kevin's code and check eQTL prediction model and    ###
### perform TWAS analysis in BEACON data                                 ###
############################################################################

rm(list=ls())


removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) >= corcutoff))]
  return(datnew)
}

library(SKAT)

### loading genotype data ###

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


### loading phenotype data ###


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


### loading summary stat data for validation ###

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

summaryfile2=paste0(summaryfolder,"EA_Bonn_autosomes.txt")
summaryfile3=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt")

#beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic")

library(data.table)
#summary stat from validation
summarydat=as.data.frame(fread(summaryfile,header=T))
summarydat1=as.data.frame(fread(summaryfile1,header=T))

summarydat2=as.data.frame(fread(summaryfile2,header=T))
summarydat3=as.data.frame(fread(summaryfile3,header=T))

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



tmp=intersect("SNP",colnames(summarydat2))
if (length(tmp)==0)
{
  colnames(summarydat2)[which(colnames(summarydat2)=="rsid")]="SNP"
  colnames(summarydat2)[which(colnames(summarydat2)=="pvalue")]="P"
  colnames(summarydat2)[which(colnames(summarydat2)=="beta")]="BETA"
  colnames(summarydat2)[which(colnames(summarydat2)=="effect-allele")]="effect_allele"
  colnames(summarydat2)[which(colnames(summarydat2)=="non-effect-allele")]="non_effect_allele"
}


tmp=intersect("SNP",colnames(summarydat3))
if (length(tmp)==0)
{
  colnames(summarydat3)[which(colnames(summarydat3)=="rsid")]="SNP"
  colnames(summarydat3)[which(colnames(summarydat3)=="pvalue")]="P"
  colnames(summarydat3)[which(colnames(summarydat3)=="beta")]="BETA"
  colnames(summarydat3)[which(colnames(summarydat3)=="effect-allele")]="effect_allele"
  colnames(summarydat3)[which(colnames(summarydat3)=="non-effect-allele")]="non_effect_allele"
}


#summary from beacon
#beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))




#summary from beacon
#beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))

##load GTEx gene expression data

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


genename = "HSP90AA1"

TWAS_SKAT <- function(genename,cvmin=1) {
result <- matrix(0,8,4)  
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
      cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=F)
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
    fit <- cv.glmnet(data.matrix(Xall),Y,nfolds=10, penalty.factor=penalty,alpha=0.5,standardize=F)
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
  if (cvmin==1) glmcoeff=as.matrix(coef(cvfit,s=lambda.min)) else glmcoeff=as.matrix(coef(cvfit,s=lambda.1se))
  
  
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
  idx1=which(sampletable$phenoBE_bca==2) #case
  idx2=which(sampletable$phenoBE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  
  result[6,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[6,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[6,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[6,4] <- out$p.value
   
    
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
  
  ##############################
  ### compare EA vs control  ###
  ##############################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoEA_bca==2) #case
  idx2=which(sampletable$phenoEA_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[7,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[7,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[7,4] <- out$p.value
  
  
  for (i in 1:ncol(Z1)) {
    Covariate2 <- cbind(Z1[,i],Covariate1)
    outfit <- glm(y~Covariate2,family="binomial")
    print(summary(outfit)$coef[2,4])
  }  
  
  
  #################################
  ### compare EA/BE vs control  ###
  #################################
  
  #print(rankMatrix(Z)[[1]])
  idx1=which(sampletable$phenoEABE_bca==2) #case
  idx2=which(sampletable$phenoEABE_bca==1)
  
  y=c(rep(1,length(idx1)),rep(0,length(idx2)))
  
  Z1=Z[c(idx1,idx2),,drop=F]
  Covariate1=as.matrix(Covariate[c(idx1,idx2),])
  tmp=rowSums(Z1)
  
  result[8,1] <- ncol(Z1)
  obj.s<-SKAT_Null_Model(y ~ Covariate1,out_type="D")
  out<-SKAT(Z1, obj.s,weights.beta=c(1,1),is_dosage=T,method="SKATO") 
  result[8,2] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights.beta=c(1,1),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[8,3] <- out$p.value
  out=SKATBinary(Z1, obj.s, weights=abs(selcoeff),kernel = "linear.weighted",r.corr=0,is_dosage=T) #1.207545e-05
  result[8,4] <- out$p.value
  
  
  
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
      
      library(survey)
      ### validation in Bonn EA ###
      
      tmp=intersect(rsid[!is.na(rsid)],summarydat2$SNP)
      valsnps1 <- selsnps[which(rsid %in% tmp)]
      rsid1 <- rsid[rsid %in% tmp]
      val1 <- summarydat2[summarydat2$SNP %in% tmp,]
      val1 <- val1[match(rsid1,val1$SNP),]
      uu <- val1$BETA
      vv <- as.numeric(val1$SE)
      ZZ <- Z1[,colnames(Z1) %in% valsnps1]
      rr <- cor(ZZ,use="pairwise.complete.obs")
      VV <- diag(vv) %*% rr %*% diag(vv)
      lamb <- eigen(VV)$values
      Q <- drop(t(uu) %*% uu)
      result[1,1] <- length(tmp)
      result[1,2] <- mean(val1$P<0.05)
      result[1,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
      result[1,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
      
      
      ### validation in Bonn EA/BE ###
      
      tmp=intersect(rsid[!is.na(rsid)],summarydat3$SNP)
      valsnps1 <- selsnps[which(rsid %in% tmp)]
      rsid1 <- rsid[rsid %in% tmp]
      val1 <- summarydat3[summarydat3$SNP %in% tmp,]
      val1 <- val1[match(rsid1,val1$SNP),]
      uu <- val1$BETA
      vv <- as.numeric(val1$SE)
      ZZ <- Z1[,colnames(Z1) %in% valsnps1]
      rr <- cor(ZZ,use="pairwise.complete.obs")
      VV <- diag(vv) %*% rr %*% diag(vv)
      lamb <- eigen(VV)$values
      Q <- drop(t(uu) %*% uu)
      result[2,1] <- length(tmp)
      result[2,2] <- mean(val1$P<0.05)
      result[2,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
      result[2,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
      
      ### validation in Bonn BE ###
      
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
        result[3,1] <- length(tmp)
        result[3,2] <- mean(val1$P<0.05)
        result[3,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
        result[3,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
        
        ### validation in Oxford BE ###
        
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
        
        result[4,1] <- length(tmp)
        result[4,2] <- mean(val2$P<0.05)
        result[4,3] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "satterthwaite")
        #sum(lamb*(1-pchisq(Q,df=1)))
        result[4,4] <- pchisqsum(Q, 1, lamb, lower.tail = FALSE,method = "saddlepoint")
        
        
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
        } 

return(result)  
        
}
        
        
        
        
        
out <- TWAS_SKAT(genename="HSP90AA1",cvmin=1)

#[,1]         [,2]         [,3]        [,4]
#[1,]   86 8.139535e-02 8.272366e-02 0.081294035
#[2,]   86 4.651163e-02 9.088799e-02 0.088920299
#[3,]   86 1.744186e-01 4.803112e-02 0.049625717
#[4,]   95 1.052632e-01 1.086846e-01 0.100721569
#[5,]   86 2.325581e-02 5.676580e-03 0.008726887
#[6,]   95 4.249421e-06 1.354467e-06 0.006065673
#[7,]   95 5.135774e-03 2.911774e-03 0.015180094
#[8,]   95 1.607161e-05 5.903575e-06 0.003862113


        
out <- TWAS_SKAT(genename="HSP90AA1",cvmin=0)


  

out <- TWAS_SKAT(genename="KXD1",cvmin=1)

#[,1]         [,2]         [,3]        [,4]
#[1,]    8 3.750000e-01 2.719933e-01 0.238072831
#[2,]    8 1.250000e-01 4.055442e-01 0.371020115
#[3,]    8 0.000000e+00 5.210786e-01 0.498698477
#[4,]    8 1.250000e-01 1.725170e-01 0.155636136
#[5,]    8 0.000000e+00 3.398899e-01 0.306474747
#[6,]    8 4.616631e-06 1.760928e-06 0.086458296
#[7,]    8 4.218945e-08 7.018217e-08 0.003723487
#[8,]    8 1.125119e-08 4.424843e-09 0.012985007



out <- TWAS_SKAT(genename="UBAC1",cvmin=1)

#[1,]   66 4.545455e-02 5.365458e-01 5.175517e-01
#[2,]   66 3.030303e-02 2.185358e-01 2.045960e-01
#[3,]   66 7.575758e-02 4.135389e-02 4.475186e-02
#[4,]   97 0.000000e+00 9.591420e-01 9.720015e-01
#[5,]   65 0.000000e+00 3.617976e-01 3.358825e-01
#[6,]  105 3.150919e-05 1.374749e-05 2.235438e-05
#[7,]  105 1.925271e-01 3.401016e-01 2.844838e-03
#[8,]  105 6.879603e-04 5.351333e-04 2.679653e-05


out <- TWAS_SKAT(genename="UBAC1",cvmin=0)

out <- TWAS_SKAT(genename="ISYNA1",cvmin=1)



out <- TWAS_SKAT(genename="ISYNA1",cvmin=1)

#[,1]         [,2]         [,3]         [,4]
#[1,]  113 9.734513e-02 9.828712e-02 9.802302e-02
#[2,]  113 9.734513e-02 8.884718e-02 8.916959e-02
#[3,]  113 5.309735e-02 6.266840e-01 6.191004e-01
#[4,]  115 1.304348e-01 2.412832e-02 2.662226e-02
#[5,]  113 8.849558e-03 5.062859e-01 4.963266e-01
#[6,]  120 1.620528e-05 4.332017e-06 5.129438e-05
#[7,]  120 1.652095e-03 7.090201e-04 7.632346e-04
#[8,]  120 5.090087e-06 1.166679e-06 6.770817e-06






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
