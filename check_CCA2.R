#!/usr/bin/env Rscript
#Use GTEx junction data

rm(list=ls())
library(data.table)
library(PMA)
library(GenomicRanges)
library(glmnet)
removeconstrows=function(dat)
{
  tmp=apply(dat,1,sd)
  idxconst=tmp==0
  # idxconst=rep(F,nrow(dat))
  # for (i in 1:nrow(dat))
  # {
  #   if (i %% 10000==0) cat(i,"..")
  #   if (var(unlist(dat[i,]))==0 | is.na(var(unlist(dat[i,])))) idxconst[i]=T
  # }
  dat=dat[!idxconst,]
}
#https://stackoverflow.com/questions/18275639/remove-highly-correlated-variables
removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
  return(datnew)
}
corquantile=function(dat=t(g_CCAX))
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  print(quantile(tmp[tmp!=0],c(0,0.25,0.5,0.75,0.9,0.99,1)))
}

#work on GTEX junction examples
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/genotype_phenotype_EC_junction.RData")
#geneexp
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
proteingenes=unique(gtex_ge_anno$Symbol[gtex_ge_anno$gene_type=="protein_coding"])
tmp=GE[rownames(GE) %in% proteingenes,]
idx=match(rownames(tmp),rownames(GE))
phenotype=tmp
phenotypepos=GE_pos[idx,]

# tmp=removeconstrows(snp)
# idx=match(rownames(tmp),rownames(snppos))
# snppos=snppos[idx,]
# snp=tmp
# #save(snp,snppos,GE,GE_pos,gtex_covariate,file="../result/genotype_phenotype_EC_junction_rmconstsnp.RData")
load("../result/genotype_phenotype_EC_junction_rmconstsnp.RData")

genes=c("MEF2B","DDX49","KXD1","CERS1","COMP","ARMC6","SGK223","TMEM161A","CRTC1","HOMER3")
genes=c("MEF2B","DDX49","KXD1","CERS1","COMP","ARMC6","TMEM161A","CRTC1","HOMER3")
idx=match(genes,phenotypepos$geneid)
gene="DDX49"
gene="CERS1"
gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
gr_snp=GRanges(seqnames = snppos$chr,ranges = IRanges(start=snppos$pos,width=1)) #snp
gwassnps=c("19:18803172_T_G", "19:18804294_T_G", "19:18817903_A_C")
#glmnet
# #predicted geneexp using crossvalidation
# fitted_cv=function(Xsel,covariateall,Y,ncv=5)
# {
#   Xall=data.matrix(cbind(Xsel,covariateall))
#   maxnumvar=floor(length(Y)*(1-1/ncv))
#   if (ncol(Xall)>maxnumvar-1) #number of covariates is greater than sample size, select subset of covariates
#   {
#     lmfit1=lm(Y~Xall)
#     lmcoeff1=summary(lmfit1)$coefficients
#     rownames(lmcoeff1)=gsub("Xall","",rownames(lmcoeff1))
#     lmleftsnp1=rownames(lmcoeff1)[rownames(lmcoeff1) %in% colnames(Xsel)] 
#     idx1=match(lmleftsnp1,rownames(lmcoeff1))
#     lmleftsnp1=lmleftsnp1[order(abs(lmcoeff1[idx1,1]),decreasing = T)]
#     idx1=match(lmleftsnp1,colnames(Xsel))
#     Xsel=Xsel[,idx1]
#     Xsel=Xsel[,1:(maxnumvar-ncol(covariateall)-1)]
#     
#     Xall=data.matrix(cbind(Xsel,covariateall))
#   }
#   
#   fitted1=rep(0,length(Y))
#   set.seed(10000)
#   permutidx=sample(1:length(Y))
#   idxs=as.integer(seq(1,length(Y)+1,length.out = ncv+1)) #boundary points of cv segments
#   
#   for (ii in 1:ncv)
#   {
#     idx_predict=rep(F,length(Y))
#     idx_predict[idxs[ii]:(idxs[ii+1]-1)]=T
#     trainfm=lm(Y[permutidx[!idx_predict]]~Xall[permutidx[!idx_predict],])
#     traincoeff=summary(trainfm)$coefficients
#     rownames(traincoeff)=gsub("Xall[permutidx[!idx_predict], ]","",rownames(traincoeff),fixed = T)
#     trainleftsnps=rownames(traincoeff)[rownames(traincoeff) %in% colnames(Xsel)]
#     numvar=length(trainleftsnps)
#     idx1=match(trainleftsnps,colnames(Xsel))
#     Xsel1=Xsel[,idx1]
#     if (numvar==1)
#     {
#       Xsel1=matrix(Xsel1,ncol=1)
#       colnames(Xsel1)=trainleftsnps
#     }
#     fitted1[permutidx[idx_predict]]=rep(traincoeff[1,1],sum(idx_predict)) #intercept term
#     idx1=match(trainleftsnps,rownames(traincoeff))
#     if (numvar>0) #to add each selected snp term
#     {
#       for (j in 1:numvar)
#       {
#         fitted1[permutidx[idx_predict]]=fitted1[permutidx[idx_predict]]+Xsel1[permutidx[idx_predict],j]*traincoeff[idx1[j],1]
#       }
#     }
#   }
#   return(fitted1)
# }

glmnetmodel=function(i,opt="1se",ncv=10,distcutoff=5e5,thesnp=snp,gr_thesnp=gr_snp) #ncv is not used here
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_thesnp,gr_pos[i])
  idx=which(tmp<distcutoff)
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  if (length(idx)>1)
  {
    X=t(thesnp[idx,])
    X0=X
    X1=t(removeconstrows(t(X0)))
    X=removehighcorr(X1,corcutoff = 0.9)

    
    set.seed(i+3000)
    cvfit=tryCatch(
      {
        cv.glmnet(data.matrix(X),Y,nlambda=100,nfolds=10,alpha=1)
      },
      error=function(e)
      {
        return(F)
      }
    )
    #plot(cvfit)
    if (is.list(cvfit))
    {
      fit=glmnet(as.matrix(X),Y,nlambda = 100, alpha=1)
      lambda_1varsel=fit$lambda[which(fit$df>0)][1] #lambda value when 1 snp variable is selected should be greater than this
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
              nummaxvar=min(nrow(X)-1,numvar)
              numvar=nummaxvar
              Xsel=Xsel[,1:nummaxvar,drop=F]
            }
            # if (numvar==1) #keep Xsel as in matrix form
            # {
            #   Xsel=matrix(Xsel,ncol=1)
            #   colnames(Xsel)=glmleftsnp
            # }
            fit1=lm(Y~Xsel) # to remove snps with NA coefficient due to colinearity
            #summary(fit1)$r.squared
            lmcoeff=summary(fit1)$coefficients
            rownames(lmcoeff)=gsub("Xsel","",rownames(lmcoeff))
            #align up coeff with Xsel
            lmleftsnp=rownames(lmcoeff)[rownames(lmcoeff) %in% colnames(Xsel)] 
            numvar=length(lmleftsnp)
            if (numvar>0)
            {
              glmflag=1
              idx1=match(lmleftsnp,rownames(lmcoeff))
              selectedsnps=paste0(rownames(lmcoeff)[idx1],collapse = "|")
              selectedsnps_coeff=paste0(lmcoeff[idx1,1],collapse = "|")
            
            }
            
          }
        }
      }
    }
  }
  return(list(glmflag=glmflag,numvar=numvar,numsnpall=length(idx),numsnp=ncol(X),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              X=X,X0=X0,X1=X1,Y=Y,fit1=fit1))
}

glmres=glmnetmodel(i=which(phenotypepos$geneid==gene),opt="min",ncv=10,distcutoff = 5e5)
summary(glmres$fit1)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          -0.2274     0.1490  -1.526  0.12865   
# Xsel19:18642448_T_A   0.6084     0.2528   2.406  0.01708 * 
#   Xsel19:19063806_G_A   0.2897     0.1007   2.877  0.00447 **
#   Xsel19:18895865_T_C   0.4587     0.2034   2.256  0.02524 * 
#   Xsel19:18916295_G_A  -0.2341     0.1007  -2.325  0.02113 * 
#   Xsel19:18618565_C_T   0.3712     0.2203   1.685  0.09363 .

glmselectedsnps=unlist(strsplit(glmres$selectedsnps,"|",fixed = T))
idx=which(rownames(phenotype)==gene)
cor(as.numeric(phenotype[idx,]),glmres$fit1$fitted.values) #[1] [1] 0.4353065
sum(glmselectedsnps %in% gwassnps) #0

##CCA
idx=which(phenotypepos$geneid==gene)
chr=phenotypepos$chr[idx]
startloc=phenotypepos$s1[idx]-5e5
startloc=max(1,startloc)
endloc=phenotypepos$s2[idx]+5e5
gr_region=GRanges(seqnames = chr,ranges = IRanges(start=startloc,end=endloc))
tmp=distance(gr_pos,gr_region)
sum(tmp==0,na.rm=T) #31 genes
idx=which(tmp==0)
CCAZ=phenotype[idx,,drop=F]
tmp=distance(gr_snp,gr_region)
sum(tmp==0,na.rm=T) #1685 snps
idx=which(tmp==0)
CCAX=snp[idx,]
tmp1=rownames(CCAX)[!rownames(CCAX) %in% colnames(glmres$X0)] #all 0 genotype
CCAX=t(glmres$X0)
Sys.time()
perm.out <- CCA.permute(t(CCAX),t(CCAZ),typex="standard",typez="standard",nperms=100,trace=F)
Sys.time()
CCAout <- CCA(t(CCAX),t(CCAZ),typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,v=perm.out$v.init)
table(CCAout$u!=0)
# FALSE  TRUE 
# 1612    73 
table(CCAout$v!=0)
# FALSE  TRUE 
# 30     1 
CCAout$u[which(rownames(CCAX) %in% glmselectedsnps)] #[1] -0.06906033  0.07292235
CCAout$v[which(rownames(CCAZ)==gene)] #[1] -0.3479262
quantile(CCAout$v)
# 0%         25%         50%         75%        100% 
# -0.36847037 -0.23609347 -0.07508655  0.00000000  0.00000000 
rownames(CCAZ)[which.max(abs(CCAout$v))] #"KIAA0892"

CCAX1=t(glmres$X) #remove high correlated snps
perm.out1 <- CCA.permute(t(CCAX1),t(CCAZ),typex="standard",typez="standard",nperms=1000,trace=F)
CCAout1 <- CCA(t(CCAX1),t(CCAZ),typex="standard",typez="standard",K=1,penaltyx=perm.out1$bestpenaltyx,penaltyz=perm.out1$bestpenaltyz,v=perm.out1$v.init)
table(CCAout1$u!=0)
# FALSE  TRUE 
# 284    28 
table(CCAout1$v!=0)
# FALSE  TRUE 
# 28     3 
#rownames(CCAX1)[CCAout1$u!=0]
CCAout1$u[which(rownames(CCAX1) %in% glmselectedsnps)] #[1]  0.00000000  0.12276873 -0.09741223  0.13678799  0.15257050  0.00000000  0.23376011  0.00000000
CCAgenes=rownames(CCAZ)[CCAout1$v!=0]
#CCAout1$v[CCAout1$v!=0]
x=matrix(CCAout1$u,nrow = 1)
x=x%*%CCAX1
idx=which(rownames(phenotype)==gene)
fit=glm(as.numeric(phenotype[idx,])~as.numeric(x))
cor(as.numeric(phenotype[idx,]),fit$fitted.values) #[1] 0.2386967
heatmap(cor(t(CCAZ)))
library(gplots)
heatmap.2(cor(t(CCAZ)), col=redgreen(75), 
          density.info="none", trace="none", dendrogram="none", 
          symm=T,symkey=T,symbreaks=T, scale="none")
tmp=cor(t(CCAZ))
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
idx=which(rownames(tmp) %in% CCAgenes)
idx1=which(rownames(tmp)==gene)
tmp1=NULL
for (i in 1:length(idx))
{
  print(rownames(tmp)[idx[i]])
  tmp1=c(tmp1,tmp[max(idx[i],idx1),min(idx[i],idx1)])
  print(tmp[max(idx[i],idx1),min(idx[i],idx1)])
}
quantile(tmp[tmp!=0])

idx=rownames(CCAZ) %in% CCAgenes
CCAZ1=rbind(CCAZ[idx,],CCAZ[!idx,])
heatmap.2(cor(t(CCAZ1)), col=redgreen(75), 
          density.info="none", trace="none", dendrogram="none", 
          symm=T,symkey=T,symbreaks=T, scale="none")


# 
# #use genotyped data------------------------
# load("/fh/fast/dai_j/CancerGenomics/EAprogression/data/TCGA_EAC_Genotype_Genexp.RData")
# idx=match(colnames(phenotype),colnames(genotypedata))
# snp1=genotypedata[,idx]
# snp1=removeconstrows(snp1)
# snp6anno=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/Tools/database/other/GenomeWideSNP_6.na35.annot.csv",header=T))
# snp6anno$`Physical Position`=as.numeric(snp6anno$`Physical Position`)
# snp1=snp1[rownames(snp1) %in% snp6anno$`Probe Set ID`[!is.na(snp6anno$`Physical Position`)],]
# idx=match(rownames(snp1),snp6anno$`Probe Set ID`)
# gr_snp1=GRanges(seqnames = snp6anno$Chromosome[idx],ranges = IRanges(start=snp6anno$`Physical Position`[idx],width=1))
# g_glmres=glmnetmodel(i=which(phenotypepos$geneid==gene),opt="min",ncv=10,distcutoff = 5e5,thesnp=snp1,gr_thesnp = gr_snp1)
# summary(g_glmres$fit1)
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)  
# # (Intercept)        0.13769    0.22695   0.607   0.5457  
# # XselSNP_A-4260408 -0.19259    0.09493  -2.029   0.0457 *
# #   XselSNP_A-4202649 -0.24283    0.13087  -1.856   0.0670 .
# g_glmselectedsnps=unlist(strsplit(g_glmres$selectedsnps,"|",fixed = T)) #"SNP_A-4260408" "SNP_A-4202649"

# #CCA--
# tmp=distance(gr_snp1,gr_region)
# sum(tmp==0,na.rm=T) #134 snps
# idx=which(tmp==0)
# g_CCAX=snp1[idx,]
# tmp1=rownames(g_CCAX)[!rownames(g_CCAX) %in% colnames(g_glmres$X0)] #all 0 genotype
# g_CCAX=t(g_glmres$X1)
# Sys.time()
# g_perm.out <- CCA.permute(t(g_CCAX),t(CCAZ),typex="standard",typez="standard",trace=F,nperms=10000)
# Sys.time()
# g_CCAout <- CCA(t(g_CCAX),t(CCAZ),typex="standard",typez="standard",K=1,penaltyx=g_perm.out$bestpenaltyx,penaltyz=g_perm.out$bestpenaltyz,v=g_perm.out$v.init)
# table(g_CCAout$u!=0)
# # FALSE  TRUE 
# # 107    13 
# rownames(g_CCAX)[which(g_CCAout$u!=0)] 
# table(g_CCAout$v!=0) 
# # FALSE  TRUE 
# # 32     1
# rownames(CCAZ)[which(g_CCAout$v!=0)] #"C19orf50"
# 
# g_CCAout1 <- CCA(t(g_CCAX),t(CCAZ),typex="standard",typez="standard",K=1,penaltyx=0.7,penaltyz=0.7,v=g_perm.out$v.init)
# table(g_CCAout1$u!=0)
# table(g_CCAout1$v!=0) 
