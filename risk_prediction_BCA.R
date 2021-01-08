###########################################
### risk prediction using BCA genotypes ###
###########################################

rm(list=ls())
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid

Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)

load(file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output.RData") 


### predicting BE ###

library(pROC)
library(glmnet)

topsnps.be <- row.names(BCAout)[BCAout[,1]<0.01]

geno <- allgenotype[row.names(allgenotype)%in% topsnps.be,sampletable$site<30]
Y <- sampletable$phenoBE_bca[sampletable$site<30]
Xall <- cbind(t(geno),Covariate[sampletable$site<30,])
Xall <- Xall[!is.na(Y),]
Y <- Y[!is.na(Y)]
Y <- (Y==2)*1
penalty=rep(1,ncol(Xall))
penalty[(nrow(geno)+1):length(penalty)]=0 #force the covariates to be included in the model
library(glmnet)
cvfit <- cv.glmnet(data.matrix(Xall),Y,family = "binomial",nfolds=5, penalty.factor=penalty,standardize=F)

plot(cvfit)

cam.y <- sampletable$phenoBE_bca[sampletable$site>=30 & !is.na(sampletable$phenoBE_bca) ]
cam.genotype <-  allgenotype[row.names(allgenotype)%in% topsnps.be,sampletable$site>=30 & !is.na(sampletable$phenoBE_bca)]

#lmcoeff=coef(cvfit,s=cvfit$lambda.min)
lset <- cvfit$lambda[cvfit$lambda>=cvfit$lambda.1se]

outroc <- rep(0,length(lset))
for (i in 1:length(lset)) {
 lmcoeff=coef(cvfit,s=lset[i])
 grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])
 roc1 <- roc(cam.y~grs)
 outroc[i] <- roc1$auc
}
plot(lset, outroc)

lmcoeff=coef(cvfit,s=cvfit$lambda.min)
grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])
roc1 <- roc(cam.y~grs)
roc1$auc
plot(roc1)

topsnps.ea <- row.names(BCAout)[BCAout[,3]<0.01]
geno <- allgenotype[row.names(allgenotype)%in% topsnps.ea,sampletable$site<30]
Y <- sampletable$phenoEA_bca[sampletable$site<30]
Xall <- cbind(t(geno),Covariate[sampletable$site<30,])
Xall <- Xall[!is.na(Y),]
Y <- Y[!is.na(Y)]
Y <- (Y==2)*1
penalty=rep(1,ncol(Xall))
penalty[(nrow(geno)+1):length(penalty)]=0 #force the covariates to be included in the model
library(glmnet)
cvfit <- cv.glmnet(data.matrix(Xall),Y,family = "binomial",nfolds=5, penalty.factor=penalty,standardize=F)

cam.y <- sampletable$phenoEA_bca[sampletable$site>=30 & !is.na(sampletable$phenoEA_bca) ]
cam.genotype <-  allgenotype[row.names(allgenotype)%in% topsnps.ea,sampletable$site>=30 & !is.na(sampletable$phenoEA_bca)]


plot(cvfit)

#lmcoeff=coef(cvfit,s=cvfit$lambda.min)
lset <- cvfit$lambda[cvfit$lambda>=cvfit$lambda.1se]

outroc <- rep(0,length(lset))
for (i in 1:length(lset)) {
  lmcoeff=coef(cvfit,s=lset[i])
  grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])
  roc1 <- roc(cam.y~grs)
  outroc[i] <- roc1$auc
}
plot(lset, outroc)





topsnps.eabe <- row.names(BCAout)[BCAout[,5]<0.01]
geno <- allgenotype[row.names(allgenotype)%in% topsnps.eabe,sampletable$site<30]
Y <- sampletable$phenoEABE_bca[sampletable$site<30]
Xall <- cbind(t(geno),Covariate[sampletable$site<30,])
Xall <- Xall[!is.na(Y),]
Y <- Y[!is.na(Y)]
Y <- (Y==2)*1
penalty=rep(1,ncol(Xall))
penalty[(nrow(geno)+1):length(penalty)]=0 #force the covariates to be included in the model
library(glmnet)
cvfit <- cv.glmnet(data.matrix(Xall),Y,family = "binomial",nfolds=5, penalty.factor=penalty,standardize=F)

plot(cvfit)

lmcoeff=coef(cvfit,s=cvfit$lambda.1se)
cam.y <- sampletable$phenoEABE_bca[sampletable$site>=30 & !is.na(sampletable$phenoEABE_bca) ]
cam.genotype <-  allgenotype[row.names(allgenotype)%in% topsnps.ea,sampletable$site>=30 & !is.na(sampletable$phenoEABE_bca)]
grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])

roc1 <- roc(cam.y~grs)
plot(roc1,ylim=c(0,1),print.auc=TRUE)


######################################
#### now consider progression risk ###
######################################

topsnps.eabe <- row.names(BCAout)[BCAout[,7]<0.001]
geno <- allgenotype[row.names(allgenotype)%in% topsnps.eabe,sampletable$site<30]

Y<- sampletable$phenoEA_bca[sampletable$site<30] 
Y[is.na(Y)] = 3
Y[Y==1] = NA


Xall <- cbind(t(geno),Covariate[sampletable$site<30,])
Xall <- Xall[!is.na(Y),]
Y <- Y[!is.na(Y)]
Y <- (Y==2)*1
penalty=rep(1,ncol(Xall))
penalty[(nrow(geno)+1):length(penalty)]=0 #force the covariates to be included in the model
cvfit <- cv.glmnet(data.matrix(Xall),Y,family = "binomial",nfolds=5, penalty.factor=penalty,standardize=F)

plot(cvfit)


#lmcoeff=coef(cvfit,s=cvfit$lambda.min)
lset <- cvfit$lambda[cvfit$lambda>=cvfit$lambda.min]
cam.y <- sampletable$phenoEA_bca[sampletable$site>=30]
cam.y[is.na(cam.y)] = 3
cam.y[cam.y==1] = NA
cam.genotype <-  allgenotype[row.names(allgenotype) %in% topsnps.eabe,sampletable$site>=30]
cam.genotype <-  cam.genotype[,!is.na(cam.y)]
cam.y <- cam.y[!is.na(cam.y)]


####
outroc <- rep(0,length(lset))
for (i in 1:length(lset)) {
  lmcoeff=coef(cvfit,s=lset[i])
  grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])
  roc1 <- roc(cam.y~grs)
  outroc[i] <- roc1$auc
}
plot(lset, outroc)
####

lmcoeff=coef(cvfit,s=cvfit$lambda.min)
grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])

roc1 <- roc(cam.y~grs)
plot(roc1,ylim=c(0,1),print.auc=TRUE)
qcut <- quantile(grs,c(0.25,0.5,0.75))
gstrata <- grs
gstrata <- ifelse(grs<qcut[1],1,gstrata)
gstrata <- ifelse(grs>=qcut[1]&grs<qcut[2],2,gstrata)
gstrata <- ifelse(grs>=qcut[2]&grs<qcut[3],3,gstrata)
gstrata <- ifelse(grs>=qcut[3],4,gstrata)

summary(glm(I(cam.y==2)~factor(gstrata),family=binomial))

summary(glm(I(cam.y==2)~factor(gstrata==4),family=binomial))

############################################
#### do GWAS and prediction among eQTLs ####
#### include interaction                ####
############################################

rm(list=ls())
library(data.table)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output_pc10.RData")
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_GERD_interaction_output.RData")
#load(file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_output.RData") 
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_cig_interaction_output.RData")
#load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_gwas_bmi_interaction_output.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_junction_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_stomach_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_blood_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_mucosa_V7_eQTL_gene_pairs.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTEx_muscularis_V7_eQTL_gene_pairs.RData")


lookuptable=fread("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt")

#jointeqtl<- jointeqtl[jointeqtl$pval_nominal < 2e-11,]

lookuptable1 <- lookuptable[lookuptable$variant_id %in% jointeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% bloodeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% stomacheqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% mucosaeqtl$variant_id,]
lookuptable1 <- lookuptable[lookuptable$variant_id %in% musculariseqtl$variant_id,]

idx <- match(musculariseqtl$variant_id,lookuptable1$variant_id)
idx <- match(mucosaeqtl$variant_id,lookuptable1$variant_id)
idx <- match(stomacheqtl$variant_id,lookuptable1$variant_id)
idx <- match(bloodeqtl$variant_id,lookuptable1$variant_id)
idx <- match(jointeqtl$variant_id,lookuptable1$variant_id)


lookuptable1 <- lookuptable1[idx,]
BCAout1 <- BCAout[bim$V2 %in%lookuptable1$rs_id_dbSNP147_GRCh37p13, ]
bim1 <- bim[bim$V2 %in%lookuptable1$rs_id_dbSNP147_GRCh37p13, ]

0.05/nrow(BCAout1)
bim1[order(BCAout1[,1]),][1:2,]
BCAout1[order(BCAout1[,1]),1][1:2]
bim1[order(BCAout1[,3]),][1:2,]
BCAout1[order(BCAout1[,3]),3][1:2]
bim1[order(BCAout1[,5]),][1:2,]
BCAout1[order(BCAout1[,5]),5][1:2]
bim1[order(BCAout1[,7]),][1:2,]
BCAout1[order(BCAout1[,7]),7][1:2]

## Joint ##
#> 0.05/nrow(BCAout1)
#[1] 1.00321e-06
#> bim1[order(BCAout1[,1]),][1:2,]
#V1         V2 V3       V4 V5 V6
#725334 19 rs10423674  0 18817903  A  B
#290435  6  rs3132545  0 31089125  A  B
#> BCAout1[order(BCAout1[,1]),1][1:2]
#rs10423674    rs3132545 
#4.239830e-06 3.679096e-05 
#> bim1[order(BCAout1[,3]),][1:2,]
#V1        V2 V3       V4 V5 V6
#303161  6 rs3904819  0 58033394  B  A
#303153  6 rs4314501  0 57858800  A  B
#> BCAout1[order(BCAout1[,3]),3][1:2]
#rs3904819    rs4314501 
#1.399651e-06 9.109568e-06 
#> bim1[order(BCAout1[,5]),][1:2,]
#V1         V2 V3        V4 V5 V6
#725334 19 rs10423674  0  18817903  A  B
#108206  2 rs16844715  0 160915106  A  B
#> BCAout1[order(BCAout1[,5]),5][1:2]
#rs10423674   rs16844715 
#2.196322e-07 2.855025e-06 
#> bim1[order(BCAout1[,7]),][1:2,]
#V1         V2 V3        V4 V5 V6
#127925  2  rs2278737  0 239079542  B  A
#498796 11 rs10769688  0   6347072  A  B
#> BCAout1[order(BCAout1[,7]),7][1:2]
#rs2278737   rs10769688 
#2.743659e-05 5.048603e-05 



## stomach ##
#> 0.05/nrow(BCAout1)
#[1] 1.146342e-06
#> bim1[order(BCAout1[,1]),][1:2,]
#V1        V2 V3       V4 V5 V6
#725325 19 rs3829671  0 18777667  A  B
#290435  6 rs3132545  0 31089125  A  B
#> BCAout1[order(BCAout1[,1]),1][1:2]
#rs3829671    rs3132545 
#5.408354e-06 3.679096e-05 
#> bim1[order(BCAout1[,3]),][1:2,]
#V1        V2 V3       V4 V5 V6
#303161  6 rs3904819  0 58033394  B  A
#303153  6 rs4314501  0 57858800  A  B
#> BCAout1[order(BCAout1[,3]),3][1:2]
#rs3904819    rs4314501 
#1.399651e-06 9.109568e-06 
#> bim1[order(BCAout1[,5]),][1:2,]
#V1         V2 V3        V4 V5 V6
#725325 19  rs3829671  0  18777667  A  B
#108206  2 rs16844715  0 160915106  A  B
#> BCAout1[order(BCAout1[,5]),5][1:2]
#rs3829671   rs16844715 
#7.487201e-07 2.855025e-06 
#> bim1[order(BCAout1[,7]),][1:2,]
#V1         V2 V3        V4 V5 V6
#416222  8 rs10875471  0 142364833  A  B
#416223  8  rs7843131  0 142366200  B  A
#> BCAout1[order(BCAout1[,7]),7][1:2]
#rs10875471    rs7843131 
#3.528438e-07 6.033983e-06 


#muscularis
#> 0.05/nrow(BCAout1)
#[1] 5.710044e-07
#> bim1[order(BCAout1[,1]),][1:2,]
#V1        V2 V3        V4 V5 V6
#62845  1  rs653737  0 236729956  B  A
#62870  1 rs4659695  0 236796979  A  B
#> BCAout1[order(BCAout1[,1]),1][1:2]
#rs653737    rs4659695 
#7.843093e-06 1.388109e-05 
#> bim1[order(BCAout1[,3]),][1:2,]
#V1        V2 V3       V4 V5 V6
#303161  6 rs3904819  0 58033394  B  A
#440021  9 rs1057713  0 96714161  B  A
#> BCAout1[order(BCAout1[,3]),3][1:2]
#rs3904819    rs1057713 
#1.399651e-06 2.701303e-06 
#> bim1[order(BCAout1[,5]),][1:2,]
#V1         V2 V3        V4 V5 V6
#108206  2 rs16844715  0 160915106  A  B
#782356 23  rs6632474  0  15387063  B  A
#> BCAout1[order(BCAout1[,5]),5][1:2]
#rs16844715    rs6632474 
#2.855025e-06 6.161993e-06 
#> bim1[order(BCAout1[,7]),][1:2,]
#V1         V2 V3        V4 V5 V6
#416222  8 rs10875471  0 142364833  A  B
#416223  8  rs7843131  0 142366200  B  A
#> BCAout1[order(BCAout1[,7]),7][1:2]
#rs10875471    rs7843131 
#3.528438e-07 6.033983e-06 


### blood ###
#> 0.05/nrow(BCAout1)
#[1] 7.387161e-07
#> bim1[order(BCAout1[,1]),][1:2,]
#V1         V2 V3       V4 V5 V6
#725330 19 rs10419226  0 18803172  A  B
#226950  5  rs2292596  0   422955  B  A
#> BCAout1[order(BCAout1[,1]),1][1:2]
#rs10419226    rs2292596 
#9.796856e-09 8.484716e-06 
#> bim1[order(BCAout1[,3]),][1:2,]
#V1         V2 V3       V4 V5 V6
#725330 19 rs10419226  0 18803172  A  B
#303161  6  rs3904819  0 58033394  B  A
#> BCAout1[order(BCAout1[,3]),3][1:2]
#rs10419226    rs3904819 
#6.496396e-07 1.399651e-06 
#> bim1[order(BCAout1[,5]),][1:2,]
#V1         V2 V3        V4 V5 V6
#725330 19 rs10419226  0  18803172  A  B
#108206  2 rs16844715  0 160915106  A  B
#> BCAout1[order(BCAout1[,5]),5][1:2]
#rs10419226   rs16844715 
#4.979043e-10 2.855025e-06 
#> bim1[order(BCAout1[,7]),][1:2,]
#V1        V2 V3       V4 V5 V6
#514988 11 rs3740629  0 68343593  A  B
#244222  5 rs4288083  0 68601031  A  B
#> BCAout1[order(BCAout1[,7]),7][1:2]
#rs3740629    rs4288083 
#5.572831e-05 6.738760e-05 



#plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10),-log(BCAout1[order(BCAout1[,1],decreasing=T),1],base=10),xlab="-log10 expected p-value",ylab="-log10 observed p-value")
#abline(0,1)

mpoint <- ceiling(quantile(1:nrow(BCAout),0.95)/200)
pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(BCAout))

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,5],decreasing=T),5],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)
bim[order(BCAout[,5]),][1:2,]
BCAout[order(BCAout[,5]),][1:2,]
#rs10419226    rs6023406 
#5.155810e-09 9.782792e-07 
plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,6],decreasing=T),6],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,10],decreasing=T),10],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,15],decreasing=T),15],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)
bim[order(BCAout[,15]),][1:2,]
BCAout[order(BCAout[,15]),15][1:2]


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,20],decreasing=T),20],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,7],decreasing=T),7],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,9],decreasing=T),9],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,11],decreasing=T),11],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,13],decreasing=T),13],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout):1)/nrow(BCAout),base=10)[pindex],-log(BCAout[order(BCAout[,15],decreasing=T),15],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)






mpoint <- ceiling(quantile(1:nrow(BCAout1),0.95)/200)
pindex <- c(1,(1:mpoint)*200,(mpoint*200+1):nrow(BCAout1))

png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/eQTL_bmi_interaction.png")
plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,1],decreasing=T),1],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value",main="eQTL in GTEx Blood: interaction with BMI",cex.main=1.5)
x1 <- -log((1)/nrow(BCAout1),base=10)
y1 <- max(-log(BCAout1[order(BCAout1[,1],decreasing=T),1],base=10))
points(x1,y1,col=2,pch=19,cex=1.5)
text(x1-0.75,y1-0.25,"rs491603",col=2,cex=1.5)
abline(0,1)
dev.off()
bim1[order(BCAout1[,5]),][1:2,]
BCAout1[order(BCAout1[,5]),][1:2,]

#> bim1[order(BCAout1[,5]),][1:2,]
#V1         V2 V3       V4 V5 V6
#13574   1   rs491603  0 36532316  A  B
#725330 19 rs10419226  0 18803172  A  B
#> BCAout1[order(BCAout1[,5]),][1:2,]
#[,1]         [,2]         [,3]       [,4]         [,5]       [,6]        [,7]
#rs491603   8.365187e-08  0.072675877 8.964474e-09 -2.1699241 2.367889e-09 0.02192083  0.03550134
#rs10419226 5.054819e-01 -0.006134448 1.164319e-01  0.4058423 1.202828e-07 0.34797211 -0.01083827
#[,8]       [,9]        [,10]        [,11]       [,12]        [,13]      [,14]
#rs491603   0.02814514 -0.9484700 0.0602803683 2.046080e-06  0.05836066 4.183309e-07 -1.7130944
#rs10419226 0.11199399  0.5187103 0.0003606238 2.319392e-01 -0.01035428 3.250944e-02  0.5169101
#[,15]     [,16]        [,17]      [,18]     [,19]      [,20]
#rs491603   4.343484e-07 0.1068461 -0.023510576 0.05367899 0.8125203 0.04002785
#rs10419226 1.311031e-08 0.4707455 -0.008023433 0.51875740 0.2066777 0.71485595



plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,3],decreasing=T),3],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,5],decreasing=T),5],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,6],decreasing=T),6],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,10],decreasing=T),10],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

bim1[order(BCAout1[,10]),][1:2,]
BCAout1[order(BCAout1[,10]),][1:2,]

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,11],decreasing=T),11],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,15],decreasing=T),15],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

bim1[order(BCAout1[,15]),][1:2,]
BCAout1[order(BCAout1[,15]),][1:2,]

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,16],decreasing=T),16],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,20],decreasing=T),20],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)




plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,11],decreasing=T),11],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,13],decreasing=T),13],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)


plot(-log((nrow(BCAout1):1)/nrow(BCAout1),base=10)[pindex],-log(BCAout1[order(BCAout1[,15],decreasing=T),15],base=10)[pindex],xlab="-log10 expected p-value",ylab="-log10 observed p-value")
abline(0,1)

#######################################
## check the validation summary sets ##
#######################################

rm(list=ls())

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"

summaryfile=paste0(summaryfolder,"BE_Bonn_autosomes.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_autosomes.txt")
summaryfile2=paste0(summaryfolder,"BE_Cambridge_autosomes.txt")
summaryfile3=paste0(summaryfolder,"EA_Bonn_autosomes.txt")
summaryfile4=paste0(summaryfolder,"BEEA_Bonn_autosomes.txt")

#beaconfile=paste0(beaconfolder,"imp_GTEx_BE_CO.assoc.logistic")

library(data.table)
#summary stat from validation
summarydat=as.data.frame(fread(summaryfile,header=T))
summarydat1=as.data.frame(fread(summaryfile1,header=T))
summarydat2=as.data.frame(fread(summaryfile2,header=T))
summarydat3=as.data.frame(fread(summaryfile3,header=T))
summarydat4=as.data.frame(fread(summaryfile4,header=T))



summaryfile=paste0(summaryfolder,"BE_Bonn_Xchromosome.txt")
summaryfile1=paste0(summaryfolder,"BE_oxford_Xchromosome.txt")

summaryfile2=paste0(summaryfolder,"EA_Bonn_Xchromosome.txt")
summaryfile3=paste0(summaryfolder,"BEEA_Bonn_Xchromosome.txt")

#summary stat from validation
summarydatx=as.data.frame(fread(summaryfile,header=T))
summarydat1x=as.data.frame(fread(summaryfile1,header=T))

summarydat2x=as.data.frame(fread(summaryfile2,header=T))
summarydat3x=as.data.frame(fread(summaryfile3,header=T))


summarydat[summarydat$SNP=="rs10875471",]
summarydat1[summarydat1$rsid=="rs10875471",]
summarydat2[summarydat2$SNP=="rs10875471",]
summarydat3[summarydat3$SNP=="rs10875471",]

summarydat[summarydat$SNP=="rs7843131",]
summarydat1[summarydat1$rsid=="rs7843131",]
summarydat2[summarydat2$SNP=="rs7843131",]
summarydat3[summarydat3$SNP=="rs7843131",]

summarydat[summarydat$SNP=="rs3904819",]
summarydat1[summarydat1$rsid=="rs3904819",]
summarydat2[summarydat2$SNP=="rs3904819",]
summarydat3[summarydat3$SNP=="rs3904819",]

summarydat[summarydat$SNP=="rs4314501",]
summarydat1[summarydat1$rsid=="rs4314501",]
summarydat2[summarydat2$SNP=="rs4314501",]
summarydat3[summarydat3$SNP=="rs4314501",]

summarydat[summarydat$SNP=="rs16844715",]
summarydat1[summarydat1$rsid=="rs16844715",]
summarydat2[summarydat2$SNP=="rs16844715",]
summarydat3[summarydat3$SNP=="rs16844715",]


summarydat[summarydat$SNP=="rs653737",]
summarydat1[summarydat1$rsid=="rs653737",]
summarydat2[summarydat2$SNP=="rs653737",]
summarydat3[summarydat3$SNP=="rs653737",]

summarydat[summarydat$SNP=="rs1057713",]
summarydat1[summarydat1$rsid=="rs1057713",]
summarydat2[summarydat2$SNP=="rs1057713",]
summarydat3[summarydat3$SNP=="rs1057713",]


summarydatx[summarydatx$SNP=="rs6632474",]
summarydat1x[summarydat1x$rsid=="rs6632474",]
summarydat2x[summarydat2x$SNP=="rs6632474",]
summarydat3x[summarydat3x$SNP=="rs6632474",]


summarydat[summarydat$SNP=="rs491603",]
summarydat1[summarydat1$rsid=="rs491603",]
summarydat2[summarydat2$SNP=="rs491603",]
summarydat3[summarydat3$SNP=="rs491603",]


summarydat[summarydat$SNP=="rs11765529",]
summarydat1[summarydat1$rsid=="rs11765529",]
summarydat2[summarydat2$SNP=="rs11765529",]
summarydat3[summarydat3$SNP=="rs11765529",]
summarydat4[summarydat4$SNP=="rs11765529",]


summarydat[summarydat$SNP=="rs7798662",]
summarydat1[summarydat1$rsid=="rs7798662",]
summarydat2[summarydat2$SNP=="rs7798662",]
summarydat3[summarydat3$SNP=="rs7798662",]
summarydat4[summarydat4$SNP=="rs7798662",]





#summary from beacon
#beacondat=as.data.frame(fread(beaconfile,header=T))

#idx=which(rownames(res_min) %in% names(fdrres[i]))
#selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))





## check whether the two SNPs in chr 2 for progression can be validated ##

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid

Covariate=sampletable[,colnames(sampletable) %in% c("sex","age","ev1_bca","ev2_bca","ev3_bca","ev4_bca")]
colnames(Covariate[which(colnames(Covariate)=="ev1_bca")])="pc1"
colnames(Covariate[which(colnames(Covariate)=="ev2_bca")])="pc2"
colnames(Covariate[which(colnames(Covariate)=="ev3_bca")])="pc3"
colnames(Covariate[which(colnames(Covariate)=="ev4_bca")])="pc4"
rownames(Covariate)=sampletable$localid
Covariate$age[is.na(Covariate$age)]=mean(Covariate$age,na.rm = T)




geno <- allgenotype[row.names(allgenotype)%in%lookuptable$rs_id_dbSNP147_GRCh37p13,sampletable$site>=30]
Y <- sampletable$phenoEA_bca[sampletable$site>=30]
Y[is.na(Y)] <- 3
Y[Y==1] <- NA

Covariate <- Covariate[sampletable$site>=30,]
Covariate <- Covariate[!is.na(Y),]

genotype <- geno[which(row.names(geno)=="rs7594872"),!is.na(Y)]
genotype <- geno[which(row.names(geno)=="rs17026352"),!is.na(Y)]
Y<- Y[!is.na(Y)]


fit1 <- glm(I(Y==2)~as.numeric(genotype)+sex+age+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=Covariate)

penalty=rep(1,ncol(Xall))
penalty[(nrow(geno)+1):length(penalty)]=0 #force the covariates to be included in the model
library(glmnet)
cvfit <- cv.glmnet(data.matrix(Xall),Y,family = "binomial",nfolds=10, penalty.factor=penalty,standardize=F)

plot(cvfit)



cam.y <- sampletable$phenoBE_bca[sampletable$site>=30 & !is.na(sampletable$phenoBE_bca) ]
cam.genotype <-  allgenotype[row.names(allgenotype)%in%lookuptable$rs_id_dbSNP147_GRCh37p13,sampletable$site>=30 & !is.na(sampletable$phenoBE_bca)]
lmcoeff=coef(cvfit,s=cvfit$lambda.min)
grs <- drop(t(cam.genotype) %*% lmcoeff[1:nrow(cam.genotype)])
roc1 <- roc(cam.y~grs)
  




######################
### try lassosum   ###
######################

install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")

library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
# We like to use dplyr for it makes codes much more readable
library(dplyr)
sum.stat <- "~/SCHARPHOME/ScharpFile/Grant/BEprogression/Height.QC.gz"
bfile <- "~/SCHARPHOME/ScharpFile/Grant/BEprogression/RiskPrediction/EUR.QC"
# Read in and process the covariates
covariate <- fread("~/SCHARPHOME/ScharpFile/Grant/BEprogression/RiskPrediction/EUR.covariate")
pcs <- fread("~/SCHARPHOME/ScharpFile/Grant/BEprogression/RiskPrediction/EUR.eigenvec")
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- as.data.frame(merge(covariate, pcs, by=c("FID", "IID")))

ld.file <- system.file("data", "Berisa.EUR.hg19.bed",package="lassosum")
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- as.data.frame(fread("~/SCHARPHOME/ScharpFile/Grant/BEprogression/RiskPrediction/EUR.height")[,c("FID", "IID", "Height")])
# Read in samples to include in the analysis
target.keep <- fread("EUR.valid.sample")[,c("FID", "IID")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Number of sample in base
size <- 253288
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Read in the LD blocks
ld <- fread(ld.file)
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
             n = size,
             sign = log(ss$OR)
)


######################
### PRS and LDpred ###
######################

rm(list=ls())
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filteredMAF_20Feb2015_t.RData")

library(readxl)
sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
sampletable <- data.frame(sampletable)
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA

allsamples=intersect(fam$V2,sampletable$localid)
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
rownames(sampletable)=sampletable$localid



load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_PRS.RData")

load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/LDpred_genotyped_PRS.RData")

names(BEprs) <- c("localid","BE.prs")
names(EAprs) <- c("localid","EA.prs")
names(BEEAprs) <- c("localid","BEEA.prs")

sampletable <- merge(sampletable,BEprs,by="localid")
sampletable <- merge(sampletable,EAprs,by="localid")
sampletable <- merge(sampletable,BEEAprs,by="localid")





fit1 <- glm(I(phenoEA_bca==2)~EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs*recurrent_HB_RF,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs*cig_smk_ever,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs*nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)


fit1 <- glm(I(phenoEA_bca==2)~EA.prs*sex,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30& sampletable$recurrent_HB_RF==1,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30& sampletable$recurrent_HB_RF==0,],y=T)
summary(fit1)


fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30 & sampletable$sex==2,],y=T)
summary(fit1)

fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30 & sampletable$sex==1,],y=T)
summary(fit1)


fit2 <- glm(I(phenoEA_bca==2)~age+sex+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable[sampletable$site<=30&sampletable$recurrent_HB_RF==1,],y=T)
summary(fit2)
roc(fit2$y~fit2$fitted.values)$auc

fit2 <- glm(I(phenoEA_bca==2)~age+sex+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable[sampletable$site<=30&sampletable$recurrent_HB_RF==0,],y=T)
summary(fit2)
roc(fit2$y~fit2$fitted.values)$auc


fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit2)
fit1 <- glm(I(phenoEA_bca==2)~EA.prs+age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit2)

summary(glm(recurrent_HB_RF~EA.prs,family=binomial,data=sampletable[sampletable$phenoEA_bca==2,]))
summary(glm(recurrent_HB_RF~EA.prs,family=binomial,data=sampletable))
summary(glm(recurrent_HB_RF~BE.prs,family=binomial,data=sampletable))

boxplot(sampletable$BE.prs~sampletable$recurrent_HB_RF)

summary(glm(cig_smk_ever~EA.prs,family=binomial,data=sampletable))



fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit2)


fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit2)

fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit2)

fit2 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit2)


fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit3)


fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc

fit3 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs+EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc


fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc

fit3 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc


fit4 <- glm(I(phenoEA_bca==2)~age+sex*EA.prs+recurrent_HB_RF*EA.prs+bmi_recent_healthy*EA.prs+cig_smk_ever*EA.prs+nsaid_ever*EA.prs,family=binomial,data=sampletable[sampletable$site<=30,],y=T)
summary(fit4)

library(pROC)
png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_ROC.png")
roc1 <- roc(fit1$y~fit1$fitted.values)
plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.5,cex.lab=1.8)
roc2 <- roc(fit3$y~fit3$fitted.values)
plot(roc2,ylim=c(0,1),print.auc=F,col=4,add=T)

roc3 <- roc(fit2$y~fit2$fitted.values)
plot(roc3,ylim=c(0,1),print.auc=F,col=4,add=T)

text(0.4,0.5,"PRS: AUC=0.55",col=3,cex=1.65)
text(0.7,0.92,"PRS+Environment: AUC=0.75",col=4,cex=1.65)
dev.off()



fit3 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc
fit3 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc


fit3 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit3)

fit1 <- glm(I(phenoEA_bca==2)~age+EA.prs+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit1)

fit2 <- glm(I(phenoEA_bca==2)~age+EA.prs+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit2)


mean(sampletable$EA.prs[!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0 & sampletable$sex==1]>q2)


mean(sampletable$EA.prs[!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& sampletable$sex==2 ]>q2)


mean(sampletable$EA.prs[!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==2 ]>q2)

mean(sampletable$EA.prs[!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==1 ]>q2)


mean(sampletable$EA.prs[!is.na(sampletable$bmi_recent_healthy)& sampletable$bmi_recent_healthy<30&!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==2 ]>q2)

mean(sampletable$EA.prs[!is.na(sampletable$bmi_recent_healthy)& sampletable$bmi_recent_healthy<30&!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==1 ]>q2)

sum(sampletable$EA.prs[!is.na(sampletable$bmi_recent_healthy)& sampletable$bmi_recent_healthy<30&!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==2 ]>q2)

sum(sampletable$EA.prs[!is.na(sampletable$bmi_recent_healthy)& sampletable$bmi_recent_healthy<30&!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0&sampletable$sex==1 ]>q2)

q1 <- quantile(sampletable$EA.prs[sampletable$EA.prs],0.2)
q2 <- quantile(sampletable$EA.prs[sampletable$EA.prs],0.8)
fit1 <- glm(I(phenoEA_bca==2)~I(EA.prs>=q2), family=binomial,data=sampletable[!is.na(sampletable$phenoEA_bca) & (sampletable$EA.prs<=q1|sampletable$EA.prs>=q2),])
summary(fit1)

q1 <- quantile(sampletable$EA.prs[sampletable$EA.prs & sampletable$sex==1],0.2)
q2 <- quantile(sampletable$EA.prs[sampletable$EA.prs & sampletable$sex==1],0.8)
fit1 <- glm(I(phenoEA_bca==2)~I(EA.prs> (-3.96)), family=binomial,data=sampletable[sampletable$sex==1 & !is.na(sampletable$phenoEA_bca) & (sampletable$EA.prs>=q1|sampletable$EA.prs<=q2),])
summary(fit1)


q1 <- quantile(sampletable$EA.prs[sampletable$EA.prs & sampletable$sex==2],0.2)
q2 <- quantile(sampletable$EA.prs[sampletable$EA.prs & sampletable$sex==2],0.8)
fit1 <- glm(I(phenoEA_bca==2)~I(EA.prs>= q2), family=binomial,data=sampletable[sampletable$sex==2 & !is.na(sampletable$phenoEA_bca) & (sampletable$EA.prs>=q2|sampletable$EA.prs<=q1),])
summary(fit1)


fit1 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==1,],y=T)
summary(fit1)
roc(fit1$y~fit1$fitted.values)$auc

fit3 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+bmi_recent_healthy+cig_smk_ever+nsaid_ever+EA.prs+EA.prs*bmi_recent_healthy,family=binomial,data=sampletable[sampletable$site<=30& sampletable$sex==2,],y=T)
summary(fit3)
roc(fit3$y~fit3$fitted.values)$auc

library(pROC)
png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_ROC_men_women.png")
roc1 <- roc(fit1$y~fit1$fitted.values)
plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction: PRS + Environment",cex.main=1.5,cex.lab=1.8)
#roc2 <- roc(fit3$y~fit3$fitted.values)
#plot(roc2,ylim=c(0,1),print.auc=F,col=2,add=T)

roc3 <- roc(fit3$y~fit3$fitted.values)
plot(roc3,ylim=c(0,1),print.auc=F,col=4,add=T)

text(0.6,0.5,"Men: AUC=0.739",col=3,cex=1.5)
text(0.8,0.9,"Women: AUC=0.815",col=4,cex=1.5)

dev.off()

fit <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF*EA.prs+bmi_recent_healthy*EA.prs+cig_smk_ever*EA.prs+nsaid_ever*EA.prs,family=binomial,data=sampletable[sampletable$site<=30,])
summary(fit)


fit <- glm(I(phenoEA_bca==2)~EA.prs*recurrent_HB_RF,family=binomial,data=sampletable)
summary(fit)

mprs <- median(sampletable$EA.prs)
fit <- glm(I(phenoEA_bca==2)~I(EA.prs>mprs)*recurrent_HB_RF,family=binomial,data=sampletable)
summary(fit)


fit <- glm(I(phenoBE_bca==2)~BE.prs*recurrent_HB_RF,family=binomial,data=sampletable)
summary(fit)
mprs <- median(sampletable$BE.prs)
fit <- glm(I(phenoBE_bca==2)~I(BE.prs>mprs)*recurrent_HB_RF,family=binomial,data=sampletable)
summary(fit)



fit <- glm(I(phenoEABE_bca==2)~BEEA.prs*recurrent_HB_RF,family=binomial,data=sampletable)
summary(fit)


fit <- glm(I(phenoBE_bca==2)~BE.prs*cig_smk_ever,family=binomial,data=sampletable)
summary(fit)


mprs <- median(sampletable$EA.prs)
fit <- glm(I(phenoEA_bca==2)~cig_smk_ever,family=binomial,data=sampletable)
summary(fit)
fit <- glm(I(phenoEA_bca==2)~I(EA.prs>mprs)*cig_smk_ever,family=binomial,data=sampletable)
summary(fit)
fit <- glm(I(phenoEA_bca==2)~I(EA.prs>mprs)+cig_smk_ever,family=binomial,data=sampletable)
summary(fit)

fit <- glm(I(phenoEA_bca==2)~EA.prs*bmi_recent_healthy,family=binomial,data=sampletable)
summary(fit)
fit <- glm(I(phenoEA_bca==2)~I(EA.prs>mprs)*bmi_recent_healthy,family=binomial,data=sampletable)
summary(fit)
fit <- glm(I(phenoEA_bca==2)~I(EA.prs>mprs)*I(bmi_recent_healthy>30),family=binomial,data=sampletable)
summary(fit)

fit <- glm(I(phenoEABE_bca==2)~BEEA.prs*bmi_recent_healthy,family=binomial,data=sampletable)
summary(fit)

fit <- glm(I(phenoBE_bca==2)~BE.prs*bmi_recent_healthy,family=binomial,data=sampletable)
summary(fit)
