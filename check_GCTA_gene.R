#check GCTA heritability

#organ="junction" #275
#organ="adipose" #393
#organ="muscularis" #385
#organ="stomach" #260
#organ="mucosa" #411
#organ="blood" #558

#prefix=paste0("GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction")
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor/") ##
#save result
gcta_opt="1grm"
heritfile=paste0(outfolder,"heritability_",gcta_opt,".txt")
heritability=read.table(heritfile,header = T)
colnames(heritability)=c("H","SE","pvalue")
#saved r2
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
all(rownames(heritability) %in% rownames(res_min))
idx=match(rownames(heritability),rownames(res_min))
heritability$r2=res_min$r2[idx]
dim(heritability)
plot(heritability$r2,heritability$H)
cor(heritability$r2,heritability$H,method = "spearman",use="pairwise.complete.obs") #0.4637717
abline(0,1,col="red")
idx=which(heritability$pvalue<0.05)
plot(heritability$r2[idx],heritability$H[idx])
cor(heritability$r2[idx],heritability$H[idx],method = "spearman",use="pairwise.complete.obs") #0.4637717
abline(0,1,col="red")
#plot(-log10(heritability$r2),-log10(heritability$H))
quantile(heritability$H)
# 0%       25%       50%       75%      100% 
# 0.0000010 0.0001105 0.0242630 0.0725275 0.9934970 
quantile(heritability$pvalue)
# 0%       25%       50%       75%      100% 
# 0.0000000 0.0026219 0.1280400 0.4968650 0.5000000
plot(-log10(heritability$H),-log10(heritability$pvalue))
sum(heritability$pvalue<0.05) #6306
sum(heritability$H>0.05) #5279
table(heritability$H>heritability$r2)
table(heritability$pvalue<0.05,heritability$r2>0.05)
#       FALSE TRUE
# FALSE  5848 3237
# TRUE   2158 4134
table(heritability$H>0.05,heritability$r2>0.05)
#       FALSE TRUE
# FALSE  6644 3467
# TRUE   1362 3904
#check diff of H and r2
idx1=which(heritability$H>0.85 &heritability$r2<0.05)
rownames(heritability)[idx1] #[1] "IFI44"   "C7orf61"
heritability[idx1,]
idx1=which(rownames(heritability) %in% c("IFI44","C7orf61"))

organs=c("junction","mucosa","stomach","muscularis","adipose","blood")
organs=c("junction","stomach")
gcta_opt="1grm"
allgenes=NULL
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor/")
  #heritfile=paste0(outfolder,"heritability_",gcta_opt,".txt")
  heritfile=paste0(outfolder,"heritability_nopeer_",gcta_opt,".txt")
  heritability=read.table(heritfile,header = T)
  colnames(heritability)=c("H","SE","pvalue")
  allgenes=unique(c(allgenes,rownames(heritability)))
}
allheritability_H=allheritability_P=data.frame(matrix(NA,nrow=length(allgenes),ncol=length(organs)),row.names = allgenes)
colnames(allheritability_H)=colnames(allheritability_P)=organs
allheritability_R2=allheritability_H
for (i in 1:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor/")
  #heritfile=paste0(outfolder,"heritability_",gcta_opt,".txt")
  heritfile=paste0(outfolder,"heritability_nopeer_",gcta_opt,".txt")
  heritability=read.table(heritfile,header = T)
  colnames(heritability)=c("H","SE","pvalue")
  idx=match(rownames(heritability),rownames(allheritability_H))
  allheritability_H[idx,i]=heritability$H
  allheritability_P[idx,i]=heritability$pvalue
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes,]
  res_min_code=res_min_code[rownames(res_min_code) %in% rownames(allheritability_R2),]
  idx=match(rownames(res_min_code),rownames(allheritability_R2))
  allheritability_R2[idx,i]=res_min_code$r2
}
plot(allheritability_H$junction,allheritability_H$stomach,xlab="Junction",ylab="Stomach",cex.axis=1.3,cex.lab=1.3)
abline(0,1,col="red")
cor(allheritability_H$junction,allheritability_H$stomach,use="pairwise.complete.obs") #0.597


dim(allheritability_H) #16640
sum(complete.cases(allheritability_H)) #13555
plot(allheritability_H$blood,allheritability_H$mucosa)
abline(0,1,col="red")
plot(allheritability_H$blood,allheritability_H$junction,ylim=c(0,1))
abline(0,1,col="red")
plot(allheritability_H$blood,allheritability_H$stomach)
abline(0,1,col="red")
plot(-log10(allheritability_P$mucosa),-log10(allheritability_P$blood))

genenames=c("IL2RB","COX7A2","FILIP1","HSP90AA1","FOXF1","LDAH","ISYNA1","UBAC1")
idx=match(genenames,rownames(allheritability_H))
allheritability_H[idx,]
#          junction   mucosa  stomach muscularis  adipose    blood
# IL2RB    0.000002 0.054629 0.000003   0.056786 0.000002 0.000004
# COX7A2   0.007528 0.003895 0.000008   0.000006 0.000007 0.003467
# FILIP1   0.000003 0.004869 0.000011   0.000003 0.063193       NA
# HSP90AA1 0.022798 0.015428 0.000002   0.023563 0.016254 0.010367
# FOXF1    0.033184 0.000005 0.000007   0.020441 0.042243       NA
# LDAH     0.345139 0.413286 0.331670   0.410278 0.282753 0.107553
# ISYNA1   0.000002 0.000002 0.000003   0.093282 0.018407 0.033009
# UBAC1    0.158246 0.238354 0.136554   0.244321 0.193690 0.011143
allheritability_R2[idx,]
#          junction      mucosa      stomach  muscularis     adipose       blood
# IL2RB    0.13664813 0.070057700 0.1129379661 0.117274831 0.066226206 0.031692339
# COX7A2   0.02723575 0.008945583 0.0293638067 0.004181465 0.006082777 0.008985393
# FILIP1   0.01114079 0.037243597 0.0329901372 0.015860827 0.037703926          NA
# HSP90AA1 0.03705637 0.002884609 0.0002164485 0.025763341 0.112081995 0.070961699
# FOXF1    0.17331328 0.075717084 0.1271439192 0.219856238 0.129063449          NA
# LDAH     0.23405357 0.185961963 0.2507573263 0.177202302 0.156241681 0.092231762
# ISYNA1   0.14106079 0.024530429 0.0123048151 0.098983742 0.040283899 0.032935734
# UBAC1    0.11455725 0.099458363 0.0582680256 0.100111848 0.128432142 0.040458952
allheritability_P[idx,]
#          junction     mucosa    stomach muscularis    adipose      blood
# IL2RB    5.0000e-01 8.5389e-02 5.0000e-01 7.0793e-02 5.0000e-01 5.0000e-01
# COX7A2   2.8462e-01 3.2178e-01 5.0000e-01 5.0000e-01 5.0000e-01 3.8062e-01
# FILIP1   5.0000e-01 2.9250e-01 5.0000e-01 5.0000e-01 1.1801e-03         NA
# HSP90AA1 1.8175e-01 1.5454e-01 5.0000e-01 1.2687e-01 2.0433e-01 1.5635e-01
# FOXF1    3.0323e-01 5.0000e-01 5.0000e-01 3.2775e-01 1.2802e-01         NA
# LDAH     3.4154e-08 1.1102e-14 4.9792e-06 2.2204e-15 1.4526e-11 3.4326e-05
# ISYNA1   5.0000e-01 5.0000e-01 5.0000e-01 1.1306e-02 1.8620e-01 6.9570e-03
# UBAC1    2.8484e-04 3.8620e-12 6.3955e-03 8.8976e-12 6.1853e-09 3.4713e-01

#check eQTL shared by stomach and junction
organs=c("junction","stomach")
organ=organs[1]
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar/")
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
res_min_code1=res_min[rownames(res_min) %in% proteingenes &res_min$numselectedsnp>0,] #15175
res_min_code1=res_min_code1[res_min_code1$r2>0.01,]
allgenes=unique(rownames(res_min_code1))
for (i in 2:length(organs))
{
  organ=organs[i]
  outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar/")
  load(paste0(outfolder,"/preidiction_michigan_model.RData"))
  res_min_code=res_min[rownames(res_min) %in% proteingenes & res_min$numselectedsnp>0,] #15630
  res_min_code=res_min_code[res_min_code$r2>0.01,]
  allgenes=intersect(allgenes,rownames(res_min_code))
}

eqtlres=data.frame(junction=rep(NA,length(allgenes)),stomach=NA,all=NA,overlap="NA")
rownames(eqtlres)=allgenes
for (i in 1:length(allgenes))
{
  #if (i %%1000==0) cat(i,'..')
  idx1=which(rownames(res_min_code1)==allgenes[i])
  junctioneqtl=unique(unlist(strsplit(res_min_code1$selectedsnps[idx1],"|",fixed=T)))
  eqtlres$junction[i]=length(junctioneqtl)
  idx2=which(rownames(res_min_code)==allgenes[i])
  stomacheqtl=unique(unlist(strsplit(res_min_code$selectedsnps[idx2],"|",fixed=T)))
  eqtlres$stomach[i]=length(stomacheqtl)
  eqtlres$all[i]=length(unique(c(junctioneqtl,stomacheqtl)))
  eqtlres$overlap[i]=length(intersect(junctioneqtl,stomacheqtl))
}
eqtlres$overlap=as.integer(eqtlres$overlap)
summary(eqtlres$junction)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   25.00   34.00   36.57   46.00  147.00
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    7.00   13.00   15.74   21.00  102.00 
summary(eqtlres$stomach)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   25.00   34.00   36.69   46.00  149.00 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00    7.00   12.00   14.77   20.00   95.00
summary(eqtlres$overlap)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   1.000   3.000   3.086   4.000  26.000 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   1.000   1.066   1.000  32.000 
nrow(res_min_code1) #14252 #6656
nrow(res_min_code) #14571 #6161
nrow(eqtlres) #13232 #3359
summary(eqtlres$overlap)

eqtlres$prop=eqtlres$overlap/eqtlres$stomach
summary(eqtlres$prop)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.03774 0.07500 0.08614 0.12121 0.75000 
hist(eqtlres$prop,col="blue",xlab="Proportion of overlap eQTL",cex.axis=1.3,cex.lab=1.3,main="",probability = T)
