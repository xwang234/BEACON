#prepare data-----------------------
outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_EAC2"
#predicted geneexp
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
predict_min=predict_min[,3:ncol(predict_min)]
predict_1se=predict_1se[,3:ncol(predict_1se)]
#read clinical table
sampletable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/bc_AT_JD1.csv",header=T,sep=",",stringsAsFactors=F)
sampletable$site[sampletable$site=="NA"]=NA
for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
geneexpsamplenames=strsplit(colnames(predict_1se),"_")
geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
  tmp=geneexpsamplenames[[x]]
  paste0(tmp[2:length(tmp)],collapse = "_")
})
colnames(predict_min)=geneexpsamplenames

#include PCs
readeigenstrat=function(eigfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pca",
                        eigsampfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_10Jan2019.pedind",
                        thesamples=geneexpsamplenames,nskip=16,opt=1)
{
  tmp=read.table(eigfile,skip=nskip,stringsAsFactors = F)
  if (opt==1)
  {
    eigsamples=read.table(eigsampfile,stringsAsFactors = F)
    eigsamples=eigsamples$V2
    idx=match(thesamples,eigsamples)
    tmp1=tmp[idx,]
    tmp1=as.data.frame(t(tmp1))
    colnames(tmp1)=thesamples
  }else #don't need to change sample names
  {
    tmp1=as.data.frame(t(tmp))
    colnames(tmp1)=thesamples
  }
  rownames(tmp1)=paste0("pc",1:nrow(tmp1))
  return(tmp1)
}
eigenstratmatrix=readeigenstrat()
allsamples=sampletable$localid[sampletable$localid %in% colnames(predict_min)]
idx=match(allsamples,colnames(predict_min))
predict_min=predict_min[,idx]
idx=match(allsamples,colnames(eigenstratmatrix))
eigenstratmatrix=eigenstratmatrix[,idx]
idx=match(allsamples,sampletable$localid)
sampletable=sampletable[idx,]
sampletable=cbind.data.frame(sampletable,t(eigenstratmatrix[1:3,]))
sampletable$sex=as.factor(sampletable$sex)

#compute 2df p-value-----------------
gene="DDX49"
idx1=which(sampletable$phenoBE_bc==2) #BE case
idx2=which(sampletable$phenoBE_bc==1) #CO
y=c(rep(1,length(idx1)),rep(0,length(idx2))) #outcome
x=as.numeric(predict_min[which(rownames(predict_min)==gene),c(idx1,idx2)]) #predicted geneexp
#test a environment factor
fit1=glm(y~x+bmi_recent_healthy+x*bmi_recent_healthy+age+sex+pc1+pc2+pc3,data=sampletable[c(idx1,idx2),],family = binomial)
summary(fit1)
# glm(formula = y ~ x + bmi_recent_healthy + x * bmi_recent_healthy + 
#       age + sex + pc1 + pc2 + pc3, family = binomial, data = sampletable[c(idx1, 
#                                                                            idx2), ])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4603  -1.1991   0.8193   1.0593   2.0314  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -2.530419   0.255198  -9.916  < 2e-16 ***
#   x                      0.812235   1.010434   0.804 0.421485    
# bmi_recent_healthy     0.076362   0.006532  11.691  < 2e-16 ***
#   age                    0.008747   0.002547   3.435 0.000593 ***
#   sex2                   0.210999   0.070875   2.977 0.002910 ** 
#   pc1                  -65.635569  11.405126  -5.755 8.67e-09 ***
#   pc2                   12.165662   3.244241   3.750 0.000177 ***
#   pc3                   -6.034269   6.425869  -0.939 0.347701    
# x:bmi_recent_healthy   0.001145   0.036028   0.032 0.974645    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 6706.8  on 4885  degrees of freedom
# Residual deviance: 6417.6  on 4877  degrees of freedom
# (578 observations deleted due to missingness)
# AIC: 6435.6
#p-value of x turns to not significant after plug in environment factor 
fit0=glm(y~bmi_recent_healthy+age+sex+pc1+pc2+pc3,data=sampletable[c(idx1,idx2),],family = binomial) #the smaller model
summary(fit0)
# glm(formula = y ~ bmi_recent_healthy + age + sex + pc1 + pc2 + 
#       pc3, family = binomial, data = sampletable[c(idx1, idx2), 
#                                                  ])
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.4991  -1.2106   0.8329   1.0628   2.0251  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -2.535275   0.253672  -9.994  < 2e-16 ***
#   bmi_recent_healthy   0.076439   0.006477  11.801  < 2e-16 ***
#   age                  0.008537   0.002539   3.363 0.000772 ***
#   sex2                 0.211992   0.070683   2.999 0.002707 ** 
#   pc1                -65.955275  11.386694  -5.792 6.94e-09 ***
#   pc2                 12.267011   3.235469   3.791 0.000150 ***
#   pc3                 -7.581304   6.407705  -1.183 0.236748    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 6706.8  on 4885  degrees of freedom
# Residual deviance: 6443.4  on 4879  degrees of freedom
# (578 observations deleted due to missingness)
# AIC: 6457.4
#2df-pvalue
p=pchisq(summary(fit0)$deviance-summary(fit1)$deviance,df=2,lower.tail=FALSE)
#[1] 2.455199e-06
