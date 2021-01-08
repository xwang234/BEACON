load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Dong23SNPs.RData")
library(CGEN)

PC <- PC[row.names(PC)%in% row.names(dat),]
PC <- PC[match(row.names(dat),row.names(PC)),]

dat <- cbind(dat,PC[,1:4])
for (i in 1:ncol(dat)) dat[which(dat[,i]==-9),i]=NA

#response
dat$BEEA=NA
dat$BEEA[which(dat$phenoBE_bc==1 | dat$phenoEA_bc==1)]=0
dat$BEEA[which(dat$phenoBE_bc==2 | dat$phenoEA_bc==2)]=1
table(dat$BEEA,useNA="ifany")
dat$BE=I(dat$phenoBE_bc==2)
dat$EA=I(dat$phenoEA_bc==2)

dat$bmi_cat <- ifelse(dat$bmi_recent_healthy<25,0,1)
dat$bmi_cat <- ifelse(dat$bmi_recent_healthy>=30,2,dat$bmi_cat)
table(dat$bmi_cat,useNA="ifany")
dat$sex=as.factor(dat$sex)

exposures=c("bmi_cat","recurrent_HB_RF","cig_smk_ever","recurrent_reflux","nsaid_ever")

#get p-values for all snps-exposure
pvalues=data.frame(matrix(NA,nrow=23,length(exposures)))
colnames(pvalues)=exposures
rownames(pvalues)=colnames(dat)[1:23]
pvaluesBE=pvaluesEA=pvalues
for (i in 1:nrow(pvalues))
{
  cat(i,'..')
  for (j in 1:length(exposures))
  {
    tmp1=additive.test(data=dat,response.var = "BEEA",snp.var =rownames(pvalues)[i],main.vars = c("age","sex"),exposure.var = exposures[j] )
    pvalues[i,j]=tmp1$pval.add #use pval.add
    tmp2=additive.test(data=dat,response.var = "BE",snp.var =rownames(pvalues)[i],main.vars = c("age","sex"),exposure.var = exposures[j] )
    pvaluesBE[i,j]=tmp2$pval.add
    if (min(table(dat[!is.na(dat$EA),rownames(pvalues)[i]],dat[!is.na(dat$EA),exposures[j]]))>=5) #count number less than 5, chi-square not applicable, for i=9 and j=4
    {
      tmp3=additive.test(data=dat,response.var = "EA",snp.var =rownames(pvalues)[i],main.vars = c("age","sex"),exposure.var = exposures[j] )
      pvaluesEA[i,j]=tmp3$pval.add
    }
    
  }
}


###########################################################

### now examine 23 SNPs interacting with 4 risk factors ###

###########################################################



library(HardyWeinberg)



dat$bmi_cat <- ifelse(dat$bmi_recent_healthy<25,0,1)

dat$bmi_cat <- ifelse(dat$bmi_recent_healthy>=30,2,dat$bmi_cat)

dat$bmi_cat[dat$bmi_recent_healthy== (-9)] <- (-9)



###

### first test the HW and swap to minor allele ###

###



phw <- rep(NA,23)



for (i in 1:23){
  
  tab1 <- table(as.numeric(dat[,i]))
  
  print(tab1)
  
  
  
  if (which.min(tab1)==1) {
    
    dat[,i] <- 2- as.numeric(dat[,i])
    
  }  
  
  
  
  HW.test <- HWChisq(table(dat[,i]))
  
  phw[i]<- HW.test$pval
  
  
  
  #p <- (2*tab1[3]+tab1[2])/sum(tab1*2)
  
  #print(2*sum(tab1)*p*(1-p))
  
  
  
} 



####################################################

### first try the individual SNPs no interaction ###

####################################################



SNPnames <- names(dat)[1:23]









out <- matrix(0,23,6)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) ,i])*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) ,])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) ,i])*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) ,])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) ,i])*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[9,4]
  
  out[i,2] <- summary(fit2)$coef[9,4]
  
  out[i,3] <- summary(fit3)$coef[9,4]
  
  
  
  fit1 <- glm(I(phenoBE_bc==2)~I(dat[dat$phenoBE_bc!= (-9) ,i]>0)*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) ,])
  
  fit2 <- glm(I(phenoEA_bc==2)~I(dat[dat$phenoEA_bc!= (-9) ,i]>0)*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) ,])
  
  fit3 <- glm(I(phenoEABE_bc==2)~I(dat[dat$phenoEABE_bc!= (-9) ,i]>0)*factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9),])
  
  
  
  out[i,4] <- summary(fit1)$coef[9,4]
  
  out[i,5] <- summary(fit2)$coef[9,4]
  
  out[i,6] <- summary(fit3)$coef[9,4]
  
  
  
} 



out1 <- out













out <- matrix(0,23,6)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i])*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[10,4]
  
  out[i,2] <- summary(fit2)$coef[10,4]
  
  out[i,3] <- summary(fit3)$coef[10,4]
  
  
  
  
  
  fit1 <- glm(I(phenoBE_bc==2)~I(dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i]>0)*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~I(dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i]>0)*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~I(dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),i]>0)*recurrent_HB_RF+factor(sex)+age + pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$recurrent_HB_RF!=(-9),])
  
  
  
  
  
  out[i,4] <- summary(fit1)$coef[10,4]
  
  out[i,5] <- summary(fit2)$coef[10,4]
  
  out[i,6] <- summary(fit3)$coef[10,4]
  
  
  
} 





out2 <- out















out <- matrix(0,23,6)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) & dat$cig_smk_ever!=(-9),i])*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) & dat$cig_smk_ever!=(-9),i])*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) & dat$cig_smk_ever!=(-9),i])*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[10,4]
  
  out[i,2] <- summary(fit2)$coef[10,4]
  
  out[i,3] <- summary(fit3)$coef[10,4]
  
  
  
  
  
  fit1 <- glm(I(phenoBE_bc==2)~I(dat[dat$phenoBE_bc!= (-9) & dat$cig_smk_ever!=(-9),i]>0)*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~I(dat[dat$phenoEA_bc!= (-9) & dat$cig_smk_ever!=(-9),i]>0)*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~I(dat[dat$phenoEABE_bc!= (-9) & dat$cig_smk_ever!=(-9),i]>0)*cig_smk_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$cig_smk_ever!=(-9),])
  
  
  
  out[i,4] <- summary(fit1)$coef[10,4]
  
  out[i,5] <- summary(fit2)$coef[10,4]
  
  out[i,6] <- summary(fit3)$coef[10,4]
  
  
  
} 





out3 <- out











out <- matrix(0,23,6)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) & dat$nsaid_ever!=(-9),i])*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) & dat$nsaid_ever!=(-9),i])*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) & dat$nsaid_ever!=(-9),i])*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[10,4]
  
  out[i,2] <- summary(fit2)$coef[10,4]
  
  out[i,3] <- summary(fit3)$coef[10,4]
  
  
  
  
  
  fit1 <- glm(I(phenoBE_bc==2)~I(dat[dat$phenoBE_bc!= (-9) & dat$nsaid_ever!=(-9),i]>0)*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~I(dat[dat$phenoEA_bc!= (-9) & dat$nsaid_ever!=(-9),i]>0)*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~I(dat[dat$phenoEABE_bc!= (-9) & dat$nsaid_ever!=(-9),i]>0)*nsaid_ever+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$nsaid_ever!=(-9),])
  
  
  
  
  
  
  
  out[i,4] <- summary(fit1)$coef[10,4]
  
  out[i,5] <- summary(fit2)$coef[10,4]
  
  out[i,6] <- summary(fit3)$coef[10,4]
  
  
  
} 



out4 <- out

















dat$bmi_standard <- dat$bmi_recent_healthy/sqrt(var(dat$bmi_recent_healthy[dat$bmi_recent_healthy != -9 & dat$phenoBE_bc==1]))

dat$bmi_standard[dat$bmi_recent_healthy== -9] <- (-9)



out <- matrix(0,23,6)

for (i in 1:23) {
  
  fit1 <- glm(I(phenoBE_bc==2)~as.numeric(dat[dat$phenoBE_bc!= (-9) & dat$bmi_standard!=(-9),i])*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~as.numeric(dat[dat$phenoEA_bc!= (-9) & dat$bmi_standard!=(-9),i])*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~as.numeric(dat[dat$phenoEABE_bc!= (-9) & dat$bmi_standard!=(-9),i])*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  
  
  out[i,1] <- summary(fit1)$coef[10,4]
  
  out[i,2] <- summary(fit2)$coef[10,4]
  
  out[i,3] <- summary(fit3)$coef[10,4]
  
  
  
  
  
  fit1 <- glm(I(phenoBE_bc==2)~I(dat[dat$phenoBE_bc!= (-9) & dat$bmi_standard!=(-9),i]>0)*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoBE_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  fit2 <- glm(I(phenoEA_bc==2)~I(dat[dat$phenoEA_bc!= (-9) & dat$bmi_standard!=(-9),i]>0)*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEA_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  fit3 <- glm(I(phenoEABE_bc==2)~I(dat[dat$phenoEABE_bc!= (-9) & dat$bmi_standard!=(-9),i]>0)*bmi_standard+factor(sex)+age+ pc1 + pc2 + pc3 + pc4,family=binomial,data=dat[dat$phenoEABE_bc!= (-9) & dat$bmi_standard!=(-9),])
  
  
  
  out[i,4] <- summary(fit1)$coef[10,4]
  
  out[i,5] <- summary(fit2)$coef[10,4]
  
  out[i,6] <- summary(fit3)$coef[10,4]
  
  
  
} 



out5 <- out

