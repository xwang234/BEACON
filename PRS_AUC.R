  rm(list=ls())
  library(readxl)
  sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
  sampletable <- data.frame(sampletable)
  for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
  idx=complete.cases(sampletable1[,c("age","sex","recurrent_HB_RF","cig_smk_ever","nsaid_ever","EA.prs","recent_bmi_healthy")])
  sampletable1 <- sampletable1[idx,]
  
  #EAfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_EA_genotyped.PRS.txt" 
  EAfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA.PRS.txt"
  #EAfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BCA_EA.PRS.txt"
  EAprs=read.table(EAfile,header=T,stringsAsFactors = F)
  names(EAprs) <- c("localid","EA.prs")
  
  sampletable1 <- merge(sampletable,EAprs,by="localid")
  sampletable1 <- sampletable1[!is.na(sampletable1$phenoEA_bca),]
  
  idx=complete.cases(sampletable1[,c("age","sex","recurrent_HB_RF","cig_smk_ever","nsaid_ever","EA.prs")])
  sampletable1 <- sampletable1[idx,]
  
  library(pROC)
  #png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_ROC.png")
  fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable1,y=T,x=T)
  summary(fit1)
  roc1 <- roc(fit1$y~fit1$fitted.values)
  roc1$auc
  
  fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable1,y=T,x=T)
  roc2 <- roc(fit2$y~fit2$fitted.values)
  print(roc2)
  #summary(fit2)
  
  mean(fit3$fitted[sampletable1$phenoEA_bca==2]>fit2$fitted[sampletable1$phenoEA_bca==2])
  mean(fit3$fitted[sampletable1$phenoEA_bca==1]<fit2$fitted[sampletable1$phenoEA_bca==1])
  
  mean(fit3$fitted[sampletable1$phenoEA_bca==2 & sampletable1$recurrent_HB_RF==0]>fit2$fitted[sampletable1$phenoEA_bca==2 & sampletable1$recurrent_HB_RF==0])
  mean(fit3$fitted[sampletable1$phenoEA_bca==2 & sampletable1$cig_smk_ever==0]>fit2$fitted[sampletable1$phenoEA_bca==2 & sampletable1$cig_smk_ever==0])
  mean(fit3$fitted[sampletable1$phenoEA_bca==2 & sampletable1$sex==2]>fit2$fitted[sampletable1$phenoEA_bca==2 & sampletable1$sex==2])
  mean(fit3$fitted[sampletable1$phenoEA_bca==1 & sampletable1$sex==2]<fit2$fitted[sampletable1$phenoEA_bca==1 & sampletable1$sex==2])
  
  fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs+ev1_bca+ev2_bca+ev3_bca+ev4_bca,family=binomial,data=sampletable1,y=T,x=T)
  fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable1,y=T,x=T)
  roc3 <- roc(fit3$y~fit3$fitted.values)
  print(roc3)
  
  fit4 <- glm(I(phenoEA_bca==2)~fit2$fitted*EA.prs,family=binomial,data=sampletable1,y=T,x=T)
  
  ## current guideline 
  crisk <- 1*(sampletable1$recurrent_HB_RF==1 & sampletable1$cig_smk_ever==1 & sampletable1$bmi_recent_healthy>25)
  
  q1 <- quantile(fit2$fitted,0.2)
  q2 <- quantile(fit2$fitted,0.4)
  q3 <- quantile(fit2$fitted,0.6)
  q4 <- quantile(fit2$fitted,0.8)
  
  qlabel <- rep(1,length(sampletable1$EA.prs))
  qlabel <- ifelse(fit2$fitted>q1 & fit2$fitted,2,qlabel)
  qlabel <- ifelse(fit2$fitted>q2 & fit2$fitted<q3,3,qlabel)
  qlabel <- ifelse(fit2$fitted>q3 & fit2$fitted<q4,4,qlabel)
  qlabel <- ifelse(fit2$fitted>q4,5,qlabel)
  
  fit4 <- glm(I(phenoEA_bca==2)~EA.prs*factor(qlabel),family=binomial,data=sampletable1,y=T,x=T)
  
  
  q1 <- quantile(sampletable1$EA.prs[sampletable1$phenoEA_bca==1],0.2)
  q2 <- quantile(sampletable1$EA.prs[sampletable1$phenoEA_bca==1],0.4)
  q3 <- quantile(sampletable1$EA.prs[sampletable1$phenoEA_bca==1],0.6)
  q4 <- quantile(sampletable1$EA.prs[sampletable1$phenoEA_bca==1],0.8)
  
  qlabel <- rep(1,length(sampletable1$EA.prs))
  qlabel <- ifelse(sampletable1$EA.prs>q1 & sampletable1$EA.prs<q2,2,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q2 & sampletable1$EA.prs<q3,3,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q3 & sampletable1$EA.prs<q4,4,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q4,5,qlabel)
  
  fit4 <- glm(I(phenoEA_bca==2)~factor(qlabel),family=binomial,data=sampletable1,y=T,x=T)
  
  
  
  
  
  
  ## pick age 50 male with GERD, smoker and obese localid=201242
  fit3$x[fit3$x[,2]==50 & fit3$x[,3]==1 & fit3$x[,4]==1 &  fit3$x[,5]==1,]
  fit3$y[fit3$x[,2]==50 & fit3$x[,3]==1 & fit3$x[,4]==1 &  fit3$x[,5]==1]
  fit3$fitted[fit3$x[,2]==50 & fit3$x[,3]==1 & fit3$x[,4]==1 &  fit3$x[,5]==1]
  
  
  fit3$fitted[fit3$x[,2]==50 & fit3$x[,3]==1 & fit3$x[,4]==0]
  
  riskcut1 <- fit2$fitted[fit2$x[,2]==50 & fit2$x[,3]==1 & fit2$x[,4]==1 &  fit2$x[,5]==1][1]
  
  
  
  riskcut1 <- median(fit2$fitted[fit2$x[,2]==50 & fit2$x[,3]==1 & fit2$x[,4]==1])
  riskcut0 <- median(fit2$fitted[fit2$x[,2]==50 & fit2$x[,3]==1 & fit2$x[,4]==0])
  riskgroup3 <- ifelse(fit3$fitted>riskcut1,1,0) + ifelse(fit3$fitted>riskcut0,1,0)
  riskgroup2 <- ifelse(fit2$fitted>riskcut1,1,0) + ifelse(fit2$fitted>riskcut0,1,0)
  
  table(fit3$y[fit3$x[,4]==0],riskgroup3[fit3$x[,4]==0])
  table(fit2$y[fit3$x[,4]==0],riskgroup2[fit3$x[,4]==0])
  
  sum(riskgroup3[fit3$y==1&fit3$x[,4]==0]>riskgroup2[fit3$y==1&fit3$x[,4]==0])
  sum(riskgroup3[fit3$y==1&fit3$x[,4]==0]<riskgroup2[fit3$y==1&fit3$x[,4]==0])
  
  
  sum(riskgroup3[fit3$y==0]>riskgroup2[fit3$y==0])
  sum(riskgroup3[fit3$y==0]<riskgroup2[fit3$y==0])
  
  
  
  table(fit3$y[fit3$x[,2]<50],fit3$fitted[fit3$x[,2]<50]>riskcut)
  
  table(fit3$y[fit3$x[,3]==1],fit3$fitted[fit3$x[,3]==1]>riskcut)
  table(fit2$y[fit2$x[,3]==1],fit2$fitted[fit2$x[,3]==1]>riskcut)
  
  table(fit3$y[fit3$x[,3]==2],fit3$fitted[fit3$x[,3]==2]>riskcut)
  table(fit2$y[fit2$x[,3]==2],fit2$fitted[fit2$x[,3]==2]>riskcut)
  
  
  fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever,family=binomial,data=sampletable1,y=T)
  roc3 <- roc(fit3$y~fit3$fitted.values)
  print(roc3)
  
  
  
  fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+EA.prs,family=binomial,data=sampletable1,y=T)
  roc3 <- roc(fit3$y~fit3$fitted.values)
  print(roc3)
  
  sampletable1$EAoutcome <- factor(1*(sampletable1$phenoEA_bca==2))
  library(rpart)
  library(randomForest)
  
  dat=sampletable1[,c("EAoutcome","age","sex","recurrent_HB_RF","cig_smk_ever","nsaid_ever","EA.prs")]
  idx=complete.cases(dat)
  
  rfit <- randomForest(EAoutcome~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs, data=sampletable1[idx,],ntree=2000,importance=TRUE,keep.forest=TRUE,keep.inbag=TRUE,na.action = na.roughfix)
  prediction_prob <- predict(rfit,type="prob")[,2]
  roc4 <- roc(sampletable1$EAoutcome[idx]~prediction_prob)
  print(roc4)
  
  
  library(e1071)
  tune.out=tune(svm ,as.matrix(dat[idx,-1]),dat$EAoutcome[idx],scale=T,  kernel ="radial", ranges =list(cost=c(0.1 ,1 ,10,20),gamma=c(0.01,0.1,0.5,1) ))
  tune.out=tune(svm ,as.matrix(dat[idx,-1]),dat$EAoutcome[idx],scale=T,  kernel ="linear", ranges =list(cost=c(0.1 ,1 ,10,20),gamma=c(0.01,0.1,0.5,1) ))
  
  summary(tune.out)
  #bestmod =tune.out$best.model
  #summary(bestmod)
  
  svmfit <- svm(EAoutcome~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs,data=dat[idx,],scale=T,kernel ="radial", cost=10,gamma = 0.01,probability=TRUE)
  svmfit <- svm(as.matrix(dat[idx,-1]),dat$EAoutcome[idx],scale=F,kernel ="linear", cost=1,gamma = 0.01,probability=TRUE)
  
  pred <- fitted(svmfit)
  table(as.numeric(dat$EAoutcome[idx]),pred)
  
  prediction_prob <- predict(svmfit,as.matrix(dat[idx,-1]), probability = TRUE)
  pred_prob <- attr(prediction_prob, "probabilities")[,2]
  roc5 <- roc(dat$EAoutcome[idx]~pred_prob)
  print(roc5)
  
  
  #remove "bmi_recent_healthy"
 
  if (sum(idx)>0)
  {
    fit3 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable1,y=T)
    #summary(fit3)
    png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_PRS_ROC.png")
    plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.2,cex.lab=1.2)
    roc2 <- roc(fit3$y~fit3$fitted.values)
    plot(roc2,ylim=c(0,1),print.auc=F,col=4,add=T)
    text(0.4,0.5,paste0("PRS: AUC=",round(roc1$auc,2)),col=3,cex=1.2)
    text(0.7,0.98,paste0("PRS+Environment: AUC=",round(roc2$auc,2)),col=4,cex=1.2)
    dev.off()
  }else
  {
    plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.2,cex.lab=1.2)
    text(0.4,0.5,paste0("PRS: AUC=",round(roc1$auc,3)),col=3,cex=1.2)
  }
  

  fit3 <- glm(I(phenoEA_bca==2)~age+recurrent_HB_RF+cig_smk_ever+nsaid_ever+EA.prs,family=binomial,data=sampletable1,y=T)
  #summary(fit3)
  png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_PRS_ROC.png")
  plot(roc1,ylim=c(0,1),print.auc=F,col=3,main="EAC risk prediction",cex.main=1.8,cex.lab=1.2)
  roc2 <- roc(fit3$y~fit3$fitted.values)
  plot(roc2,ylim=c(0,1),print.auc=F,col=4,add=T)
  text(0.4,0.55,paste0("PRS: AUC=",round(roc1$auc,2)),col=3,cex=1.4)
  text(0.7,0.95,paste0("PRS+Environment: AUC=","0.82"),col=4,cex=1.4)
  dev.off()  
  
  png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/DNN.png")
  dev.off()

  
  summary(glm(I(phenoEA_bca-1)~sex,family="binomial",data=sampletable1))
  summary(glm(I(phenoEA_bca-1)~age,family="binomial",data=sampletable1))
  summary(glm(I(phenoEA_bca-1)~recurrent_HB_RF,family="binomial",data=sampletable1))
  summary(glm(I(phenoEA_bca-1)~cig_smk_ever,family="binomial",data=sampletable1))
  summary(glm(I(phenoEA_bca-1)~bmi_recent_healthy,family="binomial",data=sampletable1))
  summary(glm(I(phenoEA_bca-1)~nsaid_ever,family="binomial",data=sampletable1))
  

  
  
  q1 <- quantile(sampletable1$EA.prs,0.2)
  q2 <- quantile(sampletable1$EA.prs,0.4)
  q3 <- quantile(sampletable1$EA.prs,0.6)
  q4 <- quantile(sampletable1$EA.prs,0.8)
  
  qlabel <- rep(1,length(sampletable1$EA.prs))
  qlabel <- ifelse(sampletable1$EA.prs>q1 & sampletable1$EA.prs<q2,2,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q2 & sampletable1$EA.prs<q3,3,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q3 & sampletable1$EA.prs<q4,4,qlabel)
  qlabel <- ifelse(sampletable1$EA.prs>q4,5,qlabel)
  
  
  fit1 <- glm(I(phenoEA_bca==2)~factor(qlabel), family=binomial,data=sampletable1[!is.na(sampletable1$phenoEA_bca),])
  summary(fit1)
  
  exp(summary(fit1)$coef[2,1])
  exp(summary(fit1)$coef[2,1]+1.96*summary(fit1)$coef[2,2])
  exp(summary(fit1)$coef[2,1]-1.96*summary(fit1)$coef[2,2])
  
  exp(summary(fit1)$coef[3,1])
  exp(summary(fit1)$coef[3,1]+1.96*summary(fit1)$coef[3,2])
  exp(summary(fit1)$coef[3,1]-1.96*summary(fit1)$coef[3,2])
  
  exp(summary(fit1)$coef[4,1])
  exp(summary(fit1)$coef[4,1]+1.96*summary(fit1)$coef[4,2])
  exp(summary(fit1)$coef[4,1]-1.96*summary(fit1)$coef[4,2])
  
  exp(summary(fit1)$coef[5,1])
  exp(summary(fit1)$coef[5,1]+1.96*summary(fit1)$coef[5,2])
  exp(summary(fit1)$coef[5,1]-1.96*summary(fit1)$coef[5,2])
  
  
  png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_EAC_RR.png")
  plot(1:4,exp(summary(fit1)$coef[2:5,1]),type="n",main="PRS association with EAC risk",xlab="PRS percentile",xaxt="n",cex.axis=1.8,ylim=c(1.0,4.2),xlim=c(0.5,4.5),ylab="relative risk vs 0-20 group",cex.main=1.5,cex.lab=1.6)
  points(1:4,exp(summary(fit1)$coef[2:5,1]),pch=19,cex=2.0,color=2)
  segments(1,exp(summary(fit1)$coef[2,1]+1.96*summary(fit1)$coef[2,2]),1,exp(summary(fit1)$coef[2,1]-1.96*summary(fit1)$coef[2,2]))
  segments(1-0.1,exp(summary(fit1)$coef[2,1]+1.96*summary(fit1)$coef[2,2]),1+0.1,exp(summary(fit1)$coef[2,1]+1.96*summary(fit1)$coef[2,2]))
  segments(1-0.1,exp(summary(fit1)$coef[2,1]-1.96*summary(fit1)$coef[2,2]),1+0.1,exp(summary(fit1)$coef[2,1]-1.96*summary(fit1)$coef[2,2]))
  
  segments(2,exp(summary(fit1)$coef[3,1]+1.96*summary(fit1)$coef[3,2]),2,exp(summary(fit1)$coef[3,1]-1.96*summary(fit1)$coef[3,2]))
  segments(2-0.1,exp(summary(fit1)$coef[3,1]+1.96*summary(fit1)$coef[3,2]),2+0.1,exp(summary(fit1)$coef[3,1]+1.96*summary(fit1)$coef[3,2]))
  segments(2-0.1,exp(summary(fit1)$coef[3,1]-1.96*summary(fit1)$coef[3,2]),2+0.1,exp(summary(fit1)$coef[3,1]-1.96*summary(fit1)$coef[3,2]))
  
  segments(3,exp(summary(fit1)$coef[4,1]+1.96*summary(fit1)$coef[4,2]),3,exp(summary(fit1)$coef[4,1]-1.96*summary(fit1)$coef[4,2]))
  segments(3-0.1,exp(summary(fit1)$coef[4,1]+1.96*summary(fit1)$coef[4,2]),3+0.1,exp(summary(fit1)$coef[4,1]+1.96*summary(fit1)$coef[4,2]))
  segments(3-0.1,exp(summary(fit1)$coef[4,1]-1.96*summary(fit1)$coef[4,2]),3+0.1,exp(summary(fit1)$coef[4,1]-1.96*summary(fit1)$coef[4,2]))
  
  segments(4,exp(summary(fit1)$coef[5,1]+1.96*summary(fit1)$coef[5,2]),4,exp(summary(fit1)$coef[5,1]-1.96*summary(fit1)$coef[5,2]))
  segments(4-0.1,exp(summary(fit1)$coef[5,1]+1.96*summary(fit1)$coef[5,2]),4+0.1,exp(summary(fit1)$coef[5,1]+1.96*summary(fit1)$coef[5,2]))
  segments(4-0.1,exp(summary(fit1)$coef[5,1]-1.96*summary(fit1)$coef[5,2]),4+0.1,exp(summary(fit1)$coef[5,1]-1.96*summary(fit1)$coef[5,2]))
  axis(side=1,at=1:4,labels=c("20-40","40-60","60-80","80-100"),cex.axis=1.5)
  dev.off()
  
  
  
  
  
  
  
  
  sum(!is.na(sampletable$bmi_recent_healthy)& sampletable$bmi_recent_healthy<30&!is.na(sampletable$recurrent_HB_RF)& sampletable$recurrent_HB_RF==0& !is.na(sampletable$cig_smk_ever)& sampletable$cig_smk_ever==0)
  sum(!is.na(sampletable$bmi_recent_healthy)&!is.na(sampletable$recurrent_HB_RF)& !is.na(sampletable$cig_smk_ever))
  
  