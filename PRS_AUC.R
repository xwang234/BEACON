
  library(readxl)
  sampletable=read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
  sampletable <- data.frame(sampletable)
  for (i in 1:ncol(sampletable)) sampletable[which(sampletable[,i]==-9),i]=NA
  EAfile="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Validation_EA.PRS.txt"
  EAprs=read.table(EAfile,header=T,stringsAsFactors = F)
  names(EAprs) <- c("localid","EA.prs")
  
  sampletable1 <- merge(sampletable,EAprs,by="localid")
  
  #png("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRS_exposure_ROC.png")
  fit1 <- glm(I(phenoEA_bca==2)~EA.prs,family=binomial,data=sampletable1,y=T,x=T)
  #summary(fit1)
  fit2 <- glm(I(phenoEA_bca==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever,family=binomial,data=sampletable1,y=T)
  roc2 <- roc(fit2$y~fit2$fitted.values)
  print(roc2)
  #summary(fit2)
  library(pROC)
  roc1 <- roc(fit1$y~fit1$fitted.values)
  #remove "bmi_recent_healthy"
  dat=sampletable1[,c("phenoEA_bca","age","sex","recurrent_HB_RF","bmi_recent_healthy","cig_smk_ever","nsaid_ever","EA.prs")]
  idx=complete.cases(dat)
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
  plot(1:4,exp(summary(fit1)$coef[2:5,1]),type="n",main="PRS association with EAC",xlab="PRS percentile",xaxt="n",cex.axis=1.8,ylim=c(1.0,3.6),xlim=c(0.5,4.5),ylab="relative risk vs 0-20 group",cex.main=1.5,cex.lab=1.6)
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
  
  