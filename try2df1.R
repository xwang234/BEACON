#prepare data-----------------------

#compute 2df p-value-----------------
compute2dfp=function(i=which(rownames(predict_min)=="DDX49"),envariable="bmi_recent_healthy",opt="min")
{
  if (opt=="min")
  {
    predict_mat=predict_min
  }else
  {
    predict_mat=predict_1se
  }
  res=data.frame(BE_p=NA,EA_p=NA,BEA_p=NA,BEEA_p=NA)
  rownames(res)=rownames(predict_mat)[i]
  idx1=which(sampletable$phenoBE_bc==2) #BE case
  idx2=which(sampletable$phenoBE_bc==1) #CO
  y=c(rep(1,length(idx1)),rep(0,length(idx2))) #outcome
  x=as.numeric(predict_mat[i,c(idx1,idx2)]) #predicted geneexp
  #test a environment factor
  fit1=glm(as.formula(paste0("y~x+",envariable,"+x*",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial)
  fit0=glm(as.formula(paste0("y~",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial) #the smaller model
  res$BE_p=pchisq(summary(fit0)$deviance-summary(fit1)$deviance,df=2,lower.tail=FALSE)
  idx1=which(sampletable$phenoEA_bc==2) #EA case
  idx2=which(sampletable$phenoEA_bc==1) #CO
  y=c(rep(1,length(idx1)),rep(0,length(idx2))) #outcome
  x=as.numeric(predict_mat[i,c(idx1,idx2)]) #predicted geneexp
  fit1=glm(as.formula(paste0("y~x+",envariable,"+x*",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial)
  fit0=glm(as.formula(paste0("y~",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial) #the smaller model
  res$EA_p=pchisq(summary(fit0)$deviance-summary(fit1)$deviance,df=2,lower.tail=FALSE)
  idx1=which(sampletable$phenoEA_bc==2) #EA case
  idx2=which(sampletable$phenoBE_bc==2) #BE case
  y=c(rep(1,length(idx1)),rep(0,length(idx2))) #outcome
  x=as.numeric(predict_mat[i,c(idx1,idx2)]) #predicted geneexp
  fit1=glm(as.formula(paste0("y~x+",envariable,"+x*",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial)
  fit0=glm(as.formula(paste0("y~",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial) #the smaller model
  res$BEA_p=pchisq(summary(fit0)$deviance-summary(fit1)$deviance,df=2,lower.tail=FALSE)
  idx1=which(sampletable$phenoEA_bc==2 | sampletable$phenoBE_bc==2) #BE EA case
  idx2=which(sampletable$phenoBE_bc==1) #CO
  y=c(rep(1,length(idx1)),rep(0,length(idx2))) #outcome
  x=as.numeric(predict_mat[i,c(idx1,idx2)]) #predicted geneexp
  fit1=glm(as.formula(paste0("y~x+",envariable,"+x*",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial)
  fit0=glm(as.formula(paste0("y~",envariable,"+age+sex+pc1+pc2+pc3+pc4")),data=sampletable[c(idx1,idx2),],family = binomial) #the smaller model
  res$BEEA_p=pchisq(summary(fit0)$deviance-summary(fit1)$deviance,df=2,lower.tail=FALSE)
  return(res)
}
res_2df_min=NULL
for (i in 1:nrow(predict_min))
{
  if (i %% 500==0) cat(i,'..')
  res_2df_min=rbind(res_2df_min,compute2dfp(i))
}
TCGA_res_2df_min=res_2df_min

mpi_compute2dfp=function(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_TCGA_PC4",
                         envariable="bmi_recent_healthy",opt="min")
{
  
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
  colnames(predict_min)=colnames(predict_1se)=geneexpsamplenames
  
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
  predict_1se=predict_1se[,idx]
  idx=match(allsamples,colnames(eigenstratmatrix))
  eigenstratmatrix=eigenstratmatrix[,idx]
  idx=match(allsamples,sampletable$localid)
  sampletable=sampletable[idx,]
  sampletable=cbind.data.frame(sampletable,t(eigenstratmatrix[1:4,]))
  sampletable$sex=as.factor(sampletable$sex)
  library(Rmpi)
  njobs=mpi.universe.size() - 1
  print(njobs)
  if (mpi.comm.size()==0)
  {
    mpi.spawn.Rslaves(nslaves=njobs,needlog = F)
  }
  mpi.bcast.Robj2slave(compute2dfp)
  mpi.bcast.Robj2slave(predict_min)
  mpi.bcast.Robj2slave(predict_1se)
  mpi.bcast.Robj2slave(sampletable)
  mpi.bcast.Robj2slave(envariable)
  if (opt=="min")
  {
    nrows=nrow(predict_min)
  }else
  {
    nrows=nrow(predict_1se)
  }
  
  nchunks=ceiling(nrows/1000)
  print(paste0("number of chunks:",nchunks))
  allres=NULL
  for (i in 1:nchunks)
  {
    cat(i,'..')
    if (i<nchunks)
    {
      seq=((i-1)*1000+1):(i*1000)
    }else
    {
      seq=((i-1)*1000+1):nrows
    }
    res=mpi.parSapply(X=seq,FUN=compute2dfp,envariable=envariable,opt=opt,job.num = njobs)
    res=t(res)
    rownames(res)=rownames(predict_min)[seq]
    allres=rbind.data.frame(allres,res)
  }
  allres=as.data.frame(allres)
  return(allres)
}

TCGA_predictmin_bmi=mpi_compute2dfp()
TCGA_predict1se_bmi=mpi_compute2dfp(opt="1se")
TCGA_predictmin_hbrf=mpi_compute2dfp(envariable="recurrent_HB_RF")
TCGA_predict1se_hbrf=mpi_compute2dfp(envariable="recurrent_HB_RF",opt="1se")
TCGA_predictmin_smk=mpi_compute2dfp(envariable="cig_smk_ever")
TCGA_predict1se_smk=mpi_compute2dfp(envariable="cig_smk_ever",opt="1se")
TCGA_predictmin_nsaid=mpi_compute2dfp(envariable="nsaid_ever")
TCGA_predict1se_nsaid=mpi_compute2dfp(envariable="nsaid_ever",opt="1se")

GTEx_predictmin_bmi=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4")
GTEx_predict1se_bmi=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",opt="1se")
GTEx_predictmin_hbrf=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="recurrent_HB_RF")
GTEx_predict1se_hbrf=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="recurrent_HB_RF",opt="1se")
GTEx_predictmin_smk=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="cig_smk_ever")
GTEx_predict1se_smk=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="cig_smk_ever",opt="1se")
GTEx_predictmin_nsaid=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="nsaid_ever")
GTEx_predict1se_nsaid=mpi_compute2dfp(outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_PC4",envariable="nsaid_ever",opt="1se")

save(TCGA_predictmin_bmi,TCGA_predict1se_bmi,TCGA_predictmin_hbrf,TCGA_predict1se_hbrf,
     TCGA_predictmin_smk,TCGA_predict1se_smk,TCGA_predictmin_nsaid,TCGA_predict1se_nsaid,
     GTEx_predictmin_bmi,GTEx_predict1se_bmi,GTEx_predictmin_hbrf,GTEx_predict1se_hbrf,
     GTEx_predictmin_smk,GTEx_predict1se_smk,GTEx_predictmin_nsaid,GTEx_predict1se_nsaid,
     file="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/try2df1_res.RData")

qqplots=function(res=TCGA_predictmin_bmi,label="bmi_recent_healthy")
{
  layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = TRUE), 
         widths=c(1,1), heights=c(0.2,1,1))
  par(mar=c(1,1,1,1))
  plot.new()
  text(x=0.5*(par("usr")[1]+par("usr")[2]),y=0.6*(par("usr")[3]+par("usr")[4]),label,cex=1.8)
  par(mar=c(2,2,1,1))
  qqplot(as.numeric(res$BE_p))
  tmp=p.adjust(as.numeric(res$BE_p),method="fdr")
  idx=which(tmp<0.05)
  if (length(idx)>0)
  {
    text(x=0.5*(par("usr")[1]+par("usr")[2]),y=0.8*(par("usr")[3]+par("usr")[4]),paste0(rownames(res)[idx],collapse = " "),cex=1.3)
  }
  qqplot(as.numeric(res$EA_p))
  tmp=p.adjust(as.numeric(res$EA_p),method="fdr")
  idx=which(tmp<0.05)
  if (length(idx)>0)
  {
    text(x=0.5*(par("usr")[1]+par("usr")[2]),y=0.8*(par("usr")[3]+par("usr")[4]),paste0(rownames(res)[idx],collapse = " "),cex=1.3)
  }
  qqplot(as.numeric(res$BEA_p))
  tmp=p.adjust(as.numeric(res$BEA_p),method="fdr")
  idx=which(tmp<0.05)
  if (length(idx)>0)
  {
    text(x=0.5*(par("usr")[1]+par("usr")[2]),y=0.8*(par("usr")[3]+par("usr")[4]),paste0(rownames(res)[idx],collapse = " "),cex=1.3)
  }
  qqplot(as.numeric(res$BEEA_p))
  tmp=p.adjust(as.numeric(res$BEEA_p),method="fdr")
  idx=which(tmp<0.05)
  if (length(idx)>0)
  {
    text(x=0.5*(par("usr")[1]+par("usr")[2]),y=0.8*(par("usr")[3]+par("usr")[4]),paste0(rownames(res)[idx],collapse = " "),cex=1.3)
  }
}
qqplots()
qqplots(res=TCGA_predictmin_hbrf,label="recurrent_HB_RF")
qqplots(res=TCGA_predictmin_smk,label="cig_smk_ever")
qqplots(res=TCGA_predictmin_nsaid,label="nsaid_ever")

qqplots(res=GTEx_predictmin_bmi)
qqplots(res=GTEx_predictmin_hbrf,label="recurrent_HB_RF")
qqplots(res=GTEx_predictmin_smk,label="cig_smk_ever")
qqplots(res=GTEx_predictmin_nsaid,label="nsaid_ever")
