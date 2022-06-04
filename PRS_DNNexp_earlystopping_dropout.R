#try to use DNN to compute PRS, use predicted geneexp
#fisrt hidden layer dropout 0.2-0.9, the rest 0.5

#load installed tensorflow
#ml TensorFlow/2.1.0-foss-2019b-Python-3.7.4

#install---
#install.packages("keras")
#install.packages("tensorflow")
#install_tensorflow(method = "conda", envname = "Tensorflow_for_R")
#conda_install(envname = "Tensorflow_for_R", packages="numpy")
#load tensorflow---
#startup::restart()
library(keras)
library(tensorflow)
reticulate::use_condaenv("Tensorflow_for_R")
reticulate::py_config()
reticulate::py_run_string('import tensorflow')
tf$constant("Hellow Tensorflow") #to see if it works

#reticulate::conda_list()
#use_condaenv("Tensorflow_for_R") #this work
#use_virtualenv("/fh/fast/dai_j/BEACON/BEACON_GRANT/code/tensorflowenv")
#locate installed tensorflow https://tensorflow.rstudio.com/installation/custom/
#Sys.setenv(RETICULATE_PYTHON="/app/software/Python/3.7.4-GCCcore-8.3.0/bin/python") #this worked
#Sys.getenv("RETICULATE_PYTHON")
tf$constant("Hellow Tensorflow")


#prepare data-----------------------------------
#split samples into training and testing----

sampletable=readxl::read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
#sampletable=xlsx::read.xlsx("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
#for EA
tmp=read.table("../result/PRSDNN_train_plinksamples.txt")
idx=match(tmp$V2,sampletable$localid)
trainsamples=tmp$V2

tmp=read.table("../result/PRSDNN_test_plinksamples.txt")
idx=match(tmp$V2,sampletable$localid)
testsamples=tmp$V2

#load expression data
#prefix="dist500K_GTEx_mucosa_Jan19" #AMOS sample uses different sampleID (not the same as those in sampletable)
prefix="dist500K_GTEx_mucosa_June11"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
#exp models
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
#predicted exp
load(paste0(outfolder,"/bca_predict_geneexp.RData"))
if ("n_totalsnp" %in% colnames(predict_min))
{
  predict_min=predict_min[,3:ncol(predict_min)]
}
if (grepl("_",colnames(predict_min)[1]))
{
  geneexpsamplenames=strsplit(colnames(predict_min),"_") #use localid
  geneexpsamplenames=sapply(1:length(geneexpsamplenames),function(x){
    tmp=geneexpsamplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
  colnames(predict_min)=geneexpsamplenames
}
tmp=rowMeans(predict_min,na.rm=T)
predict_min=predict_min[!is.na(tmp),]
#compute p-value:
pvalue=rep(NA,nrow(predict_min))
names(pvalue)=rownames(predict_min)
traindat=predict_min[,colnames(predict_min) %in% trainsamples]
idx=match(colnames(traindat),sampletable$localid)
y=as.numeric(sampletable$phenoEA_bca[idx]==2)
for (i in 1:nrow(predict_min))
{
  if (i %%1000==0) cat(i,'..')
  fit=glm(y~as.numeric(traindat[i,]),family = "binomial")
  pvalue[i]=summary(fit)$coefficients[2,4]
}

idx=match(names(pvalue),rownames(res_min))
table(res_min$r2[idx]>0.05,pvalue<0.05)

formdat=function(r2cutoff=0.05,pcutoff=0.2)
{
  comsamples=intersect(colnames(predict_min),sampletable$localid)
  allpredict=predict_min[,match(comsamples,colnames(predict_min))]
  all(rownames(allpredict) %in% rownames(res_min))
  selectgenes=intersect(rownames(res_min)[res_min$r2>r2cutoff],names(pvalue)[pvalue<pcutoff])
  trainsamples1=intersect(trainsamples,comsamples)
  traindat=allpredict[rownames(allpredict) %in% selectgenes,match(trainsamples1,colnames(allpredict))]
  idx=match(colnames(traindat),sampletable$localid)
  trainy=as.numeric(sampletable$phenoEA_bca[idx]==2)
  testsamples1=intersect(testsamples,comsamples)
  testdat=allpredict[rownames(allpredict) %in% selectgenes,match(testsamples1,colnames(allpredict))]
  all(rownames(traindat)==rownames(testdat))
  idx=match(colnames(testdat),sampletable$localid)
  testy=as.numeric(sampletable$phenoEA_bca[idx]==2)
  return(list(traindat=traindat,testdat=testdat,trainy=trainy,testy=testy))
}

dat_05_2=formdat()
dat_05_1=formdat(pcutoff=0.1)
dat_05_05=formdat(pcutoff=0.05)
dat_05=formdat(pcutoff=1)
dat_0_05=formdat(r2cutoff = 0,pcutoff = 0.05)
save(prefix,pvalue,predict_min,res_min,dat_05,dat_05_1,dat_05_2,dat_05_05,dat_0_05,file="../result/PRS_DNNexptraindat.RData")
#test run
model <- keras_model_sequential() %>% 
  layer_dense(units = 500,input_shape=c(nrow(traindat)),activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units = 250, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(1, activation = "sigmoid")
summary(model)

auc=tf$keras$metrics$AUC(name='auc')
model %>% 
  compile(
    loss = "binary_crossentropy",
    optimizer=optimizer_adam(lr = 0.001),
    metrics=auc
  )

filepath = "weights_best_dropout.hdf5"
check = callback_model_checkpoint(filepath,
                                           monitor='val_auc',  # validation AUC
                                           save_best_only=T,
                                           mode='max')
earlyStopping = callback_early_stopping(monitor="val_auc", min_delta=0,patience=3,mode='auto') 

set.seed(1000)
model %>% 
  fit(
    x = t(traindat), y = trainy,
    callbacks=list(check, earlyStopping),
    batch_size=512,
    epochs = 10,
    validation_split = 0.3,
    verbose = 3
  )

predictions <- predict(model, t(testdat))
names(predictions)=colnames(testdat)
#boxplot(predictions~testy)
plot_ROC(predictions,opt="EA")

plot_ROC=function(predictions,opt="EA")
{
  idx=match(names(predictions),sampletable$localid)
  sampletable1=sampletable[idx,]
  sampletable1$prs=predictions
  #for EA
  if (opt=="EA")
  {
    sampletable1$case=sampletable1$phenoEA_bca
  }
  if (opt=="BEEA")
  {
    sampletable1$case=sampletable1$phenoEABE_bca
  }
  sampletable1$case[sampletable1$case==-9]=NA
  fit1 <- glm(I(case==2)~prs,family=binomial,data=sampletable1,y=T)
  roc1=pROC::roc(fit1$y,fit1$fitted.values,quiet = T)
  #remove bmi_recent_healthy
  dat=sampletable1[,c("case","age","sex","recurrent_HB_RF","cig_smk_ever","nsaid_ever","prs","site")]
  idx=complete.cases(dat)
  if (sum(idx)>0)
  {
    fit3 <- glm(I(case==2)~age+sex+recurrent_HB_RF+cig_smk_ever+nsaid_ever+prs,family=binomial,data=sampletable1,y=T)
    roc3=pROC::roc(fit3$y,fit3$fitted.values,quiet = T)
    print(paste0("auc_prs=",round(roc1$auc,3)," auc_env=",round(roc3$auc,3)))
  }else
  {
    print(paste0("auc_prs=",round(roc1$auc,3)))
  }
  return(round(roc1$auc,3))
}

auc=tf$keras$metrics$AUC(name='auc')
# check = callback_model_checkpoint(filepath,
#                                   monitor='val_auc',  # validation AUC
#                                   save_best_only=T,
#                                   mode='max')
earlyStopping = callback_early_stopping(monitor="val_auc", patience=3,mode='auto') 
#3 hidden layers NN
testDNN3=function(n1=1000,n2=250,n3=50,batch_size=512,lr=0.001,dropout1=0.2,dropout2=0.2,dropout3=0.3,opt="EA",checkfile="EAexp_earlystop_dropout")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",input_shape=c(nrow(traindat))) %>% 
    layer_dropout(dropout1) %>%
    layer_dense(units = n2, activation = "relu") %>% 
    layer_dropout(dropout2) %>%
    layer_dense(units = n3, activation = "relu") %>% 
    layer_dropout(dropout3) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      #metrics = "accuracy"
      metrics=auc
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traindat), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testdat))
  names(predictions)=colnames(testdat)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  
  if (res>0.62)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,n1,"_",n2,"_",n3,"batch_size",batch_size,"_",dropout1,"_",dropout2,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

#2 hidden layers
testDNN2=function(n1=1000,n2=250,batch_size=512,lr=0.001,dropout1=0.2,dropout2=0.2,opt="EA",checkfile="EAexp_earlystop_dropout")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",input_shape=c(nrow(traindat))) %>% 
    layer_dropout(dropout1) %>% 
    layer_dense(units = n2, activation = "relu") %>% 
    layer_dropout(dropout2) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      metrics=auc
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traindat), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
 
  predictions <- predict(model, t(testdat))
  names(predictions)=colnames(testdat)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  if (res>0.62)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,n1,"_",n2,"batch_size",batch_size,"_",dropout1,"_",dropout2,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

#1 hidden layer
testDNN1=function(n1=1000,batch_size=512,lr=0.001,dropout=0.2,opt="EA",checkfile="EAexp_earlystop_dropout")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",input_shape=c(nrow(traindat))) %>% 
    layer_dropout(dropout) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      metrics=auc
      #metrics=tf$keras$metrics$AUC()
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traindat), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
 
  predictions <- predict(model, t(testdat))
  names(predictions)=colnames(testdat)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  if (res>0.62)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,n1,"batch_size",batch_size,"_",dropout,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

#Normalize test data with AGAIN mean and standart deviation of TRAINING DATA set.
library(matrixStats)
normalizedat=function(traindat,testdat)
{
  tmp1=rowMeans(traindat)
  tmp2=rowSds(as.matrix(traindat))
  dat1=(traindat-tmp1)/tmp2
  dat2=(testdat-tmp1)/tmp2
  return(list(traindat=dat1,testdat=dat2))
  
}

load("../result/PRS_DNNexptraindat.RData")
#dat_05
opt="EA"
checkfile="EAdat05_earlystop_dropout"
traindat=dat_05$traindat
testdat=dat_05$testdat
# tmp=normalizedat(traindat =dat_05$traindat,testdat = dat_05$testdat )
# traindat=tmp$traindat
# testdat=tmp$testdat
trainy=dat_05$trainy
testy=dat_05$testy
allres3=NULL

for (n1 in c(5000,4000,3000,2000,1600,1400,1000,500,200))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.6,0.2))
      {
        for (dropout2 in seq(0.2,0.6,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
          gc()
        }
      }
    }
  }
}
allres1=NULL
for (n1 in c(5000,4000,3000,2000,1600,1400,1000,800,500,300,200))
{
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.6,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout,checkfile = checkfile)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
        gc()
      }
    }
  }
}
save(allres3,file="../result/PRS_DNNexp_earlydropoutres_dat05.RData")

#dat_05_01
opt="EA"
checkfile="EAdat0501_earlystop_dropout"
# tmp=normalizedat(traindat =dat_05_1$traindat,testdat = dat_05_1$testdat )
# traindat=tmp$traindat
# testdat=tmp$testdat
traindat=dat_05_1$traindat
testdat=dat_05_1$testdat
trainy=dat_05_1$trainy
testy=dat_05_1$testy
allres3=NULL

for (n1 in c(1000,800,500,300,200,100))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.6,0.2))
      {
        for (dropout2 in seq(0.2,0.6,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(1000,800,500,300,200,100))
{
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.6,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout,checkfile = checkfile)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
    }
  }
}

save(allres3,allres1,file="../result/PRS_DNNexp_earlydropoutres_dat05_01.RData")

#dat_05_02
opt="EA"
checkfile="EAdat0502_earlystop_dropout"
# tmp=normalizedat(traindat =dat_05_2$traindat,testdat = dat_05_2$testdat )
# traindat=tmp$traindat
# testdat=tmp$testdat
traindat=dat_05_2$traindat
testdat=dat_05_2$testdat
trainy=dat_05_2$trainy
testy=dat_05_2$testy
allres3=NULL

for (n1 in c(1000,800,500,300,200,100))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.8,0.2))
      {
        for (dropout2 in seq(0.2,0.8,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(1000,800,500,300,200,100))
{
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.8,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout,checkfile = checkfile)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
    }
  }
}

save(allres3,allres1,file="../result/PRS_DNNexp_earlydropoutres_dat05_02.RData")

#dat_05_05
opt="EA"
checkfile="EAdat0505_earlystop_dropout"
# tmp=normalizedat(traindat =dat_05_05$traindat,testdat = dat_05_05$testdat )
# traindat=tmp$traindat
# testdat=tmp$testdat
traindat=dat_05_05$traindat
testdat=dat_05_05$testdat
trainy=dat_05_05$trainy
testy=dat_05_05$testy
allres3=NULL

for (n1 in c(500,400,300,200,100,80))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.8,0.2))
      {
        for (dropout2 in seq(0.2,0.8,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(500,300,200,100,80))
{
  for (batch_size in c(64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.8,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout,checkfile = checkfile)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
    }
  }
}

save(allres3,allres1,file="../result/PRS_DNNexp_earlydropoutres_dat05_05.RData")

#dat_0_05
opt="EA"
checkfile="EAdat005_earlystop_dropout"
# tmp=normalizedat(traindat =dat_0_05$traindat,testdat = dat_0_05$testdat )
# traindat=tmp$traindat
# testdat=tmp$testdat
traindat=dat_0_05$traindat
testdat=dat_0_05$testdat
trainy=dat_0_05$trainy
testy=dat_0_05$testy
allres3=NULL

for (n1 in c(1000,800,500,400,300,200,100,80))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(32,64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.8,0.2))
      {
        for (dropout2 in seq(0.2,0.8,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(1000,800,500,300,200,100,80))
{
  for (batch_size in c(32,64,256))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.8,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout,checkfile = checkfile)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
    }
  }
}

save(allres3,allres1,file="../result/PRS_DNNexp_earlydropoutres_dat0_05.RData")

# allres2=NULL
# for (n1 in c(1600,1400,1000,800))
# {
#   n2=round(n1/2)
#   for (batch_size in c(64,256,512))
#   {
#     for (lr in c(0.0001))
#     {
#       for (dropout1 in seq(0.2,0.4,0.2))
#       {
#         for (dropout2 in seq(0.2,0.4,0.2))
#         {
#           res=data.frame(opt=opt,n1=n1,n2=n2,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2)
#           res$auc=testDNN2(n1,n2,batch_size,lr,dropout1=dropout1,dropout2=dropout2,opt=opt,checkfile = checkfile)
#           allres2=rbind(allres2,res)
#         }
#       }
#     }
#   }
# }
# 

#PRS_DNN_res.RData for pvalue=0.001
#PRS_DNN_res_01.RData for pvalue=0.01,without early stopping

#EA



#to implement early stopping https://stackoverflow.com/questions/63150719/early-stopping-based-on-auc
#https://stackoverflow.com/questions/58682098/keras-callback-modelcheckpoint-doesnt-save-weights?rq=1
#https://community.rstudio.com/t/use-auc-as-metric-in-keras-for-r/84573/2
#https://cran.r-project.org/web/packages/keras/vignettes/training_callbacks.html


#check 3nn with many neurons
opt="EA"
checkfile="EA_largeNN_earlystop_p05_dropout"
allres3=NULL

for (n1 in c(10000,8000,6000,4000,3000))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.4,0.2))
      {
        for (dropout2 in seq(0.2,0.4,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}
save(allres3,file="../result/PRS_DNN_earlydropoutres_05_large3NN.RData")

#p=0.1
#check 3nn with many neurons
opt="EA"
checkfile="EA_earlystop_p1_dropout"
allres3=NULL

for (n1 in c(10000,7000,6000,5000,4000,3000))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.4,0.2))
      {
        for (dropout2 in seq(0.2,0.4,0.2))
        {
          
          res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
          res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt,checkfile = checkfile)
          allres3=rbind(allres3,res)
        }
      }
    }
  }
}
save(allres3,file="../result/PRS_DNN_earlydropoutres_p1.RData")

#work on linear model
library(glmnet)
set.seed(100)
# require(doMC)
# registerDoMC(cores=12)
cvfit=cv.glmnet(t(traingenotype), trainy, family="binomial", standardize=F,parallel = F,nfolds=10, trace.it=1)
plot(cvfit)
glmcoeff=coef(cvfit,s=cvfit$lambda.min)
idx=glmcoeff@i
snps=glmcoeff@Dimnames[[1]][idx+1]
snps=snps[!snps %in% "(Intercept)"]
snps=as.numeric(gsub("V","",snps))
snps=rownames(traingenotype)[snps]
idx=match(snps,gwas$SNP)
quantile(gwas$P[idx])
# 0%         25%         50%         75%        100% 
# 6.70900e-07 1.64000e-03 4.17550e-03 6.88925e-03 9.99700e-03 
glmnet.predictions=predict(cvfit, s=cvfit$lambda.min,newx = data.matrix(t(testgenotype)), type = "response")
names(glmnet.predictions)=colnames(testgenotype)
#boxplot(predictions~testy)
plot_ROC(glmnet.predictions,opt="EA") #p=0.01:0.565,p=0.05:0.578
plot_ROC(glmnet.predictions,opt="BEEA") #p=0.01:0.564,p=0.05:0.571

#work on svm
library(e1071)
tune.out=tune(svm ,t(traingenotype),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=c(0.1 ,1 ,10,20),gamma=c(0.01,0.1,0.5,1) ))
#summary(tune.out)
bestmod =tune.out$best.model
#summary(bestmod)
#svmfit <- svm(t(traingenotype),factor(trainy),scale=F,kernel ="radial")
svmfit <- svm(t(traingenotype),factor(trainy),scale=F,kernel ="radial", cost=bestmod$cost,gamma = bestmod$gamma)
svm.prediction=predict(svmfit, t(testgenotype))

#prediction from LDpred
tmp=read.table("../result/DNNvalidition_EA_genotyped.PRS.txt",header = T,stringsAsFactors = F)
ldpred.predictions=tmp$PRS
names(ldpred.predictions)=tmp$IID
#load prediction from NN
mymodel="../result/NNmodel/EA2_earlystop_p05_dropout_6_3000_1500_750batch_size256_0.2_0.2_0.621"
rm(reloaded_model)
reloaded_model <- load_model_tf(mymodel)
NN.predictions <- predict(reloaded_model, t(testgenotype))
names(NN.predictions)=colnames(testgenotype)
plot_ROC(NN.predictions,opt=opt)
names(testy)=colnames(testgenotype)
library(ROCR) #prediction
library(pROC) #roc,auc,plot.roc


plotroc3=function(predict1=NN.predictions,predict2=glmnet.predictions,predict3=ldpred.predictions,y=testy,main="")
{
  comsamples=intersect(names(predict1),names(predict2))
  comsamples=intersect(names(predict3),comsamples)
  predict1=predict1[match(comsamples,names(predict1))]
  predict2=predict2[match(comsamples,names(predict2))]
  predict3=predict3[match(comsamples,names(predict3))]
  y=y[match(comsamples,names(y))]
  
  predict1_=predict1[order(predict1)]
  TP1=FP1=rep(0,length(y))
  for (i in 1:length(y)) {
    TP1[i] <- mean(predict1[y==1]>=predict1_[i]) 
    FP1[i] <- mean(predict1[y==0]>=predict1_[i]) 
  }
  plot(c(FP1,0),c(TP1,0),type="l",xlim=c(0,1),ylim=c(0,1),lwd=4,xlab="1-specificity",ylab="Sensitivity",main=main,cex.lab=1.3,cex.axis=1.3,col="red",yaxs="i")
  
  #abline(0,1,lty=2)
  fit1 <- glm(I(y==1)~predict1,family=binomial,y=T)
  roc1=pROC::roc(fit1$y,fit1$fitted.values,quiet = T)
  auc1=roc1$auc
  auc1=round(as.numeric(auc1),2)
  print(paste0("auc1=",auc1))
  # pauc1=pauc(marker = predict1,status = y,fpr = 0.05)
  # print(paste0("pauc1=",round(pauc1,4)))
  #text(0.07,0.94,paste0("AUC=",auc1),cex=1.2)
  
  predict2_=predict2[order(predict2)]
  TP2=FP2=rep(0,length(y))
  for (i in 1:length(y)) {
    TP2[i] <- mean(predict2[y==1]>=predict2_[i]) 
    FP2[i] <- mean(predict2[y==0]>=predict2_[i]) 
  }
  lines(c(FP2,0),c(TP2,0),lwd=4,col="black")
  #lines(c(0.05,0.05),c(0,1),lty=2)
  fit2 <- glm(I(y==1)~predict2,family=binomial,y=T)
  roc2=pROC::roc(fit2$y,fit2$fitted.values,quiet = T)
  auc2=roc2$auc
  auc2=round(as.numeric(auc2),2)
  print(paste0("auc2=",auc2))
  # pauc2=pauc(marker = predict2,status = y,fpr = 0.05)
  # print(paste0("pauc2=",round(pauc2,4)))
  # 
  predict3_=predict3[order(predict3)]
  TP3=FP3=rep(0,length(y))
  for (i in 1:length(y)) {
    TP3[i] <- mean(predict3[y==1]>=predict3_[i]) 
    FP3[i] <- mean(predict3[y==0]>=predict3_[i]) 
  }
  lines(c(FP3,0),c(TP3,0),lwd=4,col="blue")
  fit3 <- glm(I(y==1)~predict3,family=binomial,y=T)
  roc3=pROC::roc(fit3$y,fit3$fitted.values,quiet = T)
  auc3=roc3$auc
  auc3=round(as.numeric(auc3),2)
  print(paste0("auc3=",auc3))
  # pauc3=pauc(marker = predict3,status = y,fpr = 0.05)
  # print(paste0("pauc3=",round(pauc3,4)))
  
  # auc.pvalue=roc.test(response=y,predictor1 = predict1,predictor2 = predict3)$p.value
  # print(paste0("auc.pvalue=",round(auc.pvalue,4)))
  # roc1 <- roc(y, predict1, partial.auc=c(1, 0.95), partial.auc.focus="sp")
  # roc3 <- roc(y, predict3, partial.auc=c(1, 0.95), partial.auc.focus="sp")
  # pauc.pvalue=roc.test(roc1, roc3)$p.value #0.47
  # print(paste0("pauc.pvalue=",round(pauc.pvalue,4)))
  legend=c(paste0("DNN:AUC=",round(auc1,2)),paste0("LDpred:AUC=",round(auc3,2)),paste0("Glmnet:AUC=",round(auc2,2)))
  legend("bottomright",legend=legend,col=c("red","blue","black"),lty=1,cex=1.2,bty = "n")
  # res=data.frame(auc1=auc1,pauc1=pauc1,auc2=auc2,pauc2=pauc2,auc3=auc3,pauc3=pauc3,auc.pvalue=auc.pvalue,pauc.pvalue=pauc.pvalue,stringsAsFactors = F)
  #return(res)
}
par(mar=c(6,6,2,1))
plotroc3()
