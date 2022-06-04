#try to use DNN to compute PRS
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

#NN example https://tensorflow.rstudio.com/tutorials/beginners/
mnist <- dataset_mnist()
mnist$train$x <- mnist$train$x/255
mnist$test$x <- mnist$test$x/255
model <- keras_model_sequential() %>% 
  layer_flatten(input_shape = c(28, 28)) %>% 
  layer_dense(units = 128, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(10, activation = "softmax")
summary(model)
model %>% 
  compile(
    loss = "sparse_categorical_crossentropy",
    optimizer = "adam",
    metrics = "accuracy"
  )
model %>% 
  fit(
    x = mnist$train$x, y = mnist$train$y,
    epochs = 5,
    validation_split = 0.3,
    verbose = 2
  )
predictions <- predict(model, mnist$test$x)
head(predictions, 2)

model %>% 
  evaluate(mnist$test$x, mnist$test$y, verbose = 0)
save_model_tf(object = model, filepath = "model")
reloaded_model <- load_model_tf("model")
all.equal(predict(model, mnist$test$x), predict(reloaded_model, mnist$test$x))


#prepare data-----------------------------------
#split samples into training and testing----

sampletable=readxl::read_excel("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
#sampletable=xlsx::read.xlsx("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/PLINKinputCombo_bca_07Feb2018.xls",1)
#for EA
casesamples=sampletable$localid[sampletable$phenoEA_bca==2]
contsamples=sampletable$localid[sampletable$phenoEA_bca==1]

#for BEEA
casesamples=sampletable$localid[sampletable$phenoEABE_bca==2]
contsamples=sampletable$localid[sampletable$phenoEABE_bca==1]

set.seed(1000)
idx=sample(length(casesamples),round(length(casesamples)*0.8))
traincasesamples=casesamples[idx]
idx=sample(length(contsamples),round(length(contsamples)*0.8))
traincontsamples=contsamples[idx]
testcasesamples=casesamples[!casesamples %in% traincasesamples]
testcontsamples=contsamples[!contsamples %in% traincontsamples]
trainsamples=c(traincasesamples,traincontsamples)
testsamples=c(testcasesamples,testcontsamples)
idx=match(trainsamples,sampletable$localid)
tmp=sampletable[idx,c(1,2)]
#PRSDNN_train_plinksamples.txt and PRSDNN_test_plinksamples.txt: EA
#write.table(tmp,file="../result/PRSDNN_train_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)
write.table(tmp,file="../result/PRSDNNBEEA_train_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)
idx=match(testsamples,sampletable$localid)
# tmp=sampletable[idx,c(1,2)]
# write.table(tmp,file="../result/PRSDNN_test_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)
# tmp=sampletable[idx,c(1,2)]
write.table(tmp,file="../result/PRSDNNBEEA_test_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)

#generate genotype files in bash
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
#EA
$plink  --bfile $prefix --keep ../result/PRSDNN_train_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNN_BCAtrain
$plink  --bfile $prefix --keep ../result/PRSDNN_train_plinksamples.txt \
--make-bed --out ../result/PRSDNN_BCAtrain
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
$plink  --bfile $prefix --keep ../result/PRSDNN_test_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNN_BCAtest
$plink  --bfile $prefix --keep ../result/PRSDNN_test_plinksamples.txt \
--make-bed --out ../result/PRSDNN_BCAtest

#BEEA
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
$plink  --bfile $prefix --keep ../result/PRSDNNBEEA_train_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNNBEEA_BCAtrain
$plink  --bfile $prefix --keep ../result/PRSDNNBEEA_train_plinksamples.txt \
--make-bed --out ../result/PRSDNNBEEA_BCAtrain
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
$plink  --bfile $prefix --keep ../result/PRSDNNBEEA_test_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNNBEEA_BCAtest
$plink  --bfile $prefix --keep ../result/PRSDNNBEEA_test_plinksamples.txt \
--make-bed --out ../result/PRSDNNBEEA_BCAtest

#use EA phenotype
updatefam=function(famfile="../result/PRSDNN_BCAtrain.fam",opt="EA")
{
  fam=read.table(famfile)
  idx=match(fam$V2,sampletable$localid)
  if (opt=="EA")
  {
    fam$V6=sampletable$phenoEA_bca[idx]
  }
  if (opt=="BEEA")
  {
    fam$V6=sampletable$phenoEABE_bca[idx]
  }
  if (opt=="BE")
  {
    fam$V6=sampletable$phenoBE_bca[idx]
  }
  write.table(fam,file=famfile,row.names = F,col.names = F,sep="\t",quote=F)
}
#EA
updatefam()
updatefam(famfile="../result/PRSDNN_BCAtest.fam")
#BEEA
updatefam(famfile="../result/PRSDNNBEEA_BCAtrain.fam",opt="BEEA")
updatefam(famfile="../result/PRSDNNBEEA_BCAtest.fam",opt="BEEA")

covariateforGWAS=function(famfile="../result/PRSDNN_BCAtrain.fam",prefix="../result/PRSDNN_BCAtrain")
{
  fam=read.table(famfile)
  idx=match(fam$V2,sampletable$localid)
  
  tmp=data.frame(FID=sampletable$fam[idx],
                 IID=sampletable$localid[idx],
                 pc1=sampletable$ev1_bca[idx],
                 pc2=sampletable$ev2_bca[idx],
                 pc3=sampletable$ev3_bca[idx],
                 pc4=sampletable$ev4_bca[idx],
                 age=sampletable$age[idx],
                 sex=sampletable$sex[idx])
  
  
  #for validate_twas, ind ID changed
  write.table(tmp,file=paste0(prefix,".covariate"),row.names = F,col.names = T,sep="\t",quote=F)
}
#EA
covariateforGWAS()
covariateforGWAS(famfile="../result/PRSDNN_BCAtest.fam",prefix="../result/PRSDNN_BCAtest")

#BEEA
covariateforGWAS(famfile="../result/PRSDNNBEEA_BCAtrain.fam",prefix="../result/PRSDNNBEEA_BCAtrain")
covariateforGWAS(famfile="../result/PRSDNNBEEA_BCAtest.fam",prefix="../result/PRSDNNBEEA_BCAtest")
#run gwas in bash, consider covariate
#EA
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtrain
$plink --bfile  $prefix --covar $prefix.covariate -covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out $prefix
#BEEA
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNNBEEA_BCAtrain
$plink --bfile  $prefix --covar $prefix.covariate -covar-name pc1,pc2,pc3,pc4,sex,age --logistic --hide-covar --ci 0.95 --out $prefix


#read gwas result
#EA
gwas=read.table("../result/PRSDNN_BCAtrain.assoc.logistic",header = T)
quantile(gwas$P,na.rm=T)
# 0%       25%       50%       75%      100% 
# 1.243e-06 2.457e-01 4.962e-01 7.483e-01 1.000e+00
sum(gwas$P<0.05,na.rm=T)
gwas=read.table("../result/PRSDNNBEEA_BCAtrain.assoc.logistic",header = T)
quantile(gwas$P,na.rm=T)
# 0%       25%       50%       75%      100% 
# 1.509e-07 2.399e-01 4.906e-01 7.440e-01 1.000e+00 

pcutoff=0.001
pcutoff=0.01
pcutoff=0.05
selectsnps=gwas$SNP[which(gwas$P<pcutoff)] #EA:x,9211,43584;BEEA:x,10003,46404

#read genotype
tolocalid=function(samplenames)
{
  samplenames=strsplit(samplenames,"_")
  samplenames=sapply(1:length(samplenames),function(x){
    tmp=samplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}
readgenotype=function(prefix="../result/PRSDNN_BCAtrain")
{
  tmp=data.table::fread(paste0(prefix,".traw"))
  rownames(tmp)=tmp$SNP
  idx=match(selectsnps,rownames(tmp))
  if (sum(is.na(idx))>0) warning("some snps are missing!")
  genotype=tmp[idx,7:ncol(tmp)]
  rownames(genotype)=tmp$SNP[idx]
  k <- which(is.na(genotype), arr.ind=TRUE)
  length(k)/nrow(genotype)/ncol(genotype)
  genotype[k] <- rowMeans(genotype, na.rm=TRUE)[k[,1]]
  
  return(genotype)
}

form_traintestdata=function(prefix1="../result/PRSDNN_BCAtrain",prefix2="../result/PRSDNN_BCAtest",opt="EA")
{
  traingenotype=readgenotype(prefix=prefix1)
  colnames(traingenotype)=tolocalid(samplenames = colnames(traingenotype))
  testgenotype=readgenotype(prefix=prefix2)
  colnames(testgenotype)=tolocalid(samplenames = colnames(testgenotype))
  idx1=match(colnames(traingenotype),sampletable$localid)
  idx2=match(colnames(testgenotype),sampletable$localid)
  if (opt=="EA")
  {
    trainy=as.numeric(I(sampletable$phenoEA_bca[idx1]==2))
    testy=as.numeric(I(sampletable$phenoEA_bca[idx2]==2))
  }
  if (opt=="BEEA")
  {
    trainy=as.numeric(I(sampletable$phenoEABE_bca[idx1]==2))
    testy=as.numeric(I(sampletable$phenoEABE_bca[idx2]==2))
  }
  return(list(traingenotype=traingenotype,trainy=trainy,testgenotype=testgenotype,testy=testy))
}
#EA
alldat=form_traintestdata()
#BEEA
alldat=form_traintestdata(prefix1="../result/PRSDNNBEEA_BCAtrain",prefix2="../result/PRSDNNBEEA_BCAtest",opt="BEEA")

traingenotype=alldat$traingenotype
trainy=alldat$trainy
testgenotype=alldat$testgenotype
testy=alldat$testy

#test run
model <- keras_model_sequential() %>% 
  layer_dense(units = 500, activation = "relu",input_shape=c(nrow(traingenotype)),kernel_constraint=constraint_maxnorm(3)) %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units = 250, activation = "relu",kernel_constraint=constraint_maxnorm(3)) %>% 
  layer_dropout(0.2) %>% 
  layer_dense(units = 50, activation = "relu",kernel_constraint=constraint_maxnorm(3)) %>% 
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
earlyStopping = callback_early_stopping(monitor="val_auc", min_delta=0,patience=5,mode='auto') 

set.seed(1000)
model %>% 
  fit(
    x = t(traingenotype), y = trainy,
    callbacks=list(check, earlyStopping),
    batch_size=512,
    epochs = 10,
    validation_split = 0.3,
    verbose = 3
  )

predictions <- predict(model, t(testgenotype))
names(predictions)=colnames(testgenotype)
#boxplot(predictions~testy)
plot_ROC(predictions,opt="BEEA")

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
testDNN3=function(n1=1000,n2=250,n3=50,batch_size=512,lr=0.001,dropout1=0.2,dropout2=0.2,dropout3=0.3,opt="EA",checkfile="EA_earlystop_dropout_constW")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",kernel_constraint=constraint_maxnorm(3),input_shape=c(nrow(traingenotype))) %>% 
    layer_dropout(dropout1) %>%
    layer_dense(units = n2, activation = "relu",kernel_constraint=constraint_maxnorm(3)) %>% 
    layer_dropout(dropout2) %>%
    layer_dense(units = n3, activation = "relu",kernel_constraint=constraint_maxnorm(3)) %>% 
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
      x = t(traingenotype), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  
  if (res>0.6)
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
testDNN2=function(n1=1000,n2=250,batch_size=512,lr=0.001,dropout1=0.2,dropout2=0.2,opt="EA",checkfile="EA_earlystop_dropout_constW")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",kernel_constraint=constraint_maxnorm(3),input_shape=c(nrow(traingenotype))) %>% 
    layer_dropout(dropout1) %>% 
    layer_dense(units = n2, activation = "relu",kernel_constraint=constraint_maxnorm(3)) %>% 
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
      x = t(traingenotype), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  if (res>0.6)
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
testDNN1=function(n1=1000,batch_size=512,lr=0.001,dropout=0.2,opt="EA",checkfile="EA_earlystop_dropout_constW")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_dense(units = n1, activation = "relu",kernel_constraint=constraint_maxnorm(3),input_shape=c(nrow(traingenotype))) %>% 
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
      x = t(traingenotype), y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  if (res>0.6)
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


opt="EA"
opt="BEEA"

allres3=NULL
for (n1 in c(1600,1400))
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
            res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt)
            allres3=rbind(allres3,res)
          }
      }
    }
  }
}

allres2=NULL
for (n1 in c(1600,1400))
{
  n2=round(n1/2)
  for (batch_size in c(256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.4,0.2))
      {
        for (dropout2 in seq(0.2,0.4,0.2))
        {
          res=data.frame(opt=opt,n1=n1,n2=n2,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2)
          res$auc=testDNN2(n1,n2,batch_size,lr,dropout1=dropout1,dropout2=dropout2,opt=opt)
          allres2=rbind(allres2,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(1400,1000,800,500))
{
  for (batch_size in c(256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.6,0.2))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
      
    }
  }
}
#PRS_DNN_res.RData for pvalue=0.001
#PRS_DNN_res_01.RData for pvalue=0.01,without early stopping

#EA
save(allres3,allres2,allres1,file="../result/PRS_DNN_earlydropoutconstraintWres_01.RData")
#BEEA

save(allres3,allres2,allres1,file="../result/PRS_DNNBEEA_earlydropoutconstraintWres_01.RData")

#to implement early stopping https://stackoverflow.com/questions/63150719/early-stopping-based-on-auc
#https://stackoverflow.com/questions/58682098/keras-callback-modelcheckpoint-doesnt-save-weights?rq=1
#https://community.rstudio.com/t/use-auc-as-metric-in-keras-for-r/84573/2
#https://cran.r-project.org/web/packages/keras/vignettes/training_callbacks.html

#https://machinelearningmastery.com/tensorflow-tutorial-deep-learning-with-tf-keras/

#for p=0.05
opt="EA"
opt="BEEA"
allres3=NULL
for (n1 in c(1400,1000,800,500,300))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in c(64,256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.8,0.2))
      {
        for (dropout2 in seq(0.2,0.6,0.2))
        {
            res=data.frame(opt=opt,n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2)
            res$auc=testDNN3(n1,n2,n3,batch_size,lr,dropout1=dropout1,dropout2=dropout2,dropout3=dropout2,opt=opt)
            allres3=rbind(allres3,res)

        }
        
      }
    }
  }
}

allres2=NULL
for (n1 in c(1400,1000,800,500,300))
{
  n2=round(n1/2)
  for (batch_size in c(64,256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout1 in seq(0.2,0.8,0.2))
      {
        for (dropout2 in seq(0.2,0.8,0.2))
        {
          res=data.frame(opt=opt,n1=n1,n2=n2,batch_size=batch_size,lr=lr,dropout1=dropout1,dropout2=dropout2)
          res$auc=testDNN2(n1,n2,batch_size,lr,dropout1=dropout1,dropout2=dropout2,opt=opt)
          allres2=rbind(allres2,res)
        }
      }
    }
  }
}

allres1=NULL
for (n1 in c(1400,1000,800,600,500,200))
{
  for (batch_size in c(64,256,512))
  {
    for (lr in c(0.0001))
    {
      for (dropout in seq(0.2,0.9,0.1))
      {
        res=data.frame(opt=opt,n1=n1,batch_size=batch_size,lr=lr,dropout=dropout)
        res$auc=testDNN1(n1,batch_size,lr,dropout=dropout,opt=opt)
        allres1=rbind(allres1,res)
      }
      
    }
  }
}
#EA
save(allres3,allres2,allres1,file="../result/PRS_DNN_earlydropoutconstraintWres_05.RData")
#BEEA
save(allres3,allres2,allres1,file="../result/PRS_DNNBEEA_earlydropoutconstraintWres_05.RData")

#work on linear model
library(glmnet)
set.seed(1000)
require(doMC)
registerDoMC(cores=12)
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
predictions=predict(cvfit, s=cvfit$lambda.min,newx = data.matrix(t(testgenotype)), type = "response")
names(predictions)=colnames(testgenotype)
#boxplot(predictions~testy)
plot_ROC(predictions,opt="EA") #p=0.01:0.565,p=0.05:0.578
plot_ROC(predictions,opt="BEEA") #p=0.01:0.564,p=0.05:0.571
