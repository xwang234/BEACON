#try to use CNN to compute PRS
#https://victorzhou.com/blog/intro-to-cnns-part-1/
#https://victorzhou.com/blog/keras-cnn-tutorial/
#https://adeshpande3.github.io/adeshpande3.github.io/A-Beginner's-Guide-To-Understanding-Convolutional-Neural-Networks/
#https://jjallaire.github.io/deep-learning-with-r-notebooks/notebooks/6.4-sequence-processing-with-convnets.nb.html #RNN and CNN
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
pcutoff=0.1
pcutoff=0.5
pcutoff=1

selectsnps=gwas$SNP[which(gwas$P<pcutoff)] #EA:x,9211,45384,831785;BEEA:x,10003,46404

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
  # bim=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous.bim")
  # idx=match(rownames(genotype),bim$V2)
  # sum(diff(idx)<0) #0, snps were ordered
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
# #BEEA
# alldat=form_traintestdata(prefix1="../result/PRSDNNBEEA_BCAtrain",prefix2="../result/PRSDNNBEEA_BCAtest",opt="BEEA")

traingenotype=alldat$traingenotype
trainy=alldat$trainy
testgenotype=alldat$testgenotype
testy=alldat$testy
traingenotype_1=traingenotype #pcutoff=0.1
testgenotype_1=testgenotype
traingenotype_05=traingenotype #pcutoff=0.05
testgenotype_05=testgenotype

xtrain = array(as.matrix(t(traingenotype)), dim = c(ncol(traingenotype), nrow(traingenotype), 1))
xtest = array(as.matrix(t(testgenotype)), dim = c(ncol(testgenotype), nrow(testgenotype), 1))

xtrain_1=xtrain
xtest_1=xtest
xtrain_05=xtrain
xtest_05=xtest
#test run
model <- keras_model_sequential() %>% 
  layer_conv_1d(filters = 32, kernel_size = 6,input_shape = c(nrow(traingenotype),1), activation = "relu") %>% 
  layer_max_pooling_1d() %>%
  layer_flatten() %>%
  layer_dense(units = 1, activation = "sigmoid")
summary(model)

auc=tf$keras$metrics$AUC(name='auc')
model %>%
  compile(
    loss = "binary_crossentropy",
    optimizer=optimizer_adam(lr = 0.001),
    metrics = auc
  )
# 
filepath = "weights_best_dropout.hdf5"
check = callback_model_checkpoint(filepath,
                                           monitor='val_auc',  # validation AUC
                                           save_best_only=T,
                                           mode='max')
earlyStopping = callback_early_stopping(monitor="val_auc", min_delta=0,patience=5,mode='auto')
# 
# set.seed(1000)
model %>%
  fit(
    x = xtrain, y = trainy,
    callbacks=list(check, earlyStopping),
    batch_size=512,
    epochs = 10,
    validation_split = 0.3,
    verbose = 1
  )


predictions <- predict(model, xtest)
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
earlyStopping = callback_early_stopping(monitor="val_auc", patience=3,mode='auto') 

testcnn=function(n1=32,kerneln=6,batch_size=512,lr=0.001,pool_size=2,opt="EA",checkfile="EA_CNN_earlystop")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_conv_1d(filters = n1, kernel_size = kerneln,input_shape = c(dim(xtrain)[2],1), activation = "relu") %>% 
    layer_max_pooling_1d(pool_size = pool_size) %>%
    layer_flatten() %>%
    layer_dense(units = 1, activation = "sigmoid")
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      #metrics = "accuracy"
      metrics=auc
    )
  set.seed(10000)
  model %>% 
    fit(
      x = xtrain, y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 1
    )
  
  predictions <- predict(model, xtest)
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  
  if (res>0.635)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,n1,"_kernel",kerneln,"_batch_size",batch_size,"_","lr",lr,"_","pool_size",pool_size,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

#add a dense layer
testcnn1=function(n1=32,n2=64,kerneln=6,batch_size=512,lr=0.001,pool_size=2,opt="EA",checkfile="EA_CNN_earlystop")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_conv_1d(filters = n1, kernel_size = kerneln,input_shape = c(dim(xtrain)[2],1), activation = "relu") %>% 
    layer_max_pooling_1d(pool_size = pool_size) %>%
    layer_flatten() %>%
    layer_dense(units = n2, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")
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
      x = xtrain, y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, xtest)
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  
  if (res>0.635)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,"n1_",n1,"n2_",n2,"_kernel",kerneln,"_batch_size",batch_size,"_","lr",lr,"_","pool_size",pool_size,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

testcnn2=function(n1=32,n2=16,kerneln=6,batch_size=512,lr=0.001,pool_size=2,opt="EA",checkfile="EA_CNN_earlystop")
{
  filepath = paste0("../result/weights_",checkfile,"_best.hdf5")
  check = callback_model_checkpoint(filepath,
                                    monitor='val_auc',  # validation AUC
                                    save_best_only=T,
                                    mode='max')
  model <- keras_model_sequential() %>% 
    layer_conv_1d(filters = n1, kernel_size = kerneln,input_shape = c(dim(xtrain)[2],1), activation = "relu") %>% 
    layer_max_pooling_1d(pool_size = pool_size) %>%
    layer_conv_1d(filters = n2, kernel_size = kerneln, activation = "relu") %>% 
    layer_max_pooling_1d(pool_size = pool_size) %>%
    layer_flatten() %>%
    layer_dense(units = n2, activation = "relu") %>%
    layer_dense(units = 1, activation = "sigmoid")
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
      x = xtrain, y = trainy,
      callbacks=list(check, earlyStopping),
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, xtest)
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions,opt=opt)
  
  if (res>0.635)
  {
    mymodel=paste0("../result/NNmodel/",checkfile,"covn1_",n1,"covn2_",n2,"_kernel",kerneln,"_batch_size",batch_size,"_","lr",lr,"_","pool_size",pool_size,"_",res)
    save_model_tf(object = model, filepath = mymodel)
  }
  if (file.exists(filepath))
  {
    file.remove(filepath)
  }
  return(res)
}

traingenotype0=traingenotype
testgenotype0=testgenotype
cleangenotype=function(dat=traingenotype)
{
  k <- which(dat>0 & dat<0.5, arr.ind=TRUE)
  if (nrow(k)>0) dat[k]=0
  k <- which(dat>=0.5 & dat<1.5 & dat!=1, arr.ind=TRUE)
  if (nrow(k)>0) dat[k]=1
  k <- which(dat>=1.5 & dat<2, arr.ind=TRUE)
  if (nrow(k)>0) dat[k]=2
  return(dat)
}

traingenotype=cleangenotype()
testgenotype=cleangenotype(dat=testgenotype)
# 
# save(traingenotype0,traingenotype,testgenotype0,testgenotype,trainy,testy,file="../result/PRS_DNN_EA_P05traindat.RData")

#p=0.05
opt="EA"
checkfile="EA_CNN_earlystop_p05_"

allres=NULL

for (n1 in c(4,8,16,32,64))
{
  for (kerneln in c(2,3,5,6,8,10,20))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001,0.0001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_earlydropoutres_05.RData")

#p=0.1
opt="EA"
checkfile="EA_CNN_earlystop_p1_"

allres=NULL

for (n1 in c(4,8,16,32,64,128,256,512))
{
  for (kerneln in c(2,3,5,6,8,10,20))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_earlydropoutres_1.RData")

#p=0.5
opt="EA"
checkfile="EA_CNN_earlystop_p5_"

allres=NULL

for (n1 in c(4,8,16,32,64))
{
  for (kerneln in c(2,3,5,6,8,10,20))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_earlydropoutres_5.RData")

#p=1
opt="EA"
checkfile="EA_CNN_earlystop_allsnps_"

allres=NULL

for (n1 in c(4,8,16,32,64))
{
  for (kerneln in c(2,3,5,6,8,10,20))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_earlydropoutres_allsnps.RData")

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
svmfit <- svm(t(traingenotype),factor(trainy),scale=F,kernel ="radial",probability=T)
print("svmfit done")
svm.predictions=predict(svmfit, t(testgenotype))
save(svmfit,svm.predictions,file="../result/PRSVMP1res.RData")
tune.out=tune(svm ,t(traingenotype),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=c(0.1 ,1 ,10,20),gamma=c(0.01,0.1,0.5,1) ))
#summary(tune.out)
bestmod =tune.out$best.model
#summary(bestmod)
svmfit <- svm(t(traingenotype),factor(trainy),scale=F,kernel ="radial", cost=bestmod$cost,gamma = bestmod$gamma)
#svmfit <- svm(t(traingenotype),factor(trainy),scale=F,kernel ="radial")
load("../result/PRSSVMres.RData")
svm.predictions=predict(svmfit, t(testgenotype))
plot_ROC(svm.predictions,opt="EA") #0.552
#prediction from LDpred
tmp=read.table("../result/DNNvalidition_EA_genotyped.PRS.txt",header = T,stringsAsFactors = F)
ldpred.predictions=tmp$PRS
names(ldpred.predictions)=tmp$IID

load("../result/PRSSVMres_1.RData") #PRSSVM.R
plot_ROC(svm.predictions,opt="EA") #0.552
load("../result/PRSSVMres_2.RData") #PRSSVM.R
plot_ROC(svm.predictions,opt="EA") #0.569
svm.predictions1=predict(svmfit, t(testgenotype_05),probability = T) #p=0.05
plot_ROC(svm.predictions1,opt="EA")
plot_ROC(attributes(svm.predictions1)[4]$probabilities[,1])
svm.predictions=attributes(svm.predictions1)[4]$probabilities[,1]

#load prediction from NN
mymodel="../result/NNmodel/EA2_earlystop_p05_dropout_6_3000_1500_750batch_size256_0.2_0.2_0.621"
mymodel="../result/NNmodel/EA_binary_earlystop_p05_dropout_36_8000_4000_2000batch_size256_0.4_0.2_0.631"
mymodel="../result/NNmodel/EA_CNN_earlystop_p1_dropout16_kernel3_batch_size256_0.64" #n1=16
mymodel="../result/NNmodel/EA_earlystop_p05_cnn_35_16_kernel6_batch_size256_lr0.001_pool_size2_0.649"
rm(reloaded_model)
reloaded_model <- load_model_tf(mymodel)
#MLP
#NN.predictions <- predict(reloaded_model, t(testgenotype))
#CNN
NN.predictions <- predict(reloaded_model, xtest)
names(NN.predictions)=colnames(testgenotype)
plot_ROC(NN.predictions)
names(testy)=colnames(testgenotype)
library(ROCR) #prediction
library(pROC) #roc,auc,plot.roc


plotroc3=function(predict1=NN.predictions,predict2=svm.predictions,predict3=ldpred.predictions,y=testy,main="")
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
  plot(c(FP1,0),c(TP1,0),type="l",xlim=c(0,1),ylim=c(0,1),lwd=4,xlab="1-specificity",ylab="Sensitivity",main=main,cex.lab=1.7,cex.axis=1.7,col="red",yaxs="i")
  
  #abline(0,1,lty=2)
  fit1 <- glm(I(y==1)~predict1,family=binomial,y=T)
  roc1=pROC::roc(fit1$y,fit1$fitted.values,quiet = T)
  auc1=roc1$auc
  auc1=round(as.numeric(auc1),3)
  print(paste0("auc1=",auc1))
  # pauc1=pauc(marker = predict1,status = y,fpr = 0.05)
  # print(paste0("pauc1=",round(pauc1,4)))
  #text(0.07,0.94,paste0("AUC=",auc1),cex=1.2)
  
  predict2=as.numeric(as.character(predict2))
  predict2_=predict2[order(predict2)]
  TP2=FP2=rep(0,length(y))
  for (i in 1:length(y)) {
    TP2[i] <- mean(predict2[y==1]>=predict2_[i]) 
    FP2[i] <- mean(predict2[y==0]>=predict2_[i]) 
  }
  lines(c(FP2,0),c(TP2,0),lwd=4,col="green")
  #lines(c(0.05,0.05),c(0,1),lty=2)
  fit2 <- glm(I(y==1)~predict2,family=binomial,y=T)
  roc2=pROC::roc(fit2$y,fit2$fitted.values,quiet = T)
  auc2=roc2$auc
  auc2=round(as.numeric(auc2),3)
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
  auc3=round(as.numeric(auc3),3)
  print(paste0("auc3=",auc3))
  # pauc3=pauc(marker = predict3,status = y,fpr = 0.05)
  # print(paste0("pauc3=",round(pauc3,4)))
  
  # auc.pvalue=roc.test(response=y,predictor1 = predict1,predictor2 = predict3)$p.value
  # print(paste0("auc.pvalue=",round(auc.pvalue,4)))
  # roc1 <- roc(y, predict1, partial.auc=c(1, 0.95), partial.auc.focus="sp")
  # roc3 <- roc(y, predict3, partial.auc=c(1, 0.95), partial.auc.focus="sp")
  # pauc.pvalue=roc.test(roc1, roc3)$p.value #0.47
  # print(paste0("pauc.pvalue=",round(pauc.pvalue,4)))
  legend=c(paste0("CNN : AUC=",round(auc1,3)),paste0("SVM : AUC=",round(auc2,3)),paste0("LDpred : AUC=",round(auc3,3)))
  legend("topleft",legend=legend,col=c("red","green","blue"),lty=1,cex=1.5,lwd=4,bty = "n")
  # res=data.frame(auc1=auc1,pauc1=pauc1,auc2=auc2,pauc2=pauc2,auc3=auc3,pauc3=pauc3,auc.pvalue=auc.pvalue,pauc.pvalue=pauc.pvalue,stringsAsFactors = F)
  #return(res)
}
png(filename = "../result/PRS3ROC.png",type = "cairo")
par(mar=c(6,6,2,1))
plotroc3()
dev.off()

# we used the so-called one-hot encoding (Wan
# et al. 2010) in a subset of analyses, where each genotype is
# coded as a set of three binary variables instead of a number
#tray binary genotype data
genotype_2binary=function(dat=traingenotype)
{
  tmp=length(unique(unlist(dat[,1:100])))
  if (tmp>3) stop("data include non integer")
  newdat=data.frame(matrix(NA,nrow=3*nrow(dat),ncol=ncol(dat)))
  colnames(newdat)=colnames(dat)
  rownames(newdat)[seq(1,nrow(newdat),3)]=paste0(rownames(dat),"_0")
  rownames(newdat)[seq(2,nrow(newdat),3)]=paste0(rownames(dat),"_1")
  rownames(newdat)[seq(3,nrow(newdat),3)]=paste0(rownames(dat),"_2")
  k <- which(dat==0, arr.ind=TRUE)
  tmp1=as.numeric((k[,1]-1)*3+1)
  tmp2=as.numeric(k[,2])
  tmp=as.matrix(data.frame(row=tmp1,col=tmp2))
  newdat[tmp]=0
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=0
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=0
  k <- which(dat==1, arr.ind=TRUE)
  tmp1=as.numeric((k[,1]-1)*3+1)
  tmp2=as.numeric(k[,2])
  tmp=as.matrix(data.frame(row=tmp1,col=tmp2))
  newdat[tmp]=0
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=1
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=0
  k <- which(dat==2, arr.ind=TRUE)
  tmp1=as.numeric((k[,1]-1)*3+1)
  tmp2=as.numeric(k[,2])
  tmp=as.matrix(data.frame(row=tmp1,col=tmp2))
  newdat[tmp]=0
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=0
  tmp[,1]=tmp[,1]+1
  newdat[tmp]=1
  return(newdat)
  
}
traingenotype=genotype_2binary()
testgenotype=genotype_2binary(dat=testgenotype)

#save(traingenotype,testgenotype,trainy,testy,file="../result/PRS_DNN_EA_binary_P05traindat.RData")
save(traingenotype,testgenotype,trainy,testy,file="../result/PRS_DNN_EA_binary_P1traindat.RData")
xtrain = array(as.matrix(t(traingenotype)), dim = c(ncol(traingenotype), nrow(traingenotype), 1))
xtest = array(as.matrix(t(testgenotype)), dim = c(ncol(testgenotype), nrow(testgenotype), 1))

opt="EA"
checkfile="EA_CNN_binary_earlystop_p05_"

allres=NULL

for (n1 in c(4,8,16,32,64))
{
  for (kerneln in c(3,6,9,12,15,20))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001,0.0001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_binary_earlydropoutres_05.RData")

opt="EA"
checkfile="EA_CNN_binary_earlystop_p1_"

allres=NULL

for (n1 in c(4,8,16,32,64))
{
  for (kerneln in c(3,6,9,12,15,18))
  {
    for (batch_size in c(64,256,512))
    {
      for (lr in c(0.001,0.0001))
      {
        res=data.frame(opt=opt,n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr)
        res$auc=testcnn(n1=n1,kerneln=kerneln,batch_size=batch_size,lr=lr,opt=opt,checkfile = checkfile)
        allres=rbind(allres,res)
        gc()
      }
    }
  }
}
save(allres,file="../result/PRS_CNN_binary_earlydropoutres_p1.RData")

# load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GSE81334.RData")
# BEID <- row.names(clinicaltable)[clinicaltable$type=="BE.BE"]
# BE2ID <- row.names(clinicaltable)[clinicaltable$type=="BE.EAC"]
# EACID <- row.names(clinicaltable)[clinicaltable$type=="EAC"]
# 
# 
# BEmethy  <- MEall[,names(MEall)%in% BEID]
# BEmethy2 <- MEall[,names(MEall)%in% BE2ID]
# EAmethy  <- MEall[,names(MEall)%in% EACID]
# 
# alldat=cbind(BEmethy,BEmethy2,EAmethy)
# types=c(rep("BE samples from BE patients",ncol(BEmethy)),rep("BE samples from EAC patients",ncol(BEmethy2)),rep("EAC samples from EAC patients",ncol(EAmethy)))
# allcolors=c("blue","black","red","darkorchid1","limegreen","goldenrod1","brown4","darkseagreen","darkseagreen1")
# allpch=c(15,16,17)
# load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/Yu_paired_BE_EA_output.RData")
# 
# sum(Yuout[,1] < 0.05/nrow(Yuout) |  Yuout[,2] < 0.05/nrow(Yuout))
# probes=rownames(Yuout)[Yuout[,1] < 0.05/nrow(Yuout) |  Yuout[,2] < 0.05/nrow(Yuout)]
# probes=rownames(Yuout)[Yuout[,2] < 0.05/nrow(Yuout)]
# plot.tsne=function(dat=alldat,probes=NULL,yinc=1.2,pc1=1,pc2=2,main="",opt="beta")
# {
#   if (!is.null(probes)) dat=dat[rownames(dat) %in% probes,]
#   types=factor(types,levels = unique(types))
#   row.means=rowMeans(dat,na.rm=T)
#   for (i in 1:ncol(dat))
#   {
#     idx=which(is.na(dat[,i]))
#     if (length(idx)>0)
#     {
#       dat[idx,i]=row.means[idx]
#     }
#   }
#   dat0=dat
#   if (opt!="beta")
#   {
#     dat[dat <= 0]=min(dat0[dat0!=0])/10
#     dat[dat >= 1]=max(dat0[dat0!=1])+(1-max(dat0[dat0!=1]))*0.5
#     #M values
#     dat=log(dat/(1-dat),base=2)
#   }
#   library(matrixStats)
#   #remove constant rows, var=0
#   tmp=rowSds(as.matrix(dat))
#   dat=dat[tmp>0,]
#   pcadat=prcomp(t(dat),scale = F)
#   # expl_var <- pcadat$sdev^2/sum(pcadat$sdev^2)
#   # barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
#   #         names.arg=paste0("pcadat",seq(1:50)), col="darkgreen")
#   pcadat=pcadat$x
#   library(Rtsne)
#   set.seed(1000)
#   n=as.integer(sqrt(nrow(pcadat)))
#   tsnedat = Rtsne(pcadat, pca = FALSE,perplexity = 5,theta=0,initial_dims = 500, normalize=F)
#   tsnedat=tsnedat$Y
#   rownames(tsnedat)=rownames(pcadat)
#   xmin=min(tsnedat[,pc1])
#   xmax=max(tsnedat[,pc1])
#   ymin=min(tsnedat[,pc2])
#   ymax=max(tsnedat[,pc2])*yinc
#   par(mar=c(5,5,2,1))
# 
#   pch=allpch[types]
#   library(scales)
#   
#   plot(tsnedat[,pc1],tsnedat[,pc2],col=alpha(allcolors[types],0.4),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
#        ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("TSNE",pc1),ylab=paste0("TSNE",pc2),main=main,bty='l',cex=1.3)
#   #legend("topleft",legend=c(Cancertypes),col=allcolors[1:length(Cancertypes)],pch=1,cex=1,bty = "n",ncol=3)
#   legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=allpch,cex=1.3,bty = "n",ncol=1)
# }
# 
# plot.pca=function(dat=alldat,probes=NULL,yinc=1.2,pc1=1,pc2=2,main="",prefix=NULL,opt="beta")
# {
#   types=factor(types,levels = unique(types))
#   if (!is.null(probes)) dat=dat[rownames(dat) %in% probes,]
#   row.means=rowMeans(dat,na.rm=T)
#   for (i in 1:ncol(dat))
#   {
#     idx=which(is.na(dat[,i]))
#     if (length(idx)>0)
#     {
#       dat[idx,i]=row.means[idx]
#     }
#   }
#   library(matrixStats)
#   dat0=dat
#   if (opt!="beta")
#   {
#     dat[dat <= 0]=min(dat0[dat0!=0])/10
#     dat[dat >= 1]=max(dat0[dat0!=1])+(1-max(dat0[dat0!=1]))*0.5
#     #M values
#     dat=log(dat/(1-dat),base=2)
#   }
#   #remove constant rows, var=0
#   tmp=rowSds(as.matrix(dat))
#   dat=dat[tmp>0,]
#   pcadat=prcomp(t(dat),scale = F)
#   expl_var <- pcadat$sdev^2/sum(pcadat$sdev^2)
#   # if (!is.null(prefix))
#   # {
#   #   png(paste0(resfolder,prefix,"PCs_variances.png"),width = 480, height = 480,type = "cairo")
#   #   par(mar=c(6,6,2,1))
#   #   barplot(expl_var[1:50], ylab="EXPLAINED VARIANCE",
#   #           names.arg=paste0("pcadat",seq(1:50)), col="darkgreen",las=2)
#   #   dev.off()
#   # }
#   pcadat=pcadat$x
#   xmin=min(pcadat[,pc1])
#   xmax=max(pcadat[,pc1])
#   ymin=min(pcadat[,pc2])
#   ymax=max(pcadat[,pc2])*yinc
#   pch=allpch[types]
#   plot(pcadat[,pc1],pcadat[,pc2],col=alpha(allcolors[types],0.4),pch=pch,cex.axis=1.5,cex.lab=1.5,xlim=c(xmin-0.05*(xmax-xmin),xmax+0.05*(xmax-xmin)),
#        ylim=c(ymin-0.05*(ymax-ymin),ymax+0.2*(ymax-ymin)),xlab=paste0("PC",pc1," (",round(expl_var[pc1]*100),"%)"),ylab=paste0("PC",pc2," (",round(expl_var[pc2]*100),"%)"),main=main,bty='l',cex=1.3)
#   #text(x=pcadat[,1],y=pcadat[,2],rownames(pcadat),col=colors[clinicaltable$`Platium response`],cex=1)
#   legend("topleft",legend=c(unique(as.character(types))),col=allcolors[1:length(unique(types))],pch=allpch,cex=1.3,bty = "n",ncol=1)
# }
# plot.pca(probes=probes,pc2=2)
