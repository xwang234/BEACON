#try to use DNN to compute PRS

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
    verbose = 1
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
casesamples=sampletable$localid[sampletable$phenoEA_bca==2]
contsamples=sampletable$localid[sampletable$phenoEA_bca==1]
set.seed(1000)
idx=sample(length(casesamples),round(length(casesamples)*0.7))
traincasesamples=casesamples[idx]
idx=sample(length(contsamples),round(length(contsamples)*0.7))
traincontsamples=contsamples[idx]
testcasesamples=casesamples[!casesamples %in% traincasesamples]
testcontsamples=contsamples[!contsamples %in% traincontsamples]
trainsamples=c(traincasesamples,traincontsamples)
testsamples=c(testcasesamples,testcontsamples)
idx=match(trainsamples,sampletable$localid)
tmp=sampletable[idx,c(1,2)]
write.table(tmp,file="../result/PRSDNN_train_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)
idx=match(testsamples,sampletable$localid)
tmp=sampletable[idx,c(1,2)]
write.table(tmp,file="../result/PRSDNN_test_plinksamples.txt",col.names = F,row.names = F,sep="\t",quote = F)

#generate genotype files in bash
plink=/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
$plink  --bfile $prefix --keep ../result/PRSDNN_train_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNN_BCAtrain
$plink  --bfile $prefix --keep ../result/PRSDNN_train_plinksamples.txt \
--make-bed --out ../result/PRSDNN_BCAtrain
prefix="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/bca_filtered_30Nov2018_flip_noambiguous"
$plink  --bfile $prefix --keep ../result/PRSDNN_test_plinksamples.txt \
--recode A-transpose --out ../result/PRSDNN_BCAtest
$plink  --bfile $prefix --keep ../result/PRSDNN_test_plinksamples.txt \
--make-bed --out ../result/PRSDNN_BCAtest

#use EA phenotype
updatefam=function(famfile="../result/PRSDNN_BCAtrain.fam")
{
  fam=read.table(famfile)
  idx=match(fam$V2,sampletable$localid)
  fam$V6=sampletable$phenoEA_bca[idx]
  write.table(fam,file=famfile,row.names = F,col.names = F,sep="\t",quote=F)
}
updatefam()
updatefam(famfile="../result/PRSDNN_BCAtest.fam")

#run gwas in bash, don't consider covariate
prefix=/fh/fast/dai_j/BEACON/BEACON_GRANT/result/PRSDNN_BCAtrain
$plink --bfile  $prefix --logistic --beta --hide-covar --ci 0.95 --out $prefix

#read gwas result
gwas=read.table("../result/PRSDNN_BCAtrain.assoc.logistic",header = T)
quantile(gwas$P,na.rm=T)
# 0%       25%       50%       75%      100% 
# 8.627e-08 2.217e-01 4.733e-01 7.359e-01 1.000e+00 
sum(gwas$P<0.05,na.rm=T)
pcutoff=0.001
pcutoff=0.01
selectsnps=gwas$SNP[which(gwas$P<pcutoff)] #1742,12446

#read genotype
readgenotype=function(prefix="../result/PRSDNN_BCAtrain")
{
  tmp=data.table::fread(paste0(prefix,".traw"))
  rownames(tmp)=tmp$SNP
  idx=match(selectsnps,rownames(tmp))
  if (sum(is.na(idx))>0) warning("some snps are missing!")
  genotype=tmp[idx,7:ncol(tmp)]
  k <- which(is.na(genotype), arr.ind=TRUE)
  length(k)/nrow(genotype)/ncol(genotype)
  genotype[k] <- rowMeans(genotype, na.rm=TRUE)[k[,1]]
  
  return(genotype)
}
tolocalid=function(samplenames)
{
  samplenames=strsplit(samplenames,"_")
  samplenames=sapply(1:length(samplenames),function(x){
    tmp=samplenames[[x]]
    paste0(tmp[2:length(tmp)],collapse = "_")
  })
}

traingenotype=readgenotype()
colnames(traingenotype)=tolocalid(samplenames = colnames(traingenotype))
idx=match(colnames(traingenotype),sampletable$localid)
trainy=as.numeric(I(sampletable$phenoEA_bca[idx]==2))
testgenotype=readgenotype(prefix="../result/PRSDNN_BCAtest")
colnames(testgenotype)=tolocalid(samplenames = colnames(testgenotype))
idx=match(colnames(testgenotype),sampletable$localid)
testy=as.numeric(I(sampletable$phenoEA_bca[idx]==2))
model <- keras_model_sequential() %>% 
  layer_dense(units = nrow(traingenotype),input_shape=c(nrow(traingenotype))) %>% 
  layer_dense(units = 500, activation = "relu") %>% 
  layer_dense(units = 250, activation = "relu") %>% 
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(0.2) %>% 
  layer_dense(1, activation = "sigmoid")
summary(model)

model %>% 
  compile(
    loss = "binary_crossentropy",
    optimizer=optimizer_adam(lr = 0.001),
    metrics=tf$keras$metrics$AUC()
  )

filepath = "weights_best.hdf5"
check = callback_model_checkpoint(filepath,
                                           monitor='val_auc',  # validation AUC
                                           save_best_only=T,
                                           mode='max')
earlyStopping = callback_early_stopping(monitor="val_auc", patience=0,mode='auto') 

set.seed(1000)
model %>% 
  fit(
    x = t(traingenotype), y = trainy,
    callbacks=list(check, earlyStopping),
    batch_size=512,
    epochs = 200,
    validation_split = 0.3,
    verbose = 3
  )

predictions <- predict(model, t(testgenotype))
names(predictions)=colnames(testgenotype)
#boxplot(predictions~testy)
plot_ROC(predictions)

plot_ROC=function(predictions)
{
  idx=match(names(predictions),sampletable$localid)
  sampletable1=sampletable[idx,]
  sampletable1$prs=predictions
  #for EA
  sampletable1$case=sampletable1$phenoEA_bca
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

#3 hidden layers NN
testDNN3=function(n1=1000,n2=250,n3=50,batch_size=512,lr=0.001)
{
  model <- keras_model_sequential() %>% 
    layer_dense(units = nrow(traingenotype),input_shape=c(nrow(traingenotype))) %>% 
    layer_dense(units = n1, activation = "relu") %>% 
    layer_dense(units = n2, activation = "relu") %>% 
    layer_dense(units = n3, activation = "relu") %>% 
    layer_dropout(0.2) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      metrics = "accuracy"
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traingenotype), y = trainy,
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions)
  return(res)
}

#2 hidden layers
testDNN2=function(n1=1000,n2=250,batch_size=512,lr=0.001)
{
  model <- keras_model_sequential() %>% 
    layer_dense(units = nrow(traingenotype),input_shape=c(nrow(traingenotype))) %>% 
    layer_dense(units = n1, activation = "relu") %>% 
    layer_dense(units = n2, activation = "relu") %>% 
    layer_dropout(0.2) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      metrics = "accuracy"
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traingenotype), y = trainy,
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions)
  return(res)
}

#1 hidden layer
testDNN1=function(n1=1000,batch_size=512,lr=0.001)
{
  model <- keras_model_sequential() %>% 
    layer_dense(units = nrow(traingenotype),input_shape=c(nrow(traingenotype))) %>% 
    layer_dense(units = n1, activation = "relu") %>% 
    layer_dropout(0.2) %>% 
    layer_dense(1, activation = "sigmoid")
  summary(model)
  
  model %>% 
    compile(
      loss = "binary_crossentropy",
      optimizer=optimizer_adam(lr = lr),
      metrics = "accuracy"
      #metrics=tf$keras$metrics$AUC()
    )
  set.seed(1000)
  model %>% 
    fit(
      x = t(traingenotype), y = trainy,
      batch_size=batch_size,
      epochs = 200,
      validation_split = 0.3,
      verbose = 3
    )
  
  predictions <- predict(model, t(testgenotype))
  names(predictions)=colnames(testgenotype)
  #boxplot(predictions~testy)
  res=plot_ROC(predictions)
  return(res)
}

allres3=NULL
for (n1 in c(1000,800,600,500,200,100))
{
  n2=round(n1/2)
  n3=round(n2/2)
  for (batch_size in seq(32,512,120))
  {
    for (lr in c(0.001,0.0001))
    {
      res=data.frame(n1=n1,n2=n2,n3=n3,batch_size=batch_size,lr=lr)
      res$auc=testDNN3(n1,n2,n3,batch_size,lr)
      allres3=rbind(allres3,res)
    }
  }
}

allres2=NULL
for (n1 in c(1000,800,600,500,200,100))
{
  n2=round(n1/2)
  for (batch_size in seq(32,512,120))
  {
    for (lr in c(0.001,0.0001))
    {
      res=data.frame(n1=n1,n2=n2,batch_size=batch_size,lr=lr)
      res$auc=testDNN2(n1,n2,batch_size,lr)
      allres2=rbind(allres2,res)
    }
  }
}

allres1=NULL
for (n1 in c(1000,800,600,500,200,100))
{
  for (batch_size in seq(32,512,120))
  {
    for (lr in c(0.001,0.0001))
    {
      res=data.frame(n1=n1,batch_size=batch_size,lr=lr)
      res$auc=testDNN1(n1,batch_size,lr)
      allres1=rbind(allres1,res)
    }
  }
}
#PRS_DNN_res.RData for pvalue=0.001
#PRS_DNN_res_01.RData for pvalue=0.01,without early stopping

#save(allres3,allres2,allres1,file="../result/PRS_DNN_res.RData")
save(allres3,allres2,allres1,file="../result/PRS_DNN_res_01.RData")
#to implement early stopping https://stackoverflow.com/questions/63150719/early-stopping-based-on-auc
#https://stackoverflow.com/questions/58682098/keras-callback-modelcheckpoint-doesnt-save-weights?rq=1
#https://community.rstudio.com/t/use-auc-as-metric-in-keras-for-r/84573/2
#https://cran.r-project.org/web/packages/keras/vignettes/training_callbacks.html