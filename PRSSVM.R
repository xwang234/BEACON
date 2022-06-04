#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

load("../result/PRS_DNN_EA_P1traindat.RData")

library(e1071)
opt=args[1]
print(opt)

outfile=paste0("../result/","PRSSVMP1res_",opt,".RData")
if (opt==1)
{
  tune.out=tune(svm ,t(traingenotype),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=c(0.1 ,1)))
}
if (opt==2)
{
  tune.out=tune(svm ,t(traingenotype),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=c(10,20)))
}
#summary(tune.out)
print(paste0(opt," tune done"))
save(tune.out,file=outfile)
load(outfile)
bestmod =tune.out$best.model
#summary(bestmod)
svmfit <- svm(t(traingenotype0),factor(trainy),scale=F,kernel ="radial", cost=bestmod$cost,gamma = bestmod$gamma,probability=T)
save(tune.out,svmfit,file=outfile)
svm.predictions=predict(svmfit, t(testgenotype0))
save(tune.out,svmfit,svm.predictions,file=outfile)
print(paste0(opt," done!"))