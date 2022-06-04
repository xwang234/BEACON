#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

load("../result/PRS_DNN_EA_P05traindat.RData")

library(e1071)
opt=args[1]
print(opt)

outfile=paste0("../result/","PRSSVMgammares_",opt,".RData")
if (opt==1)
{
  tune.out=tune(svm ,t(traingenotype0),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=10,gamma=c(1/10000,1/1000)))
}
if (opt==2)
{
  tune.out=tune(svm ,t(traingenotype0),factor(trainy),scale=F,  kernel ="radial", ranges =list(cost=10,gamma=c(1/100,1/10)))
}
#summary(tune.out)
print(paste0(opt," tune done"))
save(tune.out,file=outfile)

bestmod =tune.out$best.model
#summary(bestmod)
svmfit <- svm(t(traingenotype0),factor(trainy),scale=F,kernel ="radial", cost=bestmod$cost,gamma = bestmod$gamma,probability=T)
save(tune.out,svmfit,file=outfile)
svm.predictions=predict(svmfit, t(testgenotype0))
save(tune.out,svmfit,svm.predictions,file=outfile)
print(paste0(opt," done!"))