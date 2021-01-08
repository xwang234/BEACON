#!/usr/bin/env Rscript
library(data.table)

ME=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GSE72874/GSE72874-GPL13534_series_matrix.txt",skip=63,sep="\t",header = T))
rownames(ME)=ME$ID_REF
ME=ME[,-1]
for (i in 1:ncol(ME)) ME[,i]=as.numeric(ME[,i])
sum(is.na(ME))/nrow(ME)/ncol(ME) #[1] 0.0005562037
idx=grepl("^cg",rownames(ME))
table(idx) #non-CpG probes
# FALSE   TRUE 
# 1185 372376 
ME=ME[idx,]
MEall=ME

#use subset of probes
EACprobes=read.csv("../data/supp_bgw018_Copy_of_Supplementary_Table3.csv",skip = 1,header = T,stringsAsFactors = F)
BEprobes=read.csv("../data/supp_bgw018_Copy_of_Supplementary_Table4.csv",skip = 1,header = T,stringsAsFactors = F)
selectedprobes=intersect(EACprobes$TargetID,BEprobes$TargetID)
length(selectedprobes) #39191
#selectedprobes=unique(c(EACprobes$TargetID,BEprobes$TargetID))
length(selectedprobes) #63500
#use all probes
#selectedprobes=rownames(MEall)

ME=MEall[rownames(MEall) %in% selectedprobes,]


MEIDtable=read.table("../data/GSE72874/methylationID.txt",stringsAsFactors = F)
library(xlsx)
clinicaltable=read.xlsx("../data/Copy of Copy_of_Supplementary_Table2.xlsx",sheetIndex = 1,startRow = 2,header = T)
idx=match(clinicaltable$Methylation.ArrayID,MEIDtable$V3)
clinicaltable$ID=MEIDtable$V1[idx]
# The clinical features for the 250 samples, which includes 19 BE
# (11 with synchronous EAC and 8 from non-EAC patients), 125
# EAC, 85 NSE (11 squamous esophagus from healthy individuals,
#              10 squamous esophagus from patients with gastroesophageal
#              reflux disease and 64 adjacent squamous esophagus from EAC
#              or BE patients) and 21 stomach samples,
table(clinicaltable$Sample.type,clinicaltable$Freq.of.Heartburn.reflux)
#           . <Once/month Daily Monthly (few times/mo) Never Weekly (few times/wk)
# BE        7           1     6                      0     2                     3
# Control   4           1     0                      6     0                     0
# GERD      0           0     3                      3     1                     3
# Normal   50           1     7                      0     2                     4
# Stomach  21           0     0                      0     0                     0
# Tumour  111           1     3                      3     6                     1
table(clinicaltable$Sample.type,clinicaltable$BMI.Category)
#           . Healthy wt (-<25) Obese I (-<35) Obese II (-<40) Obese III (>=40) Overweight(-<30)
# BE        8                 3              0               2                1                5
# Control   2                 5              2               0                1                1
# GERD      0                 5              0               0                1                4
# Normal   51                 4              1               1                1                6
# Stomach  21                 0              0               0                0                0
# Tumour  113                 2              5               0                0                5
table(clinicaltable$Sample.type,clinicaltable$Average.Std.Alc.Drinks.Categories)
#          . High Low
# BE       4    6   9
# Control  3    4   4
# GERD     0    8   2
# Normal  21   16  27
# Stomach  0    5  16
# Tumour  34   34  57

BEEACsamples=clinicaltable$ID[clinicaltable$Sample.type %in% c("BE","Tumour")]
#nontumor squamous esophagus
#NSEsamples=clinicaltable$ID[clinicaltable$Sample.type %in% c("Control","GERD","Normal")]

clinicaltable$Array.batch=factor(clinicaltable$Array.batch,levels = c(1,2,3,4))
idx=match(BEEACsamples,colnames(ME))
BEEACME=ME[,idx]
sum(is.na(BEEACME))/nrow(BEEACME)/ncol(BEEACME)
BEEACMEall=MEall[,idx]
idx=match(BEEACsamples,clinicaltable$ID)
BEEACclinical=clinicaltable[idx,]
BEEACclinical$Gender[!BEEACclinical$Gender %in% c("F","M")]=NA
BEEACclinical$Gender=factor(BEEACclinical$Gender,levels = c("F","M"))
BEEACclinical$Age=as.numeric(BEEACclinical$Age)
table(is.na(BEEACclinical$Age))
# FALSE
#  144

BEEACclinical$Sample.type=factor(BEEACclinical$Sample.type,levels=c("BE","Tumour"))
table(BEEACclinical$BMI.Category)
# . Healthy wt (-<25)    Obese I (-<35)   Obese II (-<40)  Obese III (>=40)  Overweight(-<30) 
# 121                 5                 5                 2                 1                10 
BEEACclinical$BMI=NA
BEEACclinical$BMI[BEEACclinical$BMI.Category %in% "Healthy wt (-<25)"]="Healthy"
BEEACclinical$BMI[BEEACclinical$BMI.Category %in% "Overweight(-<30)"]="Overweight"
BEEACclinical$BMI[BEEACclinical$BMI.Category %in% c("Obese I (-<35)","Obese II (-<40)","Obese III (>=40)")]="Obese"
BEEACclinical$BMI=factor(BEEACclinical$BMI,levels = c("Healthy","Overweight","Obese"))
table(BEEACclinical$BMI)
# Healthy Overweight      Obese 
#       5         10          8 

table(BEEACclinical$Smoking.status)
# . Current smoker Current Smoker      Ex-Smoker   Never Smoker 
# 41              2              8             69             24 
BEEACclinical$Smoking=NA
BEEACclinical$Smoking[BEEACclinical$Smoking.status %in% c("Current smoker","Current Smoker")]="CSmoker"
BEEACclinical$Smoking[BEEACclinical$Smoking.status %in% "Ex-Smoker"]="ESmoker"
BEEACclinical$Smoking[BEEACclinical$Smoking.status %in% "Never Smoker"]="NSmoker"
BEEACclinical$Smoking=factor(BEEACclinical$Smoking,levels = c("NSmoker","ESmoker","CSmoker"))
table(BEEACclinical$Smoking)
# NSmoker ESmoker CSmoker 
# 24      69      10
BEEACclinical$Smoking1=as.character(BEEACclinical$Smoking)
BEEACclinical$Smoking1[as.character(BEEACclinical$Smoking) %in% c("ESmoker","CSmoker")]="Smoker"
BEEACclinical$Smoking1=factor(BEEACclinical$Smoking1,levels = c("NSmoker","Smoker"))
table(BEEACclinical$Smoking1)
# NSmoker  Smoker 
# 24      79  

table(BEEACclinical$Average.Std.Alc.Drinks.Categories)
BEEACclinical$Drink=BEEACclinical$Average.Std.Alc.Drinks.Categories
BEEACclinical$Drink[BEEACclinical$Drink=="."]=NA
BEEACclinical$Drink=factor(BEEACclinical$Drink,levels = c("Low","High"))
table(BEEACclinical$Drink)
# Low High 
# 66   40

table(BEEACclinical$Freq.of.Heartburn.reflux)
# .            <Once/month                  Daily Monthly (few times/mo)                  Never 
# 118                      2                      9                      3                      8 
# Weekly (few times/wk) 
# 4 
BEEACclinical$Freq.of.Heartburn.reflux=as.character(BEEACclinical$Freq.of.Heartburn.reflux)
BEEACclinical$reflux=BEEACclinical$Freq.of.Heartburn.reflux
BEEACclinical$reflux[BEEACclinical$reflux=="."]=NA
BEEACclinical$reflux[BEEACclinical$reflux %in% c("<Once/month","Never")]="Few"
BEEACclinical$reflux[BEEACclinical$reflux=="Monthly (few times/mo)"]="Monthly"
BEEACclinical$reflux[BEEACclinical$reflux=="Weekly (few times/wk)"]="Weekly"
BEEACclinical$reflux=factor(BEEACclinical$reflux,levels=c("Few","Monthly","Weekly","Daily"))
table(BEEACclinical$reflux)
# Few Monthly  Weekly   Daily 
# 10       3       4       9 
BEEACclinical$reflux1=as.character(BEEACclinical$reflux)
BEEACclinical$reflux1[BEEACclinical$reflux1 %in% c("Few","Monthly","Weekly")]="Other"
BEEACclinical$reflux1=factor(BEEACclinical$reflux1,levels = c("Other","Daily"))
table(BEEACclinical$reflux1)
# Other Daily 
# 17     9

BEEACclinical$reflux2=as.character(BEEACclinical$reflux)
BEEACclinical$reflux2[BEEACclinical$reflux2 %in% c("Few","Monthly")]="Other"
BEEACclinical$reflux2[BEEACclinical$reflux2 %in% c("Daily","Weekly")]="WDaily"
BEEACclinical$reflux2=factor(BEEACclinical$reflux2,levels = c("Other","WDaily"))
table(BEEACclinical$reflux2)
# Other WDaily 
# 13     13

qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

gc()
library(limma)
idx=complete.cases(BEEACclinical[,c("Smoking","Sample.type","Gender","Age")])
designmat=model.matrix(~0+Smoking+Sample.type+Gender+Age,data=BEEACclinical[idx,])
#which(!rownames(BEEACclinical) %in% rownames(designmat))
contrast.matrix=makeContrasts(SmokingCSmoker-SmokingNSmoker,SmokingCSmoker-SmokingESmoker,levels = designmat)
BEEACME1=as.matrix(BEEACME[,idx])
fit <- lmFit(BEEACME1, designmat)
fit1 <- contrasts.fit(fit, contrast.matrix)
result=eBayes(fit1,trend = T)
colnames(result$p.value)
# [1] "SmokingCSmoker - SmokingNSmoker" "SmokingCSmoker - SmokingESmoker"
#qqplot(result$p.value[,1])
#qqplot(result$p.value[,2])

gc()
#compare smoker vis non-smoker
idx=complete.cases(BEEACclinical[,c("Smoking1","Sample.type","Gender","Age")])
designmat=model.matrix(~0+Smoking1+Sample.type+Gender+Age,data=BEEACclinical[idx,])
#which(!rownames(BEEACclinical) %in% rownames(designmat))
contrast.matrix=makeContrasts(Smoking1Smoker-Smoking1NSmoker,levels = designmat)
BEEACME1=as.matrix(BEEACME[,idx])
fit <- lmFit(BEEACME1, designmat)
fit1 <- contrasts.fit(fit, contrast.matrix)
result1=eBayes(fit1,trend = T)
#qqplot(result1$p.value[,1])

#reflux
idx=complete.cases(BEEACclinical[,c("reflux1","Sample.type","Gender","Age")])
designmat=model.matrix(~0+reflux1+Sample.type+Gender+Age,data=BEEACclinical[idx,])
contrast.matrix=makeContrasts(reflux1Daily-reflux1Other,levels = designmat)
BEEACME1=as.matrix(BEEACME[,idx])
fit <- lmFit(BEEACME1, designmat)
fit1 <- contrasts.fit(fit, contrast.matrix)
result2=eBayes(fit1,trend = T)
#qqplot(result2$p.value[,1])

idx=complete.cases(BEEACclinical[,c("reflux2","Sample.type","Gender","Age")])
designmat=model.matrix(~0+reflux2+Sample.type+Gender+Age,data=BEEACclinical[idx,])
contrast.matrix=makeContrasts(reflux2WDaily-reflux2Other,levels = designmat)
BEEACME1=as.matrix(BEEACME[,idx])
fit <- lmFit(BEEACME1, designmat)
fit1 <- contrasts.fit(fit, contrast.matrix)
result3=eBayes(fit1,trend = T)
#qqplot(result3$p.value[,1])

#another way, including intercept in designmantrix, result is the same as above
# designmat=model.matrix(~reflux2+Sample.type+Gender+Age,data=BEEACclinical[idx,])
# BEEACME1=as.matrix(BEEACME[,idx])
# fit <- lmFit(BEEACME1, designmat)
# result4=eBayes(fit,trend = T)
# qqplot(result4$p.value[,2])

#drink 
idx=complete.cases(BEEACclinical[,c("Drink","Sample.type","Gender","Age")])
designmat=model.matrix(~Drink+Sample.type+Gender+Age,data=BEEACclinical[idx,])
BEEACME1=as.matrix(BEEACME[,idx])
fit <- lmFit(BEEACME1, designmat)
result5=eBayes(fit,trend = T,robust = T)
#qqplot(result5$p.value[,2])

# #remove paired EAC samples, doesn't change the result much
# tmp=BEEACclinical$Donor[duplicated(BEEACclinical$Donor)]
# idx=which(BEEACclinical$Donor %in% tmp)
# BEEACclinical$Donor=as.character(BEEACclinical$Donor)
# table(BEEACclinical$Donor[idx],as.character(BEEACclinical$Sample.type[idx]))
# #         BE Tumour
# # 40320   1      1
# # 40341   1      1
# # SOG093  1      1
# # SOG115  1      1
# # SOG132  1      1
# table(BEEACclinical$Donor[idx],as.character(BEEACclinical$reflux[idx]))
# 
# idx=BEEACclinical$Donor %in% tmp & BEEACclinical$Sample.type=="Tumour"
# idx=which(idx==F)
# BEEACME=BEEACME[,idx]
# BEEACclinical=BEEACclinical[idx,]
# all(colnames(BEEACME)==BEEACclinical$ID) #T

library(bumphunter)
library(doParallel)
registerDoParallel(cores = 15)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr)  
anno <- anno[anno$chr!="X"&anno$chr!="Y",]
anno$chr <- as.numeric(anno$chr)
anno$pos <- as.numeric(anno$pos)
#anno=anno[rownames(anno) %in% rownames(BEEACME),]
#subsetprobe
BEEACME=BEEACME[rownames(BEEACME) %in% anno$Name,]
idx=match(rownames(BEEACME),rownames(anno))
anno1=anno[idx,]
#all probes
idx=match(rownames(BEEACMEall),rownames(anno))
anno3=anno[idx,]


idx=complete.cases(BEEACclinical[,c("Smoking1","Sample.type","Gender","Age")])
designmat=model.matrix(~Smoking1+Sample.type+Gender+Age,data=BEEACclinical[idx,])
BEEACME1=as.matrix(BEEACME[,idx])
idx1=order(anno1$chr,anno1$pos)
set.seed(1000)
dmrs_smoking <- bumphunter(BEEACME1[idx1,],designmat,chr=anno1$chr[idx1],pos=anno1$pos[idx1],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(dmrs_smoking$table)
dmrs_smoking$table[1:10,]
#      chr     start       end      value      area cluster indexStart indexEnd L clusterL     p.value  fwer p.valueArea fwerArea
# 209  16   1202468   1202468 -0.2024874 0.2024874   18172      31734    31734 1        1 0.004459211 0.189  0.22974695    0.967
# 12   12  27615138  27615138  0.1855156 0.1855156   15095      26340    26340 1        1 0.014770896 0.391  0.23595157    0.968
# 110   7   8473279   8473462 -0.1507622 0.4522866    9352      16236    16239 4        4 0.022583127 0.410  0.06018212    0.708
# 87    5 156886970 156886996 -0.1432936 0.5731745    7462      12741    12744 4        5 0.032221146 0.501  0.03797410    0.579
# 71    4  54229634  54229634 -0.1749820 0.1749820    5831       9668     9668 1        1 0.031256579 0.581  0.24576949    0.969
# 216  17  47210368  47210368 -0.1736077 0.1736077   19614      34109    34109 1        1 0.034330181 0.602  0.24754169    0.970
# 109   7   1710378   1710516 -0.1555133 0.3110266    9238      16050    16051 2        5 0.047252322 0.617  0.11420938    0.858
# 161  11   1411609   1411609 -0.1707718 0.1707718   13521      23597    23597 1        1 0.042326139 0.647  0.25244491    0.970
# 67    3 183993721 183993726 -0.1403335 0.4210005    5321       8769     8771 3        6 0.063845179 0.708  0.07350234    0.761
# 194  14  74462914  74462914 -0.1655908 0.1655908   16982      29777    29777 1        1 0.060936166 0.733  0.26367523    0.973

#qqplot(dmrs_smoking$table$p.value)
BEEACMEall1=as.matrix(BEEACMEall[,idx])
idx1=order(anno3$chr,anno3$pos)
set.seed(1000)
dmrsall_smoking <- bumphunter(BEEACMEall1[idx1,],designmat,chr=anno3$chr[idx1],pos=anno3$pos[idx1],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(dmrsall_smoking$table)
dmrsall_smoking$table[1:10,]
#       chr     start       end      value      area cluster indexStart indexEnd  L clusterL      p.value  fwer p.valueArea fwerArea
# 2679  12 132865307 132865307 -0.2073801 0.2073801  119443     256325   256325  1        5 0.0002485345 0.227 0.139799144    1.000
# 663    2 113103184 113103184 -0.2073415 0.2073415   23885      50507    50507  1        1 0.0002492969 0.227 0.139878813    1.000
# 3178  17  47209827  47210822 -0.1612049 0.9672292  148483     319802   319807  6       12 0.0006914749 0.253 0.004846804    0.839
# 1916   8 121823462 121824742 -0.1298936 1.4288296   84715     182940   182951 12       13 0.0005359502 0.287 0.001204173    0.475
# 1729   7 158512380 158512671 -0.1642741 0.4928223   78681     170694   170697  4        4 0.0009537322 0.328 0.025552170    0.997
# 1007   4 111557894 111558602 -0.1537362 0.9224173   46059      95407    95412  6       10 0.0011492816 0.364 0.005568012    0.869
# 514    1 216896852 216896863 -0.1655439 0.4966316   15095      32816    32818  3        6 0.0012209449 0.398 0.025180130    0.997
# 2994  16   1202468   1202923 -0.1651237 0.4953712  136181     291582   291584  3       14 0.0012647815 0.406 0.025306685    0.997
# 2225  10 131843597 131843798 -0.1643426 0.4930277   99886     212280   212282  3        5 0.0013471181 0.415 0.025535397    0.997
# 1157   5   5025803   5025848 -0.1676570 0.3353140   49657     103038   103039  2        2 0.0015506725 0.468 0.054998936    1.000

idx=complete.cases(BEEACclinical[,c("reflux1","Sample.type","Gender","Age")])
designmat=model.matrix(~reflux1+Sample.type+Gender+Age,data=BEEACclinical[idx,])
BEEACME1=as.matrix(BEEACME[,idx])
idx1=order(anno1$chr,anno1$pos)
set.seed(1000)
dmrs_reflux <- bumphunter(BEEACME1[idx1,],designmat,chr=anno1$chr[idx1],pos=anno1$pos[idx1],pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(dmrs_reflux$table)
dmrs_reflux$table[1:10,]
#      chr     start       end     value      area cluster indexStart indexEnd L clusterL     p.value  fwer p.valueArea fwerArea
# 431   6 152957910 152957910 0.3867704 0.3867704    8973      15584    15584 1        1 0.002090800 0.143 0.212923849    0.886
# 28    1 108507078 108508207 0.2678323 2.1426585    1202       1948     1955 8        9 0.002851738 0.178 0.004089151    0.232
# 676  12  54447220  54448090 0.2561394 2.3052544   15311      26710    26718 9       14 0.002403709 0.188 0.003121978    0.208
# 374   6  17281483  17282533 0.2899471 1.7396826    7902      13502    13507 6       10 0.003363771 0.190 0.007893839    0.329
# 675  12  54403254  54403314 0.3463209 0.6926419   15300      26693    26694 2        6 0.002922854 0.192 0.080438926    0.779
# 674  12  54402699  54402888 0.3196870 0.9590611   15300      26689    26691 3        6 0.004473175 0.207 0.037335723    0.630
# 801  16  50875140  50875140 0.3617031 0.3617031   18510      32244    32244 1        1 0.005426125 0.230 0.215249332    0.886
# 656  12   4381792   4382188 0.2700836 1.8905853   14902      26017    26023 7        7 0.004128264 0.233 0.006155061    0.290
# 246   4  96470584  96471143 0.2773819 1.6642912    6020      10020    10025 6       17 0.005287450 0.252 0.009049468    0.348
# 40    1 158083159 158083475 0.2788055 1.3940277    1513       2434     2438 5        5 0.007936508 0.282 0.014450703    0.422
#qqplot(dmrs_reflux$table$p.value)
BEEACMEall1=as.matrix(BEEACMEall[,idx])
idx1=order(anno3$chr,anno3$pos)
set.seed(1000)
dmrsall_reflux <- bumphunter(BEEACMEall1[idx1,],designmat,chr=anno3$chr[idx1],pos=anno3$pos[idx1],pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(dmrsall_reflux$table)
dmrsall_reflux$table[1:10,]
#       chr     start       end     value      area cluster indexStart indexEnd  L clusterL      p.value  fwer  p.valueArea fwerArea
# 567    2 177052486 177053503 0.2683930 3.4891091   26768      56070    56082 13       14 0.0001608871 0.161 0.0003787552    0.284
# 1173   4 184826215 184827754 0.2621517 3.4079717   48319      99854    99866 13       14 0.0001899362 0.179 0.0004282875    0.305
# 1717   6 152957910 152958548 0.2843527 2.5591747   66744     144424   144432  9       16 0.0003679549 0.217 0.0015861536    0.561
# 2472  10 118032626 118034031 0.2195976 3.2939633   98701     209639   209653 15       28 0.0002487792 0.255 0.0005132002    0.343
# 2853  12  54402443  54403314 0.2889193 2.0224353  114307     245172   245178  7        8 0.0006040716 0.270 0.0036039465    0.707
# 1509   6  17281327  17282533 0.2648637 2.6486373   59473     122956   122965 10       20 0.0005448562 0.298 0.0013872792    0.535
# 642    2 238535754 238535997 0.3279577 0.6559155   30170      62651    62652  2       11 0.0008167257 0.315 0.0526741529    0.987
# 159    1 108507078 108508548 0.2622515 2.3602635    9321      20325    20333  9       10 0.0008245466 0.341 0.0021462792    0.625
# 3459  16  50875140  50875140 0.3617031 0.3617031  139548     299071   299071  1        4 0.0007638415 0.364 0.1731514365    0.997
# 2776  12   4381777   4381894 0.2628086 2.1024685  111327     238619   238626  8       49 0.0011258376 0.385 0.0031931629    0.692
#save(BEEACMEall,BEEACME,BEEACclinical,file="../data/GSE72874.RData")

#work on Grady's data:---------------------------------------------------------------
#to get the clinical info of an GEO:
#https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE104707
readGSE=function(GSEfile="../data/GSE89181/GSE89181-GPL13534_series_matrix.txt")
{
  cmd=paste0("grep -n ","Sample_title"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  studyno=unlist(strsplit(tmp1,"\t\""))
  studyno=gsub("\"","",studyno)
  studyno=studyno[2:length(studyno)]
  cmd=paste0("grep -n ","Sample_geo_accession"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  sampleid=unlist(strsplit(tmp1,"\t\""))
  sampleid=gsub("\"","",sampleid)
  sampleid=sampleid[2:length(sampleid)]
  cmd=paste0("grep -n ","ID_REF"," ",GSEfile)
  tmp1=system(cmd,intern = TRUE)
  tmp1=as.integer(unlist(strsplit(tmp1,":"))[1])
  GSE=as.data.frame(fread(GSEfile,skip=tmp1-1,header = T,sep="\t",fill=T))
  idx=grepl("^cg",GSE[,1])
  GSE=GSE[idx,]
  rownames(GSE)=GSE[,1]
  GSE=GSE[,-1]
  if(all(colnames(GSE)==sampleid)) colnames(GSE)=studyno
  # histology=NA
  # cmd=paste0("grep -n ","histology:"," ",GSEfile)
  # tmp1=system(cmd,intern = TRUE)
  # tmp1=tmp1 [length(tmp1)]
  # if (!is.na(tmp1))
  # {
  #   histology=unlist(strsplit(tmp1,"\t\""))
  #   histology=gsub("\"","",histology)
  #   histology=histology[2:length(histology)]
  #   histology=gsub("histology:","",histology)
  #   histology=gsub(" ","",histology)
  # }
  # age=NA
  # cmd=paste0("grep -n ","age:"," ",GSEfile)
  # tmp1=system(cmd,intern = TRUE)[length(tmp1)]
  # if (!is.na(tmp1))
  # {
  #   age=unlist(strsplit(tmp1,"\t\""))
  #   age=gsub("\"","",age)
  #   age=age[2:length(age)]
  #   age=gsub("age:","",age)
  #   age=gsub(" ","",age)
  # }
 return(GSE)
}
GSE89181_GPL13534=readGSE()
GSE89181_GPL18809=readGSE(GSEfile = "../data/GSE89181/GSE89181-GPL18809_series_matrix.txt")
all(rownames(GSE89181_GPL13534)==rownames(GSE89181_GPL18809)) #T
sum(colnames(GSE89181_GPL13534) %in% colnames(GSE89181_GPL18809)) #0
table_GPL13534=read.csv("../data/GSE89181/GSE89181_GPL13534.csv")
for (i in 1:ncol(table_GPL13534))
{
  idx=which(table_GPL13534[,i]=="")
  if (length(idx)>0) table_GPL13534[idx,i]=NA
}
table_GPL18809=read.csv("../data/GSE89181/GSE89181_GPL18809.csv")
for (i in 1:ncol(table_GPL18809))
{
  idx=which(table_GPL18809[,i]=="")
  if (length(idx)>0) table_GPL18809[idx,i]=NA
}
selectedsamples1=table_GPL13534$Title[table_GPL13534$Histology %in% c("BE","LGD","HGD","HGD/Ca","EAC")]
idx=match(selectedsamples1,table_GPL13534$Title)
table_GPL13534$Histology=as.character(table_GPL13534$Histology)
GSEtype1=table_GPL13534$Histology[idx]
GSEtype1[which(GSEtype1=="HGD/Ca")]="HGD"
GSEtype1=factor(GSEtype1)
GSEbmi1=table_GPL13534$Bmi[idx]
GSEbmi1[which(GSEbmi1==0)]=NA
GSEage1=table_GPL13534$Age[idx]
GSEgender1=factor(table_GPL13534$Gender[idx])
GSEsmoke1=factor(table_GPL13534$Smokecigarettes[idx])
table_GPL13534$Alcoholicrinksperweek=as.character(table_GPL13534$Alcoholicrinksperweek)
GSEalcohol1=table_GPL13534$Alcoholicrinksperweek[idx]
GSEalcohol1[which(GSEalcohol1=="Unknown")]=NA
GSEalcohol1=factor(GSEalcohol1)
run1=factor(table_GPL13534$Run[idx])

selectedsamples2=table_GPL18809$Title[table_GPL18809$Histology %in% c("BE","LGD","HGD","HGD/Ca","EAC")]
table_GPL18809$Histology=as.character(table_GPL18809$Histology)
idx=match(selectedsamples2,table_GPL18809$Title)
GSEtype2=table_GPL18809$Histology[idx]
GSEtype2[which(GSEtype2=="HGD/Ca")]="HGD"
GSEtype2=as.factor(GSEtype2)
GSEbmi2=table_GPL18809$Bmi[idx]
GSEage2=table_GPL18809$Age[idx]
GSEgender2=factor(table_GPL18809$Gender[idx])
GSEsmoke2=factor(table_GPL18809$Smokecigarettes[idx])
GSEalcohol2=factor(table_GPL18809$Alcoholicrinksperweek[idx])
run2=factor(table_GPL18809$Run[idx])
GSEselectedsamples=c(selectedsamples1,selectedsamples2)
GSEtype=unlist(list(GSEtype1,GSEtype2))
GSEage=c(GSEage1,GSEage2)
GSEbmi=c(GSEbmi1,GSEbmi2)
GSEgender=unlist(list(GSEgender1,GSEgender2))
GSEsmoke=unlist(list(GSEsmoke1,GSEsmoke2))
GSEalcohol=unlist(list(GSEalcohol1,GSEalcohol2))
GSErun=unlist(list(run1,run2))
idx1=match(selectedsamples1,colnames(GSE89181_GPL13534))
idx2=match(selectedsamples2,colnames(GSE89181_GPL18809))
GSEmeth=cbind.data.frame(GSE89181_GPL13534[,idx1],GSE89181_GPL18809[,idx2])
dim(GSEmeth)
# [1] 453444     81
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
table(rownames(GSEmeth) %in% anno$Name)
# TRUE 
# 453444 
anno$chr <- as.character(anno$chr)
anno$chr <- gsub("chr","",anno$chr)  
anno <- anno[anno$chr!="X"&anno$chr!="Y",]
anno$chr <- as.numeric(anno$chr)
anno$pos <- as.numeric(anno$pos)
anno=as.data.frame(anno)
anno=anno[order(anno$chr,anno$pos),]
table(rownames(GSEmeth) %in% anno$Name)
# FALSE   TRUE 
# 225 453219 
GSEmeth=GSEmeth[rownames(GSEmeth) %in% anno$Name,]
GSEmethall=GSEmeth
#get CpGs close to 26 SNPs
snpposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_26highrisk_SNP_POS.txt"
snppos=read.table(snpposfile,header = T,stringsAsFactors = F)
library(GenomicRanges)
idx=match(rownames(GSEmeth),anno$Name)
gr_GSEmeth=GRanges(seqnames = anno$chr[idx],ranges = IRanges(start=anno$pos[idx],width=1))
gr_snppos=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$position,width=1))
selectedCpGs=NULL
distcutoff=5e5
for (i in 1:length(gr_snppos))
{
  tmp=distance(gr_GSEmeth,gr_snppos[i])
  selectedCpGs=unique(c(selectedCpGs,rownames(GSEmeth)[which(tmp<distcutoff)]))
}
idx=match(selectedCpGs,rownames(GSEmeth))
GSEmeth=GSEmeth[idx,]
dim(GSEmeth)
# [1] 7017   81
table(GSEtype)
# BE EAC HGD LGD 
# 21  24  18  18
quantile(GSEage,na.rm=T)
# 0%   25%   50%   75%  100% 
# 40.00 62.75 71.00 77.00 93.00 
sum(is.na(GSEage)) #1
table(GSEgender)
sum(is.na(GSEgender)) #1
quantile(GSEbmi,na.rm=T)
# 0%      25%      50%      75%     100% 
# 20.63956 24.95216 29.75737 30.98900 42.76745 
table(is.na(GSEbmi))
# FALSE  TRUE 
# 42    39 
table(is.na(GSEbmi),GSEtype)
# GSEtype
# BE EAC HGD LGD
# FALSE 11   8   9  14
# TRUE  10  16   9   4
table(GSEsmoke)
# No Yes 
# 14  40 
table(is.na(GSEsmoke))
# FALSE  TRUE 
# 54    27 
table(is.na(GSEsmoke),GSEtype)
# GSEtype
# BE EAC HGD LGD
# FALSE 16  14  12  12
# TRUE   5  10   6   6
table(GSEalcohol)
# 1 to 5 drinks a week      6 to 10 drinks a week Less than one drink a week                       None 
# 14                          5                         15                         19 
table(is.na(GSEalcohol))
# FALSE  TRUE 
# 53    28
table(is.na(GSEalcohol),GSEtype)
# GSEtype
# BE EAC HGD LGD
# FALSE 16  14  12  11
# TRUE   5  10   6   7
GSEclinical=data.frame(type=GSEtype,run=GSErun,age=GSEage,gender=GSEgender,smoke=GSEsmoke,bmi=GSEbmi,alcohol=GSEalcohol)
GSEclinical$type1=as.character(GSEclinical$type)
GSEclinical$type1[GSEclinical$type1 %in% c("HGD","EAC")]="HGDEAC"
GSEclinical$type1=factor(GSEclinical$type1)
rownames(GSEclinical)=colnames(GSEmeth)
#save(GSEmethall,GSEmeth,GSEage,GSEtype,GSEgender,GSEsmoke,GSEbmi,GSEalcohol,GSEclinical,file="../data/GSE89181.RData")

#BMI
library(limma)
idx=complete.cases(GSEclinical[,c("bmi","type","gender","age")])
quantile(GSEclinical$bmi[idx])
# 0%      25%      50%      75%     100% 
# 20.63956 24.95216 29.98314 30.98900 42.76745 
designmat=model.matrix(~I(bmi>30)+type+gender+age,data=GSEclinical[idx,])
fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
result=eBayes(fit)
colnames(result$p.value)
#qqplot(result$p.value[,2])
tmp=p.adjust(result$p.value[,2],method="bonferroni")
sum(tmp<0.2) #1
min(tmp) #[1] 0.008789511
rownames(result$p.value)[which.min(tmp)] #"cg13601977"
anno$UCSC_RefGene_Name[which(anno$Name=="cg13601977")] #""
which(anno$Name=="cg13601977") #242668
unique(anno$UCSC_RefGene_Name[(242668-20):242668]) #"PHF2"
unique(anno$UCSC_RefGene_Name[242668:(242668+30)]) #"BARX1"
idx=match(rownames(result$p.value),anno$Name)
result$p.value=as.data.frame(result$p.value)
result$p.value$gene=anno$UCSC_RefGene_Name[idx]
idx=which(grepl("MSRA",result$p.value$gene))
MSRAresult=result$p.value[idx,]
MSRAresult$gene=processsemicolon(MSRAresult$gene)
qqplot(result$p.value$`I(bmi > 30)TRUE`[idx])
quantile(result$p.value$`I(bmi > 30)TRUE`[idx])
#     0%         25%         50%         75%        100% 
# 0.007928466 0.135953940 0.455220264 0.697733963 0.969974741 
tmp=p.adjust(MSRAresult$`I(bmi > 30)TRUE`,method="fdr")
quantile(tmp)
which.min(tmp)
idx=match(rownames(MSRAresult),anno$Name)
MSRAresult$enhancer=anno$Enhancer[idx]
MSRAresult$promoter=anno$Regulatory_Feature_Group[idx]
MSRAresult$group=anno$UCSC_RefGene_Group[idx]
MSRAresult$group=processsemicolon(MSRAresult$group)
MSRAresult$position=anno$pos[idx]
idx1=which(GSEclinical$bmi<=30 & !is.na(GSEclinical$gender))
idx2=which(GSEclinical$bmi>30  & !is.na(GSEclinical$gender))
MSRAresult$beta1=MSRAresult$beta0=NA
for (i in 1:nrow(MSRAresult))
{
  idx=which(rownames(GSEmeth)==rownames(MSRAresult)[i])
  MSRAresult$beta0[i]=mean(as.numeric(GSEmeth[idx,idx1]))
  MSRAresult$beta1[i]=mean(as.numeric(GSEmeth[idx,idx2]))
}
MSRAresult=MSRAresult[order(MSRAresult$position),]
plot(MSRAresult$position,-log10(MSRAresult$`I(bmi > 30)TRUE`),xlab="Position",ylab="-log10(p-value)",cex.lab=1.2,cex.axis=1.2)
idx=order(MSRAresult$`I(bmi > 30)TRUE`)
points(MSRAresult$position[idx[1]],-log10(MSRAresult$`I(bmi > 30)TRUE`[idx[1]]),col="red",pch=16)
points(MSRAresult$position[idx[2:3]],-log10(MSRAresult$`I(bmi > 30)TRUE`[idx[2:3]]),col="blue",pch=16)
# abline(v=MSRAresult$position[idx[2]])
# abline(v=MSRAresult$position[idx[3]])
idx=order(MSRAresult$`I(bmi > 30)TRUE`)
idx=which(!grepl("Body",MSRAresult$group))
points(MSRAresult$position[idx],-log10(MSRAresult$`I(bmi > 30)TRUE`[idx]),col="green")
tmp=data.frame(CpG=rownames(MSRAresult)[idx],pvalue=MSRAresult$`I(bmi > 30)TRUE`[idx],beatabmi0=MSRAresult$beta0[idx],betabmi1=MSRAresult$beta1[idx],group=MSRAresult$group[idx],enhancer=MSRAresult$enhancer[idx])
print(tmp)
tmp=data.frame(CpG=rownames(MSRAresult),pvalue=MSRAresult$`I(bmi > 30)TRUE`,beatabmi0=MSRAresult$beta0,betabmi1=MSRAresult$beta1,group=MSRAresult$group,enhancer=MSRAresult$enhancer)
print(tmp)
#include TCGA mQTL result
mqtlres_nopeer=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/mqtl_EAC_26highrisk_cn_mutation_nopeers_cis",header=T,stringsAsFactors = F)
mqtlres_nopeer$SNP=updatesnpname(oldnames = mqtlres_nopeer$SNP)
idx=match(mqtlres_nopeer$gene,anno$Name)
mqtlres_nopeer$position=anno$pos[idx]
mqtlres_nopeer=mqtlres_nopeer[mqtlres_nopeer$SNP=="rs17749155",]
mqtlres_nopeer=mqtlres_nopeer[order(mqtlres_nopeer$p.value),]
plot(mqtlres_nopeer$position,-log10(mqtlres_nopeer$p.value))
tmp=intersect(mqtlres_nopeer$gene,rownames(MSRAresult))
mqtlres_nopeer$bmipvalue=NA
idx1=match(tmp,mqtlres_nopeer$gene)
idx2=match(tmp,rownames(MSRAresult))
mqtlres_nopeer$bmipvalue[idx1]=MSRAresult$`I(bmi > 30)TRUE`[idx2]
mqtlres_nopeer$meQTLFDR=p.adjust(mqtlres_nopeer$p.value,method="fdr")
plot(-log10(mqtlres_nopeer$p.value),-log10(mqtlres_nopeer$bmipvalue),xlab='-log10(meQTL-pvalue)',ylab='-log10(BMI-pvalue)',cex.lab=1.2,cex.axis=1.2)
#use gene location
phenotypeposfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt"
phenotypepos=as.data.frame(fread(phenotypeposfile,header=T))
idx=which(phenotypepos$geneid=="MSRA")
gr_gene=GRanges(seqnames = phenotypepos$chr[idx],ranges = IRanges(start=phenotypepos$s1[idx],end=phenotypepos$s2[idx]))
idx=match(rownames(result$p.value),anno$Name)
gr_result=GRanges(seqnames = anno$chr[idx],ranges = IRanges(start=anno$pos[idx],width=1))
tmp=distance(gr_result,gr_gene)
sum(tmp==0,na.rm=T) #75
idx=which(tmp==0)
result$p.value$gene1=NA
result$p.value$gene1[idx]="MSRA" #the annotation of CpG based on gene location is consistent with that in anno 

# idx=complete.cases(GSEclinical[,c("bmi","type","gender","age","run")])
# designmat=model.matrix(~I(bmi>30)+type+run+gender+age,data=GSEclinical[idx,])
# fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
# result=eBayes(fit)
# colnames(result$p.value)
# qqplot(result$p.value[,2])

# fit <- lmFit(as.matrix(GSEmethall[,idx]), designmat)
# result=eBayes(fit)
# colnames(result$p.value)
# qqplot(result$p.value[,2])

# idx=complete.cases(GSEclinical[,c("bmi","type1","gender","age")])
# designmat=model.matrix(~I(bmi>30)+type1+gender+age,data=GSEclinical[idx,])
# fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
# result=eBayes(fit)
# colnames(result$p.value)
# qqplot(result$p.value[,2])
# tmp=p.adjust(result$p.value[,2],method="bonferroni")
# sum(tmp<0.2) #1

#Smoke
idx=complete.cases(GSEclinical[,c("smoke","type","gender","age")])
table(GSEclinical$smoke[idx])
# No Yes 
# 13  40 
designmat=model.matrix(~smoke+type+gender+age,data=GSEclinical[idx,])
fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
result=eBayes(fit)
colnames(result$p.value)
#qqplot(result$p.value[,2])
tmp=p.adjust(result$p.value[,2],method="bonferroni")
sum(tmp<0.2) #1
min(tmp) #[1] 0.02899321
rownames(result$p.value)[which.min(tmp)] #"cg05951860"
anno$UCSC_RefGene_Name[which(anno$Name=="cg05951860")] #"CTTNBP2"

#Alcohol
idx=complete.cases(GSEclinical[,c("alcohol","type","gender","age")])
table(GSEclinical$alcohol[idx])
# 1 to 5 drinks a week      6 to 10 drinks a week 
# 14                          4 
# Less than one drink a week                       None 
# 15                         19 
designmat=model.matrix(~I(alcohol!="None")+type+gender+age,data=GSEclinical[idx,])
fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
result=eBayes(fit)
#qqplot(result$p.value[,2])

designmat=model.matrix(~I(!alcohol %in% c("None","Less than one drink a week","1 to 5 drinks a week"))+type+gender+age,data=GSEclinical[idx,])
fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
result=eBayes(fit)
#qqplot(result$p.value[,2])
tmp=p.adjust(result$p.value[,2],method="bonferroni")
sum(tmp<0.05) #1
tmp[tmp<0.05]
# cg23307203 
# 0.002274389 

anno$UCSC_RefGene_Name[which(anno$Name=="cg23307203")] # "GDF15"

# #remove "Less than one drink a week"
# idx=complete.cases(GSEclinical[,c("alcohol","type","gender","age")]) & ! GSEclinical$alcohol %in% "Less than one drink a week"
# table(GSEclinical$alcohol[idx])
# # 1 to 5 drinks a week      6 to 10 drinks a week 
# # 14                          4 
# # Less than one drink a week                       None 
# # 0                         19 
# designmat=model.matrix(~I(alcohol!="None")+type+gender+age,data=GSEclinical[idx,])
# fit <- lmFit(as.matrix(GSEmeth[,idx]), designmat)
# result=eBayes(fit)
# qqplot(result$p.value[,2])

#DMR--
idx=match(rownames(GSEmethall),anno$Name)
anno2=anno[idx,]
idx=complete.cases(GSEclinical[,c("bmi","gender","age","type")])
designmat=model.matrix(~I(bmi>30)+gender+age+type,data=GSEclinical[idx,])
GSEME1=as.matrix(GSEmethall[,idx])
all(anno2$Name==rownames(GSEmethall))
idx=order(anno2$chr,anno2$pos)
set.seed(1000)
DMRs_bmi <- bumphunter(GSEME1[idx,],designmat,chr=anno2$chr[idx],pos=anno2$pos[idx],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(DMRs_bmi$table)
DMRs_bmi$table[1:10,]
#       chr     start       end      value      area cluster indexStart indexEnd  L clusterL      p.value  fwer  p.valueArea fwerArea
# 1141   2  33951647  33951647  0.3761512 0.3761512   22006      50902    50902  1        1 1.989406e-06 0.006 0.0538443984    0.999
# 5447  10  93392836  93393509  0.1987603 2.3851240  109433     250426   250437 12       16 6.820819e-05 0.122 0.0002759590    0.355
# 3955   7  27225299  27225897  0.1528767 2.4460269   79403     187328   187343 16       41 1.125435e-04 0.211 0.0002492441    0.336
# 1919   3  19752610  19752610  0.2779276 0.2779276   36023      80866    80866  1        1 1.128277e-04 0.231 0.0893206152    1.000
# 9685  11  20690628  20692113 -0.1484704 2.2270555  116012     268057   268071 15       16 1.844463e-04 0.311 0.0003791239    0.424
# 301    1  23038925  23038925  0.2589345 0.2589345    4389      10563    10563  1        1 2.708434e-04 0.380 0.1095639536    1.000
# 9427   7 133811808 133811837 -0.2072338 0.6217015   86072     201850   201852  3       11 6.201261e-04 0.421 0.0180998960    0.986
# 3947   7  27213610  27214383  0.1602638 1.7629022   79396     187277   187287 11       11 5.502127e-04 0.461 0.0009682721    0.627
# 9298   6 105628017 105628063 -0.1899832 0.9499158   72304     170031   170035  5       12 7.795628e-04 0.464 0.0064979670    0.953
# 8969   3  43020793  43021399 -0.1855921 1.1135523   37147      83244    83249  6       12 7.326696e-04 0.467 0.0041919617    0.906

idx=complete.cases(GSEclinical[,c("smoke","gender","age","type")])
designmat=model.matrix(~smoke+gender+age+type,data=GSEclinical[idx,])
GSEME1=as.matrix(GSEmethall[,idx])
idx=order(anno2$chr,anno2$pos)
set.seed(1000)
DMRs_smoke <- bumphunter(GSEME1[idx,],designmat,chr=anno2$chr[idx],pos=anno2$pos[idx],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(DMRs_smoke$table)
DMRs_smoke$table[1:10,]
#      chr     start       end      value      area cluster indexStart indexEnd  L clusterL      p.value  fwer  p.valueArea fwerArea
# 5422  20  61051032  61052060  0.1513004 3.1773077  185527     439597   439617 21       26 1.766652e-05 0.055 5.159746e-05    0.124
# 1727   5 146257347 146258615  0.1754336 3.1578041   62339     138710   138727 18       21 2.327494e-05 0.066 5.299957e-05    0.125
# 1974   6  53212618  53214258  0.1873233 2.9971726   70809     166648   166663 16       16 2.495747e-05 0.072 7.515282e-05    0.162
# 5618   1 223747670 223747670 -0.3074130 0.3074130   17210      40417    40417  1        1 2.663999e-05 0.080 6.968237e-02    1.000
# 2445   7 157484045 157484164  0.2598569 0.5197138   88213     207040   207041  2        5 3.729599e-05 0.086 2.564955e-02    0.998
# 4248  15  65669563  65670304  0.2360441 1.1802203  149082     345429   345433  5        5 4.654988e-05 0.093 3.355517e-03    0.878
# 3588  12  24715250  24716076  0.2095081 2.0950814  126575     292916   292925 10       17 5.748630e-05 0.108 4.565254e-04    0.456
# 2537   8  26724744  26724836  0.2402334 0.9609336   91118     213715   213718  4        5 5.636462e-05 0.110 5.879306e-03    0.952
# 2785   9 134151854 134151945  0.2518478 0.5036956  101962     235135   235136  2        3 5.944925e-05 0.123 2.748322e-02    1.000
# 3340  11  45687055  45687494  0.2151335 1.7210684  117059     270416   270423  8        9 7.571366e-05 0.129 9.654894e-04    0.636

idx=complete.cases(GSEclinical[,c("alcohol","gender","age","type")])
designmat=model.matrix(~I(alcohol %in% "None")+gender+age+type,data=GSEclinical[idx,])
GSEME1=as.matrix(GSEmethall[,idx])
idx=order(anno2$chr,anno2$pos)
set.seed(1000)
DMRs_alcohol <- bumphunter(GSEME1[idx,],designmat,chr=anno2$chr[idx],pos=anno2$pos[idx],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(DMRs_alcohol$table)
DMRs_alcohol$table[1:10,]
#       chr     start       end      value      area cluster indexStart indexEnd  L clusterL      p.value  fwer  p.valueArea fwerArea
# 2810   7  94537584  94538408  0.1887414 3.0198627   83397     195745   195760 16       26 7.455590e-06 0.019 5.046861e-05    0.105
# 6870  19  12876846  12877188 -0.2268638 0.9074554  176471     414354   414357  4        4 2.236677e-05 0.044 5.853499e-03    0.924
# 208    1  82266841  82267171  0.2215965 0.8863861    9247      22071    22074  4        4 3.268990e-05 0.061 6.213374e-03    0.928
# 2356   6  53212618  53214258  0.1606903 2.5710442   70809     166648   166663 16       16 3.555743e-05 0.072 1.092531e-04    0.175
# 64     1  29449845  29451471  0.1878550 1.8785496    5294      13001    13010 10       20 5.104212e-05 0.091 5.006716e-04    0.463
# 5583  19  12146163  12147014  0.1639976 2.2959663  176345     413929   413942 14       15 7.082811e-05 0.128 1.840957e-04    0.257
# 5255  17  46692422  46692859  0.1941156 1.1646935  166875     389023   389028  6       15 9.290813e-05 0.141 2.893629e-03    0.828
# 750    2 128433058 128433507  0.1916536 1.1499216   27493      62656    62661  6        9 1.106868e-04 0.151 3.004890e-03    0.836
# 2272   6  31237013  31237034  0.2086969 0.6260908   68274     155662   155664  3        4 1.106868e-04 0.152 1.442657e-02    0.981
# 1036   3   3840333   3842717  0.1665812 2.1655553   34798      78317    78329 13       14 9.204787e-05 0.153 2.529166e-04    0.319

idx=complete.cases(GSEclinical[,c("alcohol","gender","age","type")])
designmat=model.matrix(~I(alcohol %in% "6 to 10 drinks a week")+gender+age+type,data=GSEclinical[idx,])
GSEME1=as.matrix(GSEmethall[,idx])
idx=order(anno2$chr,anno2$pos)
set.seed(1000)
DMRs_alcohol1 <- bumphunter(GSEME1[idx,],designmat,chr=anno2$chr[idx],pos=anno2$pos[idx],coef=2,pickCutoff=T,nullMethod="bootstrap",B=1000,type="beta")
dim(DMRs_alcohol1$table)
DMRs_alcohol1$table[1:10,]
#       chr     start       end      value      area cluster indexStart indexEnd L clusterL      p.value  fwer p.valueArea fwerArea
# 918    5  13664584  13664584 -0.4825031 0.4825031   56147     125328   125328 1        1 0.0000383315 0.082 0.076497285    1.000
# 722   17  46671293  46671298  0.3392313 0.6784626  166860     388938   388939 2       18 0.0005769900 0.345 0.041615906    1.000
# 720   17  46671077  46671131  0.2994003 0.8982008  166860     388932   388934 3       18 0.0016724639 0.591 0.023577620    0.998
# 262    5  54516487  54516879  0.2680995 1.6085969   57545     128171   128176 6        7 0.0018828549 0.614 0.006071652    0.926
# 1080  17  14204310  14204593 -0.2558552 1.5351309  163061     379365   379370 6       18 0.0028682339 0.723 0.006844911    0.938
# 963    7  27196153  27196371 -0.2400904 1.9207234   79386     187190   187197 8       27 0.0023730946 0.725 0.003763404    0.839
# 530   11  66838934  66839191  0.2726657 1.0906627  119341     276595   276598 4        5 0.0031818029 0.743 0.015353640    0.991
# 1003   8 121775052 121775052 -0.3337487 0.3337487   95501     222929   222929 1        1 0.0034077570 0.827 0.166470254    1.000
# 337    6 136915088 136915556  0.2712611 0.8137833   73929     173379   173381 3        4 0.0051862810 0.837 0.029436288    0.999
# 724   17  46674395  46674514  0.2712268 0.8136804  166861     388952   388954 3        9 0.0051931980 0.837 0.029442628    0.999

save(dmrs_smoking,dmrsall_smoking,dmrs_reflux,dmrsall_reflux,DMRs_bmi,DMRs_smoke,DMRs_alcohol,DMRs_alcohol1,file="../result/methylation_riskfactors_DMRres.RData")

#draw manhanttan plot
source("/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/code/functions.R")
#change to FWER cutoff
draw_manhattan1=function(pvalues=NULL,fdrs=NULL,fdrthreshold=0.05,maxy=6,chrs=NULL,keepres=F,logscale=T,main=NULL,ylab=NULL,opt=1)
{
  #if fdrs is a dataframe, the first row is pvalue, second row is fdr
  if (is.null(pvalues))
  {
    pvalues=fdrs[,1]
    if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
    names(pvalues)=rownames(fdrs)
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      if (! "chr" %in% colnames(pvalues) | ! "start" %in% colnames(pvalues)) #if the dataframe have start and chr columns, it already have coordinates
      {
        tmp=rownames(pvalues)
        pvalues=pvalues[,ncol(pvalues)] #the last column is pvalues
        if(class(pvalues)=="character") pvalues=as.numeric(pvalues)
        names(pvalues)=tmp
      }
    }
  }
  if (class(pvalues)=="numeric") #without gene coordinate, needs to figure out
  {
    if (is.null(chrs))
    {
      chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
      #chrs=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
    }
    geneposition=read.table(file="/fh/fast/dai_j/CancerGenomics/Ovarian_Cancer/result/geneposition_firehose.txt",header=T,sep="\t",stringsAsFactors=F)
    #form a dataframe with gene's position and p value
    df_pvalues=data.frame(gene=names(pvalues),pvalue=pvalues)
    genetable=merge(df_pvalues,geneposition,all.x=T)
    idxkeep=which(genetable$chr %in% chrs)
    genetable=genetable[idxkeep,]
    genetable=sortgenetable(genetable)
    genenames=genetable$gene
    #form the genetable for mahanttan plot
    genetable=data.frame(chr=genetable$chr,start=genetable$start,pvalue=genetable$pvalue,stringsAsFactors=F)
    rownames(genetable)=genenames
  }else
  {
    if (class(pvalues)=="data.frame")
    {
      warning("the columns should be chr, start and pvalue")
      genetable=pvalues
      genetable=sortgenetable(genetable)
    }
  }
  
  color2<-rep(c("gray33","brown","darkgreen","aquamarine4","azure4","darkred","cyan4","chartreuse4",
                "cornflowerblue","darkblue","azure3","darkgray","cadetblue","deepskyblue4"),2)
  par(las=1, xpd=TRUE,cex.axis=1.2,cex=1,cex.lab=1.5)
  ops<-mht.control(colors=color2,yline=1,xline=2,usepos=T,srt=0,logscale=logscale)
  if (is.null(ylab))
  {
    if (logscale==T)
    {
      ylab=expression(Observed~~-log[10](italic(p)))
    }else
    {
      ylab=expression(Observed~~italic(p))
    }
  }
  
  res=mhtplot2(genetable,ops,pch=10,bg=color2,ylab=ylab)
  #res=mhtplot2(genetable,ops,pch=10,bg=color2,cex.axis=1.2,cex.lab=1.2)
  #allpos keeps the overall cooridinate
  genetable=cbind(genetable,allpos=rep(NA,nrow(genetable)))
  genetable$allpos[res$idxnoNA]=res$CM
  #lines(c(0,max(res$CM)),c(0,0),lwd=2)
  if (logscale==T)
  {
    axis(2,at=c(floor(res$y1):(ceiling(res$y2))),pos=0, cex.axis=1.2)
  }else
  {
    axis(2,pos=0, cex.axis=1.2)
  }
  xat=c(res$chrlen2)
  axis(1,pos=floor(res$y1),labels=FALSE,tick=T,at=xat)
  if (! is.null(main))
    title(main=main)
  
  #to add fdr line
  #draw the fdr cutoff line
  if (!is.null(fdrs))
  {
    dictfile="/fh/fast/dai_j/CancerGenomics/Tools/database/reference/compact/ucsc.hg19.compact.dict"
    chrlen=read.table(file=dictfile,sep="\t",skip=1,stringsAsFactors = F)
    chrlen=chrlen$V3
    chrlen=gsub("LN:","",chrlen,fixed=T)
    chrlen=as.numeric(chrlen)
    names(chrlen)=c(1:22,"X","Y")
    #available chrs length
    chrlen1=chrlen[names(chrlen) %in% genetable$chr]
    
    #change fdrs,pvalues to vector
    fdrs1=fdrs
    if (class(fdrs)=="data.frame")
    {
      tmp=rownames(fdrs)
      fdrs1=fdrs[,ncol(fdrs)] #assume fdr is in the last column
      names(fdrs1)=tmp
    }
    pvalues1=pvalues
    if (class(pvalues)=="data.frame")
    {
      tmp=rownames(pvalues)
      pvalues1=pvalues[,ncol(pvalues)] #assume pvalue is in the last column
      names(pvalues1)=tmp
    }
    
    smallfdrs=fdrs1[fdrs1<=fdrthreshold]
    if (length(smallfdrs)>0) #if have fdr<fdrthreshold
    {
      idxs=which(fdrs1 %in% smallfdrs)
      if (opt==1)
      {
        pvaluethreshold=max(pvalues1[idxs])
      }else
      {
        #
        pvaluethreshold=fdrthreshold/length(pvalues1)
      }
      segments(0,-log10(pvaluethreshold),par('usr')[2],-log10(pvaluethreshold),col="red")
      #text(sum(chrlen1)/2,-log10(pvaluethreshold)+0.5,paste0("FDR=",fdrthreshold))
      #text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FDR=",fdrthreshold),cex=1.1)
      #text(sum(chrlen1)*0.9,-log10(pvaluethreshold)+0.25,paste0("FWER=",round(max(smallfdrs),digits = 2)),cex=1.1)
      text(sum(chrlen1)*0.2,-log10(pvaluethreshold)+0.25,paste0("FWER=",round(fdrthreshold,digits = 2)),cex=1.1)
    }else
    {
      warning("no genes with small fdr were found")
      #find the gene with min fdr
      
      minfdr=round(min(fdrs1,na.rm=T),digits = 2)
      idxs=which(fdrs1==minfdr)
      
      #use the one with max p-value
      thepvalue=max(pvalues1[idxs])
      #abline(h=-log10(thepvalue),col="red")
      segments(0,-log10(thepvalue),par('usr')[2],-log10(thepvalue),col="red")
      if (class(fdrs)=="data.frame")
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FWER=",minfdr),cex=1.1)
      }else
      {
        #text(sum(chrlen1)/2,-log10(thepvalue)+0.5,paste0("FDR=",minfdr))
        text(sum(chrlen1)*0.9,-log10(thepvalue)+0.25,paste0("FWER=",minfdr),cex=1.1)
      }
    }
  }
  par(cex.axis=1,cex=1,cex.lab=1)
  if (keepres==T) return(genetable)
}

plot_manhattan=function(dat=DMRs_bmi,opt=1)
{
  dat1=data.frame(chr=dat$table$chr,start=dat$table$start,pvalue=dat$table$p.value)
  #fdrs=p.adjust(dat$table$p.value,method = "fdr")
  draw_manhattan1(pvalues=dat1,fdrs=dat$table$fwer,fdrthreshold = 0.1,opt=opt)
}
postscript("../result/manhanttan_BMI_bumphunter1.ps",
           horizontal=F,width = 12, height = 4,pointsize = 13,paper = "special")
pdf("../result/manhanttan_BMI_bumphunter1.pdf",
           width = 12, height = 4,pointsize = 13,paper = "special")
plot_manhattan(opt=2)
dev.off()
postscript("../result/manhanttan_smoke_bumphunter1.ps",
           horizontal=F,width = 12, height = 4,pointsize = 13,paper = "special")
pdf("../result/manhanttan_smoke_bumphunter1.pdf",
           width = 12, height = 4,pointsize = 13,paper = "special")
plot_manhattan(dat=DMRs_smoke,opt=1)
dev.off()

#find the genes
ucsc_refseq=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/ucsc_refseqgenes.txt",header = T,stringsAsFactors = F,comment.char = "")
genecode22=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/gencode.gene.info.v22.csv",header = T,stringsAsFactors = F)
genecode22=genecode22[genecode22$gene_type=="protein_coding" & genecode22$gene_status=="KNOWN",]
knowngenes=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/other/knownGene1.txt",header=F,stringsAsFactors = F,sep="\t")
genetable1=data.frame(genename=ucsc_refseq$name2,chr=ucsc_refseq$chrom,start=ucsc_refseq$txStart,end=ucsc_refseq$txEnd,stringsAsFactors = F)
genetable1=genetable1[!duplicated(genetable1$genename),]
genetable2=data.frame(genename=genecode22$gene_name,chr=genecode22$seqname,start=genecode22$start,end=genecode22$end,stringsAsFactors = F)
genetable2=genetable2[!duplicated(genetable2$genename),]
genetable3=data.frame(genename=knowngenes$V14,chr=knowngenes$V2,start=knowngenes$V4,end=knowngenes$V5,stringsAsFactors = F)
genetable3=genetable3[!duplicated(genetable3$genename),]
idx=!genetable2$genename %in% genetable1$genename
genetable=rbind(genetable1,genetable2[idx,])
idx=!genetable3$genename %in% genetable$genename
genetable=rbind(genetable,genetable3[idx,])
chrs=paste0("chr",c(1:22,"X","Y"))
genetable=genetable[genetable$chr %in% chrs,]
genetable$chr=gsub("chr","",genetable$chr)
genetable$chr=factor(genetable$chr,levels = c(1:22,"X","Y"))
idx=order(genetable$chr,genetable$start)
genetable=genetable[idx,]
library(GenomicRanges)
gr_genetable=GRanges(seqnames = genetable$chr,ranges = IRanges(start=genetable$start,end=genetable$end))
findgenes=function(dat=DMRs_bmi)
{
  idx=which(dat$table$fwer<0.1)
  res=data.frame(chr=dat$table$chr[idx],start=dat$table$start[idx],end=dat$table$end[idx],gene=NA,stringsAsFactors = F)
  for (i in 1:nrow(res))
  {
    gr_tmp=GRanges(seqnames = dat$table$chr[idx[i]],ranges = IRanges(start=dat$table$start[idx[i]],end=dat$table$end[idx[i]]))
    tmp=distance(gr_genetable,gr_tmp)
    idx1=which(tmp==0)
    if (length(idx1)>0)
    {
      res$gene[i]=paste0(genetable$genename[idx1],collapse=",")
      #res$numgenes[i]=length(idx1)
    }else
    {
      idx1=which.min(tmp)
      res$gene[i]=paste0(genetable$genename[idx1],"/",genetable$genename[idx1+1])
    }
  }
  return(res)
}
genesBMI=findgenes()
genesBMI$gene="MYADML"
genessmoke=findgenes(dat=DMRs_smoke)
