#https://horvath.genetics.ucla.edu/html/GeneralFramework/GBMTutorialHorvath.pdf
#ml  R/3.5.3-foss-2016b-fh1
rm(list=ls())
library(WGCNA) 
options(stringsAsFactors=F) 
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
#geneexp data
TCGAGE=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_GE_norm.txt",header=T)
rownames(TCGAGE)=TCGAGE$id
TCGAGE=TCGAGE[,-1]
colnames(TCGAGE)=gsub(".","-",colnames(TCGAGE),fixed = T)
datExpr=t(TCGAGE)
# Now we investigate soft thesholding with the power adjacency function 
powers1=c(seq(1,10,by=1),seq(12,20,by=2)) 
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers1)[[2]] 
cex1=0.7 
par(mfrow=c(1,2)) 
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab=" Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n") 
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1,cex=cex1,col="red") 
# this line corresponds to using an R^2 cut-off of h 
abline(h=0.95,col="red") 
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n") 
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red") 
## Note that at power=6, the curve has an elbow or kink, i.e. for this power the scale free topology 
beta1=6 
Connectivity=softConnectivity(datExpr,power=beta1)-1 
# Let’s create a scale free topology plot. # The black curve corresponds to scale free topology and 
# the red curve corresponds to truncated scale free topology. 
par(mfrow=c(1,1)) 
scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=F)

#Module Detection
# An important step in network analysis is module detetion.
# Here we use methods that use clustering in combination with the topological 
# overlap matrix. 
# This code allows one to restrict the analysis to the most connected genes, 
# which may speed up calculations when it comes to module detection. 
ConnectivityCut = ncol(datExpr) # number of most connected genes that will be considered  
# Incidentally, in the paper by Mischel et al (2005) we considered all 3600 #genes.  
ConnectivityRank = rank(-Connectivity) 
restConnectivity = ConnectivityRank <= ConnectivityCut 
# thus our module detection uses the following number of genes 
sum(restConnectivity) 
# Now we define the adjacency matrix for the 3600 most connected genes 
ADJ= adjacency(datExpr[,restConnectivity],power=beta1,type="signed") #try signed
ADJ= adjacency(datExpr[,restConnectivity],power=beta1) #unsigned
gc() 
# The following code computes the topological overlap matrix based on the  
# adjacency matrix. 
# TIME: This about a few minutes.... 
dissTOM=TOMdist(ADJ) 
gc() 
# Now we carry out hierarchical clustering with the TOM matrix.  
# This takes a couple of minutes. 
hierTOM = hclust(as.dist(dissTOM),method="average")
par(mfrow=c(1,1),mar=c(5,4,3,1)) 
#plot(hierTOM,labels=F)
# According to our definition, modules correspond to branches of the tree. 
# The question is what height cut-off should be used? This depends on the  
# biology. Large heigth values lead to big modules, small values lead to small  
# but tight modules.  
# In reality, the user should use different thresholds to see how robust the findings are.
# The function cutreeStatistColor colors each gene by the branches that 
# result from choosing a particular height cut-off. 
# GREY IS RESERVED to color genes that are not part of any module. 
# We only consider modules that contain at least 125 genes. 
# #option1: use static cuttree
# staticMods=cutreeStatic(hierTOM,cutHeight = 0.94, minSize = 3)
# staticColors=labels2colors(staticMods)
# #colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.95, minSize = 3) 
# length(unique(staticColors)) #142,#127
# 
# # The above should be identical to colorh1=datSummary$color1[restConnectivity] 
# par(mfrow=c(2,1),mar=c(2,4,1,1)) 
# plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub=""); 
# plotColorUnderTree(hierTOM,colors=data.frame(module=staticColors)) 
# title("Module (branch) color") 
#dynamic cuttree
# option2: Module identification using dynamic tree cut:
minModuleSize = 10
dynamicMods = cutreeDynamic(dendro = hierTOM, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
dynamicColors=dynamicMods
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(hierTOM, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# #check clustering result
# hc <- hclust(dist(USArrests))
# tmp=cutreeStaticColor(hc,cutHeight = 50, minSize = 5) 
# par(mfrow=c(2,1),mar=c(2,4,1,1)) 
# plot(hc, main="Cluster Dendrogram", xlab="", sub=""); 
# plotColorUnderTree(hc,colors=data.frame(module=tmp)) 
# title("Module (branch) color") 
# rownames(USArrests)[which(tmp=="red")]

# An alternative view of this is the so called TOM plot that is generated by the 
# function TOMplot 
# Inputs:  TOM  distance measure, hierarchical (hclust) object, color  
# Warning: for large gene sets, say more than 2000 genes  
#this will take a while. I recommend you skip this. 
#TOMplot(dissTOM , hierTOM, staticColors) 
# We also propose to use classical multi-dimensional scaling plots  
# for visualizing the network. Here we chose 2 scaling dimensions 
# This also takes about 10 minutes... 
# cmd1=cmdscale(as.dist(dissTOM),2) 
# par(mfrow=c(1,1)) 
# plot(cmd1, col=as.character(staticColors),  main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2") 

# allcolors=unique(staticMods)
# clustercolor=staticMods
#use dynamic cut result
allcolors=unique(dynamicColors)
allcolors=allcolors[allcolors!="grey"]
allcolors=allcolors[allcolors!=0]
clustercolor=dynamicColors
#clustering result matrix
clusteringres=data.frame(genes=rep(NA,length(allcolors)),numgenes=NA)
for(i in 1:length(allcolors))
{
  idx=which(clustercolor==allcolors[i])
  clusteringres$genes[i]=paste0(colnames(datExpr)[idx],collapse = ",")
  clusteringres$numgenes[i]=length(idx)
}
quantile(clusteringres$numgenes)
# 0%  25%  50%  75% 100% 
# 5    9   19   55  706 
sum(clusteringres$numgenes)




#check clusters close to highrisk snps------
library(GenomicRanges)
highrisksnpstable=read.table("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Dong23SNPs_impscore.txt",header=T,sep="\t")
gr_highrisk=GRanges(seqnames = highrisksnpstable$Chr,ranges = IRanges(start=highrisksnpstable$Position,width = 1))
clusteringres$cisgenes=NA
clusteringres$numcisgenes=NA
clusteringres$highrisksnps=NA
clusteringres$numhighrisksnps=NA
clusteringres$pairs=NA
clusteringres$numpairs=NA
GEpos=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt",header = T)
idx=match(colnames(datExpr),GEpos$geneid)
GEpos=GEpos[idx,]
gr_genes=GRanges(seqnames = GEpos$chr,ranges = IRanges(start=GEpos$s1,end=GEpos$s2))
for (i in 1:nrow(clusteringres))
{
  genes=unlist(strsplit(clusteringres$genes[i],","))
  idx=match(genes,GEpos$geneid)
  if (sum(is.na(idx))>0) warning("NA")
  gr_genes1=gr_genes[idx]
  # tmp=distanceToNearest(gr_genes1,gr_highrisk)
  # idx=which(tmp@elementMetadata$distance<5e5)
  # if (length(idx)>0)
  # {
  #   clusteringres$cisgenes[i]=paste0(unique(genes[tmp@from[idx]]),collapse = ",")
  #   clusteringres$numcisgenes[i]=length(unique(genes[tmp@from[idx]]))
  #   clusteringres$highrisksnps[i]=paste0(unique(highrisksnpstable$SNP[tmp@to[idx]]),collapse = ",")
  #   clusteringres$numhighrisksnps[i]=length(unique(highrisksnpstable$SNP[tmp@to[idx]]))
  #   tmp1=paste0(genes[tmp@from[idx]],"_",highrisksnpstable$SNP[tmp@to[idx]])
  #   clusteringres$pairs[i]=paste0(tmp1,collapse = ",")
  #   clusteringres$numpairs[i]=length(idx)
  # }
  cisgenes=highrisksnps=NULL
  for (j in 1:length(gr_highrisk))
  {
    dist1=distance(gr_genes1,gr_highrisk[j])
    idx1=which(dist1<5e5)
    if (length(idx1)>0)
    {
      cisgenes=c(cisgenes,genes[idx1])
      highrisksnps=c(highrisksnps,highrisksnpstable$SNP[j])
    }
  }
  clusteringres$cisgenes[i]=clusteringres$highrisksnps[i]=NA
  clusteringres$numcisgenes[i]=clusteringres$numhighrisksnps[i]=NA
  if (length(cisgenes)>0)
  {
    clusteringres$cisgenes[i]=paste0(cisgenes,collapse=",")
    clusteringres$numcisgenes[i]=length(cisgenes)
    clusteringres$highrisksnps[i]=paste0(highrisksnps,collapse=",")
    clusteringres$numhighrisksnps[i]=length(highrisksnps)
  }
}
sum(clusteringres$numcisgenes!=0,na.rm=T) #37
includedsnps=NULL
includedgenes=NULL
for (i in which(!is.na(clusteringres$numcisgenes)))
{
  tmp=unique(unlist(strsplit(clusteringres$highrisksnps[i],",")))
  includedsnps=unique(c(includedsnps,tmp))
  tmp=unique(unlist(strsplit(clusteringres$cisgenes[i],",")))
  includedgenes=unique(c(includedgenes,tmp))
}
length(includedsnps) #18
sum(clusteringres$numhighrisksnps,na.rm=T) #62
length(includedgenes) #82
sum(clusteringres$numcisgenes,na.rm=T) #82
idx=which(!is.na(clusteringres$cisgenes))
quantile(clusteringres$numgenes[idx])
quantile(clusteringres$numhighrisksnps[idx])
quantile(clusteringres$numgenes[idx])

#check eQTL in clusters

snpfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_EAC_23highrisk_SNP_GE.txt"
phenotypefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_GE_norm.txt"
copynumberfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_copynumber.txt"
mutationfile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_mutation_GE.txt"
covariatefile="/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_COVA_GE_PEER_clinical_used.txt"
snpdat=read.table(snpfile,header = T,stringsAsFactors = F)
phenotype=read.table(phenotypefile,header = T,stringsAsFactors = F)
cn=read.table(copynumberfile,header = T,stringsAsFactors = F)
mutation=read.table(mutationfile,header = T,stringsAsFactors = F)
covariate=read.table(covariatefile,header = T,stringsAsFactors = F)
all(colnames(snpdat)==colnames(phenotype))
all(colnames(snpdat)==colnames(cn))
all(colnames(snpdat)==colnames(mutation))
all(colnames(snpdat)==colnames(covariate))
covariate=as.data.frame(t(covariate))
colnames(covariate)=covariate[1,]
covariate=covariate[-1,]
for (i in 1:ncol(covariate)) covariate[,i]=as.numeric(as.character(covariate[,i]))
updatesnpname=function(oldnames=eqtlres$SNP)
{
  newnames=sapply(oldnames,function(x){
    unlist(strsplit(x,":"))[1]
  })
  return(newnames)
}
snpdat$id=updatesnpname(snpdat$id)
mutation$id=updatesnpname(mutation$id)
eqtlres_nopeer=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/eqtl_EAC_highrisk_cn_mutation_nopeers_cis",header=T,stringsAsFactors = F)
eqtlres_nopeer$SNP=updatesnpname(oldnames = eqtlres_nopeer$SNP)
eqtlres_nopeer=eqtlres_nopeer[order(eqtlres_nopeer$p.value),]

#get the clusters close to a highrisk snp
checkeQTL=function(snp="rs9257809")
{
  egenes=eqtlres_nopeer$gene[eqtlres_nopeer$p.value<0.9 & eqtlres_nopeer$SNP==snp]
  idx=which(grepl(snp,clusteringres$highrisksnps))
  avaiegenes=avaiegenescluster=NULL
  if (length(egenes)>0)
  {
    count=0
    for (i in 1:length(egenes))
    {
      if (sum(grepl(egenes[i],clusteringres$cisgenes[idx]))>0)
      {
        count=count+1
        avaiegenes=c(avaiegenes,egenes[i])
        for (j in 1:length(idx))
        {
          if (grepl(egenes[i],clusteringres$cisgenes[idx[j]])>0)
          {
            avaiegenescluster=c(avaiegenescluster,paste0("cluster",j))
          }
        }
      }
    }
    print (paste0(count,":",paste0(avaiegenes,collapse = ",")," in clusters"))
  }
  numclusters=length(idx)
  clustereqllist=list()
  clustereqllist[["availegenes"]]=avaiegenes
  clustereqllist[["availegenescluster"]]=avaiegenescluster
  clustereqllist[["numclusters"]]=numclusters
  clustereqllist[["egenes"]]=egenes
  nameclusters=paste0("cluster",1:length(idx))
  idxsnp=which(snpdat$id==snp)
  snpV=as.numeric(snpdat[idxsnp,2:ncol(snpdat)])
  idxmutation=which(mutation$id==snp)
  mutationV=as.numeric(mutation[idxmutation,2:ncol(mutation)])
  for (i in 1:length(idx))
  {
    genes=unlist(strsplit(clusteringres$genes[idx[i]],","))
    clustereqtl=data.frame(snp=rep(snp,length(genes)),gene=genes,egene=F,pvalue=NA)
    for (j in 1:length(genes))
    {
      if (genes[j] %in% egenes) clustereqtl$egene[j]=T
      idx1=which(phenotype$id==genes[j])
      phenotypeV=as.numeric(phenotype[idx1,2:ncol(phenotype)])
      idx1=which(cn$id==genes[j])                      
      cnV=as.numeric(cn[idx1,2:ncol(cn)])
      fm=glm(phenotypeV~snpV+cnV+mutationV+pc1+pc2+pc3+pc4,data=covariate)
      clustereqtl$pvalue[j]=summary(fm)$coefficients[2,4]
    }
    clustereqllist[[nameclusters[i]]]=clustereqtl
  }
  return(clustereqllist)
}
clustereqllist_rs9257809=checkeQTL()
qqplot=function(pvalue=NULL,main="",xlim=NULL,ylim=NULL)
{
  n=length(pvalue)
  plot(-log((1:n)/n,base=10),-log(pvalue[order(pvalue)],base=10),xlab="Expected p-value (log base 10)",
       ylab="Observed p-value (log base 10)",main=main,xlim=xlim,ylim=ylim,cex.lab=1.3,cex.axis=1.3)
  abline(0,1,lty=2)
}

#draw plot of p-value ~ association, i=1:4 for 4 egenes detected in clusters
i=4 #the first egene
egene=clustereqllist_rs9257809$availegenes[i]
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))
interestedcluster=clustereqllist_rs9257809[[clustereqllist_rs9257809$availegenescluster[i]]]
p=round(interestedcluster$pvalue[interestedcluster$gene==egene],4)
othergene=interestedcluster$gene
othergene=othergene[othergene!=egene]
idx=match(othergene,interestedcluster$gene)
xeqtl=-log10(interestedcluster$pvalue[idx])
ycorr=rep(NA,length(othergene))
idx=match(egene,rownames(TCGAGE))
egeneexp=as.numeric(TCGAGE[idx,])
for (i in 1:length(othergene))
{
  idx=match(othergene[i],rownames(TCGAGE))
  ycorr[i]=cor(egeneexp,as.numeric(TCGAGE[idx,]))
}
plot(ycorr,xeqtl,xlab="correlation",ylab="-log10(p-value)",main=paste0(egene,":p-value=",p))
abline(glm(xeqtl~ycorr),col="red")


###################################################################
## now take top 6 other genes with correlation greater than 0.85 ##
###################################################################


topcgene <- c(egene,othergene[ycorr>0.4])

idxsnp=which(snpdat$id=="rs9257809")
snpV=as.numeric(snpdat[idxsnp,2:ncol(snpdat)])

GEtopcgene <- NULL
for (j in 1:length(topcgene)) {
   idx1=which(phenotype$id==topcgene[j])
   phenotypeV=as.numeric(phenotype[idx1,2:ncol(phenotype)])
   GEtopcgene <- c(GEtopcgene,phenotypeV)
}

topSNP <- rep(snpV,times=length(topcgene))
IDx <- rep(1:length(snpV),times=length(topcgene))


topdat <- data.frame(cbind(GEtopcgene,topSNP,IDx))
names(topdat) <- c("Gene","SNP","subjID")

library(geepack)
summary(geeglm(Gene~SNP,id=subjID,data=topdat))


fit1 <- glm(topdat[,1]~topdat[,7],data=topdat)
summary(glm(topdat[,2]~topdat[,7],data=topdat))
summary(glm(topdat[,3]~topdat[,7],data=topdat))
summary(glm(topdat[,4]~topdat[,7],data=topdat))
summary(glm(topdat[,5]~topdat[,7],data=topdat))
summary(glm(topdat[,6]~topdat[,7],data=topdat))


#draw q-q plots for clusters containing egenes
par(mfrow=c(length(unique(clustereqllist_rs9257809$availegenescluster)),1))
for (clustername in unique(clustereqllist_rs9257809$availegenescluster))
{
  print(nrow(clustereqllist_rs9257809[[clustername]]))
  qqplot(clustereqllist_rs9257809[[clustername]]$pvalue)
  #to view one cluster, run:
  #View(clustereqllist_rs9257809[[clustername]])
  
}






















clustereqllist_rs10419226=checkeQTL(snp="rs10419226")


#draw plot of p-value ~ association, i=1:4 for 4 egenes detected in clusters
i=3 #the first egene
egene=clustereqllist_rs10419226$availegenes[i]
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))
interestedcluster=clustereqllist_rs10419226[[clustereqllist_rs10419226$availegenescluster[i]]]
p=round(interestedcluster$pvalue[interestedcluster$gene==egene],4)
othergene=interestedcluster$gene
othergene=othergene[othergene!=egene]
idx=match(othergene,interestedcluster$gene)
xeqtl=-log10(interestedcluster$pvalue[idx])
ycorr=rep(NA,length(othergene))
idx=match(egene,rownames(TCGAGE))
egeneexp=as.numeric(TCGAGE[idx,])
for (i in 1:length(othergene))
{
  idx=match(othergene[i],rownames(TCGAGE))
  ycorr[i]=cor(egeneexp,as.numeric(TCGAGE[idx,]))
}
plot(ycorr,xeqtl,xlab="correlation",ylab="-log10(p-value)",main=paste0(egene,":p-value=",p))
abline(glm(xeqtl~ycorr),col="red")

topcgene <- c(egene,othergene[ycorr>0.3])

idxsnp=which(snpdat$id=="rs10419226")
snpV=as.numeric(snpdat[idxsnp,2:ncol(snpdat)])

GEtopcgene <- NULL
for (j in 1:length(topcgene)) {
  idx1=which(phenotype$id==topcgene[j])
  phenotypeV=as.numeric(phenotype[idx1,2:ncol(phenotype)])
  GEtopcgene <- c(GEtopcgene,phenotypeV)
}

topSNP <- rep(snpV,times=length(topcgene))
IDx <- rep(1:length(snpV),times=length(topcgene))


topdat <- data.frame(cbind(GEtopcgene,topSNP,IDx))
names(topdat) <- c("Gene","SNP","subjID")

library(geepack)
summary(geeglm(Gene~SNP,id=subjID,data=topdat))












clustereqllist_rs2464469=checkeQTL(snp="rs2464469")


#draw plot of p-value ~ association, i=1:4 for 4 egenes detected in clusters
i=2 #the first egene
egene=clustereqllist_rs2464469$availegenes[i]
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))
interestedcluster=clustereqllist_rs2464469[[clustereqllist_rs2464469$availegenescluster[i]]]
p=round(interestedcluster$pvalue[interestedcluster$gene==egene],4)
othergene=interestedcluster$gene
othergene=othergene[othergene!=egene]
idx=match(othergene,interestedcluster$gene)
xeqtl=-log10(interestedcluster$pvalue[idx])
ycorr=rep(NA,length(othergene))
idx=match(egene,rownames(TCGAGE))
egeneexp=as.numeric(TCGAGE[idx,])
for (i in 1:length(othergene))
{
  idx=match(othergene[i],rownames(TCGAGE))
  ycorr[i]=cor(egeneexp,as.numeric(TCGAGE[idx,]))
}
plot(ycorr,xeqtl,xlab="correlation",ylab="-log10(p-value)",main=paste0(egene,":p-value=",p))
abline(glm(xeqtl~ycorr),col="red")


topcgene <- c(egene,othergene[ycorr>0.3])

idxsnp=which(snpdat$id=="rs2464469")
snpV=as.numeric(snpdat[idxsnp,2:ncol(snpdat)])

GEtopcgene <- NULL
for (j in 1:length(topcgene)) {
  idx1=which(phenotype$id==topcgene[j])
  phenotypeV=as.numeric(phenotype[idx1,2:ncol(phenotype)])
  GEtopcgene <- c(GEtopcgene,phenotypeV)
}

topSNP <- rep(snpV,times=length(topcgene))
IDx <- rep(1:length(snpV),times=length(topcgene))


topdat <- data.frame(cbind(GEtopcgene,topSNP,IDx))
names(topdat) <- c("Gene","SNP","subjID")

library(geepack)
summary(geeglm(Gene~SNP,id=subjID,data=topdat))





usnp <- unique(eqtlres_nopeer$SNP)

for (i in 1:23) {
  cat(usnp[i],"\n")
  cat(min(eqtlres_nopeer$p.value[eqtlres_nopeer$SNP==usnp[i]]),"\n")
  clustereqllist_rs2464469=checkeQTL(snp=usnp[i])
}



SNPname <- "rs1979654"

clustereqllist_rs2464469=checkeQTL(snp=SNPname)

#draw plot of p-value ~ association, i=1:4 for 4 egenes detected in clusters
i=2 #the first egene
egene=clustereqllist_rs2464469$availegenes[i]
par(mfrow=c(1,1))
par(mar=c(4,4,2,1))
interestedcluster=clustereqllist_rs2464469[[clustereqllist_rs2464469$availegenescluster[i]]]
p=round(interestedcluster$pvalue[interestedcluster$gene==egene],4)
othergene=interestedcluster$gene
othergene=othergene[othergene!=egene]
idx=match(othergene,interestedcluster$gene)
xeqtl=-log10(interestedcluster$pvalue[idx])
ycorr=rep(NA,length(othergene))
idx=match(egene,rownames(TCGAGE))
egeneexp=as.numeric(TCGAGE[idx,])
for (i in 1:length(othergene))
{
  idx=match(othergene[i],rownames(TCGAGE))
  ycorr[i]=cor(egeneexp,as.numeric(TCGAGE[idx,]))
}
plot(ycorr,xeqtl,xlab="correlation",ylab="-log10(p-value)",main=paste0(egene,":p-value=",p))
abline(glm(xeqtl~ycorr),col="red")


topcgene <- c(egene,othergene[ycorr>0.5])

idxsnp=which(snpdat$id==SNPname)
snpV=as.numeric(snpdat[idxsnp,2:ncol(snpdat)])

GEtopcgene <- NULL
for (j in 1:length(topcgene)) {
  idx1=which(phenotype$id==topcgene[j])
  phenotypeV=as.numeric(phenotype[idx1,2:ncol(phenotype)])
  GEtopcgene <- c(GEtopcgene,phenotypeV)
}

topSNP <- rep(snpV,times=length(topcgene))
IDx <- rep(1:length(snpV),times=length(topcgene))


topdat <- data.frame(cbind(GEtopcgene,topSNP,IDx))
names(topdat) <- c("Gene","SNP","subjID")

library(geepack)
summary(geeglm(Gene~SNP,id=subjID,data=topdat))



