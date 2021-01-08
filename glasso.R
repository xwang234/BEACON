#!/usr/bin/env Rscript
#https://www4.stat.ncsu.edu/~reich/BigData/code/glasso.html


rm(list=ls())
library(glasso)
TCGAGE=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_EAC_GE_norm.txt",header=T)
rownames(TCGAGE)=TCGAGE$id
TCGAGE=TCGAGE[,-1]
colnames(TCGAGE)=gsub(".","-",colnames(TCGAGE),fixed = T)
TCGAGE1=TCGAGE
#only consider cis genes
genepos=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_GE_gene_POS.txt",header=T,stringsAsFactors = F)
snppos=read.table("/fh/fast/dai_j/CancerGenomics/EAprogression/result/qtl_input/TCGA_tumors_26highrisk_SNP_POS.txt",header=T,stringsAsFactors = F)
library(GenomicRanges)
gr_genepos=GRanges(seqnames = genepos$chr,ranges = IRanges(start=genepos$s1,end=genepos$s2))
gr_snppos=GRanges(seqnames = snppos$chr,ranges=IRanges(snppos$position,width = 1))
distcutoff=1e6
selectedgenes=NULL
for (i in 1:length(gr_snppos))
{
  tmp=distance(gr_genepos,gr_snppos[i])
  idx=which(tmp<distcutoff)
  if (length(idx)>0)
  {
    selectedgenes=unique(c(selectedgenes,genepos$geneid[idx]))
  }else
  {
    print(i)
  }
}
idx=which(rownames(TCGAGE) %in% selectedgenes)
TCGAGE1=TCGAGE[idx,]
datExpr=t(TCGAGE1)

#tuning parameter using BIC
n=nrow(datExpr)
covmat=cov(datExpr)
nr=10
rho <- seq(0.0001,0.1,length=nr)
bic <- rho
for (i in 1:nr)
{
  if (i %%2==0) print(Sys.time())
  a=glasso(s=covmat,rho[i])
  p_off_d <- sum(a$wi!=0 & col(covmat)<row(covmat))
  bic[i]  <- -2*(a$loglik) + p_off_d*log(n)
}

best <- which.min(bic)
plot(rho,bic)
points(rho[best],bic[best],pch=19)

a <- glasso(covmat,rho[best])
names(a)

library(network)
P <- a$wi
A <- ifelse(P!=0 & row(P)!=col(P),1,0)
g <- network(A)
plot(g,label=1:ncol(datExpr),main="Estimated network")

#get clusters
library(igraph)
# adj=P
# for (i in 1:ncol(P))
# {
#   idx=abs(P[,i])>0.025
#   adj[idx,i]=1
#   adj[!idx,i]=0
# }
grf=graph_from_adjacency_matrix(A, mode ="undirected",diag=F)
plot(grf, edge.arrow.size=.1,vertex.label=NA)
grf_clusters=clusters(grf)
length(unique(grf_clusters$membership)) #13
table(grf_clusters$membership)
grf_egelist=get.edgelist(grf, names=TRUE)


Sys.time()
glassores=glasso(s = covmat,rho = 0.1,approx=TRUE) #it takes too long time
Sys.time()

a=glassores

#try huge
library(huge)
Sys.time()
datExpr.npn=huge.npn(datExpr)
Sys.time()
#out.npn = huge(datExpr.npn,method = "glasso", nlambda=10,lambda.min.ratio = 0.1)
out.npn = huge(datExpr.npn,nlambda=3,method="glasso",lambda.min.ratio = 0.3)
Sys.time()
out = huge(datExpr,method = "glasso", nlambda=40,lambda.min.ratio = 0.1)
Sys.time()

