
# fusion hsq, hsq p-value <0.01, hsq V3:VG(variation in genotype), V4:VP(variation in phenotype), V5:p-value
fusiontable=as.data.frame(fread("/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/weights/Esophagus_Gastroesophageal_Junction/Esophagus_Gastroesophageal_Junction.hsq"))
fusiontable$geneid=NA
fusiontable$symbol=NA
for (i in 1:nrow(fusiontable))
{
  tmp=unlist(strsplit(fusiontable$V2[i],"/"))
  fusiontable$geneid[i]=gsub("Esophagus_Gastroesophageal_Junction.","",tmp[3])
  idx=which(gtex_ge_anno$Probe_Id==fusiontable$geneid[i])
  fusiontable$symbol[i]=gtex_ge_anno$Symbol[idx]
}
all(fusiontable$symbol %in% rownames(res_min)) #T
all(fusiontable$geneid %in% gtex_ge_anno$Probe_Id) #T
sum(fusiontable$V5<0.01)
fusionfiles=list.files("/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/weights/Esophagus_Gastroesophageal_Junction/",pattern="*.RDat")
fusiongenes=rep(NA,length(fusionfiles))
for (i in 1:length(fusiongenes))
{
  tmp=unlist(strsplit(fusionfiles[i],"_"))
  tmp=gsub("Junction.","",tmp[3],fixed=T)
  fusiongenes[i]=gsub(".wgt.RDat","",tmp,fixed=T)
}
idx=match(fusiongenes,fusiontable$geneid)
fusiontable$avai=F
fusiontable$avai[idx]=T
fusiontable=rbind(fusiontable[fusiontable$avai==T,],fusiontable[fusiontable$avai==F,])
fusiontable=fusiontable[!duplicated(fusiontable$symbol),]

idx=which(fusiontable$avai==T)
quantile(fusiontable$V5[idx])
# 0%        25%        50%        75%       100% 
# 0.0000e+00 1.3445e-08 6.1550e-05 1.7485e-03 9.9930e-03 
idx=which(fusiontable$avai==F & fusiontable$V5<0.01)
load("/fh/fast/stanford_j/Xiaoyu/QTL/data/GTEx/gtex_ge_anno.RData")
gtex_ge_anno=gtex_ge_anno[gtex_ge_anno$V3=="gene",]
prefix="dist500K_GTEx_April18"
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix)
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
models1=rownames(res_min)[res_min$glmflag==1]
print(paste0("number of models in our data ",length(models1)))
#number of models in our data 22769
models2=fusiontable$symbol[fusiontable$avai]
print(paste0("number of models in fusion ",length(models2)))
#number of models in fusion 4887
print(paste0("number of common models ",length(intersect(models1,models2))))
#number of common models 3322
idx=match(models2,gtex_ge_anno$Symbol)
table(gtex_ge_anno$gene_type[idx])

tmp=intersect(fusiontable$symbol,models1)
idx=match(tmp,fusiontable$symbol)
fusiontable$avai_our=F
fusiontable$avai_our[idx]=T

#consensus hm3,not the complete list.
hm3snp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/HM3/hapmap3_r3_b36_fwd.consensus.qc.poly.map"))
allgenes=intersect(models1,models2)
idx=match(allgenes,rownames(res_min))
comres=data.frame(our_totalsnps=rep(0,length(allgenes)),our_selsnps=0,fusion_totalsnps=0,fusion_selsnps=0,com_totalsnps=0,com_selsnps=0,com_selsnps_ale=0,
               mean_cor=NA,med_cor=NA,expr_cor=NA,r2=res_min$r2[idx],hsq_p=NA)
idx=match(allgenes,fusiontable$symbol)
comres$hsq_p=fusiontable$V5[idx]
rownames(comres)=allgenes
# gene="KXD1"
# gene="UBAC1"
# 
# gene="HSP90AA1"
# gene="ISYNA1"
for (j in 1:1000)
{
 if (j %%100==0) cat(j,'..')
 gene=allgenes[j] 
 idx=which(rownames(res_min)==gene)
 totalsnps1=res_min$numtotalsnp[idx]
 selsnps1=res_min$selectedsnps[idx]
 selsnps1=unlist(strsplit(selsnps1,"|",fixed=T))
 selsnps1_=selsnps1
 for(i in 1:length(selsnps1))
 {
   tmp=unlist(strsplit(selsnps1[i],"_"))
   selsnps1_[i]=tmp[1]
 }
 coeff1=res_min$selectedsnps_coeff[idx]
 coeff1=as.numeric(unlist(strsplit(coeff1,"|",fixed=T)))
 comres$our_totalsnps[j]=totalsnps1
 comres$our_selsnps[j]=length(coeff1)
 
 #FUSION
 idx=which(fusiontable$symbol==gene)
 if (length(idx)>0) #can be found in fusion
 {
   #print(paste0("hsq-pvalue of ",gene,":",fusiontable$V5[idx]))
   if (fusiontable$avai[idx[1]]==T) #can be found in fusion significant results
   {
     
     load(paste0("/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/weights/Esophagus_Gastroesophageal_Junction/Esophagus_Gastroesophageal_Junction.",fusiontable$geneid[idx[1]],".wgt.RDat"))
     wgt.matrix=as.data.frame(wgt.matrix)
     #check HM3
     # #sum(snps$V2 %in% hm3snp$V2)
     # tmp=snps$V2[snps$V2 %in% hm3snp$V2]
     # idx=match(tmp,hm3snp$V2)
     # View(hm3snp[min(idx):max(idx),])
     totalsnps2=nrow(wgt.matrix)
     selsnps2=rownames(wgt.matrix)[wgt.matrix$enet!=0]
     
     coeff2_=wgt.matrix$enet[wgt.matrix$enet!=0]
     comres$fusion_totalsnps[j]=totalsnps2
     comres$fusion_selsnps[j]=length(coeff2_)
     if (length(coeff2_)>0) #some may select 0 SNPs
     {
       # #the info can be also found in snps matrix. No need to use biomart
       # library("biomaRt")
       # #snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
       # snpmart = useMart(biomart = "ENSEMBL_MART_SNP", 
       #                   host    = "grch37.ensembl.org", 
       #                   path    = "/biomart/martservice", 
       #                   dataset = "hsapiens_snp")
       selsnps2=data.frame(snp=selsnps2,chr=NA,pos=NA,ale1=NA,ale2=NA,stringsAsFactors = F)
       #for (i in 1:nrow(selsnps2))
       #{
       #tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart)
       #tmp=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chr_name'), filters = c('snp_filter'), values =selsnps2$snp[i],mart = snpmart)
       #tmp1=unlist(strsplit(tmp$allele,"/",fixed=T))
       
       #selsnps2$ale1[i]=tmp1[1]
       #selsnps2$ale2[i]=tmp1[2]
       #selsnps2$chr[i]=tmp$chr_name
       #selsnps2$pos[i]=tmp$chrom_start
       #}
       idx=match(selsnps2$snp,snps$V2)
       selsnps2$ale1=snps$V5[idx]
       selsnps2$ale2=snps$V6[idx]
       selsnps2$chr=snps$V1[idx]
       selsnps2$pos=snps$V4[idx]
       str11=paste0(selsnps2$chr,":",selsnps2$pos,"_",selsnps2$ale1,"_",selsnps2$ale2)
       str12=paste0(selsnps2$chr,":",selsnps2$pos,"_",selsnps2$ale2,"_",selsnps2$ale1)
       comres$com_selsnps_ale[j]=sum(str11 %in% selsnps1 | str12 %in% selsnps1) 
       selsnps2_=str11
       for(i in 1:nrow(selsnps2))
       {
         tmp=unlist(strsplit(selsnps2_[i],"_"))
         selsnps2_[i]=tmp[1]
       }
       comres$com_selsnps[j]=sum(selsnps1_ %in% selsnps2_) #0, no overlap
       
       #check correlation
       if (!exists("snp")) load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")
       idx11_=which(str11 %in% rownames(snp))
       idx12_=which(str12 %in% rownames(snp))
       idx11=which(rownames(snp) %in% str11)
       idx12=which(rownames(snp) %in% str12)
       
       dat2=snp[idx12,]
       dat2=rbind(dat2,snp[idx11,])
       coeff2=coeff2_[c(idx12_,idx11_)]
       idx=match(selsnps1,rownames(snp))
       dat1=snp[idx,]
       tmp=cor(t(dat1),t(dat2))
       # for(i in 1:ncol(tmp))
       # {
       #   print(max(as.numeric(abs(tmp[,i]))))
       # }
       #heatmap(abs(tmp))
       comres$mean_cor[j]=mean(abs(tmp))
       comres$med_cor[j]=median(abs(tmp))
       # print(paste0("total number of SNP1: ",totalsnps1))
       # print(paste0("total number of SNP2: ",totalsnps2))
       # print(paste0("total number of SNP1: ",length(coeff1)))
       # print(paste0("total number of SNP2: ",length(coeff2)))
       
       #check prop of FUSION in our data
       str11=paste0(snps$V1,":",snps$V4,"_",snps$V5,"_",snps$V6)
       str12=paste0(snps$V1,":",snps$V4,"_",snps$V6,"_",snps$V6)
       # tmp=c(str11,str12)[c(str11,str12) %in% rownames(snp)]
       # idx=match(tmp,rownames(snp))
       # tmp=snp[idx,]
       # tmp=cor(t(tmp))
       # diag(tmp)=NA
       # tmp[upper.tri(tmp)]=NA
       # hist(unlist(tmp),xlab="correlation between fusion SNPs",main=gene)
       tmp=sum(c(str11,str12) %in% rownames(snp))/totalsnps2
       #print(paste0(round(tmp,2)*100,"% of Fusion SNPs found in our data"))
       comres$com_totalsnps[j]=tmp
       expr1=as.numeric(t(dat1) %*% coeff1)
       expr2=as.numeric(t(dat2) %*% coeff2)
       comres$expr_cor[j]=cor(expr1,expr2)
       #plot(expr1,expr2,xlab="Our predicted expr",ylab="Fusion predicted expr",main=gene)
       #abline(lm(expr2~expr1),col="red")
       
     }
   }
 }
}

save(comres,fusiontable,file="../result/check_fusion_res.RData")
quantile(abs(comres$expr_cor),na.rm=T)
# 0%         25%         50%         75%        100% 
# 0.002072887 0.386580181 0.539560122 0.730014293 0.994498078 
idx=which.min(abs(comres$expr_cor))
tmp=-log10(comres$r2)
par(mar=c(5,5.5,2,1))
plot(tmp,abs(comres$expr_cor),ylab="expr correlation",xlab="-log10(R^2)")
abline(lm(abs(comres$expr_cor)~ tmp),col="red")
tmp1=-log10(comres$hsq_p)
plot(tmp1,abs(comres$expr_cor),ylab="expr correlation",xlab="-log10(hsq_p)")
abline(lm(abs(comres$expr_cor)~ tmp1),col="red")
hist(abs(comres$expr_cor),main="",xlab="expr correlation",cex.axis=1.3,cex.lab=1.3)
plot(tmp,tmp1)
#check genes only detected by fusion
fusiontable1=fusiontable[fusiontable$avai==T,]
idx=match(models2,fusiontable1$symbol)
models2_=models2[!models2 %in% models1]
idx=match(models2_,fusiontable1$symbol)
gene=fusiontable1$symbol[idx][which.min(fusiontable1$V5[idx])]
quantile(fusiontable$V5[fusiontable$avai])
# 0%        25%        50%        75%       100% 
# 0.0000e+00 1.3445e-08 6.1550e-05 1.7485e-03 9.9930e-03 
quantile(fusiontable1$V5[idx])

res_min1=res_min[rownames(res_min) %in% fusiontable$symbol,]
idx=rownames(res_min1) %in% fusiontable$symbol[fusiontable$avai==T]
hist(res_min1$r2[idx])
hist(res_min1$r2[!idx])
boxplot(res_min1$r2~idx,ylab="R2",xlab="selected by fusion")
idx=which(!rownames(res_min1) %in% fusiontable$symbol[fusiontable$avai==T])
rownames(res_min1)[idx][which.max(res_min1$r2[idx])]
