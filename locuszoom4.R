#For BEACON_Bonn data
#arguments of locuszoom are in /fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom.R

setwd("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/locuszoom/allcovar")
locuszoom="/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom"
library(data.table)
#wget -t 1220 -c --waitretry=28 --wait=33 --random-wait --referer="" --user-agent="" --limit-rate=10000k -e robots=off http://csg.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
#ml Python/2.7.15-foss-2018b #this is needed
###export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin
#export PATH=$PATH:/fh/fast/dai_j/CancerGenomics/Tools/plink #this is needed

#check if you can run locuszoom:
system("/fh/fast/dai_j/CancerGenomics/Tools/locuszoom/bin/locuszoom --help")
#if you can't run it in rstudio, use R

# Chromosome 14, DYNC1H1 and HSP90AA1
# Chromosome 6, FILIP1/COX7A2/THEM30A/SENP6/MYO6

summaryfolder="/fh/fast/dai_j/BEACON/Meta_summary_stat/"
#eqtl of ZNF641
checksnps=c("rs4760632","rs7315258","rs4760698","rs2732480")
BEEA_Bonn=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt"))
BEEA_Camb=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt"))
BEEA_BEACON=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt"))
BEEA_Bonn[BEEA_Bonn$SNP %in% checksnps,]
#               SNP position non_effect_allele effect_allele freq_non_effect_allele_cases freq_non_effect_allele_controls      OR          P     info       BETA        SE freq_effect_allele_controls    N
# 7079125 rs7315258 48572670                 T             C                    0.1521890                       0.1708827 1.14815 0.00425787 0.976040  0.1478790 0.0517348                   0.8291170 6183
# 7079243 rs4760698 48592872                 C             T                    0.1284703                       0.1435523 1.13707 0.00675178 0.986442  0.1491830 0.0550724                   0.8564480 6183
# 7079955 rs2732480 48736303                 C             A                    0.5345079                       0.5442452 1.03997 0.01143120 0.982147  0.0971242 0.0384007                   0.4557550 6183
# 7080783 rs4760632 48934248                 C             T                    0.9386092                       0.9308221 0.88007 0.08651580 0.983163 -0.1303930 0.0760719                   0.0691779 6183
BEEA_Camb[BEEA_Camb$SNP %in% checksnps,]
#         CHR       SNP position non_effect_allele effect_allele     info   all_maf cases_maf controls_maf cohort_1_hwe cases_hwe controls_hwe         P       beta        se    N
# 9310697  12 rs7315258 48572670                 T             C 0.993957 0.1583940 0.1676510    0.1533200     0.179409 0.5072860     0.261508 0.0252078 -0.1322680 0.0589706 5276
# 9310824  12 rs4760698 48592872                 C             T 0.996404 0.1342020 0.1424400    0.1296860     0.109232 0.0710822     0.493699 0.0266165 -0.1406380 0.0632695 5276
# 9311650  12 rs2732480 48736303                 C             A 1.000000 0.4245640 0.4162210    0.4291370     0.516946 0.8120090     0.327578 0.1624600 -0.0608019 0.0435517 5276
# 9312646  12 rs4760632 48934248                 C             T 0.998283 0.0688666 0.0738584    0.0661304     0.388251 1.0000000     0.262967 0.0530887  0.1618630 0.0833610 5276
BEEA_BEACON[BEEA_BEACON$SNP %in% checksnps,]
        # CHR       SNP position non_effect_allele effect_allele     info   all_maf cases_maf controls_maf cohort_1_hwe cases_hwe controls_hwe           P       beta        se     N
# 9454345  12 rs7315258 48572670                 T             C 0.998889 0.1625310 0.1721570    0.1569230     0.971584 0.6134070     0.677937 0.000911006 -0.1325120 0.0398694 10632
# 9454476  12 rs4760698 48592872                 C             T 0.996515 0.1398200 0.1497960    0.1340080     0.106453 0.0785429     0.562747 0.000217519 -0.1559750 0.0420706 10632
# 9455304  12 rs2732480 48736303                 C             A 1.000000 0.4286590 0.4174760    0.4351740     0.968412 0.8437550     0.980200 0.001105800 -0.0976552 0.0299670 10632
# 9456326  12 rs4760632 48934248                 C             T 0.998139 0.0686421 0.0738289    0.0656202     0.879275 1.0000000     0.920566 0.041624300  0.1178990 0.0577071 10632


tmp=intersect(paste0(BEEA_Bonn$SNP,"_",BEEA_Bonn$non_effect_allele,"_",BEEA_Bonn$effect_allele),paste0(BEEA_BEACON$SNP,"_",BEEA_BEACON$non_effect_allele,"_",BEEA_BEACON$effect_allele))
tmp1=BEEA_Bonn[match(tmp,paste0(BEEA_Bonn$SNP,"_",BEEA_Bonn$non_effect_allele,"_",BEEA_Bonn$effect_allele)),]
tmp2=BEEA_BEACON[match(tmp,paste0(BEEA_BEACON$SNP,"_",BEEA_BEACON$non_effect_allele,"_",BEEA_BEACON$effect_allele)),]

BEEA_BEACON=as.data.frame(fread(paste0(summaryfolder,"BEEA_BEACON_autosomes.txt"),header=T))


BE_Oxfordsummarydat=as.data.frame(fread(paste0(summaryfolder,"BE_oxford_autosomes.txt"),header=T))
tmp=intersect("SNP",colnames(BE_Oxfordsummarydat))
if (length(tmp)==0)
{
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="rsid")]="SNP"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="pvalue")]="P"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="beta")]="BETA"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="effect-allele")]="effect_allele"
  colnames(BE_Oxfordsummarydat)[which(colnames(BE_Oxfordsummarydat)=="non-effect-allele")]="non_effect_allele"
}

#add chr based on summarydata BE_Oxfordsummarydat
add_chr_tosummarydat=function(dat1=BE_Bonnsummarydat)
{
  dat2=BE_Oxfordsummarydat
  if (!"chr" %in% colnames(dat1))
  {
    dat1$chr=NA
    tmp=intersect(dat1$SNP,dat2$SNP)
    idx1=match(tmp,dat1$SNP)
    idx2=match(tmp,dat2$SNP)
    dat1$chr[idx1]=dat2$chr[idx2]
    if (is.null(dat1$position))
    {
      dat1$position=NA
      dat1$position[idx1]=dat2$position[idx2]
    }
    # idx=which(is.na(dat1$chr))
    # tmp=intersect(dat1$SNP[idx],dbsnp$V2)
    # idx1=match(tmp,dat1$SNP)
    # idx2=match(tmp,dbsnp$V2)
    # dat1$chr[idx1]=dbsnp$V1[idx2]
    print(paste0(sum(is.na(dat1$chr))/nrow(dat1),"of all snps are missing chr"))
  }
  dat1$snpname=dat1$snpname1=NA
  dat1$snpname=paste0(dat1$chr,":",dat1$position,"_",dat1$effect_allele,"_",dat1$non_effect_allele)
  dat1$snpname1=paste0(dat1$chr,":",dat1$position,"_",dat1$non_effect_allele,"_",dat1$effect_allele)
  return(dat1)
}

# library("biomaRt")
# #mart=useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
# #snpmart = useEnsembl(biomart="snp",host="grch37.ensembl.org",dataset="hsapiens_snp")
# #snpmart = useEnsembl(biomart="snp",dataset="hsapiens_snp") #hg38
# snpmart= useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# findrsid=function(selsnps=NULL,verbose=1)
# {
#   print(paste0(length(selsnps)," snprsid to be found"))
#   chr=unlist(strsplit(selsnps[1],":"))[1]
#   pos=allele1=allele2=rep(NA,length(selsnps))
#   res=data.frame(selsnps=selsnps,rsid=NA,ref=NA,alt=NA,gene=NA)
#   idx=which(is.na(res$rsid))
#   n=0
#   while(n<1 & length(idx)>0)
#   {
#     n=n+1
#     if (n %%2==0) print(n)
#     for (j in idx)
#     {
#       tmp1=unlist(strsplit(selsnps[j],":"))
#       tmp1=unlist(strsplit(tmp1[2],"_"))
#       allalleles=tmp1[2:3]
#       allele1[j]=tmp1[2] #minor allele
#       allele2[j]=tmp1[3]
#       pos[j]=as.numeric(tmp1[1])
#       res$ref[j]=allele2[j]
#       res$alt[j]=allele1[j]
#       #find rsid based on position and allles
#       #tmp=listAttributes(snpmart)
#       #sometimes biomart has connection problems
#       tmp2=tryCatch(
#         {
#           getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand','associated_gene'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = T)
#         },
#         error=function(e)
#         {
#           return(F)
#         }
#       )
#       Sys.sleep(1)
#       if (verbose==1 & class(tmp2)[1]=="logical") print(paste0(selsnps[j]," getBM can't work!"))
#       #if (class(tmp2)[1]=="logical") print(paste0(selsnps[j]," getBM can't work!"))
#       if (class(tmp2)[1]=="data.frame")
#       {
#         for (k in 1:nrow(tmp2))
#         {
#           myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
#           if (all(allalleles %in% myalleles))
#           {
#             res$rsid[j]=tmp2$refsnp_id[k]
#             res$gene[j]=tmp2$associated_gene[k]
#             break
#           }
#         }
#       }
#     }
#     idx=which(is.na(res$rsid))
#   }
#   idx=which(is.na(res$rsid))
#   if (length(idx)>0) print(paste0(length(idx)," snps not found rsid"))
#   return(res)
# }
# 
# hg38tohg19=function(snpnames=rownames(snp)[1:10])
# {
#   library(rtracklayer)
#   library(GenomicRanges)
#   chain=import.chain("/fh/fast/dai_j/CancerGenomics/Tools/database/other/hg38ToHg19.over.chain")
#   tmp=unlist(strsplit(snpnames,":"))
#   chr0=tmp[1]
#   chr=paste0("chr",tmp[1])
#   tmp=tmp[seq(2,length(tmp),2)]
#   tmp=unlist(strsplit(tmp,"_"))
#   pos=as.integer(tmp[seq(1,length(tmp),3)])
#   alt=tmp[seq(2,length(tmp),3)]
#   ref=tmp[seq(3,length(tmp),3)]
#   gr_dat=GRanges(seqnames = chr,ranges=IRanges(start=pos,width=1))
#   tmp=liftOver(gr_dat,chain)
#   newsnpnames=newpos=rep(NA,length(tmp))
#   for (i in 1:length(tmp))
#   {
#     tmp1=unlist(tmp[i])
#     if (length(tmp1)==0)
#     {
#       warning(paste0(snpnames[i]," not transformed"))
#     }else
#     {
#       if (length(tmp1)==1)
#       {
#         newpos[i]=start(tmp1)
#       }else
#       {
#         warning(paste0(snpnames[i]," transformed to ",length(tmp1)," snps"))
#         newpos[i]=start(tmp1)[1]
#       }
#     }
#   }
#   newsnpnames=paste0(chr0,":",newpos,"_",alt,"_",ref)
#   res=data.frame(snphg38=snpnames,snphg19=newsnpnames,stringsAsFactors = F)
#   return(res)
# }
# 
# find_snpname=function(selsnps,summarydat=BE_Bonnsummarydat)
# {
#   selsnps_hg19=hg38tohg19(snpnames = selsnps)
#   selsnps_hg19$selsnps=selsnps
#   selsnps_hg19$snp=NA
#   tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname)
#   idx1=match(tmp,selsnps_hg19$snphg19)
#   idx2=match(tmp,summarydat$snpname)
#   selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
#   tmp=intersect(selsnps_hg19$snphg19,summarydat$snpname1)
#   idx1=match(tmp,selsnps_hg19$snphg19)
#   idx2=match(tmp,summarydat$snpname1)
#   selsnps_hg19$snp[idx1]=summarydat$SNP[idx2]
#   return(selsnps_hg19)
# }


library(data.table)
#Meta analysis of BCA and Bonn
BE_summarydat=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BE_comsnp1.tbl"))
colnames(BE_summarydat)[which(colnames(BE_summarydat)=="MarkerName")]="SNP"
colnames(BE_summarydat)[which(colnames(BE_summarydat)=="Effect")]="BETA"
colnames(BE_summarydat)[which(colnames(BE_summarydat)=="P-value")]="P"
colnames(BE_summarydat)[which(colnames(BE_summarydat)=="Allele1")]="effect_allele"
colnames(BE_summarydat)[which(colnames(BE_summarydat)=="Allele2")]="non_effect_allele"
BE_summarydat=add_chr_tosummarydat(dat1=BE_summarydat)
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Bonn_autosomes_N.txt"))
rsid1=tmp$SNP[which(tmp$freq_non_effect_allele_cases>0.05 & tmp$freq_non_effect_allele_cases<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_Cambridge_autosomes_N.txt"))
rsid2=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BE_BEACON_autosomes_N.txt"))
rsid3=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
rsid=intersect(rsid1,rsid2)
rsid=intersect(rsid,rsid3)
BE_summarydat=BE_summarydat[BE_summarydat$SNP %in% rsid,]
BE_summarydat=BE_summarydat[!is.na(BE_summarydat$position),]

BEEA_summarydat=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BEEA_comsnp1.tbl"))
colnames(BEEA_summarydat)[which(colnames(BEEA_summarydat)=="MarkerName")]="SNP"
colnames(BEEA_summarydat)[which(colnames(BEEA_summarydat)=="Effect")]="BETA"
colnames(BEEA_summarydat)[which(colnames(BEEA_summarydat)=="P-value")]="P"
colnames(BEEA_summarydat)[which(colnames(BEEA_summarydat)=="Allele1")]="effect_allele"
colnames(BEEA_summarydat)[which(colnames(BEEA_summarydat)=="Allele2")]="non_effect_allele"
BEEA_summarydat=add_chr_tosummarydat(dat1=BEEA_summarydat)
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt"))
rsid1=tmp$SNP[which(tmp$freq_non_effect_allele_cases>0.05 & tmp$freq_non_effect_allele_cases<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt"))
rsid2=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt"))
rsid3=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
rsid=intersect(rsid1,rsid2)
rsid=intersect(rsid,rsid3)
BEEA_summarydat=BEEA_summarydat[BEEA_summarydat$SNP %in% rsid,]
BEEA_summarydat=BEEA_summarydat[!is.na(BEEA_summarydat$position),]
dim(BEEA_summarydat)

BEEA_summarydat1=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_BEEA_comsnp_new1.tbl"))
BEEA_summarydat1[BEEA_summarydat1$MarkerName %in% checksnps,]
colnames(BEEA_summarydat1)[which(colnames(BEEA_summarydat1)=="MarkerName")]="SNP"
colnames(BEEA_summarydat1)[which(colnames(BEEA_summarydat1)=="Effect")]="BETA"
colnames(BEEA_summarydat1)[which(colnames(BEEA_summarydat1)=="P-value")]="P"
colnames(BEEA_summarydat1)[which(colnames(BEEA_summarydat1)=="Allele1")]="effect_allele"
colnames(BEEA_summarydat1)[which(colnames(BEEA_summarydat1)=="Allele2")]="non_effect_allele"
BEEA_summarydat1=add_chr_tosummarydat(dat1=BEEA_summarydat1)
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Bonn_autosomes_N.txt"))
rsid1=tmp$SNP[which(tmp$freq_non_effect_allele_cases>0.05 & tmp$freq_non_effect_allele_cases<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_Cambridge_autosomes_N.txt"))
rsid2=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/BEEA_BEACON_autosomes_N.txt"))
rsid3=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
rsid=intersect(rsid1,rsid2)
rsid=intersect(rsid,rsid3)
BEEA_summarydat1=BEEA_summarydat1[BEEA_summarydat1$SNP %in% rsid,]
BEEA_summarydat1=BEEA_summarydat1[!is.na(BEEA_summarydat1$position),]
dim(BEEA_summarydat1)
BEEA_summarydat1[BEEA_summarydat1$SNP %in% checksnps,]

EA_summarydat=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/Beacon_Bonn_Cambridge_METAANALYSIS_EA_comsnp1.tbl"))
colnames(EA_summarydat)[which(colnames(EA_summarydat)=="MarkerName")]="SNP"
colnames(EA_summarydat)[which(colnames(EA_summarydat)=="Effect")]="BETA"
colnames(EA_summarydat)[which(colnames(EA_summarydat)=="P-value")]="P"
colnames(EA_summarydat)[which(colnames(EA_summarydat)=="Allele1")]="effect_allele"
colnames(EA_summarydat)[which(colnames(EA_summarydat)=="Allele2")]="non_effect_allele"
EA_summarydat=add_chr_tosummarydat(dat1=EA_summarydat)
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Bonn_autosomes_N.txt"))
rsid1=tmp$SNP[which(tmp$freq_non_effect_allele_cases>0.05 & tmp$freq_non_effect_allele_cases<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_Cambridge_autosomes_N.txt"))
rsid2=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
tmp=as.data.frame(fread("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/EA_BEACON_autosomes_N.txt"))
rsid3=tmp$SNP[which(tmp$all_maf>0.05 & tmp$all_maf<0.95)]
rsid=intersect(rsid1,rsid2)
rsid=intersect(rsid,rsid3)
EA_summarydat=EA_summarydat[EA_summarydat$SNP %in% rsid,]
EA_summarydat=EA_summarydat[!is.na(EA_summarydat$position),]
dim(EA_summarydat)

# zl_validaton1=function(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,
#                        res_min=NULL,ymax=NULL,nsnp=NULL)
# {
#   # idx=which(rownames(phenotypepos)==gene)
#   # chr=phenotypepos$chr[idx]
#   # startpos=phenotypepos$s1[idx]-5e5
#   # endpos=phenotypepos$s2[idx]+5e5
#   METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
#   if (!file.exists(METALfile))
#   {
#     #summarydat=as.data.frame(fread(summaryfile,header=T))
#     idx=which(rownames(res_min) ==gene)
#     selsnps=unique(unlist(strsplit(res_min$selectedsnps[idx],"|",fixed=T)))
#     
#     tmp=find_snpname(selsnps,summarydat=BE_Oxfordsummarydat)
#     #tmp=find_snpname(selsnps,summarydat=BE_Bonnsummarydat)
#     idx=which(is.na(tmp$snp))
#     if (length(idx)>0) #use biomart, slower method
#     {
#       pos=allele1=allele2=rep(NA,length(selsnps))
#       for (j in idx)
#       {
#         #cat(j,"\n")
#         tmp1=unlist(strsplit(tmp$selsnps[j],":"))
#         chr=tmp1[1]
#         tmp1=unlist(strsplit(tmp1[2],"_"))
#         allalleles=tmp1[2:3]
#         allele1[j]=tmp1[2] #minor allele
#         allele2[j]=tmp1[3]
#         pos[j]=as.numeric(tmp1[1])
#         #find rsid based on position and allles
#         #tmp=listAttributes(snpmart)
#         #sometimes biomart has connection problems
#         tmp2=tryCatch(
#           {
#             getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
#           },
#           error=function(e)
#           {
#             return(F)
#           }
#         )
#         
#         #tmp2=getBM(attributes=c('refsnp_id','allele','allele_1','minor_allele','chrom_start','chrom_end','chrom_strand'), filters = c('chr_name',"start","end"), values =list(chr,pos[j],pos[j]),mart = snpmart,useCache = TRUE)
#         if (class(tmp2)[1]=="logical") print(paste0(tmp$selsnps[j]," getBM can't work!"))
#         if (class(tmp2)[1]=="data.frame")
#         {
#           for (k in 1:nrow(tmp2))
#           {
#             myalleles=unlist(strsplit(tmp2$allele[k],"/",fixed=T))
#             if (all(allalleles %in% myalleles))
#             {
#               tmp$snp[j]=tmp2$refsnp_id[k]
#               break
#             }
#           }
#         }
#       }
#     }
#     selsnpstr0=rep(nrow(tmp)) #selsnps are in hg38, selsnpstr are in hg19, only use position
#     tmp$snphg19str0=NA
#     for (j in 1:nrow(tmp))
#     {
#       tmp$snphg19str0[j]=selsnpstr0[j]=unlist(strsplit(tmp$snphg19[j],"_"))[1]
#     }
#     selsnpmap=tmp
#     rsid=selsnpmap$snp
#     
#     if (sum(is.na(rsid)>0)) print("some rsid not been found!")
#     idx0=which(!rsid %in% summarydat$SNP)
#     missingrsid=rsid[!rsid %in% summarydat$SNP]
#     chr=unlist(strsplit(selsnps[1],":"))[1]
#     bcagenotype1=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc_maf001_snp/chr",chr,"_filter.traw")))
#     #bcagenotype1=bcagenotype1[,7:ncol(bcagenotype1)]
#     bcabim1=as.data.frame(fread(paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/imputation/beacondbgapcontrol_qc_hrc_maf001_snp/chr",chr,"_filter.bim")))
#     bcabim1str0=paste0(bcabim1$V1,":",bcabim1$V4)
#     summarydatstr=paste0(summarydat$chr,":",summarydat$position)
#     missingsnps=selsnpstr0[which(!selsnpstr0 %in% summarydatstr)]
#     tmp=intersect(selsnpstr0,summarydatstr)
#     if (length(missingrsid)>0)
#     {
#       #refreplacesnps1=data.frame(missingid=missingrsid,hg19pos=selsnpmap$snphg19[idx0],replaceid=NA,cor=NA,dist=NA)
#       refreplacesnps1=data.frame(snpstr0=missingsnps,replacestr0=NA,dist=NA,cor=NA)
#       for (k in 1:nrow(refreplacesnps1))
#       {
#         tmp1=as.integer(unlist(strsplit(refreplacesnps1$snpstr0[k],":"))[2])
#         idx1=which(bcabim1str0==refreplacesnps1$snpstr0[k])
#         idx_keep=which(abs(bcabim1$V4-tmp1)<=5000)
#         idx2=which(bcabim1str0[idx_keep] %in% summarydatstr)
#         if (length(idx2)>0)
#         {
#           tmp2=t(bcagenotype1[idx_keep[idx2],7:ncol(bcagenotype1)])
#           tmp3=cor(tmp2,unlist(bcagenotype1[idx1,7:ncol(bcagenotype1)]))
#           names(tmp3)=bcabim1str0[idx_keep[idx2]]
#           idx2=which.max(tmp3)
#           idx=which(summarydat$chr==chr)
#           idx3=which.min(abs(summarydat$position[idx]-tmp1))
#           refreplacesnps1$replacestr0[k]=names(tmp3)[idx2]
#           tmp4=as.integer(unlist(strsplit(names(tmp3)[idx2],":"))[2])
#           refreplacesnps1$dist[k]=abs(tmp4-tmp1)
#           refreplacesnps1$cor[k]=tmp3[idx2]
#         }
#       }
#       
#       tmp0=intersect(selsnpstr0,summarydatstr)
#       if (sum(refreplacesnps1$cor>0.6 & refreplacesnps1$dist<=50,na.rm=T)>0)
#         tmp=c(tmp0,refreplacesnps1$replacestr0[refreplacesnps1$cor>0.6 & refreplacesnps1$dist<=50])
#       
#     }
#     idx=match(tmp,summarydatstr)
#     allrsid=summarydat$SNP[idx[!is.na(idx)]]
#     
#     
#     #allrsid=c(rsid,refreplacesnps1$replaceid[which(refreplacesnps1$dist<=50)])
#     summarydat1=summarydat[summarydat$SNP %in% allrsid,]
#     
#     METALdat=data.frame(MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
#                         alt=summarydat1$effect_allele,color="red",stringsAsFactors = F)
#     
#     METALfile=paste0("./",prefix,"_",gene,"_metalfile.txt")
#     write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
#   }#else
#   {
#     METALdat=read.table(METALfile,header=T)
#   }
#   
#   #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
#   #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
#   #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
#   #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
#   #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
#   #print(cmd)
#   if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
#   #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
#   if (is.null(nsnp))
#   {
#     title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
#   }else
#   {
#     title=paste0(prefix,", ",gene,", ",nsnp," SNPs")
#   }
#   #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
#   #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
#   cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
#   
#   print(cmd)
#   system(cmd,wait=T)
#   return(0)
# }
# locusplot_validation1=function(organ="blood",genes=c("SENP6"),prefix="BEEA_Bonn",summarydat=BEEA_Bonnsummarydat,ymax=NULL,nsnp=NULL)
# {
#   prefix=paste0(organ,"_",prefix)
#   outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor_allcovar")
#   # modeldat=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/","GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction.RData")
#   # load(modeldat) #phenotypepos
#   load(paste0(outfolder,"/preidiction_michigan_model.RData")) #res_min
#   # skatres=skatres[rownames(skatres) %in% rownames(res_min)[which(res_min$r2>0.01)],]
#   
#   # load(modelres)
#   # load(skatres)
#   for (gene in genes)
#   {
#     print(gene)
#     zl_validaton1(gene=gene,summarydat = summarydat,prefix=prefix,res_min = res_min,ymax=ymax,nsnp=nsnp)
#   }
# }
# 
# locusplot_validation1(organ="adipose",genes="EXOC3",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=7,nsnp=26)
# locusplot_validation1(organ="blood",genes="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,ymax=5,nsnp=89)
# locusplot_validation1(organ="adipose",genes="KRTAP5-8",prefix="EA_BCA_Bonn",summarydat=EA_summarydat,ymax=5,nsnp=20)
# locusplot_validation1(organ="junction",genes="ZNF641",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5,nsnp=4)
# locusplot_validation1(organ="blood",genes="SENP6",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5,nsnp=48)
# locusplot_validation1(organ="stomach",genes="CFDP1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5,nsnp=5)
# locusplot_validation1(organ="junction",genes="CHST5",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,ymax=5,nsnp=20)
# locusplot_validation1(organ="blood",genes="BCAR1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,ymax=5,nsnp=328)
# 
# 
# 
# locusplot_validation1(organ="junction",genes="LDAH",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="GDF7",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="stomach",genes="YEATS2",prefix="EA_BCA_Bonn",summarydat=EA_summarydat,ymax=5)
# locusplot_validation1(organ="adipose",genes="ABCF3",prefix="EA_BCA_Bonn",summarydat=EA_summarydat,ymax=5)
# locusplot_validation1(organ="adipose",genes="BARX1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="adipose",genes="ALDH1A2",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="AQP9",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="JUND",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="SSBP4",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="ISYNA1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="blood",genes="KLHL26",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="adipose",genes="CRTC1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)
# locusplot_validation1(organ="mucosa",genes="TMEM161A",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat,ymax=5)


#this is to use meta-gwas results, not just eqtl

zl_validaton2=function(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,
                       ymax=8,nsnp=NULL)
{
  
  # idx=which(rownames(phenotypepos)==gene)
  # chr=phenotypepos$chr[idx]
  # startpos=phenotypepos$s1[idx]-5e5
  # endpos=phenotypepos$s2[idx]+5e5
  METALfile=paste0("./",prefix,"_",gene,"_metalfile1.txt")
  #if (!file.exists(METALfile))
  {
    idx=which(refgene$name2==gene)
    gr_gene=gr_gtexgene[idx]
    gr_summarydat=GRanges(seqnames = paste0("chr",summarydat$chr),ranges=IRanges(start=summarydat$position,width=1))
    tmp=distance(gr_summarydat,gr_gene)
    sum(tmp<5e5,na.rm=T)
    allrsid=summarydat$SNP[which(tmp<5e5)]
    summarydat1=summarydat[summarydat$SNP %in% allrsid,]
    mystart=min(summarydat1$position)
    myend=max(summarydat1$position)
    mychr=summarydat1$chr[1]
    METALdat=data.frame(MarkerName=summarydat1$SNP,Pvalue=summarydat1$P,ref=summarydat1$non_effect_allele,
                        alt=summarydat1$effect_allele,color="red",stringsAsFactors = F)
    idx=which(METALdat$Pvalue<1e-12)
    if (length(idx)>0)
    {
      METALdat$Pvalue[idx]=METALdat$Pvalue*1000
    }
    idx=which(METALdat$Pvalue<1e-10)
    if (length(idx)>0)
    {
      METALdat$Pvalue[idx]=METALdat$Pvalue*1000
    }
    idx=which(METALdat$Pvalue<1e-9)
    if (length(idx)>0)
    {
      METALdat$Pvalue[idx]=METALdat$Pvalue*1000
    }
    idx=which(METALdat$Pvalue<1e-8)
    if (length(idx)>0)
    {
      METALdat$Pvalue[idx]=METALdat$Pvalue*100
    }
    idx=which(METALdat$Pvalue<1e-7)
    if (length(idx)>0)
    {
      METALdat$Pvalue[idx]=METALdat$Pvalue*5
    }
    
    METALfile=paste0("./",prefix,"_",gene,"_metalfile1.txt")
    write.table(METALdat,file=METALfile,col.names=T,row.names = F,sep="\t",quote=F)
  }#else
  {
    METALdat=read.table(METALfile,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  if (is.null(nsnp))
  {
    title=paste0(gene)
  }else
  {
    title=paste0(gene)
  }
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
  cmd=paste0(locuszoom," --metal ",METALfile," --flank 10kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 100kb"," --refgene ",gene," --chr ",mychr, " --start ",mystart, " --end ",myend, " --build hg38 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}
library(GenomicRanges)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/gtexv8_ge_anno.RData")
gtexv8_ge_anno1=gtexv8_ge_anno[gtexv8_ge_anno$V3=="gene",]
idx=which(duplicated(gtexv8_ge_anno1$Symbol))
gtexv8_ge_anno1=gtexv8_ge_anno1[-idx,]
gr_gtexgene=GRanges(seqnames = gtexv8_ge_anno1$Chromosome,ranges = IRanges(start=gtexv8_ge_anno1$start,end=gtexv8_ge_anno1$end))
refgene=read.table("/fh/fast/dai_j/CancerGenomics/Tools/database/annotation/refgene_hg19",comment.char = "!",header=T)
gr_gtexgene=GRanges(seqnames = refgene$chrom,ranges = IRanges(start=refgene$txStart,end=refgene$txEnd))
locusplot_validation2=function(genes=c("SENP6"),prefix="BEEA_Bonn",summarydat=BEEA_Bonnsummarydat,ymax=8,nsnp=NULL)
{
  prefix=prefix
  for (gene in genes)
  {
    print(gene)
    zl_validaton2(gene=gene,summarydat = summarydat,prefix=prefix,ymax=ymax)
  }
}

locusplot_validation2(genes="EXOC3",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
locusplot_validation2(genes="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
locusplot_validation2(genes="KRTAP5-8",prefix="EA_BCA_Bonn",summarydat=EA_summarydat)
locusplot_validation2(genes="ZNF641",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat1)
locusplot_validation2(genes="SENP6",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
locusplot_validation2(genes="CFDP1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
locusplot_validation2(genes="CHST5",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
locusplot_validation2(genes="BCAR1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)


#mark eqtl as red, gwas as blue
zl_validaton3=function(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,
                       ymax=8,nsnp=NULL)
{
  
  # idx=which(rownames(phenotypepos)==gene)
  # chr=phenotypepos$chr[idx]
  # startpos=phenotypepos$s1[idx]-5e5
  # endpos=phenotypepos$s2[idx]+5e5
  METALfile0=paste0("./",prefix,"_",gene,"_metalfile.txt")
  METALfile1=paste0("./",prefix,"_",gene,"_metalfile1.txt")
  METALfile2=paste0("./",prefix,"_",gene,"_metalfile2.txt")
  #if (!file.exists(METALfile2))
  {
    METALdat0=read.table(METALfile0,header=T)
    METALdat1=read.table(METALfile1,header=T)
    insnps=intersect(METALdat0$MarkerName,METALdat1$MarkerName)
    addsnps=METALdat0$MarkerName[!METALdat0$MarkerName %in% METALdat1$MarkerName]
    noeqtlsnps=METALdat1$MarkerName[!METALdat1$MarkerName %in% METALdat0$MarkerName]
    if (length(addsnps)>0)
    {
      idx=match(addsnps,METALdat0$MarkerName)
      METALdat=rbind(METALdat1,METALdat0[idx,])
    }else
    {
      METALdat=METALdat1
    }
    METALdat$color="red"
    idx=match(noeqtlsnps,METALdat$MarkerName)
    METALdat$color[idx]="gray"
    #METALdat=rbind(METALdat[which(METALdat$color=="gray"),],METALdat[which(METALdat$color=="red"),])
    METALdat=rbind(METALdat[which(METALdat$color=="red"),],METALdat[which(METALdat$color=="gray"),])
    write.table(METALdat,file=METALfile2,col.names=T,row.names = F,sep="\t",quote=F)
  }#else
  {
    METALdat=read.table(METALfile2,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  if (is.null(nsnp))
  {
    title=paste0(gene)
  }else
  {
    title=paste0(gene)
  }
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
  cmd=paste0(locuszoom," --metal ",METALfile2," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9 rugAlpha=0.5 ")
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 100kb"," --refgene ",gene," --chr ",mychr, " --start ",mystart, " --end ",myend, " --build hg38 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

zl_validaton3(gene="EXOC3",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
zl_validaton3(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
zl_validaton3(gene="KRTAP5-8",prefix="EA_BCA_Bonn",summarydat=EA_summarydat)
zl_validaton3(gene="ZNF641",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat1)
zl_validaton3(gene="SENP6",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
zl_validaton3(gene="CFDP1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
zl_validaton3(gene="CHST5",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
zl_validaton3(gene="BCAR1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)

#remove some gwas snp wich is very close to eqtl
zl_validaton4=function(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat,
                       ymax=8,nsnp=NULL)
{
  
  # idx=which(rownames(phenotypepos)==gene)
  # chr=phenotypepos$chr[idx]
  # startpos=phenotypepos$s1[idx]-5e5
  # endpos=phenotypepos$s2[idx]+5e5
  METALfile0=paste0("./",prefix,"_",gene,"_metalfile.txt")
  METALfile1=paste0("./",prefix,"_",gene,"_metalfile1.txt")
  METALfile2=paste0("./",prefix,"_",gene,"_metalfile2.txt")
  #if (!file.exists(METALfile2))
  {
    METALdat0=read.table(METALfile0,header=T)
    METALdat1=read.table(METALfile1,header=T)
    insnps=intersect(METALdat0$MarkerName,METALdat1$MarkerName)
    addsnps=METALdat0$MarkerName[!METALdat0$MarkerName %in% METALdat1$MarkerName]
    noeqtlsnps=METALdat1$MarkerName[!METALdat1$MarkerName %in% METALdat0$MarkerName]
    if (length(addsnps)>0)
    {
      idx=match(addsnps,METALdat0$MarkerName)
      METALdat=rbind(METALdat1,METALdat0[idx,])
    }else
    {
      METALdat=METALdat1
    }
    METALdat$color="red"
    idx=match(noeqtlsnps,METALdat$MarkerName)
    METALdat$color[idx]="gray"
    #METALdat=rbind(METALdat[which(METALdat$color=="gray"),],METALdat[which(METALdat$color=="red"),])
    METALdat=rbind(METALdat[which(METALdat$color=="red"),],METALdat[which(METALdat$color=="gray"),])
    tmp1=METALdat[METALdat$color=="red",]
    tmp2=METALdat[METALdat$color=="gray",]
    tmp11=intersect(tmp1$MarkerName,BE_Oxfordsummarydat$SNP)
    idx=match(tmp11,BE_Oxfordsummarydat$SNP)
    tmp1$chr=BE_Oxfordsummarydat$chr[idx]
    tmp1$position=BE_Oxfordsummarydat$position[idx]
    tmp1=tmp1[tmp1$MarkerName %in% tmp11,]
    tmp22=intersect(tmp2$MarkerName,BE_Oxfordsummarydat$SNP)
    tmp2=tmp2[tmp2$MarkerName %in% tmp22,]
    idx=match(tmp2$MarkerName,BE_Oxfordsummarydat$SNP)
    tmp2$chr=BE_Oxfordsummarydat$chr[idx]
    tmp2$position=BE_Oxfordsummarydat$position[idx]
    gr_tmp1=GRanges(seqnames = tmp1$chr,ranges = IRanges(start = tmp1$position,width = 1))
    gr_tmp2=GRanges(seqnames = tmp2$chr,ranges = IRanges(start = tmp2$position,width = 1))
    tmp=distanceToNearest(gr_tmp2,gr_tmp1)
    sum(tmp@elementMetadata$distance>1000)
    idx=which(tmp@elementMetadata$distance>10000)
    #idx=which(tmp@elementMetadata$distance>2000) #BCAR1
    tmp2=tmp2[idx,]
    METALdat=rbind(tmp1,tmp2)
    write.table(METALdat,file=METALfile2,col.names=T,row.names = F,sep="\t",quote=F)
  }#else
  {
    METALdat=read.table(METALfile2,header=T)
  }
  
  #markerfile=paste0("./",prefix,"_",gene,"_markerfile.txt")
  #write.table(markerdat,file=markerfile,col.names=T,row.names = F,sep="\t",quote=F)
  #cmd=paste0("locuszoom --metal ",METALfile," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2010 --flank 500k --pvalcol Pvalue")
  #cmd=paste0("locuszoom --metal ",METALfile," --denote-markers-file ",markerfile," --chr ",chr," --start ",startpos," --end ",endpos," --build hg19 --pop EUR --source 1000G_Nov2010 --pvalcol Pvalue --plotonly")
  #cmd=paste0(locuszoom," --metal ",METALfile," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --no-ld --prefix ",prefix," colorCol=color")
  #print(cmd)
  if (is.null(ymax)) ymax=as.integer(max(-log10(METALdat$Pvalue),na.rm=T)+1)
  #title=paste0(prefix,", ",gene,", ",nrow(METALdat)," SNPs")
  if (is.null(nsnp))
  {
    title=paste0(gene)
  }else
  {
    title=paste0(gene)
  }
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title)
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, "axisTextSize=1.3 refsnpTextSize=1.3 geneFontSize=1.3 legendSize=1.3 axisSize=1.3")
  cmd=paste0(locuszoom," --metal ",METALfile2," --flank 600kb"," --refgene ",gene," --build hg19 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9 rugAlpha=0.5 ")
  #cmd=paste0(locuszoom," --metal ",METALfile," --flank 100kb"," --refgene ",gene," --chr ",mychr, " --start ",mystart, " --end ",myend, " --build hg38 --pop EUR --source 1000G_Nov2014 --pvalcol Pvalue --plotonly --prefix ",prefix," colorCol=color refsnpTextColor=transparent ","ymax=",ymax," signifLine=",-log10(0.05)," signifLineColor=red title=",title, " axisTextSize=1.4 legendSize=1.4  axisSize=1.4 xlabPos=-2.9")
  
  print(cmd)
  system(cmd,wait=T)
  return(0)
}

#change order(gray/red) and run again
zl_validaton4(gene="KRTAP5-8",prefix="EA_BCA_Bonn",summarydat=EA_summarydat)
zl_validaton4(gene="CHST5",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
zl_validaton4(gene="CFDP1",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat)
zl_validaton4(gene="HSP90AA1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)
zl_validaton4(gene="ZNF641",prefix="BEEA_BCA_Bonn",summarydat=BEEA_summarydat1)

zl_validaton4(gene="BCAR1",prefix="BE_BCA_Bonn",summarydat=BE_summarydat)