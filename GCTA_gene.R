#!/usr/bin/env Rscript
#salloc -t 6-1 --constraint=gizmok -n 51 mpirun -n 1 R --interactive
#salloc -t 6-1 --constraint=gizmok -n 41 mpirun -n 1 Rscript ./GCTA_gene.R junction 1grm &
#salloc -t 6-1 --constraint=gizmok -n 51 mpirun -n 1 Rscript ./GCTA_gene.R adipose 1grm &
#salloc -t 6-1 --constraint=gizmok -n 41 mpirun -n 1 Rscript ./GCTA_gene.R muscularis 1grm &
#salloc -t 6-1 --constraint=gizmok -n 41 mpirun -n 1 Rscript ./GCTA_gene.R stomach 1grm &
#salloc -t 6-1 --constraint=gizmok -n 41 mpirun -n 1 Rscript ./GCTA_gene.R mucosa 1grm &
#salloc -t 6-1 --constraint=gizmok -n 51 mpirun -n 1 Rscript ./GCTA_gene.R blood 1grm &
#Use GCTA to estimate heritability at gene level
args = commandArgs(trailingOnly=TRUE)
organ=args[1]
#one of mgrm,1grm,ai
gcta_opt=args[2]
plink="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/plink-1.07-x86_64/plink1.9/plink" ##
gctafolder="/fh/fast/dai_j/CancerGenomics/Tools/gcta_1.93.2beta/" ##
#organ="junction" #275
#organ="adipose" #393
#organ="muscularis" #385
#organ="stomach" #260
#organ="mucosa" #411
#organ="blood" #558
prefix=paste0("GTExV8",organ,"data_ambiguous_TPM_addcontrols_HRC_maf005_for_prediction")
outfolder=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_",organ,"_HRC_MAF005_rmhighcor/") ##
#save result
# heritfile=paste0(outfolder,"heritability_nopeer_",gcta_opt,".txt")
heritfile=paste0(outfolder,"heritability_",gcta_opt,".txt") #for used GTEx snps
heritfile=paste0(outfolder,"heritability_all_",gcta_opt,".txt") #for all GTEx snps
print(heritfile)
# oldres=NULL
# if (file.exists(heritfile))
# {
#   oldres=read.table(heritfile)
# }

genotypefolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/data/GTEx/genotype/" ##
gtexprefix="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_flip_chr"  ## this is to include all the GTEx genotype
#generate covariates
rdata=paste0("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/",prefix,".RData")
load(rdata)
dim(phenotype)
if(any(colnames(phenotype)!=rownames(covariate))) warning("check the phenotype and covariate data")
if(any(rownames(phenotype)!=rownames(phenotypepos))) warning("check the phenotype and phenotypepos data")
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
allsnp=snp
allsnppos=snppos
gr_allsnp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP ##
gr_allpos=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) ##


#generate the cov input,gender,pcr,platform
generate_cov=function(outfile=paste0(outfolder,"GCTA.covar"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=5))
  res[,1]=res[,2]=rownames(covariate)
  res[,3]=as.character(covariate$gender)
  res[res[,3]==1,3]="M"
  res[res[,3]==2,3]="F"
  res[,4]=as.character(covariate$pcr)
  res[res[,4]==0,4]="PCR0"
  res[res[,4]==1,4]="PCR1"
  res[,5]=as.character(covariate$platform)
  res[res[,5]==0,5]="Platform0"
  res[res[,5]==1,5]="Platform1"
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}

#qov,pcs+factors+age 
generate_qcov=function(outfile=paste0(outfolder,"GCTA.qcovar"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=3))
  res[,1]=res[,2]=rownames(covariate)
  res[,3]=covariate$age
  for (i in 1:4)
  {
    res=cbind(res,covariate[,paste0("pc",i)])
  }
  #no peer
  for (i in 1:15)
  {
    res=cbind(res,covariate[,paste0("factor",i)])
  }
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}

#generate fam select files to keep samples
generate_fam=function(outfile=paste0(outfolder,"GCTA.select.fam"))
{
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=2))
  res[,1]=res[,2]=rownames(covariate)
  write.table(res,file=outfile,sep=" ",row.names = F,col.names = F,quote=F)
}
generate_cov()
generate_qcov()
generate_fam()

#genename="LDAH"

readhsq=function(hsqfile=paste0(outfolder,"LDAH.hsq"),opt="mgrm")
{
  tmp1=NULL
  if (file.exists(hsqfile))
  {
    tmp=read.table(hsqfile,sep="\t",fill=T,stringsAsFactors = F,header = T)
    if (opt=="1grm")
    {
      idx1=which(tmp$Source=="V(G)/Vp")
      idx2=which(tmp$Source=="Pval")
      if (length(idx1)>0 & length(idx2)>0)
      tmp1=data.frame(H=tmp$Variance[c(idx1)],SE=tmp$SE[c(idx1)],pvalue=tmp$Variance[idx2],stringsAsFactors = F)
    }else #multiple grms
    {
      idx=which(grepl("Sum",tmp[,1]))
      if(length(idx)>0)
      tmp1=data.frame(H=tmp$Variance[idx],SE=tmp$SE[idx],pvalue=tmp$Variance[nrow(tmp)-1],stringsAsFactors = F)
    }
  }
  return(tmp1)
}

#estimate heritarbility
estimate_heri=function(genename="LDAH",distcutoff=5e5,opt="1grm")
{
  allres=NULL
  #phenotype file
  phenofile=paste0(outfolder,genename,".pheno")
  res=data.frame(matrix(NA,nrow=nrow(covariate),ncol=3))
  res[,1]=res[,2]=rownames(covariate)
  idx=which(rownames(phenotype)==genename)
  res[,3]=unlist(phenotype[idx,])
  write.table(res,file=phenofile,sep=" ",row.names = F,col.names = F,quote=F)
  
  idx=which(rownames(phenotypepos)==genename)
  tmp=distance(gr_snp,gr_pos[idx])
  idx1=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx1,]))
  idx1=idx1[tmp!=0] #remove all 0 genotypes
  if (length(idx1)>10)
  {
    bimfile=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx],".bim")
    bim=as.data.frame(data.table::fread(bimfile))
    idx2=match(snppos$pos[idx1],bim$V4)
    snpfilename=paste0(outfolder,genename,".snpname")
    tmp=data.frame(snps=bim$V2[idx2]) #for select snps
    #tmp=data.frame(snps=bim) #for all snps
    write.table(tmp,file=snpfilename,row.names = F,col.names = F,quote=F)
    famfile=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx],".fam")
    fam=read.table(famfile,stringsAsFactors = F)
    select.fam=paste0(outfolder,"GCTA.select.fam")
    
    prefix1=paste0(genotypefolder,gtexprefix,phenotypepos$chr[idx])
    genofilename=paste0(outfolder,genename)
    cmd=paste0(plink," --bfile ",prefix1," --keep ",select.fam," --extract ",snpfilename," --make-bed --out ", genofilename)
    system(cmd)
    #genotype data
    cmd=paste0(plink," --bfile ",prefix1," --keep ",select.fam," --extract ",snpfilename," --recode A-transpose --out ", genofilename)
    system(cmd)
    tmp=read.table(paste0(genofilename,".fam"))
    if(any(tmp$V1!=rownames(covariate))) warning("the fam file has different sample order!")
    
    #run gcta, use mgrm
    if (opt=="mgrm")
    {
      cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --ld-score-region 50 --thread-num 12 --out ",genofilename)
      system(cmd)
      cmd=paste0("./stratify_GCTAsnps.R ",genofilename)
      system(cmd)
      if (file.exists(paste0(genofilename,".score.ld")))
      {
        if (file.size(paste0(genofilename,"_group1.txt"))>0 & file.size(paste0(genofilename,"_group2.txt"))>0 &file.size(paste0(genofilename,"_group3.txt"))>0 & file.size(paste0(genofilename,"_group4.txt"))>0)
        {
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group1.txt --make-grm --thread-num 12 --out ",genofilename,"_group1")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group2.txt --make-grm --thread-num 12 --out ",genofilename,"_group2")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group3.txt --make-grm --thread-num 12 --out ",genofilename,"_group3")
          system(cmd)
          cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --extract ",genofilename,"_group4.txt --make-grm --thread-num 12 --out ",genofilename,"_group4")
          system(cmd)
          fileConn<-file(paste0(genofilename,"_mult_GRMs.txt"))
          writeLines(paste0(genofilename,"_group1\n",genofilename,"_group2\n",genofilename,"_group3\n",genofilename,"_group4"), fileConn)
          close(fileConn)
          cmd=paste0(gctafolder,"gcta64 --reml --mgrm ",genofilename,"_mult_GRMs.txt --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --reml-maxit 1000 --out ",genofilename)
          system(cmd)
          
          allres=readhsq(hsqfile = paste0(outfolder,genename,".hsq"))
          if (!is.null(allres)) rownames(allres)=genename
        }
      }
      cmd=paste0("rm ",outfolder,genename,".*")
      system(cmd)
      cmd=paste0("rm ",outfolder,genename,"_*.*")
      system(cmd)
    }

    #run gcta, use REML
    if (opt=="1grm")
    {
      cmd=paste0(gctafolder,"gcta64 --bfile ",genofilename," --make-grm --thread-num 12 --out ",genofilename)
      system(cmd)
      #use AI by setting reml-alg 0, it is the default value
      #--reml-alg 0
      #Specify the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML, if this option is not specified.
      #remove constrain
      #cmd=paste0(gctafolder,"gcta64 --reml --reml-alg 0 --reml-no-constrain --reml-maxit 200 --grm ",genofilename," --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --out ",genofilename)
      #use other algorithms
      #cmd=paste0(gctafolder,"gcta64 --reml --reml-alg 2 --reml-maxit 200 --grm ",genofilename," --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --out ",genofilename)
      #remove covariates
      #cmd=paste0(gctafolder,"gcta64 --reml --reml-alg 0 --reml-maxit 200 --grm ",genofilename," --pheno ",genofilename,".pheno --thread-num 12 --out ",genofilename)
      cmd=paste0(gctafolder,"gcta64 --reml --reml-alg 0 --grm ",genofilename," --pheno ",genofilename,".pheno --covar ",outfolder,"GCTA.covar --qcovar ",outfolder,"GCTA.qcovar --thread-num 12 --out ",genofilename)
      system(cmd)
      allres=readhsq(hsqfile = paste0(outfolder,genename,".hsq"),opt="1grm")
      if (!is.null(allres)) rownames(allres)=genename
    }
    cmd=paste0("rm ",outfolder,genename,".*")
    system(cmd)
  }

  return(allres)
}

allgenes=rownames(phenotype)
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/data/gtexv8_ge_anno.RData")
proteingenes=unique(gtexv8_ge_anno$Symbol[gtexv8_ge_anno$gene_type=="protein_coding"])
allgenes=allgenes[allgenes %in% proteingenes]
# 
# 
# allres=NULL
# for (j in 1:length(allgenes))
# {
#   res=estimate_heri(genename = allgenes[j],opt=gcta_opt)
#   allres=rbind(allres,res)
#   if (j %% 1000==1)
#   {
#     write.table(allres,file=heritfile,row.names = T,col.names = T, sep="\t", quote=F)
#   }
#   write.table(allres,file=heritfile,row.names = T,col.names = T, sep="\t", quote=F)
# }
# 
# print(Sys.time())
# print(paste0(organ, " done"))

# genenames=c("IL2RB","COX7A2","FILIP1","HSP90AA1","FOXF1","LDAH","ISYNA1","UBAC1")
# 
# allres3=NULL
# for (j in 1:length(genenames))
# {
#   res=estimate_heri(genename = genenames[j],opt="mgrm")
#   allres3=rbind(allres3,res)
# }
# 
# allres1=NULL
# for (j in 1:length(genenames))
# {
#   res=estimate_heri(genename = genenames[j],opt="1grm")
#   allres1=rbind(allres1,res)
# }


library(Rmpi)
njobs=mpi.universe.size() - 1
print(njobs)
mpi.spawn.Rslaves(nslaves=njobs,needlog = F)

mpi.bcast.cmd(library(data.table))
mpi.bcast.cmd(library(GenomicRanges))
mpi.bcast.Robj2slave(plink)
mpi.bcast.Robj2slave(gctafolder)
mpi.bcast.Robj2slave(outfolder)
mpi.bcast.Robj2slave(genotypefolder)
mpi.bcast.Robj2slave(gtexprefix)
mpi.bcast.Robj2slave(phenotype)
mpi.bcast.Robj2slave(phenotypepos)
divide_snp_2list=function(snpdat=snp)
{
  idxs=seq(1,nrow(snpdat),by=as.integer(nrow(snpdat)/5))
  idxs[length(idxs)]=nrow(snpdat)+1
  snplist=list()
  for (i in 1:(length(idxs)-1))
  {
    snplist[[i]]=snpdat[idxs[i]:(idxs[i+1]-1),]
  }
  return(snplist)
}

rbind_snplist=function(ii)
{
  snp=NULL
  for (i in 1:ii)
  {
    #snp=rbind(snp,snplist[[i]])
    tmpsnp=get(paste0("snp",i)) 
    snp=rbind(snp,tmpsnp)
  }
  return(snp)
}

mpi.bcast.Robj2slave(covariate)
# mpi.bcast.Robj2slave(gr_allsnp)
# mpi.bcast.Robj2slave(gr_allpos)
mpi.bcast.Robj2slave(gr_pos)
mpi.bcast.Robj2slave(readhsq)
mpi.bcast.Robj2slave(estimate_heri)

estimate_allgenes=function(opt="1grm")
{
  res=NULL
  n=njobs
  idx=match(allgenes,rownames(phenotypepos))
  for (chr in 1:22)
  {
    print(paste0("chr",chr,"----"))
    allgenes1=allgenes[phenotypepos$chr[idx]==chr]
    idx1=which(allsnppos$chr==chr)
    snppos=allsnppos[idx1,]
    snp=allsnp[idx1,]
    snplist=divide_snp_2list(snpdat=snp)
    # mpi.bcast.Robj2slave(snplist) #this doesn't work
    #send each block to slaves
    for (ii in 1:length(snplist))
    {
      #cat(ii,'..')
      tmpsnp=snplist[[ii]]
      mpi.bcast.Robj2slave(tmpsnp)
      mpi.bcast.Robj2slave(ii)
      mpi.remote.exec(assign(paste0("snp",ii),tmpsnp)) #all variables are on slaves!!!
    }
    mpi.bcast.Robj2slave(rbind_snplist)
    mpi.bcast.cmd(assign("snp",rbind_snplist(ii)))
    mpi.remote.exec(dim(snp))
    mpi.bcast.Robj2slave(snppos)
    gr_snp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP ##
    mpi.bcast.Robj2slave(gr_snp)
    rows=1:length(allgenes1)
    nchunks=ceiling(length(rows)/n)
    print(Sys.time())
    print(paste0("number of total:",nchunks))
    for (i in 1:nchunks)
    {
      cat(i,"..")
      if (i<nchunks)
      {
        seq=allgenes1[rows[((i-1)*n+1):(i*n)]]
      }else
      {
        seq=allgenes1[rows[((i-1)*n+1):length(rows)]]
      }
      tmp=mpi.parSapply(X=seq,FUN=estimate_heri,opt=opt,job.num=njobs)
      if (class(tmp)[1]=="list")
      {
        tmp1=lapply(tmp,is.null)
        tmp=tmp[!unlist(tmp1)]
      }
      #if (length(unlist(tmp)) %% 3 !=0) print (i)
      res1=matrix(unlist(tmp),ncol=3,byrow = T)
      if (class(tmp)[1]=="matrix")
      {
        rownames(res1)=seq
      }else if (class(tmp)=="list")
      {
        rownames(res1)=names(tmp)
      }else
      {
        print(paste0("somthing is wrong: ",i))
      }
      res=rbind(res,res1)
    }
  }
  return(res)
}
allres=estimate_allgenes()
write.table(allres,file=heritfile,row.names = T,col.names = T, sep="\t", quote=F)
print(paste0(organ," is done!"))
print(Sys.time())
mpi.close.Rslaves()
mpi.quit()
quit()
