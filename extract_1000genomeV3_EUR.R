thousandgtable=read.xlsx("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/20130606_sample_info.xlsx",sheetIndex = 1)

#check populations
thousandpoptable=read.table("/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/populations.txt",header = T,stringsAsFactors = F,sep="\t")
thousandgtable$superpopulation=NA
for (i in 1:nrow(thousandgtable))
{
  idx=which(thousandpoptable$Population.Code==thousandgtable$Population[i])
  thousandgtable$superpopulation[i]=thousandpoptable$Super.Population.Code[idx]
}

thousandgtable$Population.Description=as.character(thousandgtable$Population.Description)
thousandgtable$PopulationName=NA
thousandgtable$SamplingLocation=NA
for (i in 1:nrow(thousandgtable))
{
  if (thousandgtable$Population[i]=="CEU")
  {
    thousandgtable$PopulationName[i]="CEU"
    thousandgtable$SamplingLocation[i]=thousandgtable$Population.Description[i]
  }
  if (thousandgtable$Population[i]=="TSI")
  {
    thousandgtable$PopulationName[i]="Tuscan"
    thousandgtable$SamplingLocation[i]="Italy"
  }
  if (thousandgtable$Population[i]=="FIN")
  {
    thousandgtable$PopulationName[i]="Finnish"
    thousandgtable$SamplingLocation[i]="Finland"
  }
  if (thousandgtable$Population[i]=="GBR")
  {
    thousandgtable$PopulationName[i]="British"
    thousandgtable$SamplingLocation[i]="England and Scotland"
  }
  if (thousandgtable$Population[i]=="IBS")
  {
    thousandgtable$PopulationName[i]="Iberian"
    thousandgtable$SamplingLocation[i]="Spain"
  }
  if (thousandgtable$Population[i]=="MXL")
  {
    thousandgtable$PopulationName[i]="MXL"
    thousandgtable$SamplingLocation[i]=thousandgtable$Population.Description[i]
  }
  if (thousandgtable$Population[i]=="PUR")
  {
    thousandgtable$PopulationName[i]="Puerto Rican"
    thousandgtable$SamplingLocation[i]="Puerto Rico"
  }
  if (thousandgtable$Population[i]=="CLM")
  {
    thousandgtable$PopulationName[i]="Colombian"
    thousandgtable$SamplingLocation[i]="Colombia"
  }
  if (thousandgtable$Population[i]=="PEL")
  {
    thousandgtable$PopulationName[i]="Peruvian"
    thousandgtable$SamplingLocation[i]="Peru"
  }
  
  if (thousandgtable$Population[i]=="YRI")
  {
    thousandgtable$PopulationName[i]="Yoruba"
    thousandgtable$SamplingLocation[i]="Nigeria"
  }
  if (thousandgtable$Population[i]=="LWK")
  {
    thousandgtable$PopulationName[i]="Luhya"
    thousandgtable$SamplingLocation[i]="Kenya"
  }
  if (thousandgtable$Population[i]=="GWD")
  {
    thousandgtable$PopulationName[i]="Gambian"
    thousandgtable$SamplingLocation[i]="Gambia"
  }
  if (thousandgtable$Population[i]=="MSL")
  {
    thousandgtable$PopulationName[i]="Mende"
    thousandgtable$SamplingLocation[i]="Sierra Leone"
  }
  if (thousandgtable$Population[i]=="ESN")
  {
    thousandgtable$PopulationName[i]="Esan"
    thousandgtable$SamplingLocation[i]="Nigeria"
  }
  if (thousandgtable$Population[i]=="ASW")
  {
    thousandgtable$PopulationName[i]="ASW"
    thousandgtable$SamplingLocation[i]="Americans of African Ancestry in SW USA"
  }
  if (thousandgtable$Population[i]=="ACB")
  {
    thousandgtable$PopulationName[i]="African Caribbean"
    thousandgtable$SamplingLocation[i]="Barbados"
  }
}
idx=which(thousandgtable$superpopulation=="EUR")
all(thousandgtable$Sample[idx]==thousandgtable$Family.ID[idx])

thousanddir="/fh/fast/stanford_j/Xiaoyu_Oct2020/Tools/HGDP_1000G_Merge/DataBases/1000G/"
tmp=data.frame(ind=thousandgtable$Sample[idx],fam=thousandgtable$Sample[idx],gender=thousandgtable$Gender[idx],relationship=thousandgtable$Relationship[idx])
write.table(tmp[,1:2],file=paste0(thousanddir,"EURsamples.txt"),row.names = F,col.names = F,sep="\t",quote=F)

fusionfam=read.table("/fh/fast/dai_j/CancerGenomics/Tools/fusion_twas-master/LDREF/1000G.EUR.1.fam")
all(fusionfam$V1 %in% tmp$ind)
idx=match(fusionfam$V1,tmp$ind)

thousandfam=read.table(paste0(thousanddir,"ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.fam"))
all(thousandfam$V1 %in% thousandgtable$Sample )
idx=match(thousandfam$V1,thousandgtable$Sample)
table(thousandgtable$superpopulation[idx])
# AFR AMR EAS EUR SAS 
# 661 347 504 503 489 
allavailablesamples=intersect(thousandfam$V1,thousandgtable$Sample)
idx=match(fusionfam$V1,thousandgtable$Sample)
