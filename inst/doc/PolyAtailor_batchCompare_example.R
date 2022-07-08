## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = F,
  fig.width = 6, 
  fig.height = 6,
  eval=FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  # 1.Download sequence files in fastq format from NCBI
#  # 2.Align sequences and convert file formats
#  # 3.Quantitativing poly(A) tails with PolyAtailor
#  rm(list=ls())
#  library(PolyAtailor)
#  #require(devtools)
#  #install_github("BMILAB/movAPA")
#  library(movAPA)
#  library(Rsamtools)
#  #BiocManager::install("Rsamtools")
#  #Persistently prefer one function over anothe, in case
#  conflict_prefer_all("PolyAtailor", "dplyr")
#  conflict_prefer_all("PolyAtailor", "IRanges")
#  conflict_prefer("summarise", "dplyr")
#  
#  ## For full-length sequences the longRead parameter is set to T
#  #here dataset1 is subset of Dataset1 data
#  data1.path <- system.file("extdata", "./GV_algin/d1_subset.bam", package = "PolyAtailor", mustWork = TRUE)
#  D1 = tailMap(bamfile=data1.path,mcans=5,minTailLen=8,findUmi = F,longRead=T)
#  head(D1[,1:4])
#  #            read_num          chr strand     coord
#  # 1 SRR8568871.105958 NC_000001.11      - 149842185
#  # 2  SRR8568871.13378 NC_000001.11      -   8870646
#  # 3 SRR8568871.134537 NC_000001.11      +  58999617
#  # 4 SRR8568871.135358 NC_000001.11      + 236483105
#  # 5 SRR8568871.138552 NC_000001.11      - 149860081
#  # 6 SRR8568871.141140 NC_000001.11      +  23961696
#  
#  ## For NGS sequences the longread parameter is set to F
#  #here dataset2 is subset of Dataset2 data
#  data2.path <- system.file("extdata", "./GV_algin/d2_subset.bam", package = "PolyAtailor", mustWork = TRUE)
#  D2 = tailMap(bamfile=data2.path,mcans=5,minTailLen=8,findUmi = F,longRead=F)
#  head(D2[,1:4])
#  #              read_num          chr strand     coord
#  # 1  SRR5314571.1089513 NC_000001.11      -   8861006
#  # 2 SRR5314571.11256375 NC_000001.11      -  96678784
#  # 3 SRR5314571.12192933 NC_000001.11      -  96678784
#  # 4 SRR5314571.12237630 NC_000001.11      -  96678791
#  # 5 SRR5314571.12696790 NC_000001.11      - 121168671
#  # 6 SRR5314571.13095118 NC_000001.11      -  96678784
#  
#  
#  
#  # 4.Poly(A) sites determination with PolyAtailor
#  ## example
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  #BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  #BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#  bsgenome = BSgenome.Mmusculus.UCSC.mm10
#  gffFile = TxDb.Mmusculus.UCSC.mm10.knownGene
#  chrinfo = system.file("extdata", "./GV_algin/chrinfo.txt", package = "PolyAtailor", mustWork = TRUE)
#  resultpath = "./result"
#  if(!dir.exists(resultpath)){
#    dir.create(resultpath)
#  }
#  bamfile =system.file("extdata", "./GV_algin/GV1subseq.sorted.bam", package = "PolyAtailor", mustWork = TRUE)
#  conflict_prefer("reduce", "IRanges")
#  pacD1 = findAndAnnoPAs(chrinfo=chrinfo,bamfile=bamfile,
#                         resultpath=resultpath,bsgenome=bsgenome,
#                         gffFile=gffFile,sample="D1rep1",mergePAs=T,d=24)
#  #For example only, please use another set of data
#  pacD2 = ffindAndAnnoPAs(chrinfo=chrinfo,bamfile=bamfile,
#                         resultpath=resultpath,bsgenome=bsgenome,
#                         gffFile=gffFile,sample="D1rep2",mergePAs=T,d=24)

## ----eval=FALSE---------------------------------------------------------------
#  datalist = list(Batch1 = pacD1@anno, Batch2 = pacD2@anno)
#  dev.off()
#  p = batchCompare(datalist=datalist,format="upset",dimension="gene",mycolors=c("#be8ec4","#7ed321"))
#  p
#  library(VennDiagram)
#  p = batchCompare(datalist=datalist,format="veen",dimension="gene",mycolors=c("#be8ec4","#7ed321"))

## ----eval=FALSE---------------------------------------------------------------
#  #For example only
#  datalist = list(Batch1 = D1, Batch2 = D1)
#  p = batchCompare(datalist,format="bar",dimension="read",mycolors=c("#9bbfdc","#ab99c1"))
#  p

## ----eval=FALSE---------------------------------------------------------------
#  #For example only
#  datalist = list(Batch1 = D1, Batch2 = D1)
#  #BiocManager::install("PupillometryR")
#  library(PupillometryR)
#  library(ggthemes)
#  p = batchCompare(datalist,format="bar",dimension="tail",mycolors=c("#f18687","#9bbfdc"),rep=F)
#  #or
#  #batchCompare(datalist,format="bar",dimension="tail",mycolors=c("#f18687","#9bbfdc"),rep=T)
#  p

## ----eval=FALSE---------------------------------------------------------------
#  PACdslist = list(Batch1 = pacD1, Batch2 = pacD2)
#  p = batchCompare(PACdslist,dimension="PACds",findOvpPACds=T,annotateByKnownPAC=T,d=100,mycolors=c("#be8ec4","#7ed321"),rep=F)
#  head(p[['OvpPACds']])  #or head(p$OvpPACds)
#  # pacD1 pacD2 Total1 Total2 Ovp1 Ovp2 Ovp1Pct Ovp2Pct
#  # 1    1    1  31289 186626 2404 2775      8%      1%

