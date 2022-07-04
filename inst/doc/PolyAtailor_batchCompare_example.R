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
#  library(polyAtailor)
#  ## For full-length sequences the longRead parameter is set to T
#  D1 = tailMap("./dataset1.bam",mcans=5,minTailLen=8,findUmi = F,longRead=T)
#  head(D1)
#  # read_num            chr strand     coord       PAL
#  # 1 SRR8568871.100170   NC_000007.14      -  44801385  66
#  # 2 SRR8568871.100170   NC_000005.10      +  82009391  66
#  # 3 SRR8568871.100177 NW_021160006.1      -    159174  43
#  # 4 SRR8568871.100177   NC_000010.11      + 118932093  43
#  # 5 SRR8568871.100246   NC_000007.14      -  75416541  83
#  # 6 SRR8568871.100246   NC_000007.14      +  73303226  83
#  # tail
#  # 1 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 2 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 3 TTTTTTTTTTTTTTTTTTTTTTTTTGGTTTTTTTTTTTTTTTT
#  # 4 TTTTTTTTTTTTTTTTTTTTTTTTTGGTTTTTTTTTTTTTTTT
#  # 5 TTTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 6 TTTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # tailType      read_type nA        rt
#  # 1 structural two-tail-mixed 66 1.0000000
#  # 2 structural two-tail-mixed 66 1.0000000
#  # 3 structural two-tail-mixed 41 0.9534884
#  # 4 structural two-tail-mixed 41 0.9534884
#  # 5 structural two-tail-mixed 83 1.0000000
#  # 6 structural two-tail-mixed 83 1.0000000
#  ## For NGS sequences the longread parameter is set to F
#  D2 = tailMap("./dataset2.bam",mcans=5,minTailLen=8,findUmi = F,longRead=F)
#  head(D2)
#  # read_num          chr strand     coord PAL           tail
#  # 60  SRR5314571.10000235 NC_000002.12      +  32916482  14 TTTTTTTTTTCCTT
#  # 84  SRR5314571.10000325 NC_000001.11      + 179201983  11    TTTTTTTTCTT
#  # 292 SRR5314571.10001293 NC_000001.11      + 179201983  13  TTTTTTTTTTCTT
#  # 352 SRR5314571.10001518 NC_000011.10      +  66003917   8       TTTTTTTT
#  # 394 SRR5314571.10001748 NC_000015.10      +  90953981  10     TTTTTTTTAT
#  # 402  SRR5314571.1000177 NC_000002.12      +  28801419  14 TTTTTTTTTTTTTT
#  # tailType  read_type nA        rt
#  # NA  60   one-tail 12 0.8571429
#  # NA  84   one-tail 10 0.9090909
#  # NA  292  one-tail 12 0.9230769
#  # NA  352  one-tail  8 1.0000000
#  # NA  394  one-tail  9 0.9000000
#  # NA  402  one-tail 14 1.0000000
#  # 4.Poly(A) sites determination with PolyAtailor
#  ## example
#  library("BSgenome.Hsapiens.UCSC.hg38")
#  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#  bsgenome = BSgenome.Hsapiens.UCSC.hg38
#  gffFile = TxDb.Hsapiens.UCSC.hg38.knownGene
#  bamfile = "./dataset1.bam"
#  chrinfo = "./chrinfo.txt"
#  resultpath = "./result"
#  pacD1 = findAndAnnoPAs(bamfile,chrinfo,resultpath,bsgenome,gffFile,sample="D1rep1",mergePAs=T)

## ----eval=FALSE---------------------------------------------------------------
#  datalist = list(Batch1 = D1, Batch2 = D2)
#  p = batchCompare(datalist,format="upset",dimension="gene",mycolors=c("#be8ec4","#7ed321"))
#  p
#  p = batchCompare(datalist,format="veen",dimension="gene",mycolors=c("#be8ec4","#7ed321"))
#  p

## ----eval=FALSE---------------------------------------------------------------
#  p = batchCompare(datalist,dimension="read",mycolors=c("#9bbfdc","#ab99c1"))
#  p

## ----eval=FALSE---------------------------------------------------------------
#  p = batchCompare(datalist,dimension="tail",mycolors=c("#f18687","#9bbfdc"),rep=T)
#  p

## ----eval=FALSE---------------------------------------------------------------
#  PACdslist = list(Batch1 = pacD1, Batch2 = pacD2)
#  p = batchCompare(PACdslist,dimension="PACds",findOvpPACds=T,annotateByKnownPAC=T,d=100)
#  head(P[['OvpPACds']])
#  # pacD1 pacD2 Total1 Total2 Ovp1 Ovp2 Ovp1Pct Ovp2Pct
#  # 1    1    1  31289 186626 2404 2775      8%      1%

