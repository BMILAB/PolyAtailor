## -----------------------------------------------------------------------------
# library(PolyAtailor)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = F,
  fig.width = 6, 
  fig.height = 6,
  eval=FALSE
)

## ----setup--------------------------------------------------------------------
#  conflict_prefer("Position", "ggplot2")
#  conflict_prefer("filter", "dplyr")
#  conflict_prefer("mutate", "dplyr")
#  conflict_prefer("select", "dplyr")
#  conflict_prefer("rename", "dplyr")
#  conflict_prefer("summarise", "dplyr")
#  conflict_prefer("rbind", "base")
#  conflict_prefer("cbind", "base")
#  conflict_prefer("strsplit", "base")
#  conflict_prefer("count", "dplyr")
#  conflict_prefer("list", "base")
#  conflict_prefer("reduce", "IRanges")
#  conflict_prefer("geom_bar", "ggplot2")
#  conflict_prefer("first", "dplyr")
#  conflict_prefer("combine", "dplyr")
#  conflict_prefer("compose", "purrr")
#  conflict_prefer("last", "dplyr")
#  conflict_prefer("simplify", "purrr")
#  conflict_prefer("%>%", "dplyr")
#  library(polyAtailor)
#  library(movAPA)
#  #require(devtools)
#  #install_github("BMILAB/movAPA")
#  fastqfile <- system.file("extdata", "./GV_fastq/PAIso_GV1.fastq", package = "PolyAtailor", mustWork = TRUE)
#  resultpath = "./output"
#  if(!dir.exists(resultpath)){
#    dir.create(resultpath)
#  }
#  GV1tailDF<-tailScan(fastqfile,mcans=5,findUmi = F,resultpath = resultpath,samplename = "GV1",tailAnchorLen=8,minTailLen=8,realTailLen=20,maxNtail=2,mapping=F)
#  head(GV1tailDF)
#  #>read_num strand PAL tail tailType read_type nA rt
#  #SRR8798075.11991 - 9  TTTTTTTTA     structural two-tail-mixed 8  0.89 GV1
#  #SRR8798075.11991 + 24 TTTTT...ATTGT structural two-tail-mixed 19 0.79 GV1
#  #SRR8798075.12051 - 13 TTTTTTTTTTTTT structural two-tail-mixed 13 1.00 GV1
#  #SRR8798075.12051 + 57 TTTTT...TTTTG structural two-tail-mixed 56 0.98 GV1
#  #SRR8798075.12545 - 28 TTTTT...TGCTT structural two-tail-mixed 21 0.75 GV1
#  #SRR8798075.12545 + 15 TTTTT...GATCT structural two-tail-mixed 12 0.80 GV1
#  median(GV1tailDF$PAL)
#  # [1] 57

## ---- include = FALSE---------------------------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>",
#    warning = F,
#    fig.width = 6,
#    fig.height = 6,
#    eval=FALSE
#  )

## -----------------------------------------------------------------------------
#  #step1
#  #Tails and partial sequences were extracted from long reads and FASTA files were generated for alignment.
#  fastqfile <- system.file("extdata", "./GV_fastq/PAIso_GV1.fastq", package = "PolyAtailor", mustWork = TRUE)
#  library(seqRFLP)
#  faBuilderRE <- faBuilder(fastqfile,mcans=5,findUmi = F,resultpath = resultpath ,samplename = "GV1",tailAnchorLen=8,mapping=F,findTailType="both")
#  # > head(faBuilderRE[,c(1,2,3,5)])
#  #                         read_num strand PAL tail       nA
#  # 1: SRR8798075.1_unstructural_GV1      + 302 TTT...CCT 111
#  # 2: SRR8798075.2_unstructural_GV1      + 225 TTT...ACT 131
#  # 3: SRR8798075.3_unstructural_GV1      + 300 TTT...TCG 109
#  # 4: SRR8798075.4_unstructural_GV1      + 314 TTT...ACA 142
#  # 5: SRR8798075.5_unstructural_GV1      + 294 TTT...TCG 112
#  # 6:   SRR8798075.6_structural_GV1      - 260 TTT...CTT 82
#  #step2
#    ##alignment with any aligner.
#  
#  #step3
#  bamfile <- system.file("extdata", "./GV_algin/GV1subseq.sorted.bam", package = "PolyAtailor", mustWork = TRUE)
#  GV1tailMapre<-tailMap(bamfile,mcans=5,minTailLen=8,findUmi = F,longRead=T)
#  head(GV1tailMapre)
#  #           read_num                   chr strand coord PAL
#  # 1  SRR8798075.1004  ENSMUST00000009039.6      -   326  69
#  # 2 SRR8798075.10180 ENSMUST00000100052.11      -  3273   8
#  # 3 SRR8798075.10187  ENSMUST00000141115.8      +  9913  10
#  # 4 SRR8798075.10190  ENSMUST00000139787.8      -   296  40
#  # 5 SRR8798075.10378 ENSMUST00000055131.13      -  3224  63
#  # 6 SRR8798075.10430  ENSMUST00000209034.2      -  1661  57
#  #                                                                    tail
#  # 1 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 2                                                              TTTTTTTT
#  # 3                                                            TTTTTTTTTT
#  # 4                              TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 5       TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATT
#  # 6             TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  #     tailType read_type nA       rt sample
#  # 1 structural  one-tail 69 1.000000    GV1
#  # 2 structural  one-tail  8 1.000000    GV1
#  # 3 structural  one-tail 10 1.000000    GV1
#  # 4 structural  one-tail 40 1.000000    GV1
#  # 5 structural  one-tail 62 0.984127    GV1
#  # 6 structural  one-tail 57 1.000000    GV1
#  
#  #step4
#  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  data(GV1tailDF)
#  data(GV1tailMapre)
#  AnnotedTails = geneAnno(tailDF=GV1tailDF,bamdf=GV1tailMapre,GFF=TxDb.Mmusculus.UCSC.mm10.knownGene,longRead=F)
#  head(AnnotedTails)
#  #           read_num   chr strand  gene gene_type PAL
#  # 1 SRR8798075.10543  chr3      - 72007      <NA>  11
#  # 2 SRR8798075.10543  chr3      - 72007      <NA>  17
#  # 3 SRR8798075.11063 chr10      + 69412      <NA>  89
#  # 4 SRR8798075.13407 chr14      - 12159      <NA>  89
#  # 5 SRR8798075.14401 chr14      - 12159      <NA>  43
#  # 6 SRR8798075.14401 chr14      - 12159      <NA> 119
#  #   tail
#  # 1 TTTTTTTTTTT
#  # 2 TTTTTTTTTTTTTTTTT
#  # 3 TTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 4 TTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 5 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGCT
#  # 6 TTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGTGCTTTGCTTTACTTTC
#  #       tailType      read_type  nA        rt sample
#  # 1   structural two-tail-mixed  11 1.0000000    GV1
#  # 2   structural two-tail-mixed  17 1.0000000    GV1
#  # 3   structural       one-tail  89 1.0000000    GV1
#  # 4   structural       one-tail  89 1.0000000    GV1
#  # 5   structural  two-tail-same  41 0.9534884    GV1
#  # 6 unstructural  two-tail-same 111 0.9327731    GV1

## ---- include = FALSE---------------------------------------------------------
#  knitr::opts_chunk$set(
#    collapse = TRUE,
#    comment = "#>",
#    warning = F,
#    fig.width = 6,
#    fig.height = 6,
#    eval=FALSE
#  )

## -----------------------------------------------------------------------------
#  #step1
#  bamfile <- system.file("extdata", "./GV_algin/PAIso-GV1.sorted.bam", package = "PolyAtailor", mustWork = TRUE)
#  GV1tailMapre<-tailMap(bamfile,mcans=5,minTailLen=8,findUmi = F,longRead=T)
#  head(GV1tailMapre)
#  #           read_num   chr strand     coord PAL
#  # 1 SRR8798075.15955  chrM      +      5185  11
#  # 2  SRR8798075.2473 chr11      +  59211415  55
#  # 3 SRR8798075.26425  chr6      - 114079037  11
#  # 4 SRR8798075.32952  chr2      - 153826328  25
#  # 5 SRR8798075.33179  chr2      + 161287180  14
#  # 6 SRR8798075.33206  chr9      -  56438002  33
#  #                                                      tail
#  # 1                                             TTTTTTTTTTT
#  # 2 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 3                                             TTTTTTTTTTT
#  # 4                               TTTTTTTTTTTTTTTTTTTTTTTTT
#  # 5                                          TTTTTTTTTTTTAT
#  # 6                       TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAT
#  #     tailType read_type nA       rt sample
#  # 1 structural  one-tail 11 1.000000    GV1
#  # 2 structural  one-tail 55 1.000000    GV1
#  # 3 structural  one-tail 11 1.000000    GV1
#  # 4 structural  one-tail 25 1.000000    GV1
#  # 5 structural  one-tail 13 0.984127    GV1
#  # 6 structural  one-tail 32 1.000000    GV1
#  
#  #step2
#  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  data(GV1tailDF)
#  data()
#  AnnotedTails = geneAnno(tailDF=GV1tailDF,bamdf=GV1tailMapre,GFF=TxDb.Mmusculus.UCSC.mm10.knownGene,longRead=F)
#  head(AnnotedTails)
#  #           read_num   chr strand   gene gene_type PAL
#  # 1  SRR8798075.2473 chr11      +  67862      <NA>  55
#  # 2  SRR8798075.2473 chr11      +  67862      <NA> 116
#  # 3 SRR8798075.32952  chr2      -  76407      <NA>  31
#  # 4 SRR8798075.38605  chr1      +  12839      <NA>  71
#  # 5 SRR8798075.58329  chr6      + 208665      <NA>  45
#  # 6 SRR8798075.58329  chr6      + 208665      <NA>  98
#  #   tail
#  # 1 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 2 TTTTTTT...TTCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 3 TTTTTTTTTTTTTTTTTTTTTTTTTTTTGGT
#  # 4 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATGTATTA
#  # 5 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  # 6 TTTTTTT...TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
#  #       tailType      read_type  nA        rt sample
#  # 1   structural two-tail-mixed  55 1.0000000    GV1
#  # 2 unstructural two-tail-mixed 115 0.9913793    GV1
#  # 3   structural       one-tail  29 0.9354839    GV1
#  # 4   structural       one-tail  67 0.9436620    GV1
#  # 5   structural  two-tail-same  45 1.0000000    GV1
#  # 6 unstructural  two-tail-same  98 1.0000000    GV1

## -----------------------------------------------------------------------------
#  # Deciphering gff files
#  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  # Deciphering genome files
#  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#  library("BSgenome.Mmusculus.UCSC.mm10")
#  bsgenome = BSgenome.Mmusculus.UCSC.mm10
#  # Prepare the input file
#   bamfilepath = system.file("extdata", "./GV_algin/PAIso-GV1.sorted.bam", package = "PolyAtailor", mustWork = TRUE)
#  chrinfopath = system.file("extdata", "./GV_algin/chrinfo.txt", package = "PolyAtailor", mustWork = TRUE)
#  resultpath = "./"
#  # Annotated PA site
#  PAs <- findAndAnnoPAs(bamfile=bamfilepath,chrinfo=chrinfopath,resultpath=resultpath,bsgenome=bsgenome,gffFile = TxDb.Mmusculus.UCSC.mm10.knownGene,sample="GV1",mergePAs=T,d=24)
#  

## -----------------------------------------------------------------------------
#  data(GV1tailMapre)
#  library(ggthemes)
#  library(eoffice)
#  p1 <- plotPALDistribution(GV1tailMapre,"./data/figures/","global",medianPAL=T)
#  p1
#  p2 <- plotPALDistribution(AnnotedTails,"./data/figures/","gene",medianPAL=T)
#  p2

## -----------------------------------------------------------------------------
#  data(PAs)
#  p <- plotPADistribution(PAs,"./data/figures/","#9BBFDC")
#  p

## -----------------------------------------------------------------------------
#  data(PAs)
#  p <- plotGenePAnumbers(PAs,"./data/figures/","#DF7C7D")
#  p

## -----------------------------------------------------------------------------
#  library("BSgenome.Mmusculus.UCSC.mm10")
#  bsgenome = BSgenome.Mmusculus.UCSC.mm10
#  data(PAs)
#  p <- plotPASignals(PAs,"./data/figures/",bsgenome = bsgenome)
#  p

## -----------------------------------------------------------------------------
#  data(AnnotedTails)
#  library(stringi)
#  re <- nonAanalysis(AnnotedTails)
#  re$p1

## -----------------------------------------------------------------------------
#  re$p2
#  re$p3

## -----------------------------------------------------------------------------
#  data(taildf)
#  my_cutstom <- data.frame(names=c("A","C","T","G"),color=c("#3171A5","#4EAA4C","#C9C4C2","#D73D3D"))
#  library(ggmsa)
#  p <- tailViso(taildf,tailLen=100,Ntail=20,custom=my_cutstom,strand="-",faPath="D:/",showLogo=T,showReadNum= F)

## -----------------------------------------------------------------------------
#  #BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  gff <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.knownGene)
#  data(AnnotedTails)
#  files = system.file("extdata", "./output/PAs/PAs.txt", package = "PolyAtailor", mustWork = TRUE)
#  PAs <- read.table(files,header=TRUE,sep=" ")
#  diffPAL2PAgenes <- PALdsa(PAs,AnnotedTails,gff,mode="PD",SAoDMethod="ME",withViolinPlot=T,withUpsetPlot=F)
#  

