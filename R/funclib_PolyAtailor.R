# description -------------------------------------------------------------
#### This section of code contains most of polyAtailor's functional functions,
### This script is used to quantify the Poly (A) tail without alignment.
# dependence --------------------------------------------------------------
# usethis::use_package("Biostrings")
# usethis::use_package("Rsamtools")
# usethis::use_package("GenomicRanges")
# usethis::use_package("Rsamtools")
# usethis::use_package("tidyverse", type = "depends")
# usethis::use_package("magrittr")
# usethis::use_package("stringi")
# usethis::use_package("ggpubr")
# usethis::use_package("ggplot2")
# usethis::use_package("ggthemes")
# usethis::use_package("conflicted")
# usethis::use_package("dplyr")
# usethis::use_package("sqldf")
# usethis::use_package("ShortRead")
# usethis::use_package("stringr")
# usethis::use_package("esquisse")
# usethis::use_package("tidyr")
# usethis::use_package("seqRFLP")
# usethis::use_package("data.table")
# usethis::use_package("ggsci")
# usethis::use_package("gridExtra")
# usethis::use_package("plotrix")
# usethis::use_package("movAPA")
# usethis::use_package("eoffice")
# usethis::use_package("ggpubr")
# usethis::use_package("UpSetR")
# usethis::use_package("VennDiagram")
# usethis::use_package("DescTools")
# usethis::use_package("ggmsa")
# usethis::use_package("GeneAnswers")
# usethis::use_package("remotes")
# usethis::use_package("geneRal")
# usethis::use_package("GGally")
# library(Biostrings)
# library(conflicted)
# library(data.table)
# library(DescTools) #no
# library(dplyr)
# library(eoffice) #No
# library(esquisse) #No
# library(GeneAnswers) #NO
# #library(geneRal)#NO   no avialble
# library(GenomicRanges)
# library(GGally)
# library(ggmsa)#NO
# library(ggplot2)
# library(ggpubr)
# library(ggsci)
# library(ggthemes)#NO
# library(gridExtra)
# library(magrittr)
# library(movAPA)
# library(plotrix)#NO
# library(purrr)
# library(remotes)
# library(Rsamtools)
# library(seqRFLP)#No
# library(ShortRead)#No
# library(sqldf)#no
# library(stringi)
# library(stringr)
# library(tidyr)
# library(UpSetR)#no
# library(cli)
# library(VennDiagram)
# library(TxDb.Athaliana.BioMart.plantsmart28)#no
# #library(BSgenome.Athaliana.ENSEMBL.TAIR10)#no
# library(RColorBrewer)
# library(pracma)
# library(tidyverse)
# library(knitr)
# knitr::opts_chunk$set(cache = TRUE, warning = FALSE,
#                       message = FALSE, cache.lazy = FALSE)

conflict_prefer("Position", "ggplot2")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("rbind", "base")
conflict_prefer("cbind", "base")
conflict_prefer("strsplit", "base")
conflict_prefer("count", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("reduce", "IRanges")
conflict_prefer("geom_bar", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("compose", "purrr")
conflict_prefer("last", "dplyr")
conflict_prefer("simplify", "purrr")
conflict_prefer("%>%", "dplyr")

# build-in function -------------------------------------------------------
dss2df <- function(dss) {
  data.frame(width = width(dss), seq = as.character(dss), names = names(dss))
}


scan_findumi <- function(x, mcans, findUmi, lumi,
                         adapterSeq, anchorSeq, mapping, tailAnchorLen, findTailType) {
  # Initialisation
  read_nums <- c()
  strand <- c()
  umis <- c()
  tails <- c()
  tail_type <- c()
  # get the seq information first
  seq <- as.character(x[2])
  read_num <- as.character(x[3])
  pattern1 <- str_c("T{", tailAnchorLen, ",}")
  pattern2 <- str_c("A{", tailAnchorLen, ",}")
  anchor1 <- gregexpr(pattern1, seq)[[1]]
  anchor2 <- gregexpr(pattern2, seq)[[1]]
  if (findTailType == "both") {
    # return a string if no tail find
    if ((anchor1[1] == -1) & (anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #- strand
    else if ((anchor2[1] == -1) & (anchor1[1] != -1)) {
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    #+ strand
    else if ((anchor2[1] != -1) & (anchor1[1] == -1)) {
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    # all strand
    else {
      # In this case the tailSlider needs to be called twice
      tailsinfo1 <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo1$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo1$tails)
      umis <- append(umis, tailsinfo1$umis)
      tail_type <- append(tail_type, tailsinfo1$tail_type)
      # reverse
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo2 <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo2$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo2$tails)
      umis <- append(umis, tailsinfo2$umis)
      tail_type <- append(tail_type, tailsinfo2$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "A") {
    ### just find poly(A)
    # return a string if no tail find
    if ((anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #+ strand
    else if (anchor2[1] != -1) {
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "T") {
    ### just find poly(T)
    # return a string if no tail find
    if ((anchor1[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #- strand
    else if (anchor1[1] != -1) {
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$umi)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
}

scan_noumi <- function(x, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping, tailAnchorLen, findTailType) {
  # Initialisation
  read_nums <- c()
  strand <- c()
  tails <- c()
  tail_type <- c()
  # get the seq information first
  seq <- as.character(x[2])
  read_num <- as.character(x[3])
  pattern1 <- str_c("T{", tailAnchorLen, ",}")
  pattern2 <- str_c("A{", tailAnchorLen, ",}")
  anchor1 <- gregexpr(pattern1, seq)[[1]]
  anchor2 <- gregexpr(pattern2, seq)[[1]]
  if (findTailType == "both") {
    # return a string if no tail find
    if ((anchor1[1] == -1) & (anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #- strand
    else if ((anchor2[1] == -1) & (anchor1[1] != -1)) {
      # In this case there is only a negative chain, and the tailSlider is called directly
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    #+ strand
    else if ((anchor2[1] != -1) & (anchor1[1] == -1)) {
      # In this case there is only a forward chain, and the tailSlider is called directly
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    # Two-way chain
    else {
      # In this case the tailSlider needs to be called twice
      tailsinfo1 <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo1$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo1$tails)
      tail_type <- append(tail_type, tailsinfo1$tail_type)
      # reverse
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo2 <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo2$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo2$tails)
      tail_type <- append(tail_type, tailsinfo2$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "A") {
    # return a string if no tail find
    if ((anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #+ strand
    else if ((anchor2[1] != -1)) {
      # In this case there is only a forward chain, and the tailSlider is called directly
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "T") {
    # return a string if no tail find
    if ((anchor1[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #- strand
    else if ((anchor1[1] != -1)) {
      # In this case there is only a negative chain, and the tailSlider is called directly
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
}

map_findumi <- function(x, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping, tailAnchorLen, findTailType) {
  # Initialisation
  read_nums <- c()
  strand <- c()
  umis <- c()
  tails <- c()
  tail_type <- c()
  # Get seq information
  seq <- as.character(x[2])
  read_num <- as.character(x[3])
  pattern1 <- str_c("T{", tailAnchorLen, ",}")
  pattern2 <- str_c("A{", tailAnchorLen, ",}")
  anchor1 <- gregexpr(pattern1, seq)[[1]]
  anchor2 <- gregexpr(pattern2, seq)[[1]]
  if (findTailType == "both") {
    # notail
    if ((anchor1[1] == -1) & (anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #-strand
    else if ((anchor2[1] == -1) & (anchor1[1] != -1)) {
      # In this case there is only a negative chain, and the tailSlider is called directly
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    #+strand
    else if ((anchor2[1] != -1) & (anchor1[1] == -1)) {
      # In this case only the forward chain, the reverse call to tailSlider
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    # two-way chain
    else {
      # In this case the tailSlider needs to be called twice
      tailsinfo1 <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo1$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo1$tails)
      umis <- append(umis, tailsinfo1$umis)
      tail_type <- append(tail_type, tailsinfo1$tail_type)
      # reverse
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo2 <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo2$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo2$tails)
      umis <- append(umis, tailsinfo2$umis)
      tail_type <- append(tail_type, tailsinfo2$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "A") {
    # notail
    if ((anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #+strand
    else if ((anchor2[1] != -1)) {
      # In this case only the forward chain, the reverse call to tailSlider
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "T") {
    # notail
    if ((anchor1[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        umi = "not-find", tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #-strand
    else if ((anchor1[1] != -1)) {
      # In this case there is only a negative chain, and the tailSlider is called directly
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tails)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      umis <- append(umis, tailsinfo$umis)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        umi = umis, tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
}

map_noumi <- function(x, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping, tailAnchorLen, findTailType) {
  # Initialisation
  read_nums <- c()
  strand <- c()
  umis <- c()
  tails <- c()
  tail_type <- c()
  # get seq information
  seq <- as.character(x[2])
  read_num <- as.character(x[3])
  pattern1 <- str_c("T{", tailAnchorLen, ",}")
  pattern2 <- str_c("A{", tailAnchorLen, ",}")
  anchor1 <- gregexpr(pattern1, seq)[[1]]
  anchor2 <- gregexpr(pattern2, seq)[[1]]
  if (findTailType == "both") {
    # notail
    if ((anchor1[1] == -1) & (anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #-strand
    else if ((anchor2[1] == -1) & (anchor1[1] != -1)) {
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    #+strand
    else if ((anchor2[1] != -1) & (anchor1[1] == -1)) {
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
    # double chain
    else {
      tailsinfo1 <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo1$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo1$tails)
      tail_type <- append(tail_type, tailsinfo1$tail_type)
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo2 <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo2$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo2$tails)
      tail_type <- append(tail_type, tailsinfo2$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "A") {
    # notail
    if ((anchor2[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #+strand
    else if ((anchor2[1] != -1)) {
      seqReverse <- as.character(reverseComplement(DNAString(seq)))
      anchor3 <- gregexpr("T{8,}", seqReverse)[[1]]
      tailsinfo <- tailSlider(seq = seqReverse, anchors = anchor3, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "+")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
  if (findTailType == "T") {
    # notail
    if ((anchor1[1] == -1)) {
      info <- data.frame(
        read_num = read_num, strand = NA,
        tail = "not-find",
        tailType = "not-find"
      )
      return(info)
    }
    #-strand
    else if ((anchor1[1] != -1)) {
      tailsinfo <- tailSlider(seq = seq, anchors = anchor1, mcans, findUmi, lumi, adapterSeq, anchorSeq, mapping = mapping)
      tail_nums <- length(tailsinfo$tail_type)
      for (i in c(1:tail_nums)) {
        read_nums <- append(read_nums, read_num)
        strand <- append(strand, "-")
      }
      tails <- append(tails, tailsinfo$tails)
      tail_type <- append(tail_type, tailsinfo$tail_type)
      info <- data.frame(
        read_num = read_nums, strand = strand,
        tail = tails,
        tailType = tail_type
      )
      return(info)
    }
  }
}

AinB <- function(A, B, all = T) {
  if (is.factor(A)) A <- as.character(A)
  if (is.factor(B)) B <- as.character(B)
  if (!(is.vector(A) & is.vector(B))) {
    return(F)
  }
  x <- sum(A %in% B)
  if ((all & x == length(A)) | (!all & x > 0)) {
    return(T)
  } else {
    return(F)
  }
}

# tailSlider -----------------------------------------------------------
tailSlider <- function(seq, anchors, mcnas, findUmi, lumi, adapterSeq, anchorSeq, mapping) {
  # Set default parameters
  if (missing(mcnas)) {
    mcnas <- 5
  }
  # Initialisation
  # anchors=gregexpr('T{8,}',seq)
  # anchors=anchors[[1]]
  tails <- c()
  tail_type <- c()
  umis <- c()
  tailStart <- as.integer(anchors)
  currMax <- attr(anchors, "match.length")
  currMaxIdx <- tailStart + currMax
  drop <- 0
  strs <- unlist(strsplit(seq, split = ""))
  # af method
  if (!(mapping)) {
    # This is the case without the need to find umi
    if (!(findUmi)) {
      tailend <- 0
      for (j in 1:length(currMaxIdx)) {
        currMaxIdxi <- currMaxIdx[j]
        currMax1 <- currMax[j]
        curr <- currMax1
        currLen <- currMax1
        tailStartx <- tailStart[j]
        if (tailend >= tailStartx) {
          next
        }
        if (currMaxIdxi - 1 == length(strs)) {
          tail <- str_sub(seq, tailStartx, currMaxIdxi - 1)
        } else {
          for (i in currMaxIdxi:length(strs)) {
            if (strs[i] != "T") {
              curr <- curr - 1
            } else {
              curr <- curr + 1
              if (curr >= currMax1) {
                currMaxIdxi <- i
                currMax1 <- curr
              }
            }
            drop <- currMax1 - curr
            if (drop >= mcnas) {
              break
            }
          }
          tail <- substr(seq, tailStartx, currMaxIdxi - 1)
          tailend <- currMaxIdxi
        }
        tails <- append(tails, tail)
        if (anchorSeq == "miss") {
          # No structure anchor specified
          readL <- nchar(seq)
          # The 8 here should be changed to the current tail length currLen
          # if(tailStartx==1|(tailStartx+8)==readL){
          if (tailStartx == 1 | (tailStartx + currLen) == readL) {
            tail_type <- append(tail_type, "unstructural")
          } else {
            tail_type <- append(tail_type, "structural")
          }
        } else {
          # The structure anchor is specified
          anchorL <- nchar(anchorSeq)
          anchorString <- substr(seq, tailStartx - anchorL, tailStartx - 1)
          anchorType <- substr(anchorSeq, 1, 1)
          typeNum <- str_count(anchorString, anchorType)
          if (typeNum >= anchorL - 2) {
            tail_type <- append(tail_type, "structural")
          } else {
            tail_type <- append(tail_type, "unstructural")
          }
        }
      }
      tailsinfo <- base::list(tails = tails, tail_type = tail_type)
      return(tailsinfo)
    }
    # This is a case where you need to find umi
    if (findUmi) {
      tailend <- 0
      for (j in 1:length(currMaxIdx)) {
        currMaxIdxi <- currMaxIdx[j]
        currMax1 <- currMax[j]
        currLen <- currMax1
        curr <- currMax1
        tailStartx <- tailStart[j]
        if (tailend >= tailStartx) {
          next
        }
        if (currMaxIdxi - 1 == length(strs)) {
          tail <- str_sub(seq, tailStartx, currMaxIdxi - 1)
        } else {
          for (i in currMaxIdxi:length(strs)) {
            if (strs[i] != "T") {
              curr <- curr - 1
            } else {
              curr <- curr + 1
              if (curr >= currMax1) {
                currMaxIdxi <- i
                currMax1 <- curr
              }
            }
            drop <- currMax1 - curr
            if (drop >= mcnas) {
              break
            }
          }
          tail <- substr(seq, tailStartx, currMaxIdxi - 1)
          tailend <- currMaxIdxi
        }
        tails <- append(tails, tail)
        if (anchorSeq == "miss") {
          # No structure anchor specified
          ## find tail_type
          readL <- nchar(seq)
          # The 8 here should be changed to the current tail length
          # if(tailStartx==1|(tailStartx+8)==readL){
          if (tailStartx == 1 | (tailStartx + currLen) == readL) {
            tail_type <- append(tail_type, "unstructural")
          } else {
            tail_type <- append(tail_type, "structural")
          }
          ## find umi
          umi <- substr(seq, tailStartx - lumi, tailStartx - 1)
          umis <- append(umis, umi)
        } else {
          # The structure anchor is specified
          anchorL <- nchar(anchorSeq)
          anchorString <- substr(seq, tailStartx - anchorL, tailStartx - 1)
          anchorType <- substr(anchorSeq, 1, 1)
          typeNum <- str_count(anchorString, anchorType)
          if (typeNum >= anchorL - 2) {
            tail_type <- append(tail_type, "structural")
          } else {
            tail_type <- append(tail_type, "unstructural")
          }
          ## find umi
          umi <- substr(seq, tailStartx - lumi - anchorL, tailStartx - anchorL - 1)
          umis <- append(umis, umi)
        }
      }
      tailsinfo <- base::list(tails = tails, umis = umis, tail_type = tail_type)
      return(tailsinfo)
    }
  }
  # al method
  if (mapping) {
    # This is the case without the need to find umi
    if (!(findUmi)) {
      tailend <- 0
      for (j in 1:length(currMaxIdx)) {
        currMaxIdxi <- currMaxIdx[j]
        currMax1 <- currMax[j]
        curr <- currMax1
        currLen <- currMax1
        tailStartx <- tailStart[j]
        if (tailend >= tailStartx) {
          next
        }
        if (currMaxIdxi - 1 == length(strs)) {
          tail <- str_sub(seq, tailStartx, currMaxIdxi + 249)
        } else {
          for (i in currMaxIdxi:length(strs)) {
            if (strs[i] != "T") {
              curr <- curr - 1
            } else {
              curr <- curr + 1
              if (curr >= currMax1) {
                currMaxIdxi <- i
                currMax1 <- curr
              }
            }
            drop <- currMax1 - curr
            if (drop >= mcnas) {
              break
            }
          }
          tail <- substr(seq, tailStartx, currMaxIdxi + 249)
          tailend <- currMaxIdxi
        }
        tails <- append(tails, tail)
        if (anchorSeq == "miss") {
          # No structure anchor specified
          readL <- nchar(seq)
          # The 8 here should be changed to the current tail length currLen
          # if(tailStartx==1|(tailStartx+8)==readL){
          if (tailStartx == 1 | (tailStartx + currLen) == readL) {
            tail_type <- append(tail_type, "unstructural")
          } else {
            tail_type <- append(tail_type, "structural")
          }
        } else {
          # The structure anchor is specified
          anchorL <- nchar(anchorSeq)
          anchorString <- substr(seq, tailStartx - anchorL, tailStartx - 1)
          anchorType <- substr(anchorSeq, 1, 1)
          typeNum <- str_count(anchorString, anchorType)
          if (typeNum >= anchorL - 2) {
            tail_type <- append(tail_type, "structural")
          } else {
            tail_type <- append(tail_type, "unstructural")
          }
        }
      }
      tailsinfo <- base::list(tails = tails, tail_type = tail_type)
      return(tailsinfo)
    }
    # This is a case where you need to find umi
    if (findUmi) {
      tailend <- 0
      for (j in 1:length(currMaxIdx)) {
        currMaxIdxi <- currMaxIdx[j]
        currMax1 <- currMax[j]
        currLen <- currMax1
        curr <- currMax1
        tailStartx <- tailStart[j]
        if (tailend >= tailStartx) {
          next
        }
        if (currMaxIdxi - 1 == length(strs)) {
          tail <- str_sub(seq, tailStartx, currMaxIdxi + 249)
        } else {
          for (i in currMaxIdxi:length(strs)) {
            if (strs[i] != "T") {
              curr <- curr - 1
            } else {
              curr <- curr + 1
              if (curr >= currMax1) {
                currMaxIdxi <- i
                currMax1 <- curr
              }
            }
            drop <- currMax1 - curr
            if (drop >= mcnas) {
              break
            }
          }
          tail <- substr(seq, tailStartx, currMaxIdxi + 249)
          tailend <- currMaxIdxi
        }
        tails <- append(tails, tail)
        if (anchorSeq == "miss") {
          # No structure anchor specified
          ## find tail_type
          readL <- nchar(seq)
          # The 8 here should be changed to the current tail length
          # if(tailStartx==1|(tailStartx+8)==readL){
          if (tailStartx == 1 | (tailStartx + currLen) == readL) {
            tail_type <- append(tail_type, "unstructural")
          } else {
            tail_type <- append(tail_type, "structural")
          }
          ## find umi
          umi <- substr(seq, tailStartx - lumi, tailStartx - 1)
          umis <- append(umis, umi)
        } else {
          # The structure anchor is specified
          anchorL <- nchar(anchorSeq)
          anchorString <- substr(seq, tailStartx - anchorL, tailStartx - 1)
          anchorType <- substr(anchorSeq, 1, 1)
          typeNum <- str_count(anchorString, anchorType)
          if (typeNum >= anchorL - 2) {
            tail_type <- append(tail_type, "structural")
          } else {
            tail_type <- append(tail_type, "unstructural")
          }
          ## find umi
          umi <- substr(seq, tailStartx - lumi - anchorL, tailStartx - anchorL - 1)
          umis <- append(umis, umi)
        }
      }
      tailsinfo <- base::list(tails = tails, umis = umis, tail_type = tail_type)
      return(tailsinfo)
    }
  }
}


# tailFinder ------------------------------------------------------------
tailFinder <- function(fastdf, mcans, findUmi, lumi, adapterSeq, anchorSeq, resultpath, samplename, tailAnchorLen, mapping, findTailType) {
  # MCNAS:Maximum continuous non-A segment
  # Note: The read_num column in fastdf is named read_num uniformly below
  # 1.Set parameter default values
  if (missing(mcans)) {
    mcans <- 5
  }
  if (missing(tailAnchorLen)) {
    tailAnchorLen <- 8
  }
  if (missing(anchorSeq)) {
    anchorSeq <- "miss"
  }
  if (!(mapping)) {
    # 3.Loop through each reads
    # dim(fastdf)[1]
    if (findUmi) {
      testre <- apply(fastdf, 1, scan_findumi,
        findUmi = findUmi,
        lumi = lumi, mcans = mcans,
        adapterSeq = adapterSeq, anchorSeq = anchorSeq,
        mapping = mapping,
        tailAnchorLen = tailAnchorLen, findTailType = findTailType
      )
      # 4.Regularization
      testre <- data.table::rbindlist(testre, fill = TRUE)
      tailsinfo <- filter(testre, tail != "not-find")
      notail <- testre %>%
        filter(tail == "not-find") %>%
        select(read_num)
      if (dim(notail)[1] != 0) {
        notails <- fastdf %>%
          filter(read_num %in% notail$read_num)
      } else {
        notails <- c("no notail reads")
      }
      filepath <- str_c(resultpath, samplename, "_notail_reads.txt")
      write.table(notails, filepath, quote = F, row.names = F)
      tailsinfo$sample <- samplename
      tailsinfo <- tailsinfo %>%
        mutate(PAL = nchar(as.character(tail)), nA = str_count(tail, "T"), rt = as.numeric(nA / PAL)) %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      return(tailsinfo)
    } else {
      testre <- apply(fastdf, 1, scan_noumi,
        findUmi = findUmi,
        lumi = lumi, mcans = mcans,
        adapterSeq = adapterSeq, anchorSeq = anchorSeq,
        mapping = mapping,
        tailAnchorLen = tailAnchorLen, findTailType = findTailType
      )
      # 4.Regularization
      testre <- data.table::rbindlist(testre, fill = TRUE)
      tailsinfo <- filter(testre, tail != "not-find")
      notail <- testre %>%
        filter(tail == "not-find") %>%
        select(read_num)
      if (dim(notail)[1] != 0) {
        notails <- fastdf %>%
          filter(read_num %in% notail$read_num)
      } else {
        notails <- c("no notail reads")
      }
      filepath <- str_c(resultpath, samplename, "_notail_reads.txt")
      write.table(notails, filepath, quote = F, row.names = F)
      tailsinfo$sample <- samplename
      tailsinfo <- tailsinfo %>%
        mutate(PAL = nchar(as.character(tail)), nA = str_count(tail, "T"), rt = as.numeric(nA / PAL)) %>%
        select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
      return(tailsinfo)
    }
  }
  # mapping
  else {
    if (findUmi) {
      testre <- apply(fastdf, 1, map_findumi,
        findUmi = findUmi,
        lumi = lumi, mcans = mcans,
        adapterSeq = adapterSeq, anchorSeq = anchorSeq,
        mapping = mapping,
        tailAnchorLen = tailAnchorLen, findTailType = findTailType
      )
      # 4.Regularization
      testre <- data.table::rbindlist(testre, fill = TRUE)
      tailsinfo <- filter(testre, tail != "not-find")
      notail <- testre %>%
        filter(tail == "not-find") %>%
        select(read_num)
      if (dim(notail)[1] != 0) {
        notails <- fastdf %>%
          filter(read_num %in% notail$read_num)
      } else {
        notails <- c("no notail reads")
      }
      filepath <- str_c(resultpath, samplename, "_notail_reads.txt")
      write.table(notails, filepath, quote = F, row.names = F)
      tailsinfo$sample <- samplename
      tailsinfo <- tailsinfo %>%
        mutate(PAL = nchar(as.character(tail)), nA = str_count(tail, "T"), rt = as.numeric(nA / PAL)) %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      return(tailsinfo)
    } else {
      testre <- apply(fastdf, 1, map_noumi,
        findUmi = findUmi,
        lumi = lumi, mcans = mcans,
        adapterSeq = adapterSeq, anchorSeq = anchorSeq,
        mapping = mapping,
        tailAnchorLen = tailAnchorLen, findTailType = findTailType
      )
      # 4.Regularization
      testre <- data.table::rbindlist(testre, fill = TRUE)
      tailsinfo <- filter(testre, tail != "not-find")
      notail <- testre %>%
        filter(tail == "not-find") %>%
        select(read_num)
      if (dim(notail)[1] != 0) {
        notails <- fastdf %>%
          filter(read_num %in% notail$read_num)
      } else {
        notails <- c("no notail reads")
      }
      filepath <- str_c(resultpath, samplename, "_notails.txt")
      write.table(notails, filepath, quote = F, row.names = F)
      tailsinfo$sample <- samplename
      tailsinfo <- tailsinfo %>%
        mutate(PAL = nchar(as.character(tail)), nA = str_count(tail, "T"), rt = as.numeric(nA / PAL)) %>%
        select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
      return(tailsinfo)
    }
  }
}


# TailFilter --------------------------------------------------------------
tailFilter <- function(tailsinfo, findUmi, anchorSeq, minTailLen, realTailLen, maxNtail) {
  # 1.Default value setting
  if (missing(minTailLen)) {
    minTailLen <- 8
  }
  if (missing(realTailLen)) {
    realTailLen <- 30
  }
  if (missing(maxNtail)) {
    maxNtail <- 2
  }
  # 2.Determining the availability of structural inputs
  if (!(missing(anchorSeq))) {
    # If there is a structure, the structured tail is extracted first
    tailinfo1 <- tailsinfo %>%
      filter(tailType == "structural")
    # For unstructured tails
    tailinfo2 <- tailsinfo %>%
      filter(tailType == "unstructural")
    tailinfo3 <- tailinfo2 %>%
      filter(PAL >= minTailLen)
    tailsinfo <- base::rbind(tailinfo1, tailinfo3)
  } else {
    tailsinfo <- tailsinfo %>%
      filter(PAL >= minTailLen)
  }
  # 3.Trade-offs for multi-tail reads
  # The traversal method is used here
  ## First calculate the number of tails for each read
  # tailsNumber <- tailsinfo %>%
  #   dplyr::group_by(read_num) %>%
  #   dplyr::summarise(tailN = n())
  # tailsinfo1 <- tailsNumber %>%
  #   dplyr::filter(tailN<maxNtail)
  # tails1 <- tailsinfo %>%
  #   dplyr::filter(read_num%in%tailsinfo1$read_num)
  # tailsinfo2 <- tailsNumber %>%
  #   dplyr::filter(tailN>=maxNtail)
  # tails2 <- data.frame()
  # i = 0
  # for(read in as.character(tailsinfo2$read_num)){
  #   tailinfo <- tailsinfo %>%
  #     filter(read_num == read) %>%
  #     mutate(maxL = max(PAL)/5) %>%
  #     filter(!(PAL<realTailLen)|!(PAL<maxL))
  #   tails2 <- base::rbind(tails2,tailinfo)
  #   i=i+1
  #   if(i%% 10000 == 0){
  #     print(str_c("---",i," reads processed---"))
  #   }
  # }
  # tailsinfo <- base::rbind(tails1,tails2)
  # Modify here to filter
  # a.Introduce minimum length
  tailsinfo$realTailLen <- realTailLen
  # b.Introduction of group maximums
  tailsinfo <- tailsinfo %>%
    dplyr::group_by(read_num) %>%
    summarize(PALmax = max(PAL) / 5) %>%
    ungroup() %>%
    inner_join(tailsinfo)
  tailsinfo <- tailsinfo %>%
    mutate(flag1 = PAL - PALmax, flag2 = PAL - realTailLen) %>%
    filter(!(flag1 < 0 & flag2 < 0))
  if (findUmi) {
    tailsinfo <- tailsinfo %>%
      select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
  } else {
    tailsinfo <- tailsinfo %>%
      select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
  }
  return(tailsinfo)
}

# tailClassify -----------------------------------------------------------------
tailClassify <- function(tailsinfo, findUmi, maxNtail, mapping) {
  if (missing(maxNtail)) {
    maxNtail <- 2
  }
  if (mapping) {
    if (findUmi) {
      # A.First is the third screening
      tails_sum <- tailsinfo %>%
        dplyr::group_by(read_num) %>%
        dplyr::summarise(count = n())
      # 1.First filter out the total number of tails <= maxN for direct appointment
      tailsinfo1 <- tails_sum %>%
        dplyr::filter(count <= maxNtail)
      tailsinfo1 <- tailsinfo %>%
        dplyr::filter(read_num %in% as.character(tailsinfo1$read_num)) %>%
        dplyr::select(read_num, chr, strand, coord, umi, PAL, tail, tailType, nA, rt, sample)
      # 2.Filter for reads with total tails > maxN
      tailsinfo2 <- tails_sum %>%
        dplyr::filter(count > maxNtail)
      tailsinfo2 <- tailsinfo %>%
        dplyr::filter(read_num %in% as.character(tailsinfo2$read_num))
      tailsinfo2 <- tailsinfo2 %>%
        dplyr::group_by(read_num) %>%
        dplyr::summarise(tails_c = n()) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(tailsinfo2)
      tailsinfo2 <- tailsinfo2 %>%
        dplyr::group_by(read_num, tailType) %>%
        dplyr::summarise(type_c = n()) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(tailsinfo2)
      tailsinfo2$maxN <- maxNtail
      structural_num <- tailsinfo2 %>%
        dplyr::filter(tailType == "structural") %>%
        dplyr::select(read_num, type_c)
      structural_num <- structural_num[!duplicated(structural_num$read_num), ]
      structural_num <- rename(structural_num, struc_c = type_c)
      uns_num <- tailsinfo2 %>%
        filter(tailType == "unstructural" & (tails_c - type_c == 0)) %>%
        select(read_num, type_c)
      uns_num <- uns_num[!duplicated(uns_num$read_num), ]
      uns_num <- uns_num %>%
        mutate(struc_c = 0) %>%
        select(read_num, struc_c)
      structural_num <- base::rbind(uns_num, structural_num)
      tailsinfo2 <- inner_join(tailsinfo2, structural_num, by = "read_num")
      ## 2.1Take the top maxN entries that are all structrual
      sub1 <- tailsinfo2 %>%
        filter(tailType == "structural" & struc_c >= maxNtail) %>%
        group_by(read_num) %>%
        top_n(maxNtail, PAL) %>%
        ungroup() %>%
        select(read_num, chr, strand, coord, umi, PAL, tail, tailType, nA, rt, sample)
      ## 2.2Take all str's that are not all structrual and the top maxN-str_c bar in uns
      sub2 <- tailsinfo2 %>%
        filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "structural") %>%
        select(read_num, chr, strand, coord, umi, PAL, tail, tailType, nA, rt, sample)
      sub3 <- tailsinfo2 %>%
        filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "unstructural") %>%
        group_by(read_num) %>%
        mutate(topx = maxN - struc_c) %>%
        top_n(topx, PAL) %>%
        ungroup() %>%
        select(read_num, chr, strand, coord, umi, PAL, tail, tailType, nA, rt, sample)
      tailsinfo <- base::rbind(tailsinfo1, sub1, sub2, sub3)
      # B.Next is the tail classification
      ## 1.Count the number of tails in groups first
      tailsinfo <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(tail_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo)
      tailsinfo$ids <- c(1:dim(tailsinfo)[1])
      if (maxNtail == 1) {
        # Straightforward to keep with only one tail
        tailsinfo1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        # If you have multiple tails, choose the one with the longest tail
        # If there are more than one with the longest tail, choose one at random
        tailsinfo2 <- tailsinfo %>%
          filter(tail_c != 1) %>%
          group_by(read_num) %>%
          mutate(maxPAL = max(PAL)) %>%
          filter(PAL == maxPAL) %>%
          top_n(1, ids) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        # Consolidated results
        tailsinfo <- base::rbind(as.data.frame(tailsinfo1), as.data.frame(tailsinfo2))
        return(tailsinfo)
      } else if (maxNtail == 2) {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        # When the actual number of tails is two in excess: first cluster
        # according to PAL, leaving only one record per PAL; then classify all
        # remaining records according to the original method.
        tailsinfo3 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num, PAL) %>%
          top_n(1, PAL)
        tailsinfo3 <- select(tailsinfo3, -tail_c)
        tailsinfo3 <- tailsinfo3 %>%
          group_by(read_num) %>%
          summarise(tail_c = n()) %>%
          ungroup() %>%
          inner_join(tailsinfo3)
        # Only one record
        tailssub3 <- tailsinfo3 %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        # There are still two records
        tailssub41 <- tailsinfo3 %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub42 <- tailsinfo3 %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        # Greater than two records
        tailssub51 <- tailsinfo3 %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub52 <- tailsinfo3 %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub4 <- rbind(tailssub41, tailssub42)
        tailssub5 <- rbind(tailssub51, tailssub52)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsinfo <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailsinfo <- rbind(as.data.frame(tailssub3), as.data.frame(tailsinfo))
        tailsinfo <- rbind(as.data.frame(tailssub4), as.data.frame(tailsinfo))
        tailsinfo <- rbind(as.data.frame(tailssub5), as.data.frame(tailsinfo))
        return(tailsinfo)
      } else {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsa <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailssub3 <- tailsinfo %>%
          filter(!(read_num %in% as.character(tailsa$read_num))) %>%
          filter((tail_c > 2) & (tail_c <= maxNtail)) %>%
          mutate(read_type = "multi-tail") %>%
          select(read_num, chr, strand, coord, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsa), as.data.frame(tailssub3))
        return(tailsinfo)
      }
    } else {
      # A.First is the third screening
      tails_sum <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(count = n())
      # 1.First filter out the total number of tails <= maxN for direct
      # appointment
      tailsinfo1 <- tails_sum %>%
        filter(count <= maxNtail)
      tailsinfo1 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo1$read_num)) %>%
        select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
      # 2.Filter for reads with total tails > maxN
      tailsinfo2 <- tails_sum %>%
        filter(count > maxNtail)
      tailsinfo2 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo2$read_num))
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num) %>%
        summarise(tails_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num, tailType) %>%
        summarise(type_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2$maxN <- maxNtail
      structural_num <- tailsinfo2 %>%
        filter(tailType == "structural") %>%
        select(read_num, type_c)
      structural_num <- structural_num[!duplicated(structural_num$read_num), ]
      structural_num <- rename(structural_num, struc_c = type_c)
      uns_num <- tailsinfo2 %>%
        filter(tailType == "unstructural" & (tails_c - type_c == 0)) %>%
        select(read_num, type_c)
      uns_num <- uns_num[!duplicated(uns_num$read_num), ]
      uns_num <- uns_num %>%
        mutate(struc_c = 0) %>%
        select(read_num, struc_c)
      structural_num <- base::rbind(uns_num, structural_num)
      tailsinfo2 <- inner_join(tailsinfo2, structural_num, by = "read_num")
      ## 2.1Take the top maxN entries that are all structrual
      sub1 <- tailsinfo2 %>%
        filter(tailType == "structural" & type_c >= maxNtail) %>%
        group_by(read_num) %>%
        top_n(maxNtail, PAL) %>%
        ungroup() %>%
        select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
      ## 2.2Take all str's that are not all structrual and the top maxN-str_c bar
      ## in uns
      reads1 <- as.data.frame(sub1$read_num)
      reads1 <- reads1[!duplicated(reads1$`sub1$read_num`), ]
      reads2 <- as.data.frame(tailsinfo2$read_num)
      reads2 <- reads2[!duplicated(reads2$`tailsinfo2$read_num`), ]
      diff <- setdiff(reads1, reads2)
      if (!(is_empty(diff))) {
        sub2 <- tailsinfo2 %>%
          filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "structural") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        sub3 <- tailsinfo2 %>%
          filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "unstructural") %>%
          group_by(read_num) %>%
          mutate(topx = maxN - struc_c) # %>%
        top_n(topx, PAL) %>%
          ungroup() %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        tailsinfo <- base::rbind(tailsinfo1, sub1, sub2, sub3)
      } else {
        tailsinfo <- base::rbind(tailsinfo1, sub1)
      }
      tailsinfo <- as.data.frame(tailsinfo)
      # B.Next is the tail classification
      ## 1.Count the number of tails in groups first
      tailsinfo <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(tail_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo)
      if (dim(tailsinfo)[1] != 0) {
        tailsinfo$ids <- c(1:dim(tailsinfo)[1])
      } else {
        stop("something wrong!")
      }
      if (maxNtail == 1) {
        tailsinfo1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        # If you have multiple tails, choose the one with the longest tail
        # If there are more than one with the longest tail, choose one at random
        tailsinfo2 <- tailsinfo %>%
          filter(tail_c != 1) %>%
          group_by(read_num) %>%
          mutate(maxPAL = max(PAL)) %>%
          filter(PAL == maxPAL) %>%
          top_n(1, ids) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        # Consolidated results
        tailsinfo <- base::rbind(as.data.frame(tailsinfo1), as.data.frame(tailsinfo2))
      } else if (maxNtail == 2) {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        # When the actual number of tails is two in excess: first cluster
        # according to PAL, leaving only one record per PAL; then classify all
        # remaining records according to the original method.
        tailsinfo3 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num, PAL) %>%
          top_n(1, ids)
        tailsinfo3 <- select(tailsinfo3, -tail_c)
        tailsinfo3 <- tailsinfo3 %>%
          group_by(read_num) %>%
          summarise(tail_c = n()) %>%
          ungroup() %>%
          inner_join(tailsinfo3)
        # Only one record
        tailssub3 <- tailsinfo3 %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        # There are still two records
        tailssub41 <- tailsinfo3 %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub42 <- tailsinfo3 %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        # Greater than two records
        tailssub51 <- tailsinfo3 %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub52 <- tailsinfo3 %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub4 <- rbind(tailssub41, tailssub42)
        tailssub5 <- rbind(tailssub51, tailssub52)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsinfo <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailsinfo <- rbind(as.data.frame(tailssub3), as.data.frame(tailsinfo))
        tailsinfo <- rbind(as.data.frame(tailssub4), as.data.frame(tailsinfo))
        tailsinfo <- rbind(as.data.frame(tailssub5), as.data.frame(tailsinfo))
        return(tailsinfo)
      } else {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsa <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailssub3 <- tailsinfo %>%
          filter(!(read_num %in% as.character(tailsa$read_num))) %>%
          filter((tail_c > 2) & (tail_c <= maxNtail)) %>%
          mutate(read_type = "multi-tail") %>%
          select(read_num, chr, strand, coord, PAL, tail, tailType, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsa), as.data.frame(tailssub3))
        return(tailsinfo)
      }
    }
  } else {
    if (findUmi) {
      tails_sum <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(count = n())
      tailsinfo1 <- tails_sum %>%
        filter(count <= maxNtail)
      tailsinfo1 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo1$read_num)) %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      tailsinfo2 <- tails_sum %>%
        filter(count > maxNtail)
      tailsinfo2 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo2$read_num))
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num) %>%
        summarise(tails_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num, tailType) %>%
        summarise(type_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2$maxN <- maxNtail
      structural_num <- tailsinfo2 %>%
        filter(tailType == "structural") %>%
        select(read_num, type_c)
      structural_num <- structural_num[!duplicated(structural_num$read_num), ]
      structural_num <- rename(structural_num, struc_c = type_c)
      uns_num <- tailsinfo2 %>%
        filter(tailType == "unstructural" & (tails_c - type_c == 0)) %>%
        select(read_num, type_c)
      uns_num <- uns_num[!duplicated(uns_num$read_num), ]
      uns_num <- uns_num %>%
        mutate(struc_c = 0) %>%
        select(read_num, struc_c)
      structural_num <- base::rbind(uns_num, structural_num)
      tailsinfo2 <- inner_join(tailsinfo2, structural_num, by = "read_num")
      sub1 <- tailsinfo2 %>%
        filter(tailType == "structural" & struc_c >= maxNtail) %>%
        group_by(read_num) %>%
        top_n(maxNtail, PAL) %>%
        ungroup() %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      sub2 <- tailsinfo2 %>%
        filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "structural") %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      sub3 <- tailsinfo2 %>%
        filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "unstructural") %>%
        group_by(read_num) %>%
        mutate(topx = maxN - struc_c) %>%
        top_n(topx, PAL) %>%
        ungroup() %>%
        select(read_num, strand, umi, PAL, tail, tailType, nA, rt, sample)
      tailsinfo <- base::rbind(tailsinfo1, sub1, sub2, sub3)
      tailsinfo <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(tail_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo)
      tailsinfo$ids <- c(1:dim(tailsinfo)[1])
      if (maxNtail == 1) {
        tailsinfo1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo2 <- tailsinfo %>%
          filter(tail_c != 1) %>%
          group_by(read_num) %>%
          mutate(maxPAL = max(PAL)) %>%
          filter(PAL == maxPAL) %>%
          top_n(1, ids) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsinfo1), as.data.frame(tailsinfo2))
      } else if (maxNtail == 2) {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub31 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub32 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub3 <- rbind(tailssub31, tailssub32)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsinfo <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailsinfo <- rbind(as.data.frame(tailssub3), as.data.frame(tailsinfo))
        return(tailsinfo)
      } else {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsa <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailssub3 <- tailsinfo %>%
          filter(!(read_num %in% as.character(tailsa$read_num))) %>%
          filter((tail_c > 2) & (tail_c <= maxNtail)) %>%
          mutate(read_type = "multi-tail") %>%
          select(read_num, strand, umi, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsa), as.data.frame(tailssub3))
        return(tailsinfo)
      }
    } else {
      tails_sum <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(count = n())
      tailsinfo1 <- tails_sum %>%
        filter(count <= maxNtail)
      tailsinfo1 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo1$read_num)) %>%
        select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
      tailsinfo2 <- tails_sum %>%
        filter(count > maxNtail)
      tailsinfo2 <- tailsinfo %>%
        filter(read_num %in% as.character(tailsinfo2$read_num))
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num) %>%
        summarise(tails_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2 <- tailsinfo2 %>%
        group_by(read_num, tailType) %>%
        summarise(type_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo2)
      tailsinfo2$maxN <- maxNtail
      structural_num <- tailsinfo2 %>%
        filter(tailType == "structural") %>%
        select(read_num, type_c)
      structural_num <- structural_num[!duplicated(structural_num$read_num), ]
      structural_num <- rename(structural_num, struc_c = type_c)
      uns_num <- tailsinfo2 %>%
        filter(tailType == "unstructural" & (tails_c - type_c == 0)) %>%
        select(read_num, type_c)
      uns_num <- uns_num[!duplicated(uns_num$read_num), ]
      uns_num <- uns_num %>%
        mutate(struc_c = 0) %>%
        select(read_num, struc_c)
      structural_num <- base::rbind(uns_num, structural_num)
      tailsinfo2 <- inner_join(tailsinfo2, structural_num, by = "read_num")
      sub1 <- tailsinfo2 %>%
        filter(tailType == "structural" & type_c >= maxNtail) %>%
        group_by(read_num) %>%
        top_n(maxNtail, PAL) %>%
        ungroup() %>%
        select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
      reads1 <- as.data.frame(sub1$read_num)
      reads1 <- reads1[!duplicated(reads1$`sub1$read_num`), ]
      reads2 <- as.data.frame(tailsinfo2$read_num)
      reads2 <- reads2[!duplicated(reads2$`tailsinfo2$read_num`), ]
      diff <- setdiff(reads1, reads2)
      if (!(is_empty(diff))) {
        sub2 <- tailsinfo2 %>%
          filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "structural") %>%
          select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
        sub3 <- tailsinfo2 %>%
          filter(!(read_num %in% as.character(sub1$read_num)) & tailType == "unstructural") %>%
          group_by(read_num) %>%
          mutate(topx = maxN - struc_c) # %>%
        top_n(topx, PAL) %>%
          ungroup() %>%
          select(read_num, strand, PAL, tail, tailType, nA, rt, sample)
        tailsinfo <- base::rbind(tailsinfo1, sub1, sub2, sub3)
      } else {
        tailsinfo <- base::rbind(tailsinfo1, sub1)
      }
      tailsinfo <- as.data.frame(tailsinfo)
      tailsinfo <- tailsinfo %>%
        group_by(read_num) %>%
        summarise(tail_c = n()) %>%
        ungroup() %>%
        inner_join(tailsinfo)
      tailsinfo$ids <- c(1:dim(tailsinfo)[1])
      if (maxNtail == 1) {
        tailsinfo1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo2 <- tailsinfo %>%
          filter(tail_c != 1) %>%
          group_by(read_num) %>%
          mutate(maxPAL = max(PAL)) %>%
          filter(PAL == maxPAL) %>%
          top_n(1, ids) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsinfo1), as.data.frame(tailsinfo2))
      } else if (maxNtail == 2) {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub31 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub32 <- tailsinfo %>%
          filter(tail_c > 2) %>%
          group_by(read_num) %>%
          top_n(2, ids) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub3 <- rbind(tailssub31, tailssub32)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsinfo <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailsinfo <- rbind(as.data.frame(tailssub3), as.data.frame(tailsinfo))
        return(tailsinfo)
      } else {
        tailssub1 <- tailsinfo %>%
          filter(tail_c == 1) %>%
          mutate(read_type = "one-tail") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub21 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] == strand[2]) %>%
          mutate(read_type = "two-tail-same") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub22 <- tailsinfo %>%
          filter(tail_c == 2) %>%
          group_by(read_num) %>%
          filter(strand[1] != strand[2]) %>%
          mutate(read_type = "two-tail-mixed") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailssub2 <- rbind(tailssub22, tailssub21)
        tailsa <- rbind(as.data.frame(tailssub1), as.data.frame(tailssub2))
        tailssub3 <- tailsinfo %>%
          filter(!(read_num %in% as.character(tailsa$read_num))) %>%
          filter((tail_c > 2) & (tail_c <= maxNtail)) %>%
          mutate(read_type = "multi-tail") %>%
          select(read_num, strand, PAL, tail, tailType, read_type, nA, rt, sample)
        tailsinfo <- base::rbind(as.data.frame(tailsa), as.data.frame(tailssub3))
        return(tailsinfo)
      }
    }
  }
}


seq_rev <- function(char) {
  alphabets <- strsplit(char, split = "")[[1]]
  return(rev(alphabets))
}
seq_compl <- function(seq) {
  # Check if there's "U" in the sequence
  RNA <- Reduce(`|`, seq == "U")
  cmplvec <- sapply(seq, function(base) {
    # This makes DNA the default
    # As long as there's no U, the sequence is treated as DNA
    if (RNA) {
      switch(base, "A" = "U", "C" = "G", "G" = "C", "U" = "A")
    } else {
      switch(base, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
    }
  })
  return(paste(cmplvec, collapse = ""))
}
seqreverse <- function(x){
  seq <- x[6]
  seq <- seq_compl(seq_rev(seq))
  return(seq)
}
tailExtract <- function(x,minTailLen){
  seq <- x[5]
  l <- length(x)
  mcans <- x[l]
  # print(seq)
  #Initialisation
  pattern <- str_c("T{",minTailLen,",}")
  # pattern <- str_c("T{8,}")
  anchors=gregexpr(pattern = pattern,seq)
  anchors=anchors[[1]]
  tailStart=as.integer(anchors)[1]
  if(tailStart==-1){
    return("no-tail")
  }
  else{
    currMax=attr(anchors,"match.length")[1]
    currMaxIdx=tailStart+currMax
    drop=0
    curr=currMax
    strs=unlist(strsplit(seq, split=''))
    if(currMaxIdx!=length(strs)+1){
      for (i in currMaxIdx:length(strs)) {
        if (strs[i]!='T') {
          curr=curr-1
        } else {
          curr=curr+1
          if (curr>=currMax) {
            currMaxIdx=i
            currMax=curr
          }
        }
        drop=currMax-curr
        #cat(i, currMax, curr, drop, substr(s, tailStart, i), '\n')
        if (drop>=mcans) {
          break
        }
      }
      tail=substr(seq, tailStart, currMaxIdx-1)
      return(tail)
    }
    else{
      tail=substr(seq, tailStart, length(strs))
      return(tail)
      # return(seq)
    }
  }
}
ChangeChrFormat <- function(df,refPath,from,to){
  if(missing(df)){
    stop("df is indispensable!")
  }
  if(missing(from)){
    stop("from format is indispensable!")
  }
  if(missing(to)){
    stop("to format is indispensable!")
  }
  ref <- read.table(refPath,
                    header = T,sep="\t")
  if(!(from %in% colnames(ref)) | !(to %in% colnames(ref))){
    stop("from or to must be one in c('GenBank.Accn','RefSeq.Accn','UCSC.style.name')!")
  }
  ref <- select(ref,all_of(from),all_of(to))
  colnames(ref) <- str_replace(colnames(ref),from,"chr")
  df <- left_join(df,ref,by="chr")
  df <- df[!is.na(df$UCSC.style.name),]
  df <- df %>%
    select(-chr) %>%
    rename(chr = to)
  return(df)
}



# tailScan --------------------------------------------------------------------
#' @title tailScan 
#'
#' @description \code{tailScan} Returns a table containing at least read_num,
#'    tail length and tail sequence. Use a variety of methods to help you quantify 
#'   the tails in a sequence.
#'
#' @details This function quantifies the possible tails of all sequences in a
#'   FASTQ file with a non-aligned manner.You need to specify the parameters
#'   according to the structure of your sequence.We will save the found tail
#'   data and the sequence data that did not find the tail to the path you
#'   specified
#'
#' @param fastq The path of fastqfile.
#' @param mcans The maximum allowable mismatch number in the sliding window
#'   algorithm,default=5.
#' @param findUmi Boolean value.Indicates whether the sequence structure
#'   contains UMI or barcode.If it is ture, the UMI or Barcode will be extracted
#'   separately.
#' @param lumi The length of umi in reads. "lumi = 0" means there is no need to
#'   extract umi from reads.
#' @param adapterSeq character.If you enter a FASTQ file that does not
#'   remove the 3 'adapter, please provide the full sequence of adapters.
#' @param anchorSeq character.If your sequence structure has a sequence of
#'   anchor points identifying tails, enter this parameter.
#' @param samplename Specify a sample name for your data.
#' @param resultpath The path where you want to store the result data.
#' @param minTailLen Specifies the minimum tail length,default=8.
#' @param tailAnchorLen Specifies the minimum tail anchor point
#'   length,default=8.
#' @param realTailLen Specifies what you think is the true tail
#'   length,default=30.
#' @param maxNtail Specifies the maximum number of tails that should be found in
#'   a sequence,default=2.
#' @param mapping Boolean value.The default value is F.
#' @param findTailType Used to specify the type of tail to extract from the
#'   sequence. There are three values, 'A' means to search for poly(A) tails,
#'   'T' means to search for poly(T) tails or 'both' means to search at the same
#'   time. The default is 'both'.
#' @return Save the quantitative tail results table of various algorithms to the
#'   path you specify.Meanwhile, return the tail dataframe
#' @family Poly(A) Tail length quantification functions
#' @seealso [tailMap()] to quantitative tails based on sequence algin.
#' @name  tailScan
#' @examples
#' fastqfile <- "./inst/extdata/GV_fastq/PAIso_GV1.fastq"
#' GV1tailDF <- tailScan(fastqfile,
#'   mcans = 5, findUmi = F, resultpath =
#'     "./inst/extdata/output/", samplename =
#'     "GV1", tailAnchorLen = 8, minTailLen = 8, realTailLen = 20, maxNtail = 2, mapping = F
#' )
#' head(GV1tailDF)
#' @export
tailScan <- function(fastq, mcans, findUmi, lumi,
                     adapterSeq, anchorSeq, resultpath, samplename, tailAnchorLen,
                     minTailLen, realTailLen, maxNtail, mapping, mapinfo, findTailType) {
  if (!(is.logical(mapping))) {
    stop("The parameter mapping must be a Boolean!")
  } else {
    if (mapping) {
      if (missing(mapinfo)) {
        stop("The mapinfo lose")
      }
    }
  }
  if (missing(findTailType)) {
    findTailType <- "both"
  }
  if (!(findTailType %in% c("A", "T", "both"))) {
    stop("findTailType must in c('A','T','both')!")
  }
  # 0.fastdf
  rfq <- ShortRead::readFastq(fastq)
  seq <- as.data.frame(ShortRead::sread(rfq))
  seq <- dplyr::rename(seq, seq = x)
  read_num <- as.data.frame(ShortRead::id(rfq))
  read_num <- dplyr::rename(read_num, read_num = x)
  rm(rfq)
  fastdf <- base::cbind(seq, read_num)
  rm(seq)
  rm(read_num)
  fastdf <- fastdf %>%
    mutate(width = nchar(seq)) %>%
    select(width, seq, read_num)
  # change the format of read_num
  testb <- str_split(fastdf$read_num, " ")
  fastdf$read_num <- unlist(lapply(testb, function(testb) testb[[1]][1]))
  rm(testb)
  if (!(mapping)) {
    # 0.read in fastq file
    # testb <- str_split(fastdf$names," ")
    # fastdf$names <- unlist(lapply(testb, function(testb) testb[[1]][1]))
    # 1.set the default parameters
    if (missing(mcans)) {
      mcans <- 5
    }
    if (missing(minTailLen)) {
      minTailLen <- 8
    }
    if (missing(realTailLen)) {
      realTailLen <- 30
    }
    if (missing(maxNtail)) {
      maxNtail <- 2
    }
    if (missing(tailAnchorLen)) {
      tailAnchorLen <- 8
    }
    if (!(is.logical(findUmi))) {
      stop("The parameter findUmi must be a Boolean!")
    } else {
      if (findUmi) {
        if (missing(lumi) | missing(adapterSeq)) {
          stop("When the parameter findUmi is true, the parameter lumi and adapterSeq are indispensable")
        }
      }
    }
    if (missing(anchorSeq)) {
      anchorSeq <- "miss"
    }
    if (missing(resultpath)) {
      stop("resultpath lose")
    }
    # 2.find tails
    num <- dim(fastdf)[1]
    print("---Start quantifying tails---")
    print(str_c(num, " reads"))
    tailsinfo1 <- tailFinder(fastdf, mcans, findUmi, lumi, adapterSeq, anchorSeq, resultpath, samplename, tailAnchorLen, mapping = mapping, findTailType = findTailType)
    print("---Tail hunting over---")
    # 3.filter tails
    print("---Start screening tails---")
    tailsinfo2 <- tailFilter(tailsinfo1, findUmi, anchorSeq, minTailLen, realTailLen, maxNtail)
    # 4.classify tails
    print("---Start classifying tails---")
    tailsinfo3 <- tailClassify(tailsinfo2, findUmi, maxNtail, mapping)
    # 5.return
    tailsinfo3 <- dplyr::rename(tailsinfo3, polyAT = strand)
    return(tailsinfo3)
  } else {
    # 1.set the default parameters
    if (missing(mcans)) {
      mcans <- 5
    }
    if (missing(minTailLen)) {
      minTailLen <- 8
    }
    if (missing(realTailLen)) {
      realTailLen <- 30
    }
    if (missing(maxNtail)) {
      maxNtail <- 2
    }
    if (missing(tailAnchorLen)) {
      tailAnchorLen <- 8
    }
    if (!(is.logical(findUmi))) {
      stop("The parameter findUmi must be a Boolean!")
    } else {
      if (findUmi) {
        if (missing(lumi) | missing(adapterSeq)) {
          stop("When the parameter findUmi is true, the parameter lumi and adapterSeq are indispensable")
        }
      }
    }
    if (missing(anchorSeq)) {
      anchorSeq <- "miss"
    }
    if (missing(resultpath)) {
      stop("resultpath lose")
    }
    # 2.find tails
    num <- dim(mapinfo)[1]
    print("---Start quantifying tails---")
    print(str_c(num, " reads"))
    tailsinfo1 <- tailFinder(mapinfo, mcans, findUmi, lumi, adapterSeq, anchorSeq, resultpath, samplename, tailAnchorLen, mapping = mapping, findTailType = findTailType)
    print("---Tail hunting over---")
    # 3.filter tails
    print("---Start screening tails---")
    tailsinfo2 <- tailFilter(tailsinfo1, findUmi, anchorSeq, minTailLen, realTailLen, maxNtail)
    # 4.classify tails
    print("---Start classifying tails---")
    tailsinfo3 <- tailClassify(tailsinfo2, findUmi, maxNtail, mapping)
    # 5.return
    tailsinfo3 <- dplyr::rename(tailsinfo3, polyAT = strand)
    return(tailsinfo3)
  }
}









# faBuilder ---------------------------------------------------------------
#' @title  faBuilder
#'
#' @description \code{faBuilder} Tails and partial sequences were extracted from
#'   long reads and FASTA files were generated for alignment.
#'
#' @details This function is used to extract the PolyA tail and part of the
#'   sequence before the tail starting site from the FASTQ file containing
#'   Longread, and generate the FASTA file that can be input into the sequence
#'   alignment software.
#'
#' @param fastqfile The path of fastqfile.
#' @param mcans The maximum allowable mismatch number in the sliding window
#'   algorithm,default=5.
#' @param findUmi Boolean value.Indicates whether the sequence structure
#'   contains UMI or barcode.If it is ture, the UMI or Barcode will be extracted
#'   separately.
#' @param lumi The length of umi in reads. "lumi = 0" means there is no need to
#'   extract umi from reads.
#' @param adapterSeq character.If you enter a FASTQ file that does not
#'   remove the 3 'adapter, please provide the full sequence of adapters.
#' @param anchorSeq character.If your sequence structure has a sequence of
#'   anchor points identifying tails, enter this parameter.
#' @param samplename Specify a sample name for your data.
#' @param resultpath The path where you want to store the result data.
#' @param minTailLen Specifies the minimum tail length,default=8.
#' @param tailAnchorLen Specifies the minimum tail anchor point
#'   length,default=8.
#' @param mapping Boolean value.The default value is F.
#' @param findTailType Used to specify the type of tail to extract from the
#'   sequence. There are three values, 'A' means to search for poly(A) tails,
#'   'T' means to search for poly(T) tails or 'both' means to search at the same
#'   time. The default is 'both'.
#' @return Generate a FASTA file for sequence alignment and save it to the
#'   specified directory.
#' @examples
#' fastqfile <- "./inst/extdata/GV_fastq/PAIso_GV1.fastq"
#' faBuilderRE <- faBuilder(fastqfile,
#'   mcans = 5, findUmi = F, resultpath = "./inst/extdata/output/", samplename =
#'     "GV1", tailAnchorLen = 8, mapping = F
#' )
#' head(faBuilderRE)
#' @family Poly(A) Tail length quantification functions
#' @seealso [tailMap()] to quantitative tails based on sequence algin.
#' @export
faBuilder <- function(fastqfile, mcans, findUmi, lumi, adapterSeq,
                      anchorSeq, resultpath, samplename,
                      tailAnchorLen, mapping, findTailType) {
  # 0.fastdf
  rfq <- ShortRead::readFastq(fastqfile)
  seq <- as.data.frame(ShortRead::sread(rfq))
  seq <- rename(seq, seq = x)
  read_num <- as.data.frame(ShortRead::id(rfq))
  read_num <- rename(read_num, read_num = x)
  rm(rfq)
  fastdf <- base::cbind(seq, read_num)
  rm(seq)
  rm(read_num)
  fastdf <- fastdf %>%
    mutate(width = nchar(seq)) %>%
    select(width, seq, read_num)
  # Modify read_num format
  testb <- str_split(fastdf$read_num, " ")
  fastdf$read_num <- unlist(lapply(testb, function(testb) testb[[1]][1]))
  rm(testb)
  mapping <- T
  if (findUmi) {
    print("------------sub start------------")
    seqtest_fh1 <- tailFinder(fastdf,
      mcans = mcans, findUmi = findUmi, lumi = lumi,
      adapterSeq = adapterSeq, anchorSeq = anchorSeq,
      resultpath = resultpath, samplename = samplename,
      tailAnchorLen = tailAnchorLen, mapping = mapping, findTailType = findTailType
    )
    # longRead=longRead)
    seqtest_fh1 <- seqtest_fh1 %>%
      unite(read_num, read_num, umi, tailType, sample, sep = "_")
    print("------------sub end------------")
    test <- select(seqtest_fh1, read_num, tail)
    filepath <- str_c(resultpath, "subseq.fasta")
    fa <- seqRFLP::dataframe2fas(test, filepath)
    print(str_c("The sub-sequence FASTA file has been saved to the path:'", filepath, "'"))
    return(seqtest_fh1)
  } else {
    print("------------sub start------------")
    seqtest_fh1 <- tailFinder(fastdf,
      mcans = mcans, findUmi = findUmi, lumi = lumi,
      adapterSeq = adapterSeq, anchorSeq = anchorSeq,
      resultpath = resultpath, samplename = samplename,
      tailAnchorLen = tailAnchorLen, mapping = mapping, findTailType = findTailType
    )
    # longRead=longRead)
    seqtest_fh1 <- seqtest_fh1 %>%
      unite(read_num, read_num, tailType, sample, sep = "_")
    print("------------sub end------------")
    test <- select(seqtest_fh1, read_num, tail)
    filepath <- str_c(resultpath, "subseq.fasta")
    fa <- seqRFLP::dataframe2fas(test, filepath)
    print(str_c("The sub-sequence FASTA file has been saved to the path:'", filepath, "'"))
    return(seqtest_fh1)
  }
}


# tailMap ---------------------------------------------------------------
#' @title  tailMap
#'
#' @description \code{tailMap}Extract necessary comment information from BAM
#'   files and remove IP (internal priming).
#'
#' @details If your sequence is longreads, this function is
#'   used to parse BAM files, remove the IP (internal priming) from the initial
#'   tail based on the information in the cigar field in the BAM files, and
#'   finally return the dataframe with all the information of the tail.However,
#'   if your sequence type is shortreads, this function is used to extract coord
#'   annotation information from BAM files.
#'
#' @param bamfile The path of bamfile.
#' @param mcans The maximum allowable mismatch number in the sliding window
#'   algorithm,default=5.
#' @param findUmi Boolean value.Indicates whether the sequence structure
#'   contains UMI or barcode.If it is ture, the UMI or Barcode will be extracted
#'   separately.
#' @param minTailLen Specifies the minimum tail length,default=8.
#' @param maxNtail Specifies the maximum number of tails that should be found in
#'   a sequence,default=2.
#' @param mapping Boolean value. The default value is F.
#' @param longRead Boolean value. If your sequence type is longreads this
#'   parameter is T, otherwise it is F, and the default is T.
#' @param directBam Logical parameter. If you want to skip the tail_scan step
#'   and directly use tail_map to analyze the compared bam file, set this
#'   parameter to true, and the default value is true.
#' @param findTailType Used to specify the type of tail to extract from the
#'   sequence. There are three values, 'A' means to search for poly(A) tails,
#'   'T' means to search for poly(T) tails or 'both' means to search at the same
#'   time. The default is 'both'.
#'
#' @return A dataframe with all the information of the tail.Include at least
#'   read_num,tail,PAL,chr,strand,tailType,read_type and sample.
#' @examples
#' bamfile <- "./inst/extdata/GV_algin/PAIso-GV1.sorted.bam"
#' GV1tailMapre <- tailMap(bamfile, mcans = 5, minTailLen = 8, findUmi = F, longRead = T)
#' @family Poly(A) Tail length quantification functions
#' @seealso [tailScan()] to quantitative tails without sequence algin.
#' @export
tailMap <- function(bamfile, mcans, minTailLen, findUmi, maxNtail, mapping, longRead, directBam, findTailType) {
  if (missing(longRead)) {
    longRead <- T
  }
  if (missing(directBam)) {
    directBam <- T
  }
  if (missing(findTailType)) {
    findTailType <- "both"
  }
  if (!(findTailType %in% c("both", "A", "T"))) {
    stop("findTailType must in c('both','A','T')!")
  }
  pattern1 <- str_c("T{", 8, ",}")
  pattern2 <- str_c("A{", 8, ",}")
  if (longRead) {
    if (directBam) {
      if (missing(mcans)) {
        mcans <- 5
      }
      if (missing(minTailLen)) {
        minTailLen <- 8
      }
      if (missing(findUmi)) {
        stop("Please enter the parameter findUmi")
      }
      if (missing(maxNtail)) {
        maxNtail <- 2
      }
      if (missing(mapping)) {
        # 20211220change
        # stop("Please enter the parameter mapping")
        mapping <- F
      }
      # As direct extraction can lead to large discrepancies the choice was made
      # to use the extraction as a full process run after the area on the
      # comparison
      print("===Program starts execution===")
      print("---Parsing the BAM file---")
      bam <- scanBam(bamfile)
      print("---Parsing the BAM file successfully---")
      print("---Begin extracting the mapping information---")
      allcigars2 <- data.frame(
        read_num = bam[[1]][["qname"]],
        chr = bam[[1]][["rname"]],
        cigar = bam[[1]][["cigar"]],
        strand = bam[[1]][["strand"]],
        seq = as.character(bam[[1]][["seq"]]),
        coord = bam[[1]][["pos"]] + 1
      )
      allcigars2 <- allcigars2 %>%
        dplyr::mutate(
          anchor1 = str_detect(seq, pattern1),
          anchor2 = str_detect(seq, pattern2)
        )
      rm(bam)
      invisible(gc())
      print("---Mapping information extraction was successful---")
      # 1.First extract the uncompared ones to test
      # test <- filter(allcigars2,is.na(cigar))
      # 2.Extract the matched ones to allcigars2
      allcigars2 <- dplyr::filter(allcigars2, !(is.na(cigar)))
      # 3.Start identifying the exact tail
      # 3.1 Extraction of unmatched parts
      allcigars21 <- allcigars2 %>%
        dplyr::filter(strand == "-")
      allcigars21 <- allcigars21 %>%
        dplyr::mutate(
          str = str_extract(cigar, "[0-9]+S"),
          PAL = as.numeric(str_extract(str, "[0-9]+")),
          tail = str_sub(seq, 1, PAL)
        ) %>%
        dplyr::select(read_num, chr, strand, coord, PAL, tail, anchor1, anchor2)
      allcigars21 <- na.omit(allcigars21)
      allcigars22 <- allcigars2 %>%
        dplyr::filter(strand == "+")
      allcigars22 <- allcigars22 %>%
        dplyr::mutate(
          str = str_extract(cigar, "[0-9]+S"),
          PAL = as.numeric(str_extract(str, "[0-9]+")),
          tail = str_sub(seq, -PAL, -1)
        ) %>%
        dplyr::select(read_num, chr, strand, coord, PAL, tail, anchor1, anchor2)
      allcigars22 <- na.omit(allcigars22)
      allcigars3 <- base::rbind(allcigars21, allcigars22)
      allcigars3 <- na.omit(allcigars3)
      rm(allcigars2)
      rm(allcigars21)
      rm(allcigars22)
      invisible(gc())
      print("---start to extract the tail---")
      # Seeking reverse complementarity
      allcigars3 <- dplyr::filter(allcigars3, PAL >= minTailLen)
      # 3.3First all tails for reverse complementarity
      allcigars31 <- dplyr::filter(allcigars3, strand == "-")
      allcigars32 <- dplyr::filter(allcigars3, strand == "+")
      seqs <- apply(allcigars32, 1, seqreverse)
      allcigars32$tail <- seqs
      rm(seqs)
      invisible(gc())
      allcigars3 <- base::rbind(allcigars31, allcigars32)
      allcigars3 <- allcigars3 %>%
        dplyr::select(read_num, chr, strand, coord, tail, anchor1, anchor2) %>%
        dplyr::rename(seq = tail)
      allcigars3$mcans <- mcans
      seqs <- apply(allcigars3, 1, tailExtract, minTailLen)
      allcigars3$tail <- seqs
      rm(seqs)
      invisible(gc())
      # test<-allcigars3
      if (findUmi) {
        allcigars3 <- allcigars3 %>%
          dplyr::mutate(PAL = nchar(tail), nA = str_count(tail, "T"), rt = nA / PAL) %>%
          dplyr::select(read_num, chr, strand, coord, PAL, tail, nA, rt, anchor1, anchor2) %>%
          tidyr::separate(read_num, into = c("read_num", "umi", "tailType", "sample"), sep = "_")
      } else {
        allcigars3 <- allcigars3 %>%
          dplyr::mutate(PAL = nchar(tail), nA = str_count(tail, "T"), rt = nA / PAL) %>%
          dplyr::select(read_num, chr, strand, coord, PAL, tail, nA, rt, anchor1, anchor2) %>%
          tidyr::separate(read_num, into = c("read_num", "tailType", "sample"), sep = "_")
      }
      print("---Successful extraction of tail---")
      if (findTailType == "both") {
        allcigars3 <- allcigars3
      }
      if (findTailType == "A") {
        allcigars3 <- allcigars3 %>%
          dplyr::filter(anchor2 == T)
      }
      if (findTailType == "T") {
        allcigars3 <- allcigars3 %>%
          dplyr::filter(anchor1 == T)
      }
      allcigars3 <- allcigars3 %>%
        dplyr::select(-anchor1, -anchor2)
      # 3.5Tail classification
      print("---Start tail sorting tails---")
      conflict_prefer("filter", "dplyr")
      conflict_prefer("mutate", "dplyr")
      conflict_prefer("select", "dplyr")
      conflict_prefer("rename", "dplyr")
      conflict_prefer("summarise", "dplyr")
      conflict_prefer("rbind", "base")
      allcigars3 <- tailClassify(tailsinfo = allcigars3, findUmi, maxNtail, mapping = T)
      allci31 <- dplyr::filter(allcigars3, str_detect(read_type, "two-tail*"))
      allci31 <- dplyr::filter(allci31, tail != "no-tail")
      tails_false_multi <- allci31 %>%
        dplyr::group_by(read_num) %>%
        dplyr::summarise(count = n()) %>%
        dplyr::filter(count == 1)
      tails_false_multi <- tails_false_multi$read_num
      allcigars3 <- dplyr::filter(allcigars3, tail != "no-tail")
      a1 <- allcigars3 %>%
        dplyr::filter(read_num %in% tails_false_multi)
      if (dim(a1)[1] != 0) {
        a1$read_type <- "one-tail"
      }
      a2 <- allcigars3 %>%
        dplyr::filter(!(read_num %in% tails_false_multi))
      allcigars3 <- rbind(a1, a2)
      print("---Tail classification completed---")
      print("===The program is coming to an end===")
      # 4.Annotate gene information
      # test3 = annotatePAC(pac = allcigars3, aGFF = GFF, verbose = T)
      # test3 <- test3 %>%
      #   select(read_num,umi,tailType,chr,strand,gene,gene_type)
      # allcigars3 <- allcigars3 %>%
      #   select(-anchor1,-anchor2)
      return(allcigars3)
    }
    if (!directBam) {
      if (missing(mcans)) {
        mcans <- 5
      }
      if (missing(minTailLen)) {
        minTailLen <- 8
      }
      if (missing(findUmi)) {
        stop("Please enter the parameter findUmi")
      }
      if (missing(maxNtail)) {
        maxNtail <- 2
      }
      if (missing(mapping)) {
        # 20211220change
        # stop("Please enter the parameter mapping")
        mapping <- F
      }
      # As direct extraction can lead to large discrepancies the choice was made
      # to use the extraction as a full process run after the area on the
      # comparison
      print("===Program starts execution===")
      print("---Parsing the BAM file---")
      bam <- scanBam(bamfile)
      print("---Parsing the BAM file successfully---")
      print("---Begin extracting the mapping information---")
      allcigars2 <- data.frame(
        read_num = bam[[1]][["qname"]],
        chr = bam[[1]][["rname"]],
        cigar = bam[[1]][["cigar"]],
        strand = bam[[1]][["strand"]],
        seq = as.character(bam[[1]][["seq"]]),
        coord = bam[[1]][["pos"]] + 1
      )
      allcigars2 <- allcigars2 %>%
        dplyr::mutate(
          anchor1 = str_detect(seq, pattern1),
          anchor2 = str_detect(seq, pattern2)
        )
      rm(bam)
      invisible(gc())
      print("---Mapping information extraction was successful---")
      # 1.First extract the uncompared ones to test
      test <- dplyr::filter(allcigars2, is.na(cigar))
      # 2.Extract the matched ones to allcigars2
      allcigars2 <- dplyr::filter(allcigars2, !(is.na(cigar)))
      # 3.Start identifying the exact tail
      # 3.1 Extraction of unmatched parts
      allcigars21 <- allcigars2 %>%
        dplyr::filter(strand == "-")
      allcigars21 <- allcigars21 %>%
        dplyr::mutate(
          str = str_extract(cigar, "[0-9]+S$"),
          PAL = as.numeric(str_extract(str, "[0-9]+")),
          tail = str_sub(seq, -PAL, -1)
        ) %>%
        dplyr::select(read_num, chr, strand, coord, PAL, tail, anchor1, anchor2)
      allcigars22 <- allcigars2 %>%
        dplyr::filter(strand == "+")
      allcigars22 <- allcigars22 %>%
        dplyr::mutate(
          str = str_extract(cigar, "[0-9]+S"),
          PAL = as.numeric(str_extract(str, "[0-9]+")),
          tail = str_sub(seq, 1, PAL)
        ) %>%
        dplyr::select(read_num, chr, strand, coord, PAL, tail, anchor1, anchor2)
      allcigars3 <- base::rbind(allcigars21, allcigars22)
      allcigars3 <- na.omit(allcigars3)
      rm(allcigars2)
      rm(allcigars21)
      rm(allcigars22)
      invisible(gc())
      print("---start to extract the tail---")
      # Seeking reverse complementarity
      allcigars3 <- dplyr::filter(allcigars3, PAL >= minTailLen)
      # 3.3First all tails for reverse complementarity
      allcigars31 <- dplyr::filter(allcigars3, strand == "-")
      allcigars32 <- dplyr::filter(allcigars3, strand == "+")
      seqs <- apply(allcigars31, 1, seqreverse)
      allcigars31$tail <- seqs
      rm(seqs)
      invisible(gc())
      allcigars3 <- base::rbind(allcigars31, allcigars32)
      allcigars3 <- allcigars3 %>%
        dplyr::select(read_num, chr, strand, coord, tail, anchor1, anchor2) %>%
        dplyr::rename(seq = tail)
      allcigars3$mcans <- mcans
      seqs <- apply(allcigars3, 1, tailExtract, minTailLen)
      allcigars3$tail <- seqs
      rm(seqs)
      invisible(gc())
      if (findUmi) {
        allcigars3 <- allcigars3 %>%
          dplyr::mutate(PAL = nchar(tail), nA = str_count(tail, "T"), rt = nA / PAL) %>%
          dplyr::select(read_num, chr, strand, coord, PAL, tail, nA, rt, anchor1, anchor2) %>%
          tidyr::separate(read_num, into = c("read_num", "umi", "tailType", "sample"), sep = "_")
      } else {
        allcigars3 <- allcigars3 %>%
          dplyr::mutate(PAL = nchar(tail), nA = str_count(tail, "T"), rt = nA / PAL) %>%
          dplyr::select(read_num, chr, strand, coord, PAL, tail, nA, rt, anchor1, anchor2) %>%
          tidyr::separate(read_num, into = c("read_num", "tailType", "sample"), sep = "_")
      }
      # allcigars3$sample <- "flh1"
      print("---Successful extraction of tail---")
      if (findTailType == "both") {
        allcigars3 <- allcigars3
      }
      if (findTailType == "A") {
        allcigars3 <- allcigars3 %>%
          filter(anchor2 == T)
      }
      if (findTailType == "T") {
        allcigars3 <- allcigars3 %>%
          filter(anchor1 == T)
      }
      allcigars3 <- allcigars3 %>%
        select(-anchor1, -anchor2)
      # 3.5Tail classification
      print("---Start tail sorting tails---")
      allcigars3 <- tailClassify(tailsinfo = allcigars3, findUmi, maxNtail, mapping = T)
      allci31 <- dplyr::filter(allcigars3, str_detect(read_type, "two-tail*"))
      allci31 <- dplyr::filter(allci31, tail != "no-tail")
      tails_false_multi <- allci31 %>%
        dplyr::group_by(read_num) %>%
        dplyr::summarise(count = n()) %>%
        dplyr::filter(count == 1)
      tails_false_multi <- tails_false_multi$read_num
      allcigars3 <- dplyr::filter(allcigars3, tail != "no-tail")
      a1 <- allcigars3 %>%
        dplyr::filter(read_num %in% tails_false_multi)
      if (dim(a1)[1] != 0) {
        a1$read_type <- "one-tail"
      }
      a2 <- allcigars3 %>%
        dplyr::filter(!(read_num %in% tails_false_multi))
      allcigars3 <- rbind(a1, a2)
      print("---Tail classification completed---")
      print("===The program is coming to an end===")
      # 4.Annotate gene information
      # test3 = annotatePAC(pac = allcigars3, aGFF = GFF, verbose = T)
      # test3 <- test3 %>%
      #   select(read_num,umi,tailType,chr,strand,gene,gene_type)
      return(allcigars3)
    }
  } else {
    # 1.Initialization parameters
    if (missing(mcans)) {
      mcans <- 5
    }
    if (missing(minTailLen)) {
      minTailLen <- 8
    }
    if (missing(findUmi)) {
      findUmi <- T
    }
    if (missing(maxNtail)) {
      maxNtail <- 2
    }
    if (missing(mapping)) {
      mapping <- T
    }
    # 2.Extraction of necessary information
    print("===Program starts execution===")
    print("---Parsing the BAM file---")
    bam <- scanBam(bamfile)
    print("---Parsing the BAM file successfully---")
    print("---Begin extracting the mapping information---")
    allcigars2 <- data.frame(
      read_num = bam[[1]][["qname"]],
      chr = bam[[1]][["rname"]],
      strand = bam[[1]][["strand"]],
      coord = bam[[1]][["pos"]] + 1
    )
    rm(bam)
    print("---Mapping information extraction was successful---")
    return(allcigars2)
  }
}



# geneAnno function
# geneAnno ----------------------------------------------------------------
#' @title geneAnno
#'
#' @description \code{geneAnno}Add genetic information to the tail after
#'   alignment.
#'
#' @details This function annotates the gene information, including the gene ID
#'   and gene type, using the tail of the matched gene from the GFF file.
#'
#' @param tailDF The tailScan dataframe.This parameter is required if the
#'   sequence type is shortreads, and omitted if the sequence type is longreads.
#' @param refPath The path of Reference conversion table.
#' @param bamdf The output dataframe of the function tailMap. Include at least
#'   read_num,chr,strand,coord.
#' @param GFF GFF file, can be GFF/GTF/GFF3, or TXDB comment package.It is
#'   recommended to use the TXDB annotation package, but be careful about the
#'   correspondence of the reference genome version.
#' @param longRead Boolean value.If your sequence type is longreads this
#'   parameter is T, otherwise it is F, and the default is T.
#'
#' @return Added tail table of gene annotation information.Include at least
#'   read_num,tail,PAL,chr,strand,gene,gene_type,tailType,read_type and sample.
#' @examples
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' load("./data/output/GV1tailDF.RData")
#' load("./data/output/GV1tailMapre.RData")
#' AnnotedTails <- geneAnno(tailDF = GV1tailDF, bamdf = GV1tailMapre, GFF = TxDb.Mmusculus.UCSC.mm10.knownGene, longRead = F)
#' @export
geneAnno <- function(tailDF, refPath, bamdf, GFF, longRead) {
  if (missing(GFF)) {
    stop("Parameter GFF is missing,Please check!")
  }
  if (missing(bamdf)) {
    stop("Parameter bamdf is missing,Please check!")
  }
  if (missing(longRead)) {
    longRead <- T
  }
  if (missing(refPath)) {
    refPath <- "./data/GRCH38_referenc_report.txt"
  }
  ## change the chr name format
  GFF <- parseGenomeAnnotation(GFF)
  if (str_detect(bamdf$chr[1], "^NC_")) {
    if (str_detect(GFF[["anno.need"]][["seqnames"]][1], "^chr")) {
      print("changing the format of chr name")
      bamdf <- ChangeChrFormat(bamdf, refPath = refPath, from = "RefSeq.Accn", to = "UCSC.style.name")
    } else {
      stop("Chr names of GFF and PAC not the same, please use the ChangeChrFormat() function first!")
    }
  }
  if (longRead) {
    if (missing(tailDF)) {
      tailDF <- " "
    }
    test3 <- annotatePAC(pac = bamdf, aGFF = GFF, verbose = T)
    test3 <- test3 %>%
      dplyr::filter(gene != "character(0)" & gene != " " & !is.na(gene))
    return(test3)
  } else {
    if (missing(tailDF)) {
      stop("Parameter tailDF is missing,Please check!")
    }
    # 3.Annotate gene
    # bamdf$chr <- str_replace(bamdf$chr,"ChrC","Pt")
    # bamdf$chr <- str_replace(bamdf$chr,"ChrM","Mt")
    # bamdf$chr <- str_replace(bamdf$chr,"Chr1","1")
    # bamdf$chr <- str_replace(bamdf$chr,"Chr2","2")
    # bamdf$chr <- str_replace(bamdf$chr,"Chr3","3")
    # bamdf$chr <- str_replace(bamdf$chr,"Chr4","4")
    # bamdf$chr <- str_replace(bamdf$chr,"Chr5","5")
    test3 <- annotatePAC(pac = bamdf, aGFF = GFF, verbose = T)
    test3 <- test3 %>%
      select(read_num, chr, strand, gene, gene_type)
    tailDF <- tailDF %>%
      select(-strand)
    if (dim(test3)[1] <= dim(tailDF)[1]) {
      test4 <- left_join(test3, tailDF, by = "read_num")
    } else {
      test4 <- left_join(tailDF, test3, by = "read_num")
    }
    test4 <- test4 %>%
      dplyr::filter(gene != "character(0)" & gene != " " & !is.na(gene))
    return(test4)
  }
}


# findPAs function
### This script is used to calculate the PA site.
### There are two ways to calculate PA sites: counting and not counting.
### The count method not only returns information about the PA sites but also
### counts how many reads have the same PA sites.
# findPAs -----------------------------------------------------------
#' @title findPAs
#'
#' @description Each chromosome was traversed and PA sites were identified by
#'   a-rich at the end. Note that the BAM file needs to be indexed.
#'   
#' @details The findPAs function takes the information from the BAM file and the
#'   CHRinfo file that you input, calculates the PA site for each read based on
#'   the A-rich in the sequence, and consolidates and saves the PA sites
#'   information to the result path that you specify.
#'   
#' @param chrinfo The path of chrinfo file.
#' @param bamfile The path of bam file.
#' @param resultpath The path of result file.
#' @param count Boolean parameter. If count=T, the PA list after the count will
#'   be returned. There will be one more column of count in the PA list. If
#'   count=F, all PA site information is returned without counting.
#' @return Save the calculated PA site results to the result path that you specify.
#' @examples
#' bamfilepath <- "./data/GV_algin/PAIso-GV1.sorted.bam"
#' chrinfopath <- "./data/GV_algin/chrinfo.txt"
#' resultpath <- "./data/output/PAs"
#' test <- findPAs(chrinfo = chrinfopath, bamfile = bamfilepath, resultpath = resultpath)
#' @export
findPAs <- function(chrinfo, bamfile, resultpath, count) {
  result <- data.frame()
  chr.g <- readChrInfo(chrinfo)
  what <- c("qname", "rname", "strand", "pos", "cigar", "seq")
  for (i in c(1:length(chr.g))) {
    print(chr.g[i])
    if (count) {
      resultson <- findChrTails(bamfile, which = chr.g[i], what, count = T)
    }
    if (!count) {
      resultson <- findChrTails(bamfile, which = chr.g[i], what, count = F)
    }
    result <- base::rbind(result, as.data.frame(resultson))
  }
  if (count) {
    resultpath <- str_c(resultpath, "/PAscount.txt")
  }
  if (!count) {
    resultpath <- str_c(resultpath, "/PAs.txt")
  }

  write.table(result, file = resultpath, quote = FALSE, row.names = F)

  return(result)
}

# readChrInfo function
# readChrInfo -------------------------------------------------------------
#' @title readChrInfo
#'
#' @description \code{readChrInfo} returns the S4 format data of all the chromosome information you input.
#'
#' @details  You need to input the TXT text of the HEAD information of the Sam file after
#' the alignment, and then the function will calculate and return the name of
#' the chromosome, the range on the reference genome, the direction of the
#' sequence, etc
#'
#' @param file The path of TXT text of the HEAD information of the Sam file.
#' @return The S4 format data of all the chromosome information.
#' @examples
#' chrinfopath <- "./data/GV_algin/chrinfo.txt"
#' test <- readChrInfo(chrinfopath)
#' @export
readChrInfo <- function(file) {
  chr_length <- read.delim(file = file)
  colnames(chr_length) <- c("flag", "chr", "coord")
  chr_length$chr <- gsub("SN:", "", chr_length$chr)
  chr_length$coord <- as.numeric(gsub("LN:", "", chr_length$coord))
  chr_length <- subset(chr_length, chr != "MT")
  chr.g <- with(chr_length, GenomicRanges::GRanges(seqnames = chr, ranges = IRanges(start = 1, end = coord)))
  return(chr.g)
}


# findChrTails function
# findChrTails ------------------------------------------------------------
#' @title findChrTails
#'
#' @description Each chromosome was traversed and PA sites were identified by
#'   a-rich at the end. Note that the BAM file needs to be indexed.
#'   
#' @param bamfile The path of bam file.
#' @param which Refer to ScanBamParam parameters for Rsamtools for explanation.
#' @param what Refer to ScanBamParam parameters for Rsamtools for explanation.
#' @param count Boolean parameter. If count=T, the PA list after the count will
#'   be returned. There will be one more column of count in the PA list. If
#'   count=F, all PA site information is returned without counting.
#' @return Returns a PA list annotated with chromosome information, chain
#'   orientation, and reference genomic coordinates.
#' @examples
#' bamfilepath <- "./data/GV_algin/PAIso-GV1.sorted.bam"
#' which <- chr.g[i]
#' what <- c("qname", "rname", "strand", "pos", "cigar", "seq")
#' PAlist <- findChrTails(bamfile = bamfilepath, which, what, count = T)
#' @export
findChrTails <- function(bamfile, which, what, count) {
  library(dplyr)
  param <- Rsamtools::ScanBamParam(what = what, which = which)
  gal1 <- GenomicAlignments::readGAlignments(bamfile, use.names = TRUE, param = param)
  s_1 <- (grepl("[0-9.*]+M[0-9.*]+S$", gal1@cigar) & as.vector(gal1@strand) == "+")
  s_2 <- (grepl("^[0-9.*]+S[0-9.*]+M", gal1@cigar) & as.vector(gal1@strand) == "-")

  bam1 <- gal1[s_1]
  bam2 <- gal1[s_2]

  bam1 <- bam1[grepl("A{10,}", bam1@elementMetadata@listData$seq)]
  bam2 <- bam2[grepl("T{10,}", bam2@elementMetadata@listData$seq)]

  final_bam1 <- data.frame(read_num = as.vector(names(bam1)), chr = as.vector(seqnames(bam1)), strand = as.vector(strand(bam1)), coord = end(bam1))
  final_bam2 <- data.frame(read_num = as.vector(names(bam2)), chr = as.vector(seqnames(bam2)), strand = as.vector(strand(bam2)), coord = start(bam2))

  bam <- base::rbind(final_bam1, final_bam2)
  bam <- dplyr::group_by(bam, read_num, chr, strand, coord)

  if (count) {
    bam <- dplyr::group_by(bam, chr, strand, coord) %>% dplyr::summarise(count = n())
  }

  return(bam)
}


# findAandAnnoPAs
# findAndAnnoPAs ------------------------------------------------------------
#' @title findAndAnnoPAs
#'
#' @description Find and annotate PA sites.
#' 
#' @param chrinfo The path of chrinfo file.
#' @param bamfile The path of bam file.
#' @param resultpath The path of result file.
#' @param bsgenome BSgenome package for your species, please refer to movAPA
#'   package for details.
#' @param gffFile Genome annotation stored in a GFF/GTF file or a TXDB R object
#'   can be used for annotating PACs. Please refer to movAPA for details.
#' @param sample The sample name.
#' @param mergePAs TRUE/FALSE. If TRUE, then will mergePACds groups nearby PACs
#'   from single/multiple PACdataset objects.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @return The function returns data in PACds format with the annotation
#'   information for all PA sites. In addition, we will also generate a
#'   "pac_data.txt" file to store the PA annotation table under the result path
#'   you specified, and a "pac_data.coldata.txt" file to store the additional
#'   comment information.And an "ACTGpdf.PDF" file, which holds a legend of the
#'   base composition around the PACs, has two images representing the base
#'   composition of the plus and minus chains.
#' @examples
#' # Deciphering gff files
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' # Deciphering genome files
#' BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10
#' gffFile <- TxDb.Mmusculus.UCSC.mm10.knownGene
#' bamfile <- system.file("extdata", "./GV_algin/GV1subseq.sorted.bam", package = "PolyAtailor", mustWork = TRUE)
#' chrinfo <- system.file("extdata", "./GV_algin/chrinfo.txt", package = "PolyAtailor", mustWork = TRUE)
#' resultpath <- "./data/output/PAs"
#' # Annotated PA site
#' PAs <- findAndAnnoPAs(bamfile = bamfilepath, chrinfo = chrinfopath, resultpath = resultpath, bsgenome = bsgenome, gffFile = gffFile, sample = "GV1", mergePAs = T, d = 24)
#' @export
findAndAnnoPAs <- function(chrinfo, bamfile, resultpath, bsgenome, gffFile, sample, mergePAs, d) {
  # 1.Determine the PA site
  print("findind PAs...")
  PAtab <- findPAs(bamfile = bamfile, chrinfo = chrinfo, resultpath = resultpath, count = F)
  print("findind PAs finished")
  PAtab$sample <- sample
  gff <- parseGenomeAnnotation(gffFile)
  ## unify the chr names
  if (str_detect(PAtab$chr[1], "^NC_") || str_detect(PAtab$chr[1], "^NT_")) {
    if (str_detect(gff[["anno.need"]][["seqnames"]][1], "^chr")) {
      print("changing the format of chr name")
      PAtab <- ChangeChrFormat(PAtab, from = "RefSeq.Accn", to = "UCSC.style.name")
    } else {
      stop("Chr names of gff and PAC not the same, please use the ChangeChrFormat() function first!")
    }
  }
  # 2.Parse annotation file
  PACds <- movAPA::readPACds(PAtab)
  # 3.Preprocessing of PAC data
  # 1)Remove internal priming artifacts
  PACdsIP <- removePACdsIP(PACds, bsgenome, chrCheck = F, returnBoth = TRUE, up = -10, dn = 10, conA = 6, sepA = 7)
  PACds <- PACdsIP$real
  # 2)Group nearby cleavage sites
  if (missing(mergePAs)) {
    print("The parameter mergePAs will use the default value TRUE!")
  }
  if (mergePAs) {
    if (missing(d)) {
      print("The d parameter will use the default value 24!")
      PACds <- mergePACds(PACds, d = 24)
    } else {
      PACds <- mergePACds(PACds, d = d)
    }
  }

  # 3)Normalization with "TMM"
  PACds <- normalizePACds(PACds, method = "TMM")
  # # # 4)ext3UTR
  # ext3UTRPACds(PACds, ext3UTRlen=1000)
  # 4.Annotate PACs
  # #统一chrname
  # PACds@anno[["chr"]] <- str_remove(PACds@anno[["chr"]],"chr")
  PACds <- annotatePAC(pac = PACds, aGFF = gff, verbose = T)
  # 5.ext3UTR
  PACds <- ext3UTRPACds(PACds, ext3UTRlen = 1000)
  pac_data_path <- str_c(resultpath, "pac_data.txt", sep = "")
  pac_coldata_path <- str_c(resultpath, "pac_data.coldata.txt", sep = "")
  writePACds(PACds, file = pac_data_path, colDataFile = pac_coldata_path)
  # 6.Base compostions and k-grams
  resultpdf <- str_c(resultpath, "ACTGpdf.pdf", sep = "")
  pdf(resultpdf, width = 15, height = 15, onefile = T)
  # #再次修改chrname
  # PACds@anno[["chr"]]<-paste("chr",PACds@anno[["chr"]],sep="")
  faFiles0 <- faFromPACds(PACds, bsgenome,
    what = "updn", fapre = "updn",
    up = -300, dn = 100, byGrp = c("ftr")
  )



  # faFiles4=c("updn.3UTR.fa", "updn.cds.fa", "updn.intergenic.fa", "updn.intron.fa")
  p0 <- plotATCGforFAfile(faFiles0, ofreq = FALSE, opdf = TRUE, refPos = 301, mergePlots = TRUE, filepre = "ATCG_all")
  # topptx(filename =str_c(resultpath,"plot.pptx"))
  faFiles <- faFromPACds(PACds, bsgenome,
    what = "updn", fapre = "updn",
    up = -300, dn = 100, byGrp = c("ftr", "strand")
  )

  temp.name <- names(table(paste0(PACds@anno$ftr, PACds@anno$strand)))
  temp.name <- paste0("updn.", temp.name, ".fa")
  faFiles <- faFiles[which(faFiles %in% temp.name)]
  ## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
  # faFiles2=c("updn.3UTR+.fa", "updn.cds+.fa", "updn.intergenic+.fa", "updn.intron+.fa")
  faFiles2 <- faFiles[grep("\\+", faFiles)]
  p1 <- plotATCGforFAfile(faFiles2, ofreq = FALSE, opdf = TRUE, refPos = 301, mergePlots = TRUE, filepre = "ATCG_Forward")
  ## Plot single nucleotide profiles using the extracted sequences and merge all plots into one.
  faFiles3 <- faFiles[grep("\\-", faFiles)]
  p2 <- plotATCGforFAfile(faFiles3, ofreq = FALSE, opdf = TRUE, refPos = 301, mergePlots = TRUE, filepre = "ATCG_reverse")
  # print(p0)
  # print(p1)
  # print(p2)
  # dev.off()
  # dev.off()
  # dev.off()
  return(PACds)
}

# Polyatailor contains abundant analysis and visualization tools of poly(A) tail
# base composition, including statistics of the number distribution of reads in
# different non-A base combinations, analysis of non-A base content in polyAtail
# of different lengths, and direct visual comparative analysis of intrested
# tails.
### The following functions are used for non-A base analysis
# ccf --------------------------------------------------
#' @title ccf
#'
#' @description Convert chromosome name format.
#' 
#' @details For various reasons, we sometimes need to convert the format of
#'   chromosome names in files. This function can help users convert the format.
#'   
#' @param ref Reference table of chromosome information from NCBI database.
#' @param PAdf Files that need to convert chromosome names, dataframe format.
#' @param fromFormat Primitive chromosome nomenclature.
#' @param toFormat Target format of chromosome name.
#' @return This function returns a data frame with the chromosome name changed.
#' @examples
#' ref <- read.table("./data/HELA1/GRCH38p13_sequences_report.txt",
#'   header = T, sep = "\t", blank.lines.skip = T
#' )
#' PAdf <- read.table("***")
#' result <- ccf(ref, PAdf, "RefSeq.Accn", "UCSC.style.name")
#' @export
ccf <- function(ref, PAdf, fromFormat, toFormat) {
  ref <- ref %>%
    select(fromFormat, toFormat) %>%
    rename(chr = fromFormat)
  PAdf <- left_join(PAdf, ref)
  PAdf <- PAdf %>%
    select(-chr) %>%
    rename(chr = toFormat)
  PAdf <- na.omit(PAdf)
  return(PAdf)
}

# Read a PACdataset -------------------------------------------------------
#' @title readPACds
#'
#' @description readPACds reads PAC counts and sample annotation into a PACdataset.
#'
#' @usage readPACds(pacFile, colDataFile, noIntergenic=TRUE, PAname='PA')
#' @param pacFile a file name or a data frame. If it is a file, it should have header, with at least (chr, strand, coord) columns.
#' This file could have other columns, including gff cols (gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end) and user-defined sample columns.
#' If there are at least one non-numeric columns other than above gff cols, then all remaining columns are considered as annotation columns.
#' If all remaining columns are numeric, then they are all treated as sample columns.
#' Use annotatePAC() first if need genome annotation of coordinates.
#' @param colDataFile a file name or a data frame. If it is a file, then it is an annotation file of samples with header,
#' rownames are samples (must be all in pacFile), columns names are sample groups.
#' There could be single or multiple columns to define the groups of samples.
#' When colDataFile=NULL, then readPACds will automately retreive sample columns and gff columns (if any) from pacFile.
#' If there is no sample columns, then will set colData as a data frame with 1 column (=group) and 1 row (=tag), and element=group1.
#' If pacfile or colDataFile is a character, then it is a file name, so readPACds will read data from file.
#' @param noIntergenic TRUE/FALSE. If TRUE, then will remove PACs in intergenic (ftr='^inter')
#' @param PAname specify how to set the name (rowname) of PACs.
#' PAname=PA, the PA name is set as 'gene:PAN'; PAname=coord, then 'gene:coord'.
#' @return A PACdataset object, with @anno being a data frame with at least three columns chr/strand/coord. If there is no sample column, then will add one sample named tag in group1.
#' @examples
#' data(PACds)
#' ## read simple PACfile that only has columns chr/strand/coord
#' pacFile <- PACds@anno[, c("chr", "strand", "coord")]
#' colDataFile <- NULL
#' p <- readPACds(pacFile, colDataFile)
#'
#' ## read PACfile that has columns chr/strand/coord and sample columns
#' pacFile <- PACds@anno[, c("chr", "strand", "coord")]
#' pacFile <- cbind(pacFile, PACds@counts[, c("anther1", "embryo1", "anther2")])
#' colDataFile <- NULL
#' p <- readPACds(pacFile, colDataFile)
#' p@colData
#' head(p@counts)
#'
#' ## read PACfile that has columns chr/strand/coord, sample columns, and gff cols like gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end
#' pacFile <- PACds@anno
#' pacFile <- cbind(pacFile, PACds@counts[, c("anther1", "embryo1", "anther2")])
#' colDataFile <- NULL
#' p <- readPACds(pacFile, colDataFile)
#' p@colData
#' head(p@counts)
#' head(p@anno)
#'
#' ##read from data frame of PACfile and colDataFile
#' pacFile <- PACds@anno
#' smps <- c("anther1", "embryo1", "anther2")
#' pacFile <- cbind(pacFile, PACds@counts[, smps])
#' colDataFile <- as.data.frame(matrix(c("group1", "group2", "group1"), ncol = 1, dimnames = list(smps, "group")))
#' p <- readPACds(pacFile, colDataFile)
#' p@colData
#' head(p@counts)
#' head(p@anno)
#' ## read from file names of PACfile and colDataFile
#' write.table(pacFile, file = "pacFile", row.names = FALSE)
#' write.table(colDataFile, file = "colDataFile", row.names = TRUE)
#' p <- readPACds(
#'   pacFile = "pacFile",
#'   colDataFile = "colDataFile", noIntergenic = TRUE, PAname = "PA"
#' )
#'
#' @name readPACds
#' @family PACdataset functions
#' @export
readPACds <- function(pacFile, colDataFile = NULL, noIntergenic = TRUE, PAname = "PA") {
  if (!(PAname %in% c("PA", "coord"))) {
    stop("PAname must be PA or coord")
  }

  if (is.character(pacFile)) {
    d <- read.table(pacFile, header = T)
  } else {
    d <- pacFile
  }

  gffcols <- c("gene", "gene_type", "ftr", "ftr_start", "ftr_end", "UPA_start", "UPA_end")

  if (sum(!(c("chr", "strand", "coord") %in% colnames(d))) != 0) {
    stop("chr,strand,coord not all in header of pacfile")
  }

  cat(nrow(d), "PACs\n")

  if ("ftr" %in% colnames(d)) {
    if (noIntergenic) {
      d <- d[getNonItgFtrId(d$ftr), ]
      cat(nrow(d), "No intergenic PACs\n")
    }
    # WB's data 20190620
    d$ftr[d$ftr == "three_prime_UTR"] <- "3UTR"
    d$ftr[d$ftr == "five_prime_UTR"] <- "5UTR"
  }

  # d[d=='unkown']=NA

  if ("ftr_start" %in% colnames(d)) {
    if (!is.numeric(d$ftr_start)) {
      d$ftr_start <- as.numeric(d$ftr_start)
    }
    if (!is.numeric(d$ftr_end)) {
      d$ftr_end <- as.numeric(d$ftr_end)
    }
  }

  if ("gene" %in% colnames(d)) {
    idx <- which(d$gene == "NULL")
    if (length(idx) > 0) {
      cat(length(idx), "gene name is NULL, change to chrStrand")
      d$gene[idx] <- paste0(d$chr[idx], d$strand[idx])
    }

    # order 5' to 3'
    d1 <- d[d$strand == "+", ]
    d1 <- d1[order(d1$gene, d1$coord, decreasing = FALSE), ]
    d2 <- d[d$strand == "-", ]
    d2 <- d2[order(d2$gene, d2$coord, decreasing = TRUE), ]
    d <- rbind(d1, d2)

    if (PAname == "coord") {
      paid <- paste0(d$gene, ":", d$coord)
    } else if (PAname == "PA") {
      rg <- rle(d$gene)
      paid <- paste0(d$gene, ":PA", unlist(lapply(rg$lengths, seq)))
    }
  } else {
    paid <- paste0("PA", 1:nrow(d))
  }

  # something new
  # rownames(d)=paid
  rownames(d) <- make.names(paid, unique = TRUE)

  allcols <- colnames(d)


  # group defination of sample columns
  if (!is.null(colDataFile)) {
    if (is.vector(colDataFile)) {
      colData <- read.table(colDataFile, colClasses = "character")
    } else {
      colData <- colDataFile
    }

    if (sum(rownames(colData) %in% colnames(d)) != nrow(colData)) {
      stop("rownames of annofile not all in columns of pacfile")
    }
  } else {
    # remove gffcolid and chr/strand/coord, other columns are sample columns, and the sample group is 'group', groupname is 'group1'
    smpcols <- which(allcols %in% c(gffcols, "chr", "strand", "coord"))
    smpcols <- allcols[-smpcols]
    # no tag columns, then add tag=1 columns
    if (length(smpcols) == 0) {
      smpcols <- "tag"
      d$tag <- 1
    } else { # one columns is chr, then they are all annotations
      for (i in smpcols) {
        if (!(is.numeric(d[, i]))) {
          smpcols <- "tag"
          d$tag <- 1
          break
        }
      }
    }
    colData <- as.data.frame(matrix(rep("group1", length(smpcols)), ncol = 1, dimnames = list(smpcols, "group")))
  }

  for (i in 1:ncol(colData)) {
    colData[, i] <- factor(colData[, i])
  }
  colData <- colData[order(rownames(colData)), , drop = F]
  anno <- d[, -which(colnames(d) %in% rownames(colData))]
  anno[anno == "unkown"] <- NA

  counts <- d[, rownames(colData), drop = F]
  PACds <- new("PACdataset", counts = counts, colData = colData, anno = anno)
  return(PACds)
}

# Merge multiple PACdatasets ----------------------------------------------
#' @title mergePACds
#'
#' @description  mergePACds groups nearby PACs from single/multiple PACdataset objects.
#'
#' @details  This function is particularlly useful for grouping nearby cleavage sites into PACs.
#' It is also useful When you have multiple PA or PAC files, each file is from one sample.
#' Then you need to merge these PACds into one PACds for DE or other analyses.
#' But after grouping and/or merging, you may need call annotatePAC to annotate the merged PACs by a GFF annotation.
#' 
#' @usage mergePACds(PACdsList, d=24)
#' @param PACdsList a PACdataset, or a list of multiple PACdataset objects. The PACds@anno should have columns chr/strand/coord.
#' If there is no colData in PACds, then the sample label will be set as groupN.
#' If PACdsList is a PACdataset, then will treat it as PA and group nearby PAs into PACs.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @return A merged PACdataset. The counts slot stores counts of merged samples.
#' If sample names from different PACdataset objects are duplicated, then will be set as .x,.y.
#' The colData slot stores the merge sample annotation from the first column of each @colData.
#' The anno slot contains these columns: chr, strand, coord, tottag, UPA_start, UPA_end, nPA, maxtag.
#' @examples
#' ## Group PA into PACs
#' data(PACds)
#' PACds@counts <- rbind(PACds@counts, PACds@counts)
#' PACds@anno <- rbind(PACds@anno, PACds@anno)
#' ds <- mergePACds(PACds, d = 24)
#' ## merge two PACds
#' ds1 <- makeExamplePACds()
#' ds2 <- makeExamplePACds()
#' ds <- mergePACds(list(ds1, ds2), d = 24)
#' @name mergePACds
#' @family PACdataset functions
#' @seealso [annotatePAC()] to annotate a PACdataset; [rbind()] to combine multiple PACdatasets of the same format.
#' @export
mergePACds <- function(PACdsList, d = 24) {

  # library(GenomicRanges, verbose = FALSE)
  # library(dplyr, verbose = FALSE)

  if (class(PACdsList) == "PACdataset") PACdsList <- list(PACdsList)

  for (i in 1:length(PACdsList)) {
    if (!AinB(c("chr", "strand", "coord"), colnames(PACdsList[[i]]@anno))) stop("chr/strand/coord not in PACdsList@anno")
    if (nrow(PACdsList[[i]]@colData) == 0) {
      PACdsList[[i]]@colData <- as.data.frame(matrix(rep(
        paste0("group", i),
        ncol(PACdsList[[i]]@counts)
      ),
      ncol = 1,
      dimnames = list(colnames(PACdsList[[i]]@counts), "group")
      ))
    }
  }

  allpa <- PACds2PAdf(PACdsList[[1]])
  # return(allpa)
  if (length(PACdsList) >= 2) {
    for (j in 2:length(PACdsList)) {
      pa2 <- PACds2PAdf(PACdsList[[j]])
      allpa <- merge(allpa, pa2, all = T, by.x = c("chr", "strand", "coord"), by.y = c("chr", "strand", "coord"))
    }
  }
  allpa[is.na(allpa)] <- 0
  invisible(gc())


  gr <- with(allpa, GRanges(
    seqnames = chr,
    ranges = IRanges(
      start = coord,
      end = coord
    ),
    strand = strand
  ))

  # resize width as d+1
  gr <- resize(gr, width = d + 1, fix = "start", use.names = TRUE, ignore.strand = FALSE)

  itv <- IRanges::reduce(gr, drop.empty.ranges = TRUE)

  cat("Group PA to PACs\n")
  ov <- findOverlaps(gr, itv,
    maxgap = -1L, minoverlap = 1L,
    type = c("any"), select = "all",
    ignore.strand = FALSE
  )
  ov <- as.data.frame(ov)
  allpa$idx <- 1:nrow(allpa)
  allpa <- merge(allpa, ov, by.x = "idx", by.y = "queryHits")

  allpa$idx <- NULL
  smpcols <- colnames(allpa)[!(colnames(allpa) %in% c("chr", "strand", "coord", "subjectHits"))]

  # sum tag per interval
  cat("count tot tag for each sample within each PAC\n")
  allpa$tottag <- rowSums(allpa[, smpcols, drop = F])
  byItv <- dplyr::group_by(allpa, subjectHits)
  dots <- sapply(smpcols, function(x) substitute(sum(x), list(x = as.name(x))))
  dots[["tottag"]] <- substitute(sum(x), list(x = as.name("tottag")))
  pac <- do.call(dplyr::summarise, c(list(.data = byItv), dots))

  cat("Annotate the range of each PAC\n")
  # get interval range...
  pacAnno <- byItv %>% dplyr::summarise(nPA = n(), UPA_start = min(coord), UPA_end = max(coord), maxtag = max(tottag), coord = coord[which.max(tottag)], chr = chr[1], strand = strand[1])

  pac <- merge(pac, pacAnno, by.x = "subjectHits", by.y = "subjectHits")

  # pacds
  smpcols <- smpcols[smpcols != "tottag"]
  counts <- pac[, smpcols, drop = F]
  anno <- pac[, c("chr", "strand", "coord", "tottag", "UPA_start", "UPA_end", "nPA", "maxtag")]

  colData <- as.data.frame(matrix(unlist(lapply(PACdsList, function(ds) {
      return(as.character(ds@colData[, 1]))
    })),
    ncol = 1, dimnames = list(colnames(counts), "group")
  ))
  d <- new("PACdataset", counts = counts, colData = colData, anno = anno)
  return(d)
}


# get3UTRAPAds ------------------------------------------------
#' @title get3UTRAPAds
#'
#' @description  get3UTRAPAds subset a PACdataset to get all PACs of genes with 3'UTR APA sites.
#'
#'  @details  Genes with 3'UTR APA sites are genes with multiple PACs in its 3'UTR.
#' If the column toStop is not in PACds@anno, then will calculate the 3'UTR length first.
#' First filter by paramters like avgPACtag, etc., and then by choose2PA.
#'
#' @usage get3UTRAPAds(pacds, sortPA=TRUE, choose2PA=NULL, avgPACtag=0, avgGeneTag=0, clearPAT=0)
#' @param sort If TRUE, then order PACs by the respective 3'UTR length in each gene.
#' @param avgPACtag if >0, then to filter by PAC tag num, >=avgPACtag, see subsetPACds().
#' @param avgGeneTag if >0, then to filter by PAC tag num, >=avgGeneTag, see subsetPACds().
#' @param clearPAT if >0, then to filter by clearPAT, see subsetPACds().
#' @param choose2PA specify whether and how to choose two PACs when there are >2 PACs. The value can be NULL (use all PACs), PD (choose only proximal and distal sites),
#' farest (choose two PACs that are with the longest distance), most (choose two PACs with the most abundance).
#' @return A PACdataset with only 3'UTR APA sites. If there is no result, return an empty PACdataset with 0-row @anno and @counts.
#' @examples
#' data(PACds)
#' ## Get all 3'UTR APA
#' pp1 <- get3UTRAPAds(PACds, sortPA = TRUE)
#' summary(pp1)
#' is3UTRAPAds(pp1)
#' is3UTRAPAds(PACds)
#' ## Get proximal distal PA
#' pp1 <- get3UTRAPAds(PACds, sortPA = TRUE, choose2PA = "most")
#' summary(pp1)
#' is3UTRAPAds(pp1)
#' is3UTRAPAds(PACds)
#' head(pp1@anno)
#' pp1@anno[pp1@anno$gene == "PH02Gene00018", ]
#' @name get3UTRAPAds
#' @family PACdataset functions
#' @export
get3UTRAPAds <- function(pacds, sortPA = TRUE, choose2PA = NULL, avgPACtag = 0, avgGeneTag = 0, clearPAT = 0) {
  if (!is.null(choose2PA)) {
    sortPA <- TRUE
    choose2PA <- toupper(choose2PA)
    if (!(choose2PA %in% c("PD", "MOST"))) stop("choose2PA must be PD or MOST\n")
  }

  # filter by tagnum
  pacds <- subsetPACds(pacds, avgPACtag = avgPACtag, avgGeneTag = avgGeneTag, noIntergenic = TRUE, clearPAT = clearPAT, verbose = FALSE)

  # filter 3UTR APA
  if (!is3UTRAPAds(pacds)) {
    pacds@anno <- pacds@anno[which(pacds@anno$ftr == "3UTR"), ]
    genes <- unique(pacds@anno$gene[duplicated(pacds@anno$gene)])
    if (length(genes) == 0) { # no 3utr APA
      pacds@anno <- pacds@anno[character(0), ]
      pacds@counts <- pacds@counts[character(0), ]
      return(pacds)
    }
    pacds@anno <- pacds@anno[pacds@anno$gene %in% genes, ]
    if (!("toStop" %in% colnames(pacds@anno))) {
      if ("three_UTR_length" %in% colnames(pacds@anno)) {
        pacds@anno$toStop <- pacds@anno$three_UTR_length
      } else {
        pacds@anno$toStop <- pacds@anno$coord
        id <- grep("\\+", pacds@anno$strand)
        pacds@anno$toStop[id] <- pacds@anno$coord[id] - pacds@anno$ftr_start[id] + 1
        id <- grep("\\-", pacds@anno$strand)
        pacds@anno$toStop[id] <- pacds@anno$ftr_end[id] - pacds@anno$coord[id] + 1
      }
    }
    if (min(pacds@anno$toStop) < 0) stop("get3UTRAPAds: error toStop (<0), please check coord/ftr_start/ftr_end\n")
    data <- pacds@counts[rownames(pacds@anno), ]
    pacds@counts <- subset(pacds@counts, rownames(pacds@counts) %in% rownames(pacds@anno))
    pacds@counts$tag <- data
  }

  # sort by 3UTR len
  if (sortPA) {
    pacds <- pacds[order(pacds@anno$gene, pacds@anno$toStop, decreasing = FALSE)]
  }

  if (!is.null(choose2PA)) {
    rs <- rowSums(pacds@counts)
    d <- pacds@anno[, c("gene", "toStop")]
    d$tag <- rs
    d$PA <- rownames(pacds@anno)

    .choose2 <- function(x, method) {
      if (nrow(x) == 2) {
        return(x$PA)
      }
      if (method == "PD") {
        return(x$PA[c(1, nrow(x))])
      } else if (method == "MOST") {
        x <- x[order(x$tag, decreasing = TRUE), ]
        return(x$PA[c(1:2)])
      }
    }

    sp <- split(d, d$gene, drop = T)
    PAs <- unlist(lapply(sp, .choose2, method = choose2PA))
    # x=matrix(unlist(PAs), ncol=2, byrow = TRUE)
    # PAs=as.data.frame(cbind(gene=names(PAs), x))

    # PAs=plyr::ddply(d,.(gene), .choose2, method=choose2PA)
    # PAs=unlist(PAs[,2:3])

    pacds@anno <- pacds@anno[PAs, ]
    data <- pacds@counts[PAs, ]
    pacds@counts <- subset(pacds@counts, rownames(pacds@counts) %in% rownames(pacds@anno))
    pacds@counts$tag <- data
    if (sortPA) {
      pacds <- pacds[order(pacds@anno$gene, pacds@anno$toStop, decreasing = FALSE)]
    }
  }
  return(pacds)
}



# function PACds2PAdf
PACds2PAdf <- function(PACds) {
  if (!AinB(c("chr", "strand", "coord"), colnames(PACds@anno))) stop("chr/strand/coord not in PACds@anno")
  df <- cbind(PACds@anno[, c("chr", "strand", "coord")], PACds@counts)
  return(df)
}

# internal function -------------------------------------------------------
# function1 for apply
PAtag <- function(x, y, d) {
  coord <- x["coord"]
  gene1 <- x["gene"]
  z <- y %>%
    filter(gene == gene1)
  middle <- (z[1, ]$coord + z[2, ]$coord) / 2
  if (coord >= middle) {
    return("PA1")
  } else {
    return("PA2")
  }
}
# function for DSA
DSA <- function(x, y) {
  library(DescTools)
  genex <- x["gene"]
  PA1_PALs <- y %>%
    filter(gene == genex & PAID == "PA1") %>%
    select(PAL)
  PA1_PALs <- as.integer(PA1_PALs$PAL)
  PA2_PALs <- y %>%
    filter(gene == genex & PAID == "PA2") %>%
    select(PAL)
  PA2_PALs <- as.integer(PA2_PALs$PAL)
  # ks-test
  ks_re <- ks.test(PA1_PALs, PA2_PALs)$p.value
  ks_re <- round(ks_re, 5)
  # wilcox-test
  wilcox_re <- wilcox.test(PA1_PALs, PA2_PALs)$p.value
  wilcox_re <- round(wilcox_re, 5)
  # moses-test
  if (length(PA1_PALs) >= 2) {
    moses_re <- MosesTest(PA1_PALs, PA2_PALs)$p.value
    moses_re <- round(moses_re, 5)
  } else {
    moses_re <- "lose"
  }
  # Mann-Whitney U test
  MU_re <- wilcox.test(PA1_PALs, PA2_PALs, alternative = "greater", exact = F)$p.value
  MU_re <- round(MU_re, 5)
  res <- c(ks_re, wilcox_re, moses_re, MU_re)
  res[is.na(res)] <- "lose"
  res[res == "NaN"] <- "lose"
  res <- str_c(res[1], res[2], res[3], res[4], sep = "_")
  return(res)
}
myFun2 <- function(x, y) {
  library(DescTools)
  genex <- x["gene"]
  PA1_PALs <- y %>%
    filter(gene == genex & PAID == "distal") %>%
    select(PAL)
  PA1_PALs <- as.integer(PA1_PALs$PAL)
  PA2_PALs <- y %>%
    filter(gene == genex & PAID == "proximal") %>%
    select(PAL)
  PA2_PALs <- as.integer(PA2_PALs$PAL)
  # ks-test
  ks_re <- ks.test(PA1_PALs, PA2_PALs)$p.value
  ks_re <- round(ks_re, 5)
  # wilcox-test
  wilcox_re <- wilcox.test(PA1_PALs, PA2_PALs)$p.value
  wilcox_re <- round(wilcox_re, 5)
  # moses-test
  if (length(PA1_PALs) >= 2) {
    moses_re <- MosesTest(PA1_PALs, PA2_PALs)$p.value
    moses_re <- round(moses_re, 5)
  } else {
    moses_re <- "lose"
  }
  # Mann-Whitney U test
  MU_re <- wilcox.test(PA1_PALs, PA2_PALs, alternative = "greater", exact = F)$p.value
  MU_re <- round(MU_re, 5)
  res <- c(ks_re, wilcox_re, moses_re, MU_re)
  res[is.na(res)] <- "lose"
  res[res == "NaN"] <- "lose"
  res <- str_c(res[1], res[2], res[3], res[4], sep = "_")
  return(res)
}
# function for tag dpPAs
dpPAtag <- function(x, y) {
  coord <- x["coord"]
  gene1 <- x["gene"]
  z <- y %>%
    filter(gene == gene1)
  if (!is.na(z[1, ]$coord)) {
    if (z[1, ]$coord == "+") {
      middle <- (z[1, ]$coord + z[2, ]$coord) / 2
      if (coord >= middle) {
        return("distal")
      } else {
        return("proximal")
      }
    } else {
      middle <- (z[1, ]$coord + z[2, ]$coord) / 2
      if (coord <= middle) {
        return("distal")
      } else {
        return("proximal")
      }
    }
  } else {
    return("lose")
  }
}

# diffPAL2PAgene -----------------------------------------------------------------
#' @title diffPAL2PAgene
#'
#' @description Analysis of significant difference in tail length of genes with
#'   two PAs.
#'   
#' @details This function was used to analyze the data for genes with two PA
#'   sites showing significant differences in tail length between the two PA
#'   sites.The function will use four difference significance test methods to
#'   find all genes with significant difference in tail length in the provided
#'   data and make a UpSet graph.
#'   
#' @param PAdf A dataframe of PA infomation.
#' @param PALdf A dataframe that contains the length information of all reads
#'   tails.
#' @param gff Genome annotation stored in a GFF/GTF file or a TXDB R object can
#'   be used for annotating PACs. Please refer to movAPA for details.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @return The function returns a dataframe of result. The result data frame
#'   contains 7 column contents, the first column is geneID, and the second and
#'   third columns are the median tail lengths under different conditions. The
#'   remaining four columns are p-values calculated by different difference
#'   significance test methods, each column corresponds to a difference
#'   significance test method, in which "lose" or "NAN" means that for some
#'   reason (probably too little data) this method has not been able to
#'   calculate the exact p-values.
#' @examples
#' gff <- parseGenomeAnnotation(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' diffPAL2PAgenes <- diffPALgene(PAdf1, al1, gff)
#' @export
diffPAL2PAgene <- function(PAdf, PALdf, gff, d) {
  # check the paraments
  if (missing(PAdf) | missing(PALdf)) {
    stop("Missing important parameters: PAdf or PALdf!")
  }
  if (missing(gff)) {
    stop("Missing the annotations file!please check!")
  }
  if (missing(d)) {
    d <- 24
  }
  if (!is.integer(d) & !is.numeric(d)) {
    stop("Parameter d must be an integer!")
  }
  # read data
  PAdf <- na.omit(PAdf)
  PACds <- readPACds(pacFile = PAdf, colDataFile = NULL)
  PACds@counts <- rbind(PACds@counts, PACds@counts)
  PACds@anno <- rbind(PACds@anno, PACds@anno)
  # merge PAs to PACs
  ds <- mergePACds(PACds, d = d)
  # Annotate genetic information
  ds <- annotatePAC(ds, gff)
  ds <- cbind(ds@anno, ds@counts)
  # filter the genes having 2 PAs
  ds1 <- ds %>%
    select(gene, nPA) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(PAcounts = n()) %>%
    dplyr::filter(PAcounts == 2 & gene != "character(0)") %>%
    select(gene)
  ## cobvert to geneSymbol
  # if(str_detect(ds1$gene[1],"^[0-9]*$") || str_detect(ds1$gene[1],"^ENSMUSG[0-9]*")){
  #   library(annotate)
  #   library(org.Hs.eg.db)
  #   ds1$symbol <- getSYMBOL(ds1$gene, data='org.Hs.eg.db')
  # }
  genes_2_PAS <- ds1$gene
  # Annotate raw data genetic information
  PAdf2 <- annotatePAC(PAdf, gff)
  PAdf2 <- PAdf2 %>%
    filter(!is.na(gene))
  if (str_detect(PAdf2$read_num[1], "_")) {
    PAs <- PAdf2 %>%
      tidyr::separate(read_num, c("read_num", "umi", "structure_type"), sep = "_") %>%
      dplyr::filter(gene %in% genes_2_PAS) %>%
      dplyr::select(read_num, gene, coord)
  } else {
    PAs <- PAdf2 %>%
      dplyr::filter(gene %in% genes_2_PAS) %>%
      dplyr::select(read_num, gene, coord)
  }
  PACs <- ds %>%
    dplyr::filter(gene %in% genes_2_PAS) %>%
    dplyr::select(gene, UPA_start, UPA_end, coord) %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(PAID = c("PA1", "PA2")) %>%
    dplyr::ungroup()
  PAs$PAID <- apply(PAs, 1, PAtag, PACs, d)
  PAs <- PAs %>%
    select(-coord)
  # PALdf
  PALdf <- PALdf %>%
    filter(gene %in% genes_2_PAS) %>%
    select(read_num, PAL)
  re <- left_join(PAs, PALdf)
  re <- na.omit(re)
  re <- re %>%
    group_by(gene, PAID) %>%
    dplyr::summarise(counts = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(counts = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(counts == 2)
  genes_down <- re$gene
  re <- PAs %>%
    left_join(PALdf) %>%
    na.omit() %>%
    filter(gene %in% genes_down)
  result <- re %>%
    group_by(gene, PAID) %>%
    dplyr::summarise(midPAL = median(PAL)) %>%
    tidyr::pivot_wider(names_from = PAID, values_from = midPAL) %>%
    dplyr::rename(PA1_medianPAL = PA1, PA2_medianPAL = PA2)
  # plot upset
  result$DSAre <- apply(result, 1, DSA, re)
  result <- tidyr::separate(result, DSAre, c("KS_re", "wilcox_re", "moses_re", "MU_re"), sep = "_")
  return(result)
}

# diffPALdpPA-------------------------------------------------------------------
#' @title diffPALdpPA
#'
#' @description Analyzing whether there is a significant difference
#'   in polyA tail length between the proximal and distal PA sites of genes with
#'   two PACs.
#'   
#' @details This function analyzes whether there is a significant difference
#'   in polyA tail length between the proximal and distal PA sites of genes with
#'   two PACs.
#' @param PAdf A dataframe of PA infomation.
#' @param PALdf A dataframe that contains the length information of all reads
#'   tails.
#' @param gff Genome annotation stored in a GFF/GTF file or a TXDB R object can
#'   be used for annotating PACs. Please refer to movAPA for details.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @return The function returns a dataframe of result. The result data frame
#'   contains 7 column contents, the first column is geneID, and the second and
#'   third columns are the median tail lengths under different conditions. The
#'   remaining four columns are p-values calculated by different difference
#'   significance test methods, each column corresponds to a difference
#'   significance test method, in which "lose" or "NAN" means that for some
#'   reason (probably too little data) this method has not been able to
#'   calculate the exact p-values.
#' @examples
#' gff <- parseGenomeAnnotation(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' DSA_result <- diffPALdpPA(PAdf1, al1, gff)
#' @export
diffPALdpPA <- function(PAdf, PALdf, gff, d) {
  # check the paraments
  if (missing(PAdf) | missing(PALdf)) {
    stop("Missing important parameters: PAdf or PALdf!")
  }
  if (missing(gff)) {
    stop("Missing the annotations file!please check!")
  }
  if (missing(d)) {
    d <- 24
  }
  if (!is.integer(d) & !is.numeric(d)) {
    stop("Parameter d must be an integer!")
  }
  # input
  PACds <- readPACds(pacFile = PAdf, colDataFile = NULL)
  PACds@counts <- rbind(PACds@counts, PACds@counts)
  PACds@anno <- rbind(PACds@anno, PACds@anno)
  ds <- mergePACds(PACds, d = d)
  ds <- annotatePAC(ds, gff)
  ds <- get3UTRAPAds(ds, sortPA = TRUE, choose2PA = "PD")
  ds <- cbind(ds@anno, ds@counts)
  ds1 <- ds %>%
    dplyr::select(gene, nPA) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(PAcounts = n()) %>%
    dplyr::filter(PAcounts == 2 & gene != "character(0)") %>%
    dplyr::select(gene)
  # if(str_detect(ds1$gene[1],"^[0-9]*$") || str_detect(ds1$gene[1],"^ENSMUSG[0-9]*")){
  #   library(annotate)
  #   #library(org.Mm.eg.db)
  #   ds1$symbol <- getSYMBOL(ds1$gene, data=gff)
  # }
  genes_2_PAS <- ds1$gene
  PAdf2 <- annotatePAC(PAdf, gff)
  # PAdf2 <- PAdf2 %>%
  #   dplyr::filter(!is.na(gene))
  PAdf2 <- PAdf2[which(!is.na(PAdf2$gene)), ]
  if (str_detect(PAdf2$read_num[1], "_")) {
    PAs <- PAdf2 %>%
      tidyr::separate(read_num, c("read_num", "umi", "structure_type"), sep = "_") %>%
      filter(gene %in% genes_2_PAS) %>%
      select(read_num, gene, coord)
  } else {
    if (any(duplicated(names(PAdf2)))) {
      PAdf2 <- PAdf2[, -(which(duplicated(names(PAdf2))))]
    }
    PAs <- PAdf2 %>%
      filter(gene %in% genes_2_PAS) %>%
      select(read_num, gene, coord)
  }
  PACs <- ds %>%
    filter(gene %in% genes_2_PAS) %>%
    select(gene, strand, coord, UPA_start, UPA_end) %>%
    group_by(gene) %>%
    mutate(PAID = c("PA1", "PA2")) %>%
    ungroup()
  PAs_new <- data.frame()
  for (genes in genes_2_PAS) {
    sub_PAs <- filter(PAs, gene == genes)
    sub_PAs$PAID <- apply(sub_PAs, 1, dpPAtag, PACs)
    PAs_new <- rbind(PAs_new, sub_PAs)
  }
  PAs <- PAs_new
  rm(PAs_new)
  rm(sub_PAs)
  invisible(gc())
  PAs <- PAs %>%
    select(-coord)
  al1 <- PALdf %>%
    filter(gene %in% genes_2_PAS) %>%
    select(read_num, PAL)
  re <- invisible(left_join(PAs, al1))
  re <- na.omit(re)
  re <- re %>%
    dplyr::group_by(gene, PAID) %>%
    dplyr::summarise(counts = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(counts = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(counts == 2)
  genes_down <- re$gene
  re <- PAs %>%
    left_join(al1) %>%
    na.omit() %>%
    filter(gene %in% genes_down)
  result <- re %>%
    dplyr::group_by(gene, PAID) %>%
    dplyr::summarise(midPAL = median(PAL)) %>%
    tidyr::pivot_wider(names_from = PAID, values_from = midPAL) %>%
    dplyr::rename(distal_PAL = distal, proximal_PAL = proximal)
  result$DSAre <- apply(result, 1, myFun2, re)
  result <- tidyr::separate(result, DSAre, c("KS_re", "wilcox_re", "moses_re", "MU_re"), sep = "_")
  reback <- list(re1 = re, re2 = result)
  return(reback)
}



getNonItgFtrId <- function(ftr) {
  id <- grep("^inter", ftr)
  if (length(id) > 0) {
    return((1:length(ftr))[-id])
  } else {
    return(1:length(ftr))
  }
}


# If a PACds is 3UTRAPA
# Judge by only ftr=3UTR, have toStop column, and APA
# will not order PACds
is3UTRAPAds <- function(pacds) {
  if (nrow(pacds@anno) == 0) {
    return(FALSE)
  }
  if (sum(pacds@anno$ftr != "3UTR") > 0) {
    return(FALSE)
  }
  if (!("toStop" %in% colnames(pacds@anno))) {
    return(FALSE)
  }
  if (min(table(pacds@anno$gene)) <= 1) {
    return(FALSE)
  } # not APA
  return(TRUE)
}
