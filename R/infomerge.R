# description -------------------------------------------------------------
# This part of the code is used to consolidate data.

# infomerge1 ---------------------------------------------------------------
#' infomerge1
#'
#' @description Merging PA information and tail information.
#' @details Although we have identified the PA site and tail information for
#'   each sequence, the two tables are currently independent. We need to use
#'   the common column (reads_num) of the two tables to merge the data to
#'   produce the final data table.
#' @param tailinfo A dataframe that contains the tail information of all reads.
#' @param PAinfo A dataframe that contains the PA information of all reads.
#' @return This function returns a dataframe that contains both PA information
#'   (including PA gene, type based, chromosome, and so on) and tail information.
#' @usage
#' tailinfo = "the/path/of/tailinfo.txt"
#' PAinfo = "the/path/of/PAinfo.txt"
#' finaldata <- infomerge1(tailinfo,PAinfo)
#' @export
infomerge1 <- function(tailinfo,PAinfo){
  # 1. Build indexes
  df.tmp <- PAinfo[,-1]
  len <- length(rownames(df.tmp))
  c <- paste("PA", 1:len, sep = "")
  c <- data.frame(PAID = c)
  PAinfo <- cbind(c,df.tmp)
  # 2. find overlap
  PAlist <- c()
  genelist <- c()
  coord <- tail$coord
  i <- 0
  j <- 1
  l <- dim(PAinfo)[1]
  for(co in coord){
    i = i + 1
    if (i%%1000 == 0){
      print(i)
    }
    for(j in c(1:l)){
      start = PAinfo[j,]$UPA_start
      end = PAinfo[j,]$UPA_start
      if (start <= co & co <= end){
        PAlist <- append(PAlist,PAinfo[j,]$PAID)
        genelist <- append(genelist,PAinfo[j,]$gene)
        break
      }
    }
    continue
  }
  tailinfo$PAID <- PAlist
  tailinfo$gene <- genelist
}

# infomerge ---------------------------------------------------------------
#' infomerge
#'
#' @description Merging PA information and tail information.
#' @details Although we have identified the PA site and tail information for
#'   each sequence, the two tables are currently independent. We need to use
#'   the common column (reads_num) of the two tables to merge the data to
#'   produce the final data table.
#' @param tailinfo A dataframe that contains the tail information of all reads.
#' @param PAinfo A dataframe that contains the PA information of all reads.
#' @return This function returns a dataframe that contains both PA information
#'   (including PA gene, type based, chromosome, and so on) and tail information.
#' @usage
#' tailinfo = "the/path/of/tailinfo.txt"
#' PAinfo = "the/path/of/PAinfo.txt"
#' finaldata <- infomerge(tailinfo,PAinfo)
#' @export
infomerge <- function(tailinfo,PAinfo){
  # 1. Build indexes
  df.tmp <- PAinfo[,-1]
  len <- length(rownames(df.tmp))
  c <- paste("PA", 1:len, sep = "")
  c <- data.frame(PAID = c)
  PAinfo <- cbind(c,df.tmp)
  # 2. find overlap
  PAlist <- c()
  genelist <- c()
  coord <- tail$coord
  i <- 0
  j <- 1
  l <- dim(PAinfo)[1]
  for(co in coord){
    i = i + 1
    if (i%%1000 == 0){
      print(i)
    }
    for(j in c(1:l)){
      start = PAinfo[j,]$UPA_start
      end = PAinfo[j,]$UPA_start
      if (start <= co & co <= end){
        PAlist <- append(PAlist,PAinfo[j,]$PAID)
        genelist <- append(genelist,PAinfo[j,]$gene)
        break
      }
    }
    continue
  }
  tailinfo$PAID <- PAlist
  tailinfo$gene <- genelist
}

# PASTailMerge ---------------------------------------------------------------
#' PASTailMerge
#'
#' @description Merging PA information and tail information.
#' @details Although we have identified the PA site and tail information for
#'   each sequence, the two tables are currently independent. We need to use
#'   the common column (reads_num) of the two tables to merge the data to
#'   produce the final data table.
#' @param tailinfo A dataframe that contains the tail information of all reads.
#' @param PAinfo A dataframe that contains the PA information of all reads.
#' @return This function returns a dataframe that contains both PA information
#'   (including PA gene, type based, chromosome, and so on) and tail information.
#' @importFrom utils read.table
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @usage
#' tailinfo = "the/path/of/tailinfo.txt"
#' PAinfo = "the/path/of/PAinfo.txt"
#' finaldata <- PASTailMerge(tailinfo,PAinfo)
#' @export
PASTailMerge <- function(tailinfo,PAinfo){
  # 1. Build indexes
  taildf <- read.table(tailinfo,header=T,sep="\t")
  PAdf <- read.table(PAinfo,header = T,sep="\t")
  PAdf <- PAdf %>%
    separate(read_num,into = c("read_num","umi","tailType"),sep = "_") %>%
    select(-chr,-strand,-coord,-tailType,-sample,-gene_type)
  finaldf <- left_join(PAdf,taildf,by="read_num")
  return(finaldf)
}