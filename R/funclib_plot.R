# description -------------------------------------------------------------
# This section of code contains most of polyAtailor's visualization capabilities.


# innner functions --------------------------------------------------------
filter_lims <- function(x){
  l <- boxplot.stats(x)$stats[1]
  u <- boxplot.stats(x)$stats[5]
  
  for (i in 1:length(x)){
    x[i] <- ifelse(x[i]>l & x[i]<u, x[i], NA)
  }
  return(x)
}

# plotPALDistribution ---------------------------------------------------------------
#' plotPALDistribution
#'
#' @description Plot the PAL distribution of gene, umi or global.
#' @details Using the PoltPAL family of functions, the tail length distribution
#'   of the gene or UMI of different sample will be plotted, as well as the
#'   global tail length distribution.
#' @param taildf A dataframe that contains the tail information of all reads,
#'   Contains at least three columns: PAL, gene, umi, and sample.
#' @param resultpath The specified path to save barplot.
#' @param dType The tail length affiliation relationship that you want to draw,
#'   which includes three options: "global", "gene", "umi". "global" A stands
#'   for uniform treatment of all tail lengths. "gene" stands for tail length
#'   plotted according to genetic statistics. "umi" stands for tail length
#'   plotted according to umi statistics
#' @param medianPAL A Boolean value. T means Annotating peak. F means the
#'   opposite.
#' @return A density curve with peak annotations.
#' @examples
#' load("./data/output/AnnotatedTails.RData")
#' p1 <- plotPALDistribution(AnnotedTails,"./inst/extdata/figures/","global",medianPAL=T)
#' p2 <- plotPALDistribution(AnnotedTails,"./inst/extdata/figures/","gene",medianPAL=T)
#' @export
plotPALDistribution <- function(taildf,resultpath,dType,medianPAL){
  library(RColorBrewer)
  library(pracma)
  if(missing(taildf)){
    stop("Tail data must be provided")
  }
  if(missing(resultpath)){
    stop("You must provide a path to save the results")
  }
  if(missing(dType)){
    stop("The DTYPE must be specified")
  }
  else{
    if(!(dType %in% c("global","gene","umi"))){
      stop("dType must comes from c('global','gene','umi')")
    }
  }
  if(missing(medianPAL)){
    medianPAL = T
  }
  if(dType=="global" & medianPAL){
    data <- taildf %>%
      select(PAL,sample) %>%
      group_by(sample,PAL) %>%
      summarise(count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$PAL <- density(sub.data$PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=PAL,colour=sample))+
      xlab("PAL(nt)")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$PAL[i],linetype = "longdash",color=colors[i],size = 0.7)+
        # annotate(geom = "rect", xmin = 10, xmax = 18, ymin = 0.05, ymax = 0.07, alpha = .1, fill="navy")+
        annotate(geom = "text", fontface = "bold", color=colors[i],
                 x = (peak.table$PAL[i])+5,y=0.000+(i-1)*0.0003,
                 label = round((peak.table$PAL[i]),2), size=3)
    }
    plot(p)
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
  if(dType=="global" & !medianPAL){
    data <- taildf %>%
      select(PAL,sample) %>%
      group_by(sample,PAL) %>%
      summarise(count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$PAL <- density(sub.data$PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=PAL,colour=sample))+
      xlab("PAL(nt))")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$PAL[i],linetype = "longdash",color=colors[i],size = 0.7)
    }
    plot(p)
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
  if(dType=="gene" & medianPAL){
    data <- taildf %>%
      select(gene,PAL,sample) %>%
      group_by(sample,gene) %>%
      summarise(md_PAL = median(PAL),count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             md_PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$md_PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$md_PAL <- density(sub.data$md_PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=md_PAL,colour=sample))+
      xlab("PAL(nt) per GENE")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$md_PAL[i],linetype = "longdash",color=colors[i],size = 0.7)+
        # annotate(geom = "rect", xmin = 10, xmax = 18, ymin = 0.05, ymax = 0.07, alpha = .1, fill="navy")+
        annotate(geom = "text", fontface = "bold", color=colors[i],
                 x = (peak.table$md_PAL[i])+5,y=0.000+(i-1)*0.0003,
                 label = round((peak.table$md_PAL[i]),2), size=3)
    }
    plot(p)
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
  if(dType=="gene" & !medianPAL){
    data <- taildf %>%
      select(gene,PAL,sample) %>%
      group_by(sample,gene) %>%
      summarise(md_PAL = median(PAL),count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             md_PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$md_PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$md_PAL <- density(sub.data$md_PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=md_PAL,colour=sample))+
      xlab("PAL(nt) per GENE")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$md_PAL[i],linetype = "longdash",color=colors[i],size = 0.7)
    }
    p
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
  if(dType=="umi" & medianPAL){
    data <- taildf %>%
      select(umi,PAL,sample) %>%
      group_by(sample,umi) %>%
      summarise(md_PAL = median(PAL),count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             md_PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$md_PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$md_PAL <- density(sub.data$md_PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=md_PAL,colour=sample))+
      xlab("PAL(nt) per UMI")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$md_PAL[i],linetype = "longdash",color=colors[i],size = 0.7)+
        # annotate(geom = "rect", xmin = 10, xmax = 18, ymin = 0.05, ymax = 0.07, alpha = .1, fill="navy")+
        annotate(geom = "text", fontface = "bold", color=colors[i],
                 x = (peak.table$md_PAL[i])+5,y=0.000+(i-1)*0.0003,
                 label = round((peak.table$md_PAL[i]),2), size=3)
    }
    plot(p)
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
  if(dType=="umi" & !medianPAL){
    data <- taildf %>%
      select(umi,PAL,sample) %>%
      group_by(sample,umi) %>%
      summarise(md_PAL = median(PAL),count=n())
    sample.names <- unique(data$sample)
    peak.table <- data.frame(sample=NULL,
                             umi=NULL,
                             md_PAL=NULL,
                             filename=NULL,
                             count=NULL)
    for(i in 1:length(sample.names)){
      sub.data <- dplyr::filter(.data =data, sample==sample.names[i])
      index <- which.max(density(sub.data$md_PAL)$y)
      
      peak.id <-as.data.frame(sub.data[findpeaks(sub.data$count, npeaks=1, threshold=0, sortstr=TRUE)[,2],])
      peak.id$md_PAL <- density(sub.data$md_PAL)$x[index]
      peak.table  <- base::rbind(peak.table,peak.id)
    }
    #geom_density
    data$sample  <- factor(data$sample,levels = sample.names)
    p <- ggplot(data,aes(x=md_PAL,colour=sample))+
      xlab("PAL(nt) per UMI")+
      geom_density(size=1.1)+
      theme_base()+
      # scale_color_manual(values=mycolors)+
      scale_color_brewer(palette = "Dark2")+
      theme(legend.position="top")
    colors <- brewer.pal(8,"Dark2")
    for(i in 1:nrow(peak.table)){
      p <- p+geom_vline(xintercept = peak.table$md_PAL[i],linetype = "longdash",color=colors[i],size = 0.7)
    }
    plot(p)
    fil = str_c(resultpath,"PALDistribution.pptx")
    topptx(filename = fil)
    return(p)
  }
}


# plotPADistribution ---------------------------------------------------------------
#' plotPADistribution
#'
#' @description Plot the gene region distribution of PA sites.
#' @details Plot the gene region distribution of PA sites, including "3UTR",
#'   "5UTR", "CDS", "exon", "intergenic", "intron". Finally, we will save the PA
#'   locus distribution barpllot in PPTX format to the path you specified.
#' @param PAdf A dataframe that contains the PA information of all reads,
#'   Contains at least three columns: PAID, gene, ftr.
#' @param resultpath The specified path to save barplot.
#' @param colorOfBar The color of barplot, as input in hexadecimal form,
#'   defaults to "#196BE1".
#' @return A barplot.
#' @examples
#' load("./data/output/PAs.RData")
#' p <- plotPADistribution(PAs,"./inst/extdata/figures/","#9BBFDC")
#' @export
plotPADistribution <- function(PAdf,resultpath,colorOfBar){
  if(missing(colorOfBar)){
    colorOfBar = "#196BE1"
  }
  if(missing(resultpath)){
    stop("The result path must be specified")
  }
  if(missing(PAdf)){
    stop("PA information must be entered")
  }
  data1 <- PAdf %>%
    select(gene,ftr) %>%
    filter(!is.na(gene)&gene!=("character(0)"))
  data2 <- data1 %>%
    group_by(ftr) %>%
    summarise(count=n())
  data2 <- mutate(data2,rt=count/sum(count)*100)
  p <- data2 %>%
    ggplot() +
    aes(x = reorder(ftr,-rt),y=rt) +
    geom_bar(stat = 'identity',position = 'dodge',color=colorOfBar,fill=colorOfBar) +
    scale_fill_manual(values=c(colorOfBar))+
    scale_color_manual(values=c(colorOfBar))+
    labs(x = "", y = "proportion(%)") +
    geom_text(aes(label = count), position = position_dodge(0.9),vjust=-0.5)+
    ggthemes::theme_base()+
    theme(legend.position="top",
          axis.text.x  = element_text(angle=30, vjust=0.5))
  plot(p)
  fil = str_c(resultpath,"barPlotOfPASdistribution.pptx")
  topptx(filename = fil)
  return(p)
}





# plotGenePAnumbers ---------------------------------------------------------------
#' plotGenePAnumbers
#'
#' @description Plot the gene frequency distribution with different number of PA sites.
#' @details Plot the gene frequency distribution with different number of PA
#'   sites. Genes were divided into five groups, each containing the number of
#'   PA sites 1, 2, 3, 4, 5 and greater than 5, respectively.
#' @param PAdf A dataframe that contains the PA information of all reads,
#'   Contains at least three columns: PAID, gene, ftr.
#' @param resultpath The specified path to save barplot.
#' @param colorOfBar The color of barplot, as input in hexadecimal form,
#'   defaults to "#196BE1".
#' @return A barplot.
#' @examples
#' load("./data/output/PAs.RData")
#' p <- plotGenePAnumbers(PAs,"./inst/extdata/figures/","#DF7C7D")
#' @export
plotGenePAnumbers <- function(PAdf,resultpath,colorOfBar){
  if(missing(colorOfBar)){
    colorOfBar = "#B681BD"
  }
  if(missing(PAdf)){
    stop("PA information must be entered")
  }
  if(missing(resultpath)){
    stop("The result path must be specified")
  }
  if("PAID" %in% colnames(PAdf)){
    data1 <- PAdf %>%
      select(PAID,gene,ftr) %>%
      filter(!is.na(gene)&gene!=("character(0)"))
    data9D<-data1 %>%
      group_by(gene,PAID)%>%
      summarise(counts=n()) %>%
      group_by(gene) %>%
      summarise(PAs=n())
  }
  else{
    data1 <- PAdf %>%
      select(gene,ftr) %>%
      filter(!is.na(gene)&gene!=("character(0)"))
    data9D<-data1 %>%
      group_by(gene)%>%
      summarise(PAs=n())
  }
  #data9D <- table(data9D$PAs)
  data9D$type.new <- data9D$PAs
  data9D$type.new[which(data9D$type.new>5)] <- ">5"
  data9D$type.new <- factor(data9D$type.new,levels=c("1","2","3","4","5",">5"))
  
  data9D.fre <- data.frame(fre=names(table(data9D$type.new)),number=as.integer(table(data9D$type.new)) )
  data9D.fre$fre<-factor(data9D.fre$fre,levels = unique(data9D.fre$fre))
  
  #table(data9D$type.new)
  library(RColorBrewer)
  p <- data9D.fre %>%
    ggplot() +
    aes(x = factor(fre),y=number,fill="#196BE1") +
    geom_bar(stat = 'identity',position = 'dodge')+# +scale_fill_brewer(palette="Dark2")+
    # scale_fill_manual(values=c("#196BE1"))+
    #scale_color_manual(values=c("#B681BD"))+
    labs(x = "Number of PA sites", y = "Gene Frequency") +
    geom_text(aes(label = number), position = position_dodge(0.9),vjust=-0.5)+
    ggthemes::theme_base()+
    theme(legend.position="none",
          axis.text.x  = element_text(angle=0, vjust=0.5))
  plot(p)
  fil = str_c(resultpath,"PAnumbersOfGene.pptx")
  topptx(filename = fil)
  return(p)
}




# plotPASignals ---------------------------------------------------------------
#' plotPASignals
#'
#' @description Draw the frequency distribution diagram of PA signal.
#' @details plot the probability distribution of the occurrence of
#'   user-specified or default PA signal within the first 50bp base range of PA
#'   site.
#' @param PAdf A dataframe that contains the PA information of all reads,
#'   Contains at least three columns: PAID, gene, ftr.
#' @param resultpath The specified path to save barplot.
#' @param colorOfBar The color of barplot, as input in hexadecimal form,
#'   defaults to "#196BE1".
#' @param bsgenome Reference genomes of the species of interest.
#' @param signals The set of PA signal sequences of interest, such as
#'   signals=c("AATAAA","ATTAAA","AAGAAA","AATATA").
#' @return A barplot.
#' @examples
#' library("BSgenome.Mmusculus.UCSC.mm10")
#' bsgenome = BSgenome.Mmusculus.UCSC.mm10
#' load("./data/output/PAs.RData")
#' p <- plotPASignals(PAs,"./inst/extdata/figures/",bsgenome = bsgenome)
#' @export
plotPASignals <- function(PAdf,resultpath,bsgenome,signals,colorOfBar){
  if(missing(colorOfBar)){
    colorOfBar = "#7ED321"
  }
  if(missing(PAdf)){
    stop("PA information must be entered")
  }
  if(missing(signals)){
    print("No PA signal is specified. Statistics are performed according to the default values.")
    signals=c("AATAAA","ATTAAA","AAGAAA","AAAAAG","AATACA","TATAAA","ACTAAA","AGTAAA","GATAAA",
              "AATATA")
  }
  if(missing(bsgenome)){
    stop("Reference genomes must be specified")
  }
  if(missing(resultpath)){
    stop("The result path must be specified")
  }
  library(movAPA)
  PAdf1 <- PAdf %>%
    dplyr::select(chr,strand,coord,ftr)
  PACds <- readPACds(PAdf1)
  priority=c(1,2,rep(3, length(signals)-2))
  PACdsHUMAN=annotateByPAS(PACds, bsgenome, grams=signals,
                           priority=priority, from=-50, to=-1, label='HU')
  test2 <- table(PACdsHUMAN@anno$HU_gram)
  data9F <- as.data.frame(test2)
  all_counts = sum(data9F$Freq)
  data9F <- data9F %>%
    mutate(rt = (Freq / all_counts) *100)
  p <- data9F %>%
    ggplot() +
    aes(x = reorder(Var1,-rt),y=rt) +
    geom_bar(stat = 'identity',position = 'dodge',color="#57B355",fill="#57B355") +
    scale_fill_manual(values=c("#57B355"))+
    scale_color_manual(values=c("#57B355"))+
    labs(x = "", y = "ratio(%)") +
    geom_text(aes(label = Freq), position = position_dodge(0.9),vjust=-0.5)+
    ggthemes::theme_base()+
    theme(legend.position="top",
          axis.text.x  = element_text(angle=30, vjust=0.5))
  plot(p)
  fil = str_c(resultpath,"barPlotOfPASignaldistribution.pptx")
  topptx(filename = fil)
  return(p)
}

# nonAanalysis ---------------------------------------------------
#' nonAanalysis
#'
#' @description Quantify and Statistics the nonA base in tail.
#' @details After users input the tail information table, we will quantify the
#'   non-A base content in each tail, and count the number of reads with non-A
#'   base tail in each sample and the proportion of non-A base content in each
#'   sample.
#' @param tailinfo A tail information table containing at least PAL, sample
#'   name, tail sequence, chain direction, and gene infomation.
#' @param justOneTail Boolean value. Indicates whether to analyze only data of
#'   type one-tail, which defaults to T.
#' @return A large list containing quantitative data tables and statistical
#'   graphs.The two figures respectively represent the column chart of the
#'   number and proportion of reads with non-A bases in each sample and the
#'   specific occurrence frequency of non-A bases in each sample.
#' @examples
#' load("./data/AnnotatedTails.RData")
#' re <- nonAanalysis(AnnotedTails)
#' re$p1
#' re$p2
#' re$p3
#' nonAre <- re$nonAinfo
#' @export
nonAanalysis <- function(tailinfo,justOneTail){
  if(missing(justOneTail)){
    justOneTail <- T
  }
  if(!is.logical(justOneTail)){
    stop("The justOneMail parameter must be a logical value!")
  }
  ##First extract the useful columns of the entered dataframe
  resultlist <- list()
  inputlist <- tailinfo %>%
    select(read_num, sample, strand, tail, PAL)
  ##Find the reverse complementation of all the reverse chains:polyA input
  tailinfo <- tailinfo %>%
    dplyr::mutate(nA=str_count(tail,"A"),rA=nA/PAL)
  subset1 <- tailinfo %>%
    dplyr::filter(rA<=0.5)
  subset2 <-tailinfo %>%
    dplyr::filter(rA>0.5)
  for(i in c(1:dim(subset1)[1])){
    #Build the mapping
    y <- c("T"="A", "A"="T", "G"="C" , "C"="G")
    #Calculate the reverse complementary chain
    x <- unlist(strsplit(subset1$tail[i], "", fixed=TRUE))
    x <- y[x]
    names(x) = NULL
    subset1$tail[i] <- str_c(x,collapse = "")
  }
  inputlist <- rbind(subset1,subset2)
  ##Calculate nonA bases
  inputlist <- inputlist %>%
    mutate(Acounts = stri_count(tail,fixed = "A"),
           Ccounts = stri_count(tail,fixed = "C"),
           Gcounts = stri_count(tail,fixed = "G"),
           Tcounts = stri_count(tail,fixed = "T"),
           nonA_number = PAL-Acounts,
           nonA_rate = nonA_number/PAL
    )
  inputlist <- inputlist %>%
    dplyr::filter(rt>=0.9)
  if(justOneTail){
    inputlist <- inputlist %>%
      dplyr::filter(read_type == "one-tail")
  }
  resultlist$nonAinfo <- inputlist
  ##upset plot
  #1.A simple statistical chart counted the number of reads containing nonA Base tail
  ##count the reads number of all kinds
  input <- c(
    "C"=length(which(inputlist[,"Ccounts"]!=0 & inputlist[,"Gcounts"]==0 & inputlist[,"Tcounts"]==0)),
    "G"=length(which(inputlist[,"Ccounts"]==0 & inputlist[,"Gcounts"]!=0 & inputlist[,"Tcounts"]==0)),
    "T"=length(which(inputlist[,"Ccounts"]==0 & inputlist[,"Gcounts"]==0 & inputlist[,"Tcounts"]!=0)),
    "C&G"=length(which(inputlist[,"Ccounts"]!=0 & inputlist[,"Gcounts"]!=0 & inputlist[,"Tcounts"]==0)),
    "C&T"=length(which(inputlist[,"Ccounts"]!=0 & inputlist[,"Gcounts"]==0 & inputlist[,"Tcounts"]!=0)),
    "G&T"=length(which(inputlist[,"Ccounts"]==0 & inputlist[,"Gcounts"]!=0 & inputlist[,"Tcounts"]!=0)),
    "C&G&T"=length(which(inputlist[,"Ccounts"]!=0 & inputlist[,"Gcounts"]!=0 & inputlist[,"Tcounts"]!=0))
  )
  library(UpSetR)
  data <- UpSetR::fromExpression(input)
  p1 <- upset(data,
              mainbar.y.label = "reads counts",
              sets.x.label = "reads counts",
              queries = list(list(query=intersects, params=list("C", "G"), color="#F18687", active=T),
                             list(query=intersects, params=list("C", "T"), color="#92D050", active=T),
                             list(query=intersects, params=list("G", "T"), color="#9BBFDC", active=T),
                             List(query=intersects, params=list("C","G","T"),color="#F6BF18",active=T)
              )
  )
  # p1
  resultlist$p1 <- p1
  #2.The line chart1
  base=c("<12","12-24","24-36","36-48","48-60","60-72","72-84","84-96",">=96")
  C_c=c()
  G_c=c()
  T_c=c()
  nonA_c=c()
  all=c()
  for (i in c(1:9)) {
    k = (i-1)*12
    start = 0 + k
    end = 12 + k
    if(i!=9){
      inputlist1 <- inputlist %>%
        dplyr::filter(PAL >= start & PAL < end)
      allc = dim(inputlist1)[1]
      Cc=length(which(inputlist1[,"Ccounts"]!=0 & inputlist1[,"Gcounts"]==0 & inputlist1[,"Tcounts"]==0))
      Gc=length(which(inputlist1[,"Ccounts"]==0 & inputlist1[,"Gcounts"]!=0 & inputlist1[,"Tcounts"]==0))
      Tc=length(which(inputlist1[,"Ccounts"]==0 & inputlist1[,"Gcounts"]==0 & inputlist1[,"Tcounts"]!=0))
      nonAc=length(which(inputlist1[,"Ccounts"]!=0 | inputlist1[,"Gcounts"]!=0 | inputlist1[,"Tcounts"]!=0))
      C_c[i] = Cc
      G_c[i] = Gc
      T_c[i] = Tc
      nonA_c[i] = nonAc
      all[i] = allc
    }
    else{
      inputlist1 <- inputlist %>%
        dplyr::filter(PAL >= start)
      allc = dim(inputlist1)[1]
      Cc=length(which(inputlist1[,"Ccounts"]!=0 & inputlist1[,"Gcounts"]==0 & inputlist1[,"Tcounts"]==0))
      Gc=length(which(inputlist1[,"Ccounts"]==0 & inputlist1[,"Gcounts"]!=0 & inputlist1[,"Tcounts"]==0))
      Tc=length(which(inputlist1[,"Ccounts"]==0 & inputlist1[,"Gcounts"]==0 & inputlist1[,"Tcounts"]!=0))
      nonAc=length(which(inputlist1[,"Ccounts"]!=0 | inputlist1[,"Gcounts"]!=0 | inputlist1[,"Tcounts"]!=0))
      C_c[i] = Cc
      G_c[i] = Gc
      T_c[i] = Tc
      nonA_c[i] = nonAc
      all[i] = allc
    }
  }
  data <- data.frame(base=base,C_c=C_c,G_c=G_c,T_c=T_c,nonA_c=nonA_c,all=all)
  data <- data %>%
    mutate(Crt = C_c / all * 100,
           Grt = G_c / all * 100,
           Trt = T_c / all * 100,
           nonArt = nonA_c / all * 100)
  d1 <- data %>%
    select(base,Crt) %>%
    mutate(type="C") %>%
    rename(rt=Crt)
  d2 <- data %>%
    select(base,Grt) %>%
    mutate(type="G") %>%
    rename(rt=Grt)
  d3 <- data %>%
    select(base,Trt) %>%
    mutate(type="T") %>%
    rename(rt=Trt)
  d4 <- data %>%
    select(base,nonArt) %>%
    mutate(type="nonA") %>%
    rename(rt=nonArt)
  data10B <- rbind(d1,d2,d3,d4)
  data10B$base <- factor(data10B$base,levels=c("<12","12-24","24-36","36-48","48-60","60-72","72-84","84-96",">=96"))
  p2 <- ggplot(data10B, aes(x=factor(base), y=rt, colour=type,group=type)) +
    geom_line(size=1.5)+
    geom_point(size=2.5)+
    ggthemes::theme_base()+
    scale_fill_manual(values=c("#7ED321","#F28D8E","#AA87C2","#9BBFDC"))+
    scale_color_manual(values=c("#7ED321","#F28D8E","#AA87C2","#9BBFDC"))+
    labs(x = "PAL(nt)", y = "reads with non_A base/total one-tail reads(%)") +
    theme(legend.position="top",axis.text.x  = element_text(angle=30, vjust=0.5))
  # +
  #   scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5))
  # p2
  resultlist$p2 <- p2
  #3.The line chart2
  C_c=c()
  G_c=c()
  T_c=c()
  nonA_c=c()
  all=c()
  for (i in c(1:9)) {
    k = (i-1)*12
    start = 0 + k
    end = 12 + k
    if(i!=9){
      inputlist1 <- inputlist %>%
        dplyr::filter(PAL >= start & PAL < end)
      allc = sum(inputlist1$PAL)
      Cc=sum(inputlist1$Ccounts)
      Gc=sum(inputlist1$Gcounts)
      Tc=sum(inputlist1$Tcounts)
      nonAc=sum(inputlist1$nonA_number)
      C_c[i] = Cc
      G_c[i] = Gc
      T_c[i] = Tc
      nonA_c[i] = nonAc
      all[i] = allc
    }
    else{
      inputlist1 <- inputlist %>%
        dplyr::filter(PAL >= start)
      allc = sum(inputlist1$PAL)
      Cc=sum(inputlist1$Ccounts)
      Gc=sum(inputlist1$Gcounts)
      Tc=sum(inputlist1$Tcounts)
      nonAc=sum(inputlist1$nonA_number)
      C_c[i] = Cc
      G_c[i] = Gc
      T_c[i] = Tc
      nonA_c[i] = nonAc
      all[i] = allc
    }
  }
  data10C <- data.frame(base=base,C_c=C_c,G_c=G_c,T_c=T_c,nonA_c=nonA_c,all=all)
  data10C <- data10C %>%
    mutate(Crt = C_c / all * 1000,
           Grt = G_c / all * 1000,
           Trt = T_c / all * 1000,
           nonArt = nonA_c / all * 1000)
  d1 <- data10C %>%
    select(base,Crt) %>%
    mutate(type="C") %>%
    rename(rt=Crt)
  d2 <- data10C %>%
    select(base,Grt) %>%
    mutate(type="G") %>%
    rename(rt=Grt)
  d3 <- data10C %>%
    select(base,Trt) %>%
    mutate(type="T") %>%
    rename(rt=Trt)
  data10C <- rbind(d1,d2,d3)
  data10C$base <- factor(data10C$base,levels=c("<12","12-24","24-36","36-48","48-60","60-72","72-84","84-96",">=96"))
  # data10B$type <- factor(data10B$type,levels=c("non_A","C","G","T"))
  p3 <- ggplot(data10C, aes(x=factor(base), y=rt, colour=type,group=type)) +
    geom_line(size=1.5)+
    geom_point(size=2.5)+
    ggthemes::theme_base()+
    scale_fill_manual(values=c("#7ED321","#F28D8E","#9BBFDC"))+
    scale_color_manual(values=c("#7ED321","#F28D8E","#9BBFDC"))+
    labs(x = "PAL(nt)", y = "Frequency of non-A residues (\u2030)") +
    theme(legend.position="top",axis.text.x  = element_text(angle=30, vjust=0.5))
  # +
  #   scale_y_continuous(limits=c(0,10), breaks=seq(0,10,1))
  resultlist$p3 <- p3
  return(resultlist)
}



# tailViso ---------------------------------------------------
#' tailViso
#'
#' @description Visualize tails using heat maps and logos.
#' @details Visualize the tail input by the user so that the base composition of
#'   different tails can be compared intuitively.
#' @param tails A dataframe. The set of tails that the user wants to visualize
#'   contains at least one set of tail sequences and a unique identifier for
#'   each tail.
#' @param tailLen Maximum tail length to draw.
#' @param Ntail Number of tails to draw.
#' @param custom A dataframe. Custom base color matching information, one column
#'   is the base name, the other column is the corresponding color scheme.
#' @param strand The direction of the tail chain, + represents the polyA tail
#'   and - represents the polyT tail.
#' @param faPath Tail sequence file path.
#' @param showLogo Logical value indicating whether to draw the LOGO diagram,
#'   The default is T.
#' @param showReadNum Logical value indicating whether to display the sequence
#'   identifier, default to T.
#' @return Return the picture.
#' @examples
#' load(./data/taildf.RData)
#' my_cutstom <- data.frame(names=c("A","C","T","G"),color=c("#3171A5","#4EAA4C","#C9C4C2","#D73D3D"))
#' p <- tailViso(taildf,tailLen=100,Ntail=20,custom=my_cutstom,strand="-",faPath="D:/",showLogo=T,showReadNum= F)
#' @export
tailViso <- function(taildf,tailLen,Ntail,custom,strand,faPath,showLogo,showReadNum){
  ##check the param
  if(missing(showReadNum)){
    showReadNum=F
  }
  if(missing(showLogo)){
    showLogo <- T
  }
  if(!is.logical(showLogo)){
    stop("ShowLogo must be a logical value")
  }
  if(!is.logical(showReadNum)){
    stop("showReadNum must be a logical value")
  }
  if(missing(taildf)){
    stop("taildf lossing")
  }
  if(missing(faPath)){
    stop("faPath lossing")
  }
  if(!(strand %in% c("+","-"))){
    stop("The STRAND parameter must be + or -")
  }
  if(missing(custom)){
    my_cutstom <- data.frame(names=c("A","C","T","G"),color=c("#3171A5","#4EAA4C","#C9C4C2","#D73D3D"))
  }
  else{
    my_cutstom <- custom
  }
  if(missing(Ntail)){
    stop("Ntail lossing")
  }
  if(missing(tailLen)){
    stop("tailLen lossing")
  }
  ##plot
  taildf <- select(taildf,read_num,tail,PAL)
  if(Ntail>120){
    print("Too many tails were input, only 120 tails were extracted for processing")
    Ntail = 120
  }
  if(Ntail>dim(taildf)[1]){
    Ntail = dim(taildf)[1]
  }
  if(strand == "+"){
    seqs <- apply(taildf, 1, seqreverse)
    taildf$tail <- seqs
  }
  taildf <- taildf[sample(nrow(taildf), Ntail), ]
  maxTailL <- max(taildf$PAL)
  if(maxTailL>=tailLen){
    print("The input tail is too long, only the length you specify is cut for processing")
    maxTailL = tailLen
  }
  taildf <- select(taildf,read_num,tail)
  taildf$tail <- str_sub(taildf$tail,1,maxTailL)
  taildf$tail <- str_pad(taildf$tail,maxTailL,side = "right",pad = "-")
  filepath=str_c(faPath,"subseq.fasta")
  fa <- dataframe2fas(as.data.frame(taildf),filepath)
  if(showLogo){
    p <- ggmsa(filepath, start=1,end=80,font = T,custom_color=my_cutstom,
               char_width=0.7,seq_name=showReadNum) +
      geom_seqlogo("helvetical",top=T,adaptive=T,custom_color=my_cutstom)
  }
  else{
    p <- ggmsa(filepath, start=1,end=80,font = T,custom_color=my_cutstom,
               char_width=0.7,seq_name=showReadNum)
  }
  plot(p)
  return(p)
}

# DSA_ViolinPlot ----------------------------------------------------------
#' DSA_ViolinPlot
#'
#' @description Draw a violin diagram of the result of differential tail length
#'   analysis.
#' @details The tail length difference analysis results of the four difference
#'   significance test methods are displayed through the violin diagram. The
#'   tail length distribution under different conditions is represented by
#'   different colors to facilitate users to view the differences.
#' @param DSA_result The output of PALdsa(). 
#' @param SAoDMethod one of "KS", "MWU", "ME", "Wilcox" and "ALL".
#' @return A ViolinPlot.
#' @examples
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' gff <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' load("./data/AnnotedTails.RData")
#' PAs <- read.table("./inst/extdata/output/PAs/PAs.txt",header=T,sep=" ")
#' DSA_result <- PALdsa(PAs,AnnotedTails,gff,mode="PD",SAoDMethod="ME",withViolinPlot=T,withUpsetPlot=F)
#' p <- DSA_ViolinPlot(DSA_result,SAoDMethod="KS")
#' @export
DSA_ViolinPlot <- function(DSA_result,SAoDMethod){
  #data input
  #conflict_prefer("filter", "dplyr")
  KS_re <- dplyr::filter(DSA_result$re2,KS_re<=0.05)[,1]
  KS_count <- str_c("K-S test\n",dim(KS_re)[1]," genes")
  KS_re <- DSA_result$re1 %>%
    dplyr::filter(gene %in% KS_re$gene) %>%
    dplyr::select(PAID,PAL) %>%
    dplyr::mutate(method=KS_count)
  Wilcox_re <- filter(DSA_result$re2,wilcox_re<=0.05)[,1]
  Wilcox_count <- str_c("Wilcox test\n",dim(Wilcox_re)[1]," genes")
  Wilcox_re <- DSA_result$re1 %>%
    dplyr::filter(gene %in% Wilcox_re$gene) %>%
    dplyr::select(PAID,PAL) %>%
    dplyr::mutate(method=Wilcox_count)
  moses_re <- dplyr::filter(DSA_result$re2,moses_re<=0.05)[,1]
  moses_count <- str_c("MosesExtreme reaction\n",dim(moses_re)[1]," genes")
  moses_re <- DSA_result$re1 %>%
    dplyr::filter(gene %in% moses_re$gene) %>%
    dplyr::select(PAID,PAL) %>%
    dplyr::mutate(method=moses_count)
  MU_re <- dplyr::filter(DSA_result$re2,MU_re<=0.05)[,1]
  MU_count <- str_c("M-W U test\n",dim(MU_re)[1]," genes")
  MU_re <- DSA_result$re1 %>%
    dplyr::filter(gene %in% MU_re$gene) %>%
    dplyr::select(PAID,PAL) %>%
    dplyr::mutate(method=MU_count)
  if(SAoDMethod == "ALL"){
    plot_data <- rbind(KS_re,Wilcox_re,moses_re,MU_re)
    P <- ggplot(plot_data, aes(x=method, y=PAL, fill=PAID)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.2,position=position_dodge(0.9),color="black")+
      scale_fill_manual(values = c("#ef857d","#00b0f0"))+
      scale_color_manual(values = c("#ef857d","#00b0f0"))+
      theme_base()+
      labs(x = "", y = "PAL per 3'UTR isform") +
      theme(legend.position="top",
            axis.text.y=element_text(family="Arial",size=10,face="plain"),
            axis.title.y=element_text(family="Arial",size = 12,face="plain"),
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
            legend.text=element_text(face="plain", family="Arial", colour="black",
                                     size=10),
            legend.title=element_text(face="plain", family="Arial", colour="black",
                                      size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(P)
  }
  else if(SAoDMethod == "KS"){
    plot_data <- KS_re
    P <- ggplot(plot_data, aes(x=method, y=PAL, fill=PAID)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.2,position=position_dodge(0.9),color="black")+
      scale_fill_manual(values = c("#ef857d","#00b0f0"))+
      scale_color_manual(values = c("#ef857d","#00b0f0"))+
      theme_base()+
      labs(x = "", y = "PAL per 3'UTR isform") +
      theme(legend.position="top",
            axis.text.y=element_text(family="Arial",size=10,face="plain"),
            axis.title.y=element_text(family="Arial",size = 12,face="plain"),
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
            legend.text=element_text(face="plain", family="Arial", colour="black",
                                     size=10),
            legend.title=element_text(face="plain", family="Arial", colour="black",
                                      size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(P)
  }
  else if(SAoDMethod == "MWU"){
    plot_data <- MU_re
    P <- ggplot(plot_data, aes(x=method, y=PAL, fill=PAID)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.2,position=position_dodge(0.9),color="black")+
      scale_fill_manual(values = c("#ef857d","#00b0f0"))+
      scale_color_manual(values = c("#ef857d","#00b0f0"))+
      theme_base()+
      labs(x = "", y = "PAL per 3'UTR isform") +
      theme(legend.position="top",
            axis.text.y=element_text(family="Arial",size=10,face="plain"),
            axis.title.y=element_text(family="Arial",size = 12,face="plain"),
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
            legend.text=element_text(face="plain", family="Arial", colour="black",
                                     size=10),
            legend.title=element_text(face="plain", family="Arial", colour="black",
                                      size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(P)
  }
  else if(SAoDMethod == "Wilcox"){
    plot_data <- Wilcox_re
    P <- ggplot(plot_data, aes(x=method, y=PAL, fill=PAID)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.2,position=position_dodge(0.9),color="black")+
      scale_fill_manual(values = c("#ef857d","#00b0f0"))+
      scale_color_manual(values = c("#ef857d","#00b0f0"))+
      theme_base()+
      labs(x = "", y = "PAL per 3'UTR isform") +
      theme(legend.position="top",
            axis.text.y=element_text(family="Arial",size=10,face="plain"),
            axis.title.y=element_text(family="Arial",size = 12,face="plain"),
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
            legend.text=element_text(face="plain", family="Arial", colour="black",
                                     size=10),
            legend.title=element_text(face="plain", family="Arial", colour="black",
                                      size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(P)
  }
  else if(SAoDMethod == "ME"){
    plot_data <- moses_re
    P <- ggplot(plot_data, aes(x=method, y=PAL, fill=PAID)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.2,position=position_dodge(0.9),color="black")+
      scale_fill_manual(values = c("#ef857d","#00b0f0"))+
      scale_color_manual(values = c("#ef857d","#00b0f0"))+
      theme_base()+
      labs(x = "", y = "PAL per 3'UTR isform") +
      theme(legend.position="top",
            axis.text.y=element_text(family="Arial",size=10,face="plain"),
            axis.title.y=element_text(family="Arial",size = 12,face="plain"),
            panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
            legend.text=element_text(face="plain", family="Arial", colour="black",
                                     size=10),
            legend.title=element_text(face="plain", family="Arial", colour="black",
                                      size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    return(P)
  }
  else{
    stop("SAoDMethod must be one of c('KS', 'MWU', 'ME', 'Wilcox', 'ALL')!")
  }
}


# DSA_UpsetPlot ----------------------------------------------------------
#' DSA_UpsetPlot
#'
#' @description Draw a Upset diagram of the result of differential tail length
#'   analysis.
#' @details Overlapping analysis of the results of differential tail length
#'   analysis of the four difference significance testing methods can help users
#'   quickly determine which results of differential tail length analysis are
#'   highly confident. The consensus results of the four methods are marked in
#'   red to represent highly confident tail length genes.
#' @param DSA_result The output of PALdsa(). 
#' @return A UpsetPlot.
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' gff <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' data(AnnotedTails)
#' files = system.file("extdata", "./output/PAs/PAs.txt", package = "PolyAtailor", mustWork = TRUE)
#' PAs <- read.table(files,header=TRUE,sep=" ")
#' DSA_result <- PALdsa(PAdf=PAs,PALdf=AnnotedTails,d=24,gff=gff,mode="PD",SAoDMethod="ME",withViolinPlot=T,withUpsetPlot=F)
#' p <- DSA_UpsetPlot(DSA_result)
#' @export
DSA_UpsetPlot <- function(DSA_result){
  plot_data <- DSA_result[,-c(2,3)]
  plot_data[plot_data == "lose"] <- "0"
  plot_data[is.na(plot_data)] <- "0"
  plot_data[plot_data == "NaN"] <- "0"
  plot_data_col1 <- as.character(plot_data$gene)
  plot_data <- plot_data[,-1]
  plot_data <- as.data.frame(lapply(plot_data,as.numeric))
  plot_data[plot_data < 0.05] <- 0.05
  plot_data[plot_data > 0.05] <- 0
  plot_data[plot_data == 0.05] <- 1
  plot_data$gene <- plot_data_col1
  plot_data <- plot_data %>%
    rename(KS=KS_re,Wilcox=wilcox_re,Moses=moses_re,MU=MU_re)
  p <- upset(plot_data, sets = c("KS", "Wilcox", "Moses", "MU"),
             mb.ratio = c(0.6, 0.4), order.by = "freq",
             queries = list(list(query=intersects,
                                 params=list("KS", "Wilcox", "Moses", "MU"),
                                 color="#DD3737", active=T)),
             sets.bar.color = c("#b681bd","#377eb8","#7cc47a","#e06868"),
             nsets = 7, number.angles = 0, point.size = 6, line.size = 1,
             mainbar.y.label = "overlap gene counts",
             sets.x.label = "gene counts", text.scale = c(1, 1,1))
  return(p)
}



# PALdsa ------------------------------------------------------------------
#' PALdsa
#'
#' @description Analysis of significant difference of tail length under
#'   different conditions.
#' @details The poly(A) tail length of genes with two PACs was analyzed for
#'   significant difference.It is feasible to directly analyze whether there is
#'   a significant difference in tail length between the two PACs, or to
#'   annotate PACs as the proximal PA site and the distal PA site before
#'   conducting the significance analysis of tail length difference.
#' @param PAdf A dataframe of PA infomation without merge.
#' @param PALdf A dataframe that contains the length information of all reads
#'   tails.
#' @param gff Genome annotation stored in a GFF/GTF file or a TXDB R object can
#'   be used for annotating PACs. Please refer to movAPA for details.
#' @param d distance to group nearby PACds, default is 24 nt.
#' @param mode The PAL comparison mode, "PD" refers to comparing the tail
#'   length difference between the proximal and distal PA sites of each gene
#'   with two PACs."2PA" represents a direct comparison of tail length
#'   differences between the two PA loci of genes with two PACs, default for
#'   "PD".
#' @param withViolinPlot Logical value variable, Whether to draw a violin
#'   diagram, see function DSA_ViolinPlot() for details, default is true.
#' @param withUpsetPlot Logical value variable, Whether to draw a upset
#'   diagram, see function DSA_UpsetPlot() for details, default is true.
#' @param SAoDMethod one of "KS", "MWU", "ME", "Wilcox" and "ALL".
#' @return By default, this function returns a list of two graphs and a
#'   dataframe, or only a dataframe if no drawing operation has been performed.
#'   The result dataframe contains 7 column contents, the first column is
#'   geneID, and the second and third columns are the median tail lengths under
#'   different conditions. The remaining four columns are p-values calculated by
#'   different difference significance test methods, each column corresponds to
#'   a difference significance test method, in which "lose" or "NAN" means that
#'   for some reason (probably too little data) this method has not been able to
#'   calculate the exact p-values.
#' @examples
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' gff <- parseGenomeAnnotation(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' load("./data/AnnotedTails.RData")
#' PAs <- read.table("./inst/extdata/output/PAs/PAs.txt",header=T,sep=" ")
#' diffPAL2PAgenes <- PALdsa(PAs,AnnotedTails,gff,mode="PD",SAoDMethod="ME",withViolinPlot=T,withUpsetPlot=F)
#' @export
PALdsa <- function(PAdf,PALdf,gff,d,mode,withViolinPlot,withUpsetPlot,SAoDMethod){
  #check the paraments
  if(missing(PAdf) | missing(PALdf)){
    stop("Missing important parameters: PAdf or PALdf!")
  }
  if(missing(gff)){
    stop("Missing the annotations file!please check!")
  }
  if(missing(d)){
    d = 24
  }
  if(missing(SAoDMethod)){
    SAoDMethod = "ALL"
  }
  if(!is.integer(d) & !is.numeric(d)){
    stop("Parameter d must be an integer!")
  }
  if(missing(mode)){
    mode = "PD"
  }
  if(!(mode %in% c("PD","2PA"))){
    stop("The mode parameter must be one of c('PD','2PA')!")
  }
  if(missing(withViolinPlot)){
    withViolinPlot = T
  }
  if(missing(withUpsetPlot)){
    withUpsetPlot = T
  }
  #functions
  if(mode == "PD"){
    all_re = diffPALdpPA(PAdf,PALdf,gff,d=d)
  }
  if(mode == "2PA"){
    all_re = diffPAL2PAgene(PAdf,PALdf,gff,d=d)
  }
  re <- list()
  re$DSAresult  = all_re$re2
  if(withViolinPlot){
    p1 = DSA_ViolinPlot(all_re,SAoDMethod)
    re$ViolinPlot  = p1
    plot(p1)
  }
  if(withUpsetPlot){
    p2 = DSA_UpsetPlot(all_re$re2)
    re$UpsetPlot  = p2
  }
  return(re)
}

# batchCompare -----------------------------------------------------------
#' batchCompare
#'
#' @description The batchCompare function can be used to compare the calculation
#'   results of multiple different batches of data.
#' @details The batchCompare function provides the ability to quickly compare
#'   and visualize tail quantification results from multiple different batches
#'   of data. Users can use this function to compare computational results from
#'   multiple different datasets, such as one from NGS sequencing technology and
#'   another from PACBIO sequencing technology. batchCompare provides three
#'   comparable dimensions, namely "gene", "tail", "PACds" or "read".
#' @param datalist Multiple sets of analysis results stored in list format.
#' @param format Plot format, "upset" or "veen" or "bar".
#' @param dimension The dimension to compare, "gene", "tail", "PACds" or "read".
#' @param mycolors Drawing color card, must be specified, no default value. The
#'   color number is the same as the number of data sets.
#' @param rep Logical variable. When plotting the comparison of poly(A) tail
#'   length distribution for different batches of data, the parameter 'rep' can
#'   be used to specify whether to display data from different biological
#'   replicates separately. If the parameter 'rep' is true, the data of
#'   different repetitions are displayed separately, otherwise they are not
#'   separated.
#'  @param findOvpPACds Logical variable. Indicates whether to output
#'    overlapping PACds, defaults to 'False'.
#'  @param annotateByKnownPAC Logical variable. annotateByKnownPAC returns
#'    overlapping status and min distance with given PACdatasets.
#' @return Return the comparison result graph.
#' @examples
#' batchCompare(datalist,"upset","read",mycolors=c("#b681bd","#377eb8"))
#' @export
batchCompare <- function(datalist,format,dimension,mycolors,rep,findOvpPACds,annotateByKnownPAC,d){
  if(missing(datalist)){
    stop("missing input data!")
  }
  if(dimension != "tail"){
    if(missing(format)){
      print("The drawing method is not specified, the default selection is the upset drawing.")
      format = "upset"
    }
  }
  if(missing(dimension)){
    print("No comparison dimension is specified, gene is selected by default.")
  }
  if(missing(mycolors)){
    stop("Color card must be specified!")
  }
  if(dimension != "tail"){
    if(!(format %in% c("upset","veen","bar"))){
      stop("Format must in c('upset','veen','bar')!")
    }
  }
  if(!(dimension %in% c("gene","tail","read","PACds"))){
    stop("Dimension must in c('gene','tail','PACds','read')!")
  }
  if(dimension == "tail" & missing(rep)){
    rep = F
  }
  if(dimension=="PACds"){
    if(missing(findOvpPACds)){
      findOvpPACds=F
    }
    if(missing(annotateByKnownPAC)){
      annotateByKnownPAC=F
    }
    if(missing(d)){
      d=0
    }
  }
  #STEP1:initialize data
  N = length(datalist)
  dataNumber = 1:N
  #STEP2:plot
  if(dimension=="gene"){
    dataToPlot = {}
    batchList = {}
    for(i in dataNumber){
      batchLab = paste("Batch",i,sep = "")
      batchList = append(batchList,batchLab)
      if(length(datalist[[1]][["a"]])){
        stop("Please check the input data!")
      }
      else{
        dataToPlot[[batchLab]] = na.omit(datalist[[i]][["gene"]])
      }
    }
    if(format=="upset"){
      p<- UpSetR::upset(UpSetR::fromList(dataToPlot),order.by = "freq",
                queries = list(list(query=UpSetR::intersects,
                                    params=batchList,
                                    color="#DD3737", active=T)),
                sets.bar.color = mycolors)
    }
    if(format=="veen"){
      venn.diagram(x =dataToPlot,
                   filename = "veen_plot.png",
                   height = 800, width= 800,
                   resolution =300,
                   imagetype="png",
                   lwd=0.6,
                   fill=mycolors,
                   cex=0.45,
                   cat.cex=0.45)
      print("The Venn diagram has been saved to the current location, named 'veen_plot.png'.")
    }
  }
  else if(dimension=="read"){
    dataToPlot = data.frame()
    if(length(datalist[[1]][["read_type"]])==0 | length(datalist[[1]][["tailType"]])==0){
      stop("Please check the input data!")
    }
    else{
      for(i in dataNumber){
        data = as.data.frame(table(datalist[[i]]$read_type))
        data$group = paste("Batch",i,sep = "")
        data$rt = data$Freq/(sum(data$Freq))
        data1 = as.data.frame(table(datalist[[i]]$tailType))
        #data1$group = paste("Batch",i,sep = "")
        #data1$rt = data1$Freq/(sum(data1$Freq))
        #dataToPlot = base::rbind(dataToPlot,data,data1)
      }
    }
    dataToPlot$type <- factor(dataToPlot$Var1,levels=c("one-tail","two-tail-same","two-tail-mixed","structural","unstructural"))
    dataToPlot <- dataToPlot %>% dplyr::mutate(Percentage=round(rt*100,2))
    p<-ggplot(dataToPlot,                                   
              aes(x = type,
                  y = Percentage,
                  fill = group)) +
      ggplot2::geom_bar(stat = "identity",
               position = "dodge",
               alpha=0.8)+
      scale_colour_manual(breaks = dataToPlot$group,values=c("#9bbfdc","#ab99c1"))+
      scale_fill_manual(values=c("#9bbfdc","#ab99c1"))+
      xlab("") +
      ylab("Percentage") +
      ggthemes::theme_base()+
      theme(axis.text.x = element_text(size = 10,color ="black"),
            axis.text.y = element_text(size = 10,color ="black"),
            legend.text=element_text(size=10),
            legend.position="right",
            axis.title=element_text(size=12,color ="black"),
            legend.title=element_text(size=12,color ="black"),
            plot.title = element_text(size = 12,color ="black"))+
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  }
  else if(dimension=="tail" & rep == T){
    dataToPlot = data.frame()
    for(i in dataNumber){
      data <- datalist[[i]] %>% dplyr::select(PAL,rep)
      data$sample <- paste("Batch",i,sep = "")
      dataToPlot <- base::rbind(dataToPlot,data)
    }
    p <- dataToPlot %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(PAL2 = PAL) %>%
      ggplot(aes(x=rep,y=PAL2))+
      geom_flat_violin(aes(fill=sample),position=position_nudge(x=.1),trim = FALSE,)+
      scale_fill_manual(values=mycolors)+
      scale_color_manual (values=mycolors)+
      scale_x_discrete(expand=c(0.8, 0))+
      geom_boxplot(width=.07)+
      coord_flip()+
      ggthemes::theme_base()+
      ylab('Poly(A) tial length')+
      xlab('')+
      facet_wrap(.~sample, scales = "free")
  }
  else if(dimension=="tail" & rep == F){
    dataToPlot = data.frame()
    for(i in dataNumber){
      data <- datalist[[i]] %>% dplyr::select(PAL)
      data$sample <- paste("Batch",i,sep = "")
      dataToPlot = base::rbind(dataToPlot,data)
    }
    p <- dataToPlot %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(PAL2 = PAL) %>%
      ggplot(aes(x="",y=PAL2))+
      PupillometryR::geom_flat_violin(aes(fill=sample),position=position_nudge(x=.1),trim = FALSE,)+
      scale_fill_manual(values=mycolors)+
      scale_color_manual (values=mycolors)+
      scale_x_discrete(expand=c(0.8, 0))+
      geom_boxplot(width=.07)+
      coord_flip()+
      ggthemes::theme_base()+
      ylab('Poly(A) tial length')+
      xlab('')+
      facet_wrap(.~sample, scales = "free")
  }
  else if(dimension=="PACds"){
    p = {}
    p[["OvpPACds"]] = countOvpPACds(datalist[[1]], datalist[[2]], mode1='point', mode2='point', d=d, pats1=1, pats2=1)
    if(findOvpPACds){
      # p[["OvpPACdsList"]] = findOvpPACds(datalist[[1]], datalist[[2]], mode1='point', mode2='point', d=d, pats1=1, pats2=1)
      p[["OvpPACdsList"]] = findOvpPACds(datalist[[1]], datalist[[2]], d=d)
    }
    if(annotateByKnownPAC){
      #p[["OvpPACdsAnnotate"]] = annotateByKnownPAC(datalist[[1]], datalist[[2]], mode1='point', mode2='point', d=d, pats1=1, pats2=1)
      p[["OvpPACdsAnnotate"]] = annotateByKnownPAC(datalist[[1]], datalist[[2]], d=d,labels=1)
    }
  }
  else{
    stop("Dimension must in c('gene','tail','PACds','read')!")
  }
  return(p)
}