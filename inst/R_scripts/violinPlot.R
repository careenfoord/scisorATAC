#'
#'
`%>%` <- magrittr::`%>%`
violinPlot <- function(){

  CellTypeFolder <- getwd()
  setwd(CellTypeFolder)
  celltypes <- c()
  Dataframe1 <- data.frame(PercentSig= numeric(0), source= numeric(0))
  for (dir in list.dirs()[-1]) {
    setwd(dir)
    #cat("dir:", dir,"\n")
    if (file.exists("Percent_Sig_Downsampled.tab")){
      cat("dir:", dir,"\n")
      assign(paste0("Data"), read.table("Percent_Sig_Downsampled.tab"))

      setwd(CellTypeFolder)
      temp <- stringr::str_split(dir, "/")
      celltype <- temp[[1]][2]
      Data$source <- celltype
      Dataframe1 <- rbind(Dataframe1,Data)
      celltypes <- append(celltypes,celltype)
    } else {
      setwd(CellTypeFolder)
    }
  }


  pdf("Downsampling_ViolinPlot.pdf")

  a <- ggplot2::ggplot(Dataframe1,ggplot2::aes(x=source, y=V1, color=source)) + ggplot2::geom_violin() +
    ggplot2::theme_classic() + ggplot2::xlab("Cell Type")+ ggplot2::ylab("Sampled % Exons Significant") + ggplot2::ggtitle("Sampled % Exons Significant By Group") +
    ggplot2::theme(legend.position="none",axis.text= ggplot2::element_text(size=26,color = 'black', angle = 45, hjust = 1) ,axis.title= ggplot2::element_text(size=26), plot.title = ggplot2::element_text(size=26,hjust = 0.5))

  plot(a)
  dev.off()

  #Dataframe1$source <- as.factor(Dataframe1$source)
  #stat.test <- Dataframe1 %>% rstatix::wilcox_test(V1 ~ source)
  #stat.test

  #write.table(file = "Downsampling_wilcoxtest.tsv", x = stat.test, quote = FALSE, sep = "\t", row.names = FALSE)
}

violinPlot()
