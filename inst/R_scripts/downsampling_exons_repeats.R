#' @export
#descriptor
#check dplyr

`%>%` <- magrittr::`%>%`
args <- commandArgs(trailingOnly=TRUE)
Num_Exons_Selected <- args[1]
Num_Repeats <- args[2]
example <-args[3]

if (example == TRUE){
  set.seed(7)
  cat("setting example seed", "\n")
}else {
  rm(.Random.seed, envir=globalenv())
}

downsampling_exons_repeats <- function(Num_Exons_Selected,Num_Repeats){



  #Num_Exons_Selected <- args[2]
  Num_Exons_Selected <- as.numeric(Num_Exons_Selected)
  #Num_Repeats <- args[3]
  Num_Repeats <- as.numeric(Num_Repeats)

  CellTypeFolder <- getwd()
  setwd(CellTypeFolder)
  for (dir in list.dirs()[-1]) {
    setwd(dir)
    #cat("dir:", dir,"\n")
    if (file.exists("Sampling_DPSI_Table.tab")){
      #cat("dir:", dir,"\n")
      table1 <- read.table("Sampling_DPSI_Table.tab")
      table1 <- table1 %>%
        tidyr::separate(V1, c("Chr", "start","end","Gene","Strand"), "_")
      BF = 0.05/Num_Exons_Selected
      cat("Num_Exons_Selected",Num_Exons_Selected,"BF:",BF,"\n")

      pval <- c()
      percent_sig <- c()

      #selecting 1 exon/gene with MAX DPSI
      #table1_1exonpergene <- table1 %>% group_by(Gene) %>% slice_max(order_by = abs(V11), n = 1)
      #table1_1exonpergene <- table1_1exonpergene %>% group_by(Gene) %>% slice_max(order_by = abs(V10), n = 1)
      #table1_1exonpergene <- table1_1exonpergene %>% group_by(Gene) %>% filter(n() <= 1)

      #selecting 1 exon/gene with Random DPSI
      table1_1exonpergene <- table1 %>% dplyr::group_by(Gene) %>% dplyr::sample_n(size = 1)

      #check gene # == exon #
      table1_1exonpergene_counts <- table1_1exonpergene %>% dplyr::group_by(Gene) %>% dplyr::summarise(n=dplyr::n())
      cat("genes:",dim(table1_1exonpergene_counts)[1],"exons:",sum(table1_1exonpergene_counts$n), "\n")

      my_range <- 1:Num_Repeats
      for (i in my_range) {
        pval <- c()

        table_sampled <- table1_1exonpergene[sample(nrow(table1_1exonpergene), Num_Exons_Selected),]
        for (row in 1:nrow(table_sampled)) {
          include1=as.numeric(table_sampled[row,10])
          exclude1=as.numeric(table_sampled[row,11])
          include2=as.numeric(table_sampled[row,12])
          exclude2=as.numeric(table_sampled[row,13])

          dataframe <- data.frame(c(include1,exclude1),c(include2,exclude2))
          fisher <- fisher.test(dataframe)
          pval <- append(pval, fisher[["p.value"]])
        }


        table_sampled$Pval <- pval
        num_sig_trial <- dim(table_sampled[table_sampled$Pval <=BF,])[1]
        percent_sig_trial <- (num_sig_trial/Num_Exons_Selected)*100

        percent_sig <- append(percent_sig, percent_sig_trial)
      }

      percent_sig <- as.data.frame(percent_sig)
      write.table(file = "Percent_Sig_Downsampled.tab",x=percent_sig, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      setwd(CellTypeFolder)
    } else {
      #cat("Sampling_DPSI_Table.tab not in folder", dir, "\n")
      setwd(CellTypeFolder)
    }
  }
}

downsampling_exons_repeats(Num_Exons_Selected,Num_Repeats)
