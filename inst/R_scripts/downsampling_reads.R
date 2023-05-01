

args <- commandArgs(trailingOnly=TRUE)
Num_Downsampled_Reads <- args[1]
example <- args[2]

if (example == TRUE){
  set.seed(7)
  cat("setting example seed","\n")
}else {
  rm(.Random.seed, envir=globalenv())
}

downsampling_reads <-  function(Num_Downsampled_Reads){

  Num_Downsampled_Reads <- as.numeric(Num_Downsampled_Reads)

  ## current path is ../ from celltype subfolders

  CellTypeFolder <- getwd()
  setwd(CellTypeFolder)
  for (dir in list.dirs()[-1]) {
    setwd(dir)
    #cat("dir:", dir,"\n")

    if (file.exists("cases_vs_controls.counts.passed-chi-sq-crit.tab.gz")){
      cat("dir:", dir,"\n")
      table <- read.table("cases_vs_controls.counts.passed-chi-sq-crit.tab.gz")

      sample1_included = c()
      sample1_excluded = c()
      sample2_included = c()
      sample2_excluded = c()

      table <- table[table$V2 + table$V3 >=Num_Downsampled_Reads, ]
      table <- table[table$V4 + table$V5 >=Num_Downsampled_Reads, ]
      for (row in 1:nrow(table)) {

        value1=table[row,2]
        value2=table[row,3]
        value3=table[row,4]
        value4=table[row,5]

        list1e <- rep(1, value1)
        list1i <- rep(0, value2)
        list1 <- c(list1e,list1i)

        if (length(list1) <Num_Downsampled_Reads){
          #sample2 <- sample(list1, 10, replace=TRUE)
          sample1 <-sample(list1, 0, replace=TRUE)
        }else {
          sample1 <- sample(list1, Num_Downsampled_Reads)

        }
        e1 <- length(which(sample1 ==1))
        i1 <- length(which(sample1 ==0))
        sample1_included <- append(sample1_included, e1)
        sample1_excluded <- append(sample1_excluded, i1)

        list2e <- rep(1, value3)
        list2i <- rep(0, value4)
        list2 <- c(list2e,list2i)
        if (length(list2)<Num_Downsampled_Reads){
          #sample2 <- sample(list2, 10, replace=TRUE)
          sample2 <-sample(list2, 0, replace=TRUE)
        }else {
          sample2 <- sample(list2, Num_Downsampled_Reads)
        }
        e2 <- length(which(sample2 ==1))
        i2 <- length(which(sample2 ==0))
        sample2_included <- append(sample2_included, e2)
        sample2_excluded <- append(sample2_excluded, i2)
      }

      table$sample1_included <- sample1_included
      table$sample1_excluded <- sample1_excluded
      table$sample2_included <- sample2_included
      table$sample2_excluded <- sample2_excluded

      table$DPSI_og <- (table$V2/(table$V2+table$V3)) - (table$V4/(table$V4+table$V5))
      table$DPSI_new <- (table$sample1_included/(table$sample1_included+table$sample1_excluded)) - (table$sample2_included/(table$sample2_included+table$sample2_excluded))

      pdf("CorrelationPlot.pdf")
      plot(table$DPSI_new, table$DPSI_og)
      dev.off()

      cor(table$DPSI_new, table$DPSI_og)
      cat("Corr:",cor(table$DPSI_new, table$DPSI_og), "\n")

      write.table(file = "Sampling_DPSI_Table.tab", x = table, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
      setwd(CellTypeFolder)
    } else {
      #cat("cases_vs_controls.counts.passed-chi-sq-crit.tab.gz not in folder", dir, "\n")
      setwd(CellTypeFolder)
    }
  }
}

downsampling_reads(Num_Downsampled_Reads)
