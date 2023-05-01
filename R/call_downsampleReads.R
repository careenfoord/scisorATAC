#' Downsamples number of reads per exon in CasesVControls output: cases_vs_controls.counts.passed-chi-sq-crit.tab.gz
#'
#' This function converts input temperatures in Fahrenheit to Celsius.
#' @param Num_Downsampled_Reads Num of minimum reads per exon and how many reads will be downsampled.
#'
#' @return Sampling_DPSI_Table.tab
#' @return Correlation plot
#'
#'
#' @export
#'
#'
############ Downsample exons to having x number of reads

downsampleReads <- function(Num_Downsampled_Reads = 10, example=FALSE) {
  file2 <- system.file("R_scripts", "downsampling_reads.R", package = "scisorATAC")
  string2 <- paste("Rscript", file2, Num_Downsampled_Reads, example)
  system(string2)
}
