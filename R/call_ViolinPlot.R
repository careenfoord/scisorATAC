#' Downsamples number of reads, exons, and then repeats to generate a distribution of significant events from CasesVControls output: cases_vs_controls.counts.passed-chi-sq-crit.tab.gz
#'
#' This function converts input temperatures in Fahrenheit to Celsius.

#' @return violin plot from output of downsampleExonsAndIterate.R
#' @export

############ Downsample exons to having x number of reads

ViolinPlot <- function() {
  file4 <- system.file("R_scripts", "violinPlot.R", package = "scisorATAC")
  string4 <- paste("Rscript", file4)
  system(string4)
}
