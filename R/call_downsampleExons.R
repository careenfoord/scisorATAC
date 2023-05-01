#' Downsamples number exons and then repeats to generate a distribution of significant events
#'
#'
#' .
#' @param Num_Exons_Selected Number of exons to select to calculate % significant
#' @param Num_Repeats Number of times to loop over randomly selecting exons
#' @return Percent_Sig_Downsampled.tab
#' @export
#'
#'
############ Downsample exons to having x number of reads


downsampleExonsAndIterate <- function(Num_Exons_Selected,Num_Repeats,example = FALSE) {
  file3 <- system.file("R_scripts", "downsampling_exons_repeats.R", package = "scisorATAC")
  string3 <- paste("Rscript", file3, Num_Exons_Selected,Num_Repeats, example)
  system(string3)
}
