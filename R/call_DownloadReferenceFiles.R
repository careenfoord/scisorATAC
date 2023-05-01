#' Exons are also outputted into cases_vs_controls.counts.pVal.FDR.LOR.tab.gz
#'
#' @return downloaded example files
#'
#' @export
############ Downsample exons to having x number of reads

DownloadRefs <- function() {
  file1 <- system.file("bash", "DownloadReferenceFiles.sh", package = "scisorATAC")
  string1 <- paste("sh", file1)
  system(string1)
}
