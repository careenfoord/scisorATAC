#'Calling differential accessible peaks by comparing two conditions of one specific cell type ####################
#'
#'
#' @param ATACobj_path the object with chromatin assay created.
#' @param AssayName The assay names of the chromatin assay, default = "ATAC"
#' @param condition.query  the condition in query.
#' @param celltypeA the cell type A for comparison.
#' @param celltypeB the cell type B for comparison.
#' @param cellnum the number of cells to be subsampled, default = 500.
#' @param peaknum the number of peaks to be subsampled, default = 5000.
#' @param MinCellRatio Only test peaks that are detected in a minimum fraction of MinCellRatio cells in either of the two conditions, the default setting is 0.02.
#' @param random.repeats Random subsampling times to be performed, default = 10.
#' @param outputDir Directory path to the output filee, which will be created if not existed.
#' @param savePeakRobj The peaks called by MACS2 for each subsampling will be stored as assay ‘peaks’. It will be saved as Robj for downstream analysis, default = FALSE.
#' @param harmony will perform RunHarmony if set to be TRUE, default = FALSE.
#' @param MACS2_path path to MACS2
#' 
#' @return downsampled ATAC cells
#' @export
#'

DAPeaks_ByCelltype <- function(ATACobj_path, AssayName = "ATAC", condition.query, celltypeA, celltypeB, cellnum = 500, peaknum = 5000, MinCellRatio = 0.02, random.repeats = 10, harmony = FALSE, outputDir, savePeakRobj = FALSE, MACS2_path) {
  file6 <- system.file("R_scripts", "DAPeaks_ByCelltype.R", package = "scisorATAC")
  string6 <- paste("Rscript", file6, ATACobj_path, AssayName,condition.query,celltypeA,celltypeB,cellnum,peaknum,MinCellRatio,random.repeats,harmony, outputDir,savePeakRobj,MACS2_path)
  system(string6)
}
