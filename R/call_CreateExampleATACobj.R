#' create an example seurat object includes the chromatin assay with the fragment files and barcode metrics, condition and cell type per cell will be assigned.
#' the object wil be names as 'combined' and can be applied as input for the random subsampling
#'
#' .
#'
#' @param outDir the directory for the output file
#' @param harmony will perform RunHarmony if set to be TRUE, default = FALSE.
#' @return example atac object
#' @export
#'
#'
############ Downsample exons to having x number of reads


Create_Example_ATACobj <- function(example.data.path, outDir="OutputDir",harmony = FALSE) {
  file7 <- system.file("R_scripts", "Create_Example_ATACobj.new.R", package = "scisorATAC")
  string7 <- paste("Rscript", file7, example.data.path, outDir, harmony)
  system(string7)
}
