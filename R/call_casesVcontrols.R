#' This function performs the casesVcontrols exon analysis on 2 AllInfo.gz files.
#'
#'
#' Exons must pass chi sqaure crit. to be tested. Those exons are outputted into cases_vs_controls.counts.passed-chi-sq-crit.tab.gz with counts returned.
#' Exons are also outputted into cases_vs_controls.counts.pVal.FDR.LOR.tab.gz
#'
#' @param caseList file with path to "case" AllInfo.gz
#' @param controlList file with path to "control" AllInfo.gz
#' @param chrom_file file with list of chromosomes desired for testing
#' @param numThreads number of threads to be used
#' @param annotation_path path to gencode formatted annotation
#' @param ci_low min percent spliced inclusion considered; recommended starting value = 0.05
#' @param ci_upper max percent spliced inclusion considered; recommended starting value = 0.95
#' @param min_reads minimum number of reads for sum of 2 allInfos for a given exon
#' @param zipping_function command for unzipping files; Linux = "zcat"
#' @param OL_fraction the fraction of the reads for a given position must be either inclusion or exclusion; recommended starting value = 0.8
#' @param CellTypeFile file with list of all celltypes to be considered
#' @param OutputDir name of OutputDir. Default is OutputDir
#'
#' @return cases_vs_controls.counts.passed-chi-sq-crit.tab.gz : exon_gene  includedcases excludedcases includedcontrols excludedcontrols
#' @return cases_vs_controls.counts.pVal.FDR.LOR.tab.gz
#' @return all.altExons.matrix.tab.gz, all.altExons.tab.gz, all.internalExons.tab.gz, CASES.sampleBoth.inc.exc.tab.gz, CONTROLS.sampleBoth.inc.exc.tab.gz
#' @return LargeDPSI_exons above |DPSI| > 0.5 , MediumDPSI_exons above |DPSI| > 0.25, sampleBoth.inc.exc.tab. per cell types
#'
#'
#'
#' @export
############ 

casesVcontrols <- function(caseList,controlList,chrom_file,numThreads = 10,annotation_path,ci_low = 0.05,ci_high = 0.95,min_reads = 10,zipping_function = "zcat",OL_fraction = 0.8,CellTypeFile, OutputDir="OutputDir") {
  file1 <- system.file("bash", "casesVcontrols.sh", package = "scisorATAC")
  string1 <- paste("sh", file1, caseList,controlList,chrom_file,numThreads,annotation_path,ci_low,ci_high,min_reads,zipping_function,OL_fraction,CellTypeFile,OutputDir)
  system(string1)
}
