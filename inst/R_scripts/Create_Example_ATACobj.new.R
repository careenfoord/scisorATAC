args <- commandArgs(trailingOnly=TRUE)
example.data.path <-args[1]
outDir <- args[2]
harmony <- args[3]

CreateExampleATACobj <- function(example.data.path,outDir, harmony)
{
  library(Signac)
  library(Seurat)
### create an output directory #######
dir.create(outDir)
##### load annotation granges file ########
load(paste0(example.data.path,"/Macaca.mulatta.Mmul.10.104.granges.chr.annotation.NCBI.Robj"))

#### import the peaks.bed file as query features #######
all.atac.peaks.bed <- read.table(paste0(example.data.path,file = "/combined.7K.cells.ATAC.peaks.bed"), sep = "\t", header = T)
all.atac.peaks.bed <- all.atac.peaks.bed[,1:3]
colnames(all.atac.peaks.bed) <-  c("chr", "start", "end")
all.atac.peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(all.atac.peaks.bed)
combined.peaks <- Signac::reduce(x = all.atac.peaks.gr, drop.empty.ranges = FALSE)

###### four sample names ########
sample.names <- c("M1_PFC","M1_VIS","M2_PFC","M2_VIS")

############################ create chromatin assay for 4 samples ########################################

for (i in 1:length(sample.names))
{
  ##### set path to the fragment files and barcode metrics #########
  barcode_metrics.path <- paste0(example.data.path,"/",sample.names[i],c("_per_barcode_metrics.csv"))
  fragpath <- paste0(example.data.path,"/",sample.names[i],c("_atac_fragments.tsv.gz"))

  ## import the barcode_metrics ####
  metadata <- read.csv(
    file = barcode_metrics.path,
    header = TRUE,
    row.names = 1
  )

  options(warn=-1)

  # Create a Fragment object with the fragments file ####3
  frags <- Signac::CreateFragmentObject(
    path = fragpath
  )

  options(warn=0)

  ## Construct a feature x cell matrix from the fragments object #####
  counts <- Signac::FeatureMatrix(
    fragments = frags,
    features = combined.peaks
  )

  ### create chromatin assay ###
  chrom_assay <- Signac::CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = fragpath,
    ranges = combined.peaks,
    annotation = chr.annotation.ncbi,
  )

  ## create the object with the chromatin assay ########
  s <- Seurat::CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    meta.data = metadata,
    project = sample.names[i]
  )

  ##### calculate percent.mt and save in s
  s[["percent.mt"]] <- Seurat::PercentageFeatureSet(s, pattern = "^MT")

  ###### Computing QC Metrics #######
  Seurat::DefaultAssay(s) <- "ATAC"
  s <- Signac::NucleosomeSignal(s)
  s <- Signac::TSSEnrichment(s)

  #### assign the chromatin assay object for each sample #########
  assign(sample.names[i], s)
  ### clear the tmp variables #####
  rm(s,chrom_assay,frags,counts,barcode_metrics.path,fragpath,metadata)
}

##### assign condition information for each sample ######
M1_PFC$condition <- 'PFC'
M2_PFC$condition <- 'PFC'
M1_VIS$condition <- 'VIS'
M2_VIS$condition <- 'VIS'

###### merge all datasets, adding a sample ID to make sure cell names are unique ######
combined <- merge(
  x = M1_PFC,
  y = list(M2_PFC, M1_VIS, M2_VIS),
  add.cell.ids = c("M1PFC", "M2PFC", "M1VIS", "M2VIS")
)

### remove sample objects before merge #####
rm(list = sample.names)

######### Normalization and linear dimensional reduction #####
Seurat::DefaultAssay(combined) <- "ATAC"
combined <- Signac::RunTFIDF(combined)
combined <- Signac::FindTopFeatures(combined, min.cutoff = 'q0')
combined <- Signac::RunSVD(combined)

##### Run harmony for merged data if harmony == TRUE  #########
if (harmony == TRUE)
{
  combined.harmony <- harmony::RunHarmony(
    object = combined,
    group.by = 'orig.ident',
    reduction = 'lsi',
    assay.use = 'ATAC',
    project.dim = FALSE
  )
}


### import the cell assignment file #####
cta <- read.table(paste0(example.data.path,"/barcode2celltype.txt"), sep = "\t", header = T)

#### assign the cell type for each cell by barcode  ##########
combined$barcode <- rownames(combined@meta.data)
combined$celltype <- combined$barcode
combined$celltype <- plyr::mapvalues(x = combined$celltype, from = as.vector(cta$barcode), to = as.vector(cta$celltype))


#### check data structure ####
combined[['ATAC']]

print("summary of sample VS celltype")
table(combined$orig.ident, combined$celltype)

print("summary of condition VS celltype")
table(combined$condition, combined$celltype)

#### save the object ######
print("save the seurat object in outDir")
save(combined, file = paste0(outDir,"/combined.7K.ATAC.Robj"))

}

CreateExampleATACobj(example.data.path, outDir, harmony)
