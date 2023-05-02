
args <- commandArgs(trailingOnly=TRUE)
ATACobj_path <- args[1]
AssayName <- args[2]
celltype.query <- args[3]
conditionA <- args[4]
conditionB <- args[5]
cellnum <- args[6]
peaknum <- args[7]
MinCellRatio <- args[8]
random.repeats <- args[9]
harmony <- args[10]
outputDir <- args[11]
savePeakRobj <- args[12]
MACS2_path <- args[13]

DAPeaksByCondition <- function(ATACobj_path, AssayName, celltype.query, conditionA, conditionB , cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj,MACS2_path)
{
  ATACobj=readRDS(ATACobj_path)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  dir.create(outputDir)
  Seurat::DefaultAssay(ATACobj) <- AssayName
  ###### set the random subsampling name #####
  rand.times <- c(1:random.repeats)
  rand.version <- paste0(c("rand.V"),rand.times)
############ subset the query cell type ###############
  ATACobj <- subset(ATACobj, subset = celltype == celltype.query)
 ######### subset the Robj of specific conditionA and celltype B for comparison ##########
 ATACobj.A <- subset(ATACobj,subset = condition == conditionA)
 ATACobj.B <- subset(ATACobj,subset = condition == conditionB)

conditionA.total <-  ncol(ATACobj.A)
conditionB.total <-  ncol(ATACobj.B)

###### print out cell numbers of query celltypes
print(paste0(conditionA," ",conditionA.total))
print(paste0(conditionB," ",conditionB.total))

rowname.listA <- as.data.frame(rownames(ATACobj.A@meta.data))
rowname.listB <- as.data.frame(rownames(ATACobj.B@meta.data))
#### remove tmp files
rm(ATACobj.A,ATACobj.B)

######if number of celltype A and celltype B exceed the random sampling size cellnum ######
# if (conditionA.total>cellnum & conditionB.total>cellnum)
 #{

   for (k in 1:random.repeats)
   {
 ######## get random subsampled cells with a specified number cellnum ########
rand.cellA <- c(rowname.listA[sample(nrow(rowname.listA), cellnum), ])
rand.cellB <- c(rowname.listB[sample(nrow(rowname.listB), cellnum), ])

### subsample cells and subset the signac obj ####
 rand.cells <- as.data.frame(c(rand.cellA, rand.cellB))
 colnames(rand.cells) <- c("cell.barcode")
ATACobj.rand <- subset(ATACobj, cells = rand.cells$cell.barcode)

Seurat::DefaultAssay(ATACobj.rand) <- AssayName

#### calling peaks with MACS2
macs.peaks.rand <- Signac::CallPeaks(
  object = ATACobj.rand,
  group.by = "condition",
  macs2.path = MACS2_path
)

macs.peaks.rand.chr <- GenomeInfoDb::keepStandardChromosomes(macs.peaks.rand, pruning.mode = "coarse")

# create a new assay using the MACS2 peak set and add it to the Seurat object
macs2.counts.ATACobj.rand <- Signac::FeatureMatrix(
  fragments = Fragments(ATACobj.rand),
  features = macs.peaks.rand.chr,
  cells = colnames(ATACobj.rand)
)


ATACobj.rand[["peaks"]] <- Signac::CreateChromatinAssay(
  counts = macs2.counts.ATACobj.rand,
  annotation = annotation.gr
)


### run harmony if set harmony = TRUE
if (harmony == TRUE)
{
ATACobj.rand <- harmony::RunHarmony(
  object = ATACobj.rand,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)
}

######## save signac obj ##################3
if (savePeakRobj == TRUE)
{
  save(ATACobj.rand, file = paste0(outputDir,"/",rand.version[k],"_",celltype.query,"_",conditionA,".VS.",conditionB,"_Signac.Robj"))
}
#### export annotation information of all subsampled peaks ####
rand.gr <- Signac::granges(ATACobj.rand)
#### save granges file
save(rand.gr,file = paste0(outputDir,"/",rand.version[k],"_",celltype.query,"_",conditionA,".VS.",conditionB,"_all.peaks.granges.Robj"))

###### subsampling random peaks #######
rand.peaks <- as.data.frame(rand.gr)
rand.peaks$peak.region <- paste0(rand.peaks$seqnames, "-", rand.peaks$start,"-", rand.peaks$end)
rand.subsample.peaks <- rand.peaks[sample(nrow(rand.peaks), peaknum), ]
write.csv(rand.subsample.peaks, file = paste0(outputDir,"/",rand.version[k],"_",celltype.query,"_",conditionA,".VS.",conditionB,"_subsampled.peaks.gr.csv"),row.names= F,quote = F)

###### find differential accessible peaks among random subsampled peaks by comparing conditionA and conditionB ##########
print(paste0("Find differential accessible peaks between.",conditionA," and ",conditionB))
ATACobj.rand$celltype.prefix <- paste0(ATACobj.rand$condition, c("."),ATACobj.rand$celltype)
Seurat::Idents(ATACobj.rand) <- ATACobj.rand$celltype.prefix
ident_A <- paste0(conditionA,c("."),celltype.query)
ident_B <- paste0(conditionB,c("."),celltype.query)

#### features set to be the subsampled peaks, only peaks detected in MinCellRatio cells are considered for test, the test method is set to be logistic regression #####
rand.da.peaks <- Seurat::FindMarkers(
  object = ATACobj.rand,
  ident.1 = ident_A,
  ident.2 = ident_B,
  features = rand.subsample.peaks$peak.region,
  min.pct = MinCellRatio,
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 0
)

rand.da.peaks$peakID <- rownames(rand.da.peaks)
#### save the DA peaks test result ########
write.csv(rand.da.peaks, file = paste0(outputDir,"/",rand.version[k],"_",celltype.query,"_",conditionA,".VS.",conditionB,".DA.peaks.csv"),row.names= F,quote = F)
rm(ATACobj.rand,rand.da.peaks,rand.peaks,rand.subsample.peaks,rand.gr,rand.cellA,rand.cellB,rand.cells)
## remove seed ####
rm(.Random.seed, envir=globalenv())
print(paste0(Sys.time()," Random subsampling V", k ," is done"))
   }

   print(paste0(Sys.time(), "All random subsamplings are done"))

   ######## perform stats of subsampling results ########
   print("Starting stats of tested peaks")

   ### create a stats table
   stats.group.name <- c("celltype","condition.comparison","sig.peaks","tested.peaks","sig.peak.pct")
   peaks.stats <- matrix(NA, nrow = as.numeric(random.repeats), ncol = 5)
   colnames(peaks.stats) <- stats.group.name

   for (i in 1:random.repeats)
   {
   dif.peaks.path <- paste0(outputDir,"/",rand.version[i],"_",celltype.query,"_",conditionA,".VS.",conditionB,".DA.peaks.csv")
   dif.peaks <- read.csv(dif.peaks.path, sep = ",", header = T,na = "NA")

   ##### get significant peaks by filtering for p_val_adj < 0.05
   dif.sig.peaks <- dif.peaks[which(dif.peaks$p_val_adj < 0.05),]

   #col1: query celltype, col2: conditions for comparison, col3: No.of significant peaks, col4: No.of all tested peaks, col5: pct of significant peak
   peaks.stats[i,1] <- celltype.query
   peaks.stats[i,2] <- paste0(conditionA,".VS.",conditionB)
   peaks.stats[i,3] <- nrow(dif.sig.peaks)
   peaks.stats[i,4] <- nrow(dif.peaks)
   peaks.stats[i,5] <- 100*nrow(dif.sig.peaks)/nrow(dif.peaks)
 }

##### print out the stats
write.table(peaks.stats, file = paste0(outputDir,"/","Sig.Peak.Stats_",random.repeats,"random.subsampling_",celltype.query,"_",conditionA,".VS.",conditionB,".txt"), sep = "\t" ,col.names = TRUE, row.names = FALSE, quote = FALSE)

#}

####### if given subsampling size is larger than the cell number of either one of the query celltypes
#if (conditionA.total < cellnum | conditionB.total < cellnum)
#{
#  print("subsampling size is larger than total cell number")
#}

}

DAPeaksByCondition(ATACobj_path, AssayName, celltype.query, conditionA, conditionB, cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj,MACS2_path)

