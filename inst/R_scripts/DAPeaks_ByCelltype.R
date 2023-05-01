args <- commandArgs(trailingOnly=TRUE)
ATACobj <- args[1]
annotation.gr <- args[2]
AssayName <- args[3]
condition.query <- args[4]
celltypeA <- args[5]
celltypeB <- args[6]
cellnum <- args[7]
peaknum <- args[8]
MinCellRatio <- args[9]
random.repeats <- args[10]
harmony <- args[11]
outputDir <- args[12]
savePeakRobj <- args[13]

DAPeaksByCelltype <- function(ATACobj, annotation.gr, AssayName, condition.query, celltypeA, celltypeB, cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj)
{
  library(Signac)
  library(Seurat)
  dir.create(outputDir)
  ###### set the random subsampling name #####
  rand.times <- c(1:random.repeats)
  rand.version <- paste0(c("rand.V"),rand.times)

  Seurat::DefaultAssay(ATACobj) <- AssayName
  ##### subset the query condition ########
  ATACobj <- subset(ATACobj, subset = condition == condition.query)
 ######### subset the Robj of specific celltype A and celltype B for comparison ##########
 ATACobj.A <- subset(ATACobj,subset = celltype == celltypeA)
 ATACobj.B <- subset(ATACobj,subset = celltype == celltypeB)

celltypeA.total <-  ncol(ATACobj.A)
celltypeB.total <-  ncol(ATACobj.B)

###### print out cell numbers of query celltypes
print(paste0(celltypeA," ",celltypeA.total))
print(paste0(celltypeB," ",celltypeB.total))

rowname.listA <- as.data.frame(rownames(ATACobj.A@meta.data))
rowname.listB <- as.data.frame(rownames(ATACobj.B@meta.data))
#### remove tmp files
rm(ATACobj.A,ATACobj.B)

######if number of celltype A and celltype B exceed the random sampling size cellnum ######
 if (celltypeA.total>cellnum & celltypeB.total>cellnum)
 {
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
  group.by = "celltype"
)

macs.peaks.rand.chr <- GenomeInfoDb::keepStandardChromosomes(macs.peaks.rand, pruning.mode = "coarse")

# create a new assay using the MACS2 peak set and add it to the Seurat object
macs2.counts.ATACobj.rand <- Signac::FeatureMatrix(
  fragments = Fragments(ATACobj.rand),
  features = macs.peaks.rand.chr,
  cells = colnames(ATACobj.rand)
)

#### if annotation file is supplied ########
if (file.exists("annotation.gr") == TRUE)
{
  ATACobj.rand[["peaks"]] <- Signac::CreateChromatinAssay(
    counts = macs2.counts.ATACobj.rand,
    annotation = annotation.gr
  )
}

#### if annotation file is NULL ########
if (file.exists("annotation.gr") == FALSE)
{
  ATACobj.rand[["peaks"]] <- Signac::CreateChromatinAssay(
    counts = macs2.counts.ATACobj.rand
  )
}

### run harmony if harmony = TRUE
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
##########################################################
##### Normalization and linear dimensional reduction #####
Seurat::DefaultAssay(ATACobj.rand) <- "peaks"
ATACobj.rand <- Signac::FindTopFeatures(ATACobj.rand, min.cutoff = 5)
ATACobj.rand <- Signac::RunTFIDF(ATACobj.rand)
ATACobj.rand <- Signac::RunSVD(ATACobj.rand)

######## save signac obj ##################3
if (savePeakRobj == TRUE)
{
  save(ATACobj.rand, file = paste0(outputDir,"/",rand.version[k],"_",condition.query,"_",celltypeA,".VS.",celltypeB,"_Signac.Robj"))
}

#### export annotation information of all subsampled peaks ####
rand.gr <- Signac::granges(ATACobj.rand)
#### save granges file
save(rand.gr,file = paste0(outputDir,"/",rand.version[k],"_",condition.query,"_",celltypeA,".VS.",celltypeB,"_all.peaks.granges.Robj"))

###### subsampling random peaks #######
rand.peaks <- as.data.frame(rand.gr)
rand.peaks$peak.region <- paste0(rand.peaks$seqnames, "-", rand.peaks$start,"-", rand.peaks$end)
rand.subsample.peaks <- rand.peaks[sample(nrow(rand.peaks), peaknum), ]
write.csv(rand.subsample.peaks, file = paste0(outputDir,"/",rand.version[k],"_",condition.query,"_",celltypeA,".VS.",celltypeB,".subsampled.peaks.gr.csv"),row.names= F,quote = F)

###### find differential accessible peaks among random subsampled peaks by comparing celltypeA and celltypeB ##########
print(paste0("Find differential accessible peaks between.",celltypeA," and ",celltypeB))

Seurat::Idents(ATACobj.rand) <- ATACobj.rand$celltype
#### features set to be the subsampled peaks, only peaks detected in MinCellRatio cells are considered for test, the test method is set to be logistic regression #####
rand.da.peaks <- Seurat::FindMarkers(
  object = ATACobj.rand,
  ident.1 = celltypeA,
  ident.2 = celltypeB,
  features = rand.subsample.peaks$peak.region,
  min.pct = MinCellRatio,
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 0
)

rand.da.peaks$peakID <- rownames(rand.da.peaks)
#### save the DA peaks test result ########
write.csv(rand.da.peaks, file = paste0(outputDir,"/",rand.version[k],"_",condition.query,"_",celltypeA,".VS.",celltypeB,".DA.peaks.csv"),row.names= F,quote = F)
rm(ATACobj.rand,rand.da.peaks,rand.gr,rand.subsample.peaks,rand.peaks,rand.cellA,rand.cellB,rand.cells)
## remove seed ####
rm(.Random.seed, envir=globalenv())
print(paste0("Random subsampling V", k ," is done"))
   }

print(paste0("All Random subsamplings are done"))
######## perform stats of subsampling results ########
print("Starting stats of tested peaks")

### create a stats table
stats.group.name <- c("condition","celltype.comparison","sig.peaks","tested.peaks","sig.peak.pct")
peaks.stats <- data.frame(matrix(NA, nrow = random.repeats, ncol = length(stats.group.name)))
colnames(peaks.stats) <- stats.group.name

for (i in 1:random.repeats)
{
  dif.peaks.path <- paste0(outputDir,"/",rand.version[i],"_",condition.query,"_",celltypeA,".VS.",celltypeB,".DA.peaks.csv")
  dif.peaks <- read.csv(dif.peaks.path, sep = ",", header = T,na = "NA")

  ##### get significant peaks by filtering for p_val_adj < 0.05
  dif.sig.peaks <- dif.peaks[which(dif.peaks$p_val_adj < 0.05),]

  #col1: query celltype, col2: conditions for comparison, col3: No.of significant peaks, col4: No.of all tested peaks, col5: pct of significant peak
  peaks.stats[i,1] <- condition.query
  peaks.stats[i,2] <- paste0(celltypeA,".VS.",celltypeB)
  peaks.stats[i,3] <- nrow(dif.sig.peaks)
  peaks.stats[i,4] <- nrow(dif.peaks)
  peaks.stats[i,5] <- 100*nrow(dif.sig.peaks)/nrow(dif.peaks)
}

##### print out the stats
write.table(peaks.stats, file = paste0(outputDir,"/","Sig.Peak.Stats_",random.repeats,"random.subsamplings_",condition.query,"_",celltypeA,".VS.",celltypeB,".txt"),sep = "\t"  ,col.names = TRUE, row.names = FALSE, quote = FALSE)

}

####### if given subsampling size is larger than the cell number of either one of the query celltypes
if (celltypeA.total < cellnum | celltypeB.total < cellnum)
{
  print("subsampling size is larger than the total cell number")
}


}

DAPeaksByCelltype(ATACobj, annotation.gr, AssayName, condition.query, celltypeA, celltypeB, cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj)


