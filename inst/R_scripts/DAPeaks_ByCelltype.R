args <- commandArgs(trailingOnly=TRUE)
ATACobj_path <- args[1]
AssayName <- args[2]
condition.query <- args[3]
celltypeA <- args[4]
celltypeB <- args[5]
cellnum <- args[6]
peaknum <- args[7]
MinCellRatio <- args[8]
random.repeats <- args[9]
harmony <- args[10]
outputDir <- args[11]
savePeakRobj <- args[12]
MACS2_path <- args[13]

DAPeaksByCelltype <- function(ATACobj_path,AssayName,condition.query, celltypeA, celltypeB, cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj, MACS2_path)
{
  ATACobj=readRDS(ATACobj_path)

  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  dir.create(outputDir)
  ###### set the random subsampling name #####
  rand.times <- c(1:random.repeats)
  rand.version <- paste0(c("rand.V"),rand.times)

  DefaultAssay(ATACobj) <- AssayName
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
 #if (celltypeA.total>cellnum & celltypeB.total>cellnum)
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
  group.by = "celltype",
	macs2.path = MACS2_path
)

macs.peaks.rand.chr <- GenomeInfoDb::keepStandardChromosomes(macs.peaks.rand, pruning.mode = "coarse")

# create a new assay using the MACS2 peak set and add it to the Seurat object
macs2.counts.ATACobj.rand <- Signac::FeatureMatrix(
  fragments = Fragments(ATACobj.rand),
  features = macs.peaks.rand.chr,
  cells = colnames(ATACobj.rand)
)

peak.assay <- Signac::CreateChromatinAssay(counts = macs2.counts.ATACobj.rand,annotation = annotation.hg38)

Dimnames <- peak.assay@data@Dimnames
peak.cells <- Dimnames[[2]]
ATACobj.rand <- subset(ATACobj.rand, cells = peak.cells)

	   
#### if annotation file is supplied ########
 ATACobj.rand[["peaks"]] <- Signac::CreateChromatinAssay(
    counts = macs2.counts.ATACobj.rand
  )

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
Seurat::DefaultAssay(ATACobj.rand) <- "peaks"
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
peaks.stats <- matrix(NA, nrow = as.numeric(random.repeats), ncol = 5)
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

#}

####### if given subsampling size is larger than the cell number of either one of the query celltypes
#if (celltypeA.total < cellnum | celltypeB.total < cellnum)
#{
#  print("subsampling size is larger than the total cell number")
#}


}

DAPeaksByCelltype(ATACobj_path, AssayName, condition.query, celltypeA, celltypeB, cellnum, peaknum, MinCellRatio, random.repeats, harmony, outputDir, savePeakRobj, MACS2_path)


