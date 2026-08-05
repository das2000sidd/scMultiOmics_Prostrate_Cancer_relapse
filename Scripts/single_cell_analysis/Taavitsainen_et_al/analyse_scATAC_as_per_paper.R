setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


library(Seurat)
library(ggplot2)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)
library(dplyr)


setMethod(
  "seqinfo",
  signature(x = "ChromatinAssay"),
  function(x) {
    GenomeInfoDb::seqinfo(x@ranges)
  }
)


#dmso_count <- readMM("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/matrix.mtx.gz")
#dmso_peak <- read.table("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/peaks.bed.gz",sep="\t", header = F)

#dmso_barcodes <- readLines("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/barcodes.tsv.gz")
load_scATAC_sample <- function(matrix_path,
                               peaks_path,
                               barcodes_path,
                               fragments_path,
                               sample_id = "sample",
                               genome = "hg38") {
  
  # libraries
  library(Matrix)
  library(GenomicRanges)
  library(Signac)
  library(Seurat)
  
  # -------------------------
  # 1. Load count matrix
  # -------------------------
  counts <- readMM(matrix_path)
  
  # -------------------------
  # 2. Load peaks
  # -------------------------
  peaks <- read.table(peaks_path, sep = "\t", header = FALSE)
  
  granges <- GRanges(
    seqnames = peaks$V1,
    ranges = IRanges(start = peaks$V2, end = peaks$V3)
  )
  
  # -------------------------
  # 3. Load barcodes
  # -------------------------
  barcodes <- readLines(barcodes_path)
  colnames(counts) <- barcodes
  
  # -------------------------
  # 4. Create ChromatinAssay
  # -------------------------
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    ranges = granges,
    genome = genome,
    fragments = fragments_path
  )
  
  # -------------------------
  # 5. Create Seurat object
  # -------------------------
  seurat_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  # -------------------------
  # 6. Add metadata
  # -------------------------
  seurat_obj$Sample_Id <- sample_id
  
  return(seurat_obj)
}

dmso_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/fragments.tsv.gz",
  sample_id = "DMSO",
  genome = "hg38"
)

enz48_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/fragments.tsv.gz",
  sample_id = "DMSO",
  genome = "hg38"
)

resa_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/fragments.tsv.gz",
  sample_id = "DMSO",
  genome = "hg38"
)

resb_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/fragments.tsv.gz",
  sample_id = "DMSO",
  genome = "hg38"
)


# Compute QC metrics
# Nucleosome signal
compute_QC_metric <- function(seu_obj, annotations, blacklist_regions, genome = "hg38"){
  
  # fix genome metadata
  gr <- seu_obj[["peaks"]]@ranges
  genome(gr) <- genome
  seu_obj[["peaks"]]@ranges <- gr
  
  # fix annotation genome
  genome(annotations) <- genome
  Annotation(seu_obj) <- annotations
  
  # old QC functions
  # seu_obj <- ATACqc(seu_obj, verbose = TRUE)
  seu_obj <- NucleosomeSignal(seu_obj)
  seu_obj <- TSSEnrichment(seu_obj)
  #seu_obj <- CountFragments(seu_obj)
  # required for pct_reads_in_peaks
  #seu_obj <- CountFragments(seu_obj)
  #seu_obj <-  FractionCountsInRegion(seu_obj)
  
  # derived metrics (only if columns exist)
  if(all(c("peak_region_fragments", "passed_filters") %in% colnames(seu_obj@meta.data))){
    
    seu_obj$pct_reads_in_peaks <- 
      seu_obj$peak_region_fragments / seu_obj$passed_filters * 100
  }
  
  seu_obj$blacklist_ratio <- FractionCountsInRegion(
    object = seu_obj, 
    assay = "peaks",
    regions = blacklist_regions
  )
  
  return(seu_obj)
}


library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)

# -----------------------------
# Rebuild annotations
# -----------------------------

annotations <- GetGRangesFromEnsDb(ensdb_v98)

# Convert chromosome naming safely
seqlevelsStyle(annotations) <- "UCSC"

# Remove problematic genome metadata first
seqinfo(annotations) <- Seqinfo(
  seqnames = seqlevels(annotations),
  genome = "hg38"
)

# Set genome
genome(annotations) <- "hg38"


# -----------------------------
# Fix ChromatinAssay genome
# -----------------------------

# Check current seqlevels
print(seqlevels(granges(dmso_s))[1:10])
print(seqlevels(annotations)[1:10])


# Extract ranges
assign_annotation <- function(seu_obj, annotations){
  
  peak_ranges <- granges(seu_obj)
  
  
  # Replace genome metadata of assay
  seqinfo(peak_ranges) <- Seqinfo(
    seqnames = seqlevels(peak_ranges),
    genome = "hg38"
  )
  
  # Put corrected ranges back into ChromatinAssay
  DefaultAssay(seu_obj) <- "peaks"
  
  seu_obj[["peaks"]]@ranges <- peak_ranges
  
  # Verify
  genome(granges(seu_obj))
  
  # -----------------------------
  # Verify before assignment
  # -----------------------------
  
  print(genome(annotations))
  print(genome(granges(seu_obj)))
  
  print(all(seqlevels(granges(seu_obj)) %in% seqlevels(annotations)))
  
  
  # -----------------------------
  # Assign annotation
  # -----------------------------
  
  Annotation(seu_obj) <- annotations
  
  return(seu_obj)
  
}


dmso_s <- assign_annotation(dmso_s, annotations)
enz48_s <- assign_annotation(enz48_s, annotations)
resa_s <- assign_annotation(resa_s, annotations)
resb_s <- assign_annotation(resb_s, annotations)


frag_path_dmso <- "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/fragments.tsv.gz"
frag_path_resa <- "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/fragments.tsv.gz"
frag_path_resb <- "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/fragments.tsv.gz"
frag_path_enz48 <- "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/fragments.tsv.gz"


#Fragments(dmso_s) <- CreateFragmentObject(
#  path = frag_path_dmso,
#  cells = colnames(dmso_s)
#)

dmso_s <- compute_QC_metric(dmso_s, annotations, blacklist_regions)
enz48_s <- compute_QC_metric(enz48_s, annotations, blacklist_regions)
resa_s <- compute_QC_metric(resa_s, annotations, blacklist_regions)
resb_s <- compute_QC_metric(resb_s, annotations, blacklist_regions)


DensityScatter(dmso_s, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(enz48_s, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(resa_s, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
DensityScatter(resb_s, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(dmso_s, features = "TSS.enrichment")
VlnPlot(enz48_s, features = "TSS.enrichment")
VlnPlot(resa_s, features = "TSS.enrichment")
VlnPlot(resb_s, features = "TSS.enrichment")


dmso_s <- subset(
  x = dmso_s,
  subset = nCount_peaks > 3000 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

#hist(dmso_s_keep$nCount_peaks)
#hist(dmso_s_keep_f$nCount_peaks)

enz48_s <- subset(
  x = enz48_s,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 20000 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

#hist(enz48_s_keep$nCount_peaks)
#hist(enz48_s_keep_f$nCount_peaks)

resa_s <- subset(
  x = resa_s,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 20000 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)

#hist(resa_s_keep$nCount_peaks)
#hist(resa_s_keep_f$nCount_peaks)


resb_s <- subset(
  x = resb_s,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 25000 &
    blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 4
)


preprocess_data <- function(seurat_atac){
  seurat_atac <- RunTFIDF(seurat_atac)
  seurat_atac <- FindTopFeatures(
    seurat_atac,
    min.cutoff = "q50"   # top 50% most accessible/variable peaks
  )
  seurat_atac <- RunSVD(seurat_atac)
  seurat_atac <- RunUMAP(
    seurat_atac,
    reduction = "lsi",
    dims = 2:30
  )
  return(seurat_atac)
  
}

dmso_s <- preprocess_data(dmso_s)
enz48_s <- preprocess_data(enz48_s)
resa_s <- preprocess_data(resa_s)
resb_s <- preprocess_data(resb_s)

compute.gene.act <- function(seu_obj){
  
  gene.activities <- GeneActivity(
    seu_obj,
    extend.upstream = 2000
  )
  
  seu_obj[["ACTIVITY"]] <- CreateAssayObject(
    counts = gene.activities
  )
  
  DefaultAssay(seu_obj) <- "ACTIVITY"
  
  seu_obj <- NormalizeData(
    seu_obj, verbose=TRUE
  )
  
  return(seu_obj)
  
} 

dmso_s <- compute.gene.act(dmso_s)
enz48_s <- compute.gene.act(enz48_s)
resa_s <- compute.gene.act(resa_s)
resb_s <- compute.gene.act(resb_s)


dmso_s <- FindNeighbors(dmso_s, reduction = "lsi",dims = 2:30, verbose = TRUE)
enz48_s <- FindNeighbors(enz48_s, reduction = "lsi",dims = 2:30, verbose = TRUE)
resa_s <- FindNeighbors(resa_s, reduction = "lsi",dims = 2:30, verbose = TRUE)
resb_s <- FindNeighbors(resb_s, reduction = "lsi",dims = 2:30, verbose = TRUE)

dmso_s <- FindClusters( dmso_s, resolution = c(0.2, 0.5, 0.8), algorithm = 4, leiden_method = "igraph", verbose = TRUE)
enz48_s <- FindClusters( enz48_s, resolution = c(0.2, 0.5, 0.8), algorithm = 4, leiden_method = "igraph", verbose = TRUE)
resa_s <- FindClusters( resa_s, resolution = c(0.2, 0.5, 0.8), algorithm = 4, leiden_method = "igraph", verbose = TRUE)
resb_s <- FindClusters( resb_s, resolution = c(0.2, 0.5, 0.8), algorithm = 4, leiden_method = "igraph", verbose = TRUE, graph.name = "peaks_snn")


DimPlot(
  dmso_s,
  reduction = "umap",
  group.by = c("peaks_snn_res.0.8"),
  split.by = c("Sample_Id"),
  label = FALSE
)


DimPlot(
  enz48_s,
  reduction = "umap",
  group.by = c("peaks_snn_res.0.8"),
  split.by = c("Sample_Id"),
  label = FALSE
)

DimPlot(
  resa_s,
  reduction = "umap",
  group.by = c("peaks_snn_res.0.5"),
  split.by = c("Sample_Id"),
  label = FALSE
)

DimPlot(
  resb_s,
  reduction = "umap",
  group.by = c("peaks_snn_res.0.8"),
  split.by = c("Sample_Id"),
  label = FALSE
)

## Now assigning clusters from best resolution of scATAC
Idents(dmso_s) <- dmso_s@meta.data$peaks_snn_res.0.8
Idents(enz48_s) <- enz48_s@meta.data$peaks_snn_res.0.8
Idents(resa_s) <- resa_s@meta.data$peaks_snn_res.0.5
Idents(resb_s) <- resb_s@meta.data$peaks_snn_res.0.8



rna_clustered <- readRDS("Batch_corrected_clustered_scRNAseq_object.rds")

DimPlot(
  rna_clustered,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.8"),
  split.by = c("Sample_id"),
  label = FALSE
)

dmso_rna <- subset(rna_clustered, subset = Sample_id == "DMSO")
enz48_rna <- subset(rna_clustered, subset = Sample_id == "ENZ48")
resa_rna <- subset(rna_clustered, subset = Sample_id == "RESA")
resb_rna <- subset(rna_clustered, subset = Sample_id == "RESB")


DimPlot(
  dmso_rna,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.6"),
  label = FALSE
) ## 0.6 looks okay

DimPlot(
  enz48_rna,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.6"),
  label = FALSE
) ## 0.6 looks okay

DimPlot(
  resa_rna,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.6"),
  label = FALSE
) ## 0.6 looks okay


DimPlot(
  resb_rna,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.6"),
  label = FALSE
) ## 0.6 looks okay

## Now assigning clusters from best resolution of scATAC
Idents(dmso_rna) <- dmso_rna@meta.data$RNA_snn_res.0.6
Idents(enz48_rna) <- enz48_rna@meta.data$RNA_snn_res.0.6
Idents(resa_rna) <- resa_rna@meta.data$RNA_snn_res.0.6
Idents(resb_rna) <- resb_rna@meta.data$RNA_snn_res.0.6


library(clustree)
clustree(
  resb_rna@meta.data,
  prefix = "RNA_snn_res."
)

annotate.atac.with.rna.cluster <- function(seu_rna, seu_atac){
  
  anchors.atac <- FindTransferAnchors(
    reference = seu_rna,
    query = seu_atac,
    reference.assay = "RNA",
    query.assay = "ACTIVITY",
    reduction = "cca",
    dims = 1:30,
    verbose = TRUE
  )
  
  predictions.atac <- TransferData(
    anchorset = anchors.atac,
    refdata = Idents(seu_rna),
    weight.reduction = "cca.l2",
    dims = 1:30,
    verbose = TRUE
  )
  
  
  #seu_atac@meta.data$Predicted_cluster <- predictions.atac
  
  return(predictions.atac)

}

#anchors.dmso <- FindTransferAnchors(
#  reference = dmso_rna,
#  query = dmso_s,
#  reference.assay = "RNA",
#  query.assay = "ACTIVITY",
#  reduction = "cca",
#  dims = 1:30,
#  verbose = TRUE
#)

#predictions.dmso.atac <- TransferData(
#  anchorset = anchors.dmso,
#  refdata = Idents(dmso_rna),
#  weight.reduction = "cca.l2",
#  dims = 1:30,
#  verbose = TRUE
#  )

dmso_s_pred_from_rna <- annotate.atac.with.rna.cluster(dmso_rna,dmso_s) ## original res is 0.8

stopifnot(rownames(dmso_s_pred_from_rna)==rownames(dmso_s@meta.data))


dmso_s@meta.data$Predicted_cluster <- dmso_s_pred_from_rna$predicted.id


table(
  dmso_s_annotated@meta.data$peaks_snn_res.0.8,
  dmso_s@meta.data$Predicted_cluster
)


df.atac.rna.clusters <- as.data.frame(table(
  dmso_s_annotated@meta.data$peaks_snn_res.0.8,
  dmso_s@meta.data$Predicted_cluster
))

colnames(df.atac.rna.clusters) <- c("ATAC", "RNA", "Overlap")

df.atac.rna.clusters.grouped.by.atac.c <- df.atac.rna.clusters %>%
  group_by(ATAC) %>%
  mutate(
    proportion = Overlap / sum(Overlap)
  )

library(tidyr)

df.atac.rna.clusters.grouped.by.atac.c <- as.data.frame(df.atac.rna.clusters.grouped.by.atac.c)

heatmap_matrix <- df.atac.rna.clusters.grouped.by.atac.c %>%
  dplyr::select(ATAC, RNA, proportion) %>%
  pivot_wider(
    names_from = RNA,
    values_from = proportion,
    values_fill = 0
  )


hm <- as.matrix(heatmap_matrix[, -1])

rownames(hm) <- paste0(
  "scATAC_cluster ",
  heatmap_matrix$ATAC
)

colnames(hm) <- paste0(
  "scRNA_cluster ",
  colnames(hm)
)

row_annotation <- data.frame(
  Modality = rep("scATAC", nrow(hm))
)

rownames(row_annotation) <- rownames(hm)



col_annotation <- data.frame(
  Type = "scRNA"
)

rownames(col_annotation) <- colnames(hm)

p <- pheatmap(
  hm,
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f"
)

pdf(
  "Heatmap_of_proportion_of_DMSO_scATAC_predicted_to_belong_to_scRNA_clusters.pdf",
  width = 8,
  height = 5
)

grid::grid.draw(p$gtable)

dev.off()



