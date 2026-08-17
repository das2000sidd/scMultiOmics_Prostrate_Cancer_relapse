setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


library(Seurat)
library(ggplot2)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(Matrix)
library(dplyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
library(ChIPpeakAnno)


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
  sample_id = "ENZ48",
  genome = "hg38"
)

resa_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/fragments.tsv.gz",
  sample_id = "RESA",
  genome = "hg38"
)

resb_s <- load_scATAC_sample(
  matrix_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/matrix.mtx.gz",
  peaks_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/peaks.bed.gz",
  barcodes_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/barcodes.tsv.gz",
  fragments_path = "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/fragments.tsv.gz",
  sample_id = "RESB",
  genome = "hg38"
)


frags.dmso <- fread(
  "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/DMSO/fragments.tsv.gz",
  header = FALSE
)

frags.enz48 <- fread(
  "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/ENZ48/fragments.tsv.gz",
  header = FALSE
)

frags.resa <- fread(
  "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESA/fragments.tsv.gz",
  header = FALSE
)

frags.resb <- fread(
  "/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168667_RAW_ATAC/RESB/fragments.tsv.gz",
  header = FALSE
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

annotations <- keepStandardChromosomes(
  annotations,
  pruning.mode = "coarse"
);

seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# add the gene information to the object
add_genome_to_peaks <- function(seu_obj){
  
  peaks <- granges(seu_obj[["peaks"]])
  
  genome(peaks) <- "hg38"
  
  seu_obj[["peaks"]]@ranges <- peaks
  
  genome(granges(seu_obj[["peaks"]]))
  
  Annotation(seu_obj[["peaks"]]) <- annotations
  
  return(seu_obj)
  
}

dmso_s <- add_genome_to_peaks(dmso_s)
enz48_s <- add_genome_to_peaks(enz48_s)
resa_s <- add_genome_to_peaks(resa_s)
resb_s <- add_genome_to_peaks(resb_s)


# compute nucleosome signal score per cell
dmso_s <- NucleosomeSignal(object = dmso_s, verbose = TRUE)
enz48_s <- NucleosomeSignal(object = enz48_s, verbose = TRUE)
resa_s <- NucleosomeSignal(object = resa_s, verbose = TRUE)
resb_s <- NucleosomeSignal(object = resb_s, verbose = TRUE)

# compute TSS enrichment score per cell
dmso_s <- TSSEnrichment(object = dmso_s, verbose = TRUE)
enz48_s <- TSSEnrichment(object = enz48_s, verbose = TRUE)
resa_s <- TSSEnrichment(object = resa_s, verbose = TRUE)
resb_s <- TSSEnrichment(object = resb_s, verbose = TRUE)

## total number of peaks in encode blacklisted regions
library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")

blacklist_regions <- ah[['AH107305']]

dmso_s$blacklist_ratio <- FractionCountsInRegion(
  object = dmso_s, 
  assay = 'peaks',
  regions = blacklist_regions
)

enz48_s$blacklist_ratio <- FractionCountsInRegion(
  object = enz48_s, 
  assay = 'peaks',
  regions = blacklist_regions
)

resa_s$blacklist_ratio <- FractionCountsInRegion(
  object = resa_s, 
  assay = 'peaks',
  regions = blacklist_regions
)

resb_s$blacklist_ratio <- FractionCountsInRegion(
  object = resb_s, 
  assay = 'peaks',
  regions = blacklist_regions
)

hist(resb_s$nCount_peaks)
hist(resb_s$nucleosome_signal)
hist(resb_s$TSS.enrichment)
hist(resb_s$blacklist_ratio)


hist(resb_s_f$nCount_peaks)
hist(resb_s_f$nucleosome_signal)
hist(resb_s_f$TSS.enrichment)
hist(resb_s_f$blacklist_ratio)


dmso_s_f <- subset(
  dmso_s,
  subset =
    nCount_peaks > 2000 & nCount_peaks < 25000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    blacklist_ratio < 0.01
)


enz48_s_f <- subset(
  enz48_s,
  subset =
    nCount_peaks > 1000 & nCount_peaks < 20000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    blacklist_ratio < 0.01
)


resa_s_f <- subset(
  resa_s,
  subset =
    nCount_peaks > 1000 & nCount_peaks < 20000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    blacklist_ratio < 0.01
)

resb_s_f <- subset(
  resb_s,
  subset =
    nCount_peaks > 1000 & nCount_peaks < 25000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    blacklist_ratio < 0.01
)

#enz48_s_f$Sample_Id <- "ENZ48"
#resa_s_f$Sample_Id <- "RESA"
#resb_s_f$Sample_Id <- "RESB"


combined.atac <- merge(
  dmso_s_f,
  y = list(enz48_s_f, resa_s_f, resb_s_f),
  add.cell.ids = c("DMSO", "ENZ48", "RESA", "RESB")
)


DefaultAssay(combined.atac) <- "peaks"

combined.atac <- RunTFIDF(combined.atac, verbose = TRUE)

combined.atac <- FindTopFeatures(combined.atac,min.cutoff = "q50", verbose = TRUE)

combined.atac <- RunSVD(
  combined.atac,
  features = VariableFeatures(combined.atac),
  verbose = TRUE
)

DepthCor(combined.atac)

combined.atac <- FindNeighbors(
  combined.atac,
  reduction = "lsi",
  dims = 2:30,
  verbose = TRUE
)

combined.atac <- FindClusters(
  combined.atac,
  graph.name = "peaks_snn",
  resolution = c(0.2, 0.4, 0.6, 0.8),
  algorithm = 4,
  verbose = TRUE,
  method = "igraph"
)

combined.atac <- RunUMAP(
  combined.atac,
  reduction = "lsi",
  dims = 2:30,
  verbose = TRUE
)

DimPlot(object = combined.atac, label = TRUE, group.by = "Sample_Id") +
  NoLegend() +
  ggtitle('UMAP')

run_lsi_each_sample <- function(seu_obj){
  
  DefaultAssay(seu_obj) <- "peaks"
  
  seu_obj <- RunTFIDF(seu_obj)
  seu_obj <- FindTopFeatures(seu_obj, min.cutoff = "q50")
  seu_obj <- RunSVD(seu_obj)
  return(seu_obj)
}

dmso_s_f <- run_lsi_each_sample(dmso_s_f)
enz48_s_f <- run_lsi_each_sample(enz48_s_f)
resa_s_f <- run_lsi_each_sample(resa_s_f)
resb_s_f <- run_lsi_each_sample(resb_s_f)

dmso_s_f <- RenameCells(dmso_s_f, add.cell.id = "DMSO")
enz48_s_f <- RenameCells(enz48_s_f, add.cell.id = "ENZ48")
resa_s_f <- RenameCells(resa_s_f, add.cell.id = "RESA")
resb_s_f <- RenameCells(resb_s_f, add.cell.id = "RESB")

atac.obj.list <- SplitObject(combined.atac, split.by = "Sample_Id")


## Hrmony integration does not work. Using Seurat integration method
integration.anchors <- FindIntegrationAnchors(
  object.list = atac.obj.list,
  anchor.features = rownames(combined.atac),
  reduction = "rlsi",
  dims = 2:30,
  verbose = TRUE
)

rm(dmso_s, enz48_s, resa_s, resb_s)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined.atac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 2:30,
  k.weight = 5
)

Reductions(integrated)

dim(integrated)

dim(Embeddings(integrated, "integrated_lsi")) # 12450    29


head(integrated[,1:5])

integrated[["integrated_lsi"]]

integrated@reductions$integrated_lsi <- CreateDimReducObject(
  embeddings = integrated,
  key = "integratedLSI_",
  assay = DefaultAssay(integrated)
)


integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:29, k.param=13, verbose = TRUE)


integrated <- FindClusters(object = integrated, verbose = TRUE, algorithm = 4, resolution = c(0.2,0.4,0.6,0.8), leiden_method = "igraph")

integrated <- RunUMAP(
  integrated,
  reduction = "integrated_lsi",
  dims = 2:29,
  verbose = TRUE
)


DimPlot(object = integrated, label = TRUE, group.by = "Sample_Id") +
  NoLegend() +
  ggtitle('UMAP') ## overcorrection, whole thing is a blob

DimPlot(object = integrated, label = TRUE, group.by = "ident") +
  NoLegend() +
  ggtitle('UMAP') ## no clear cluster separation



## Use Harmony as mentioned in paper
library(harmony)

combined.atac <- RunHarmony(
  object = combined.atac,
  group.by.vars = "Sample_Id",
  reduction.use = "lsi",
  dims.use = 2:30,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)

combined.atac <- RunUMAP(
  combined.atac,
  reduction = "harmony",
  dims = 1:29,
  verbose = TRUE
)

DimPlot(object = combined.atac, label = FALSE, group.by = "Sample_Id") +
  NoLegend() +
  ggtitle('UMAP') ## overcorrection, whole thing is a blob


combined.atac <- FindNeighbors(
  combined.atac,
  reduction = "harmony",
  dims = 1:29,
  verbose = TRUE
)

combined.atac <- FindClusters(
  combined.atac,
  algorithm = 4,
  resolution = c(0.2, 0.4, 0.6,0.8),
  method = "igraph"
)

Idents(combined.atac) <- combined.atac@meta.data$peaks_snn_res.0.6

## 0.8 overclusters it, 0.6 maybe okay
DimPlot(object = combined.atac, label = TRUE) +
  NoLegend() +
  ggtitle('UMAP')

pdf("scATAC_UMAP_after_clustering_as_per_paper.pdf", width = 9)
DimPlot(object = combined.atac, label = TRUE, split.by = "Sample_Id") +
  NoLegend() +
  ggtitle('UMAP')
dev.off()

saveRDS(
  combined.atac,
  file = "combined_atac_harmony_integrated.rds"
)


## Overlap the peaks
peaks.dmso <- dmso_s_f@assays$peaks@ranges
peaks.enz48 <- enz48_s_f@assays$peaks@ranges
peaks.resa <- resa_s_f@assays$peaks@ranges
peaks.resb <- resb_s_f@assays$peaks@ranges


ol <- findOverlapsOfPeaks(peaks.dmso,
                          peaks.enz48,
                          peaks.resa,
                          peaks.resb,
                          ignore.strand = FALSE,
                          connectedPeaks="keepAll")

## venn diagram to show the overlaps
pdf("Overlap_scATAC_peaks_4_treatments_venn_diagram.pdf", width = 10)
makeVennDiagram(ol, connectedPeaks = "keepAll", by=c("region"),NameOfPeaks = c("DMSO", "ENZ48", "RESA","RESB"), cex = 0.8)
dev.off()

overlap.venn$vennCounts


overlap.venn.matrix <- as.matrix(overlap.venn$vennCounts)

overlap.venn.matrix <- as.data.frame(overlap.venn.matrix)

#upset_df <- upset_df[
#  rowSums(upset_df[, c(
#    "peaks.dmso",
#    "peaks.enz48",
#    "peaks.resa",
#    "peaks.resb"
#  )]) > 0,
#]

#upset_expanded <- upset_df[
#  rep(seq_len(nrow(upset_df)), upset_df$Counts),
#  c("peaks.dmso", "peaks.enz48", "peaks.resa", "peaks.resb")
#]

#upset_expanded <- as.data.frame(upset_expanded)

#upset(
#  upset_expanded,
#  sets = c(
#    "peaks.dmso",
#    "peaks.enz48",
#    "peaks.resa",
#    "peaks.resb"
#  ),
#  order.by = "freq",
#  decreasing = TRUE,
#  keep.order = TRUE
#)



# Build peak-by-sample overlap matrix
#ov1 <- countOverlaps(consensus.peaks, peaks.dmso) > 0
#ov2 <- countOverlaps(consensus.peaks, peaks.enz48) > 0
#ov3 <- countOverlaps(consensus.peaks, peaks.resa) > 0
#ov4 <- countOverlaps(consensus.peaks, peaks.resb) > 0

#peak_df <- data.frame(
#  DMSO = ov1,
#  ENZ48 = ov2,
#  RESA = ov3,
#  RESB = ov4
#)

#library(UpSetR)

#peak_df_copy <- peak_df

#peak_df_copy$DMSO <- ifelse(peak_df$DMSO == "FALSE", 0,1)
#peak_df_copy$ENZ48 <- ifelse(peak_df$ENZ48 == "FALSE", 0,1)
#peak_df_copy$RESA <- ifelse(peak_df$RESA == "FALSE", 0,1)
#peak_df_copy$RESB <- ifelse(peak_df$RESB == "FALSE", 0,1)


#pdf("Overlap_scATAC_peaks_4_treatments.pdf", width = 13)

#UpSetR::upset(
#  peak_df_copy,
#  sets = c("DMSO","ENZ48","RESA","RESB"),
#  order.by = "freq",
#  decreasing = TRUE,
#  keep.order = TRUE
#)

#dev.off()

## plot fxn of cells by treatment per cluster
combined.meta <- combined.atac@meta.data

combined.meta$peaks_snn_res.0.6 <- as.integer(as.character(combined.meta$peaks_snn_res.0.6))

combined.meta$seurat_clusters <- combined.meta$peaks_snn_res.0.6

combined.meta$Col_1 <- 1

## Using 0.6 resolution for determining cluster fxn
cluster.by.treatment <- as.data.frame(
  table(
    combined.meta$seurat_clusters,
    combined.meta$Sample_Id
  )
)

colnames(cluster.by.treatment) <- c(
  "seurat_clusters",
  "Treatment",
  "n_cells"
)

cell.by.treatment <- as.data.frame(
  table(
    combined.meta$Sample_Id
  )
)

cluster.by.treatment <- left_join(cluster.by.treatment, cell.by.treatment, by = c("Treatment" = "Var1"))

cluster.by.treatment$Fxn_cells <- cluster.by.treatment$n_cells / cluster.by.treatment$Freq



dmso <- cluster.by.treatment[cluster.by.treatment$Treatment == "DMSO,"]
  
  
 
library(ggplot2)

pdf("Fraction_of_cells_per_cluster_by_treatment_for_scATAC_seq.pdf", width = 5, height = 5)

ggplot(
  cluster.by.treatment,
  aes(x = Treatment, y = Fxn_cells, fill = Treatment)
) +
  geom_col() +
  facet_wrap(~ seurat_clusters, nrow = 1) +
  labs(
    y = "Fraction of cells",
    fill = "Treatment"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  scale_fill_manual(
    values = c(
      DMSO = "grey50",
      ENZ48 = "steelblue",
      RESA = "orange",
      RESB = "firebrick"
    ),
    breaks = c("DMSO", "ENZ48", "RESA", "RESB")
  )

dev.off()


## plot enrichment of 4 treatments at various gene groups
house_keeping_genes <- read.csv(file="Housekeeping_Genes_list.csv",header = T,stringsAsFactors = F)
androgen_genes <- read.table(file="HALLMARK_ANDROGEN_RESPONSE.v2026.1.Hs.gmt",header = T,stringsAsFactors = F,sep="\t")
myc_genes <- read.csv(file="MYC_target_genes.csv",header = T,stringsAsFactors = F)

androgen_genes <- colnames(androgen_genes)
androgen_genes <- androgen_genes[c(3:length(androgen_genes))]


plot_TSS_enrcihment_roi <- function(atac_seurat, roi_gene_ranges,roi_name){
  
  atac_seurat <- TSSEnrichment(
    atac_seurat, 
    fast = FALSE, 
    tss.positions = roi_gene_ranges
  )
  
  tss <- atac_seurat[["peaks"]]@positionEnrichment$TSS
  
  # collapse into TSS-centered profile
  sample <- atac_seurat$Sample_Id
  
  group_means <- t(
    sapply(unique(sample), function(s) {
      Matrix::colMeans(tss[sample == s, , drop = FALSE])
    })
  )
  
  ## get proper TSS scaling
  position <- as.numeric(colnames(tss))
  
  df <- data.frame(
    position = position,
    t(group_means)
  )
  
  library(reshape2)
  df_long <- melt(df, id.vars = "position")
  
  
  
  ## plot
  library(ggplot2)
  
  pdf(paste0("ATAC_TSS_enrichment_all_regions_",roi_name,".pdf"), width = 13)
  p <- ggplot(df_long, aes(position, value, color = variable)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.2, linewidth = 1.2) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_classic(base_size = 14) +
    xlab("Distance from TSS (bp)") +
    ylab("Accessibility") +
    scale_color_manual(values = c(
      DMSO = "#2c7bb6",
      ENZ48 = "#fdae61",
      RESA  = "#abd9e9",
      RESB  = "#d7191c"
    )) + ylim(0,5)
  
  print(p)
  
  dev.off()
  
}

length(intersect(house_keeping_genes$Housekeeping_genes,androgen_genes))

## enrichment at 5000 random gene promoters
all_gene_annot <- Annotation(combined.atac)


hkg_gene_annot <- all_gene_annot[all_gene_annot$gene_name %in% house_keeping_genes$Housekeeping_genes]
androgen_gene_annot <- all_gene_annot[all_gene_annot$gene_name %in% androgen_genes]
myc_gene_annot <- all_gene_annot[all_gene_annot$gene_name %in% myc_genes$Gene]



promoters_hkg <- promoters(
  hkg_gene_annot,
  upstream = 1000,
  downstream = 1000
)

promoters_androgen <- promoters(
  androgen_gene_annot,
  upstream = 1000,
  downstream = 1000
)

promoters_myc <- promoters(
  myc_gene_annot,
  upstream = 1000,
  downstream = 1000
)

plot_TSS_enrcihment_roi(combined.atac, promoters_hkg, "housekeeping_genes_test")
plot_TSS_enrcihment_roi(combined.atac, promoters_androgen, "androgen_genes")
plot_TSS_enrcihment_roi(combined.atac, promoters_myc, "myc_genes")



random.genes.5000 <- sample(unique(all_gene_annot$gene_id),size=5000)


random_gene_annot <- all_gene_annot[all_gene_annot$gene_id %in% random.genes.5000]


promoters_random_genes <- GenomicRanges::promoters(
  random_gene_annot,
  upstream = 1000,
  downstream = 1000
)

plot_TSS_enrcihment_roi(combined_ATAC, promoters_random_genes, "random_genes")




'''
Set of differentially accessible chromatin regions (DARs) in each scATAC-seq cluster (compared to all other
clusters), referred to as marker differentially accessible regions.
'''

table(Idents(combined.atac))


DefaultAssay(combined.atac) <- "peaks"

## Set of DARs in each scATAC-seq cluster in each sample, compared to all other clusters in the sample.
marker_DARs <- FindAllMarkers(
  object = combined.atac,
  assay = "peaks",
  test.use = "LR",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  only.pos = TRUE,
  verbose = TRUE
)

dmso.from.combined <- subset(combined.atac, subset = Sample_Id == "DMSO")
resa.from.combined <- subset(combined.atac, subset = Sample_Id == "RESA")
resb.from.combined <- subset(combined.atac, subset = Sample_Id == "RESB")
enz.from.combined <- subset(combined.atac, subset = Sample_Id == "ENZ48")


find_DARs_and_annotate <- function(
    seurat_obj,
    test.use = "LR",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    only.pos = TRUE,
    verbose = TRUE
) {
  
  # Make sure peak assay is used
  DefaultAssay(seurat_obj) <- "peaks"
  
  # Find marker DARs: each cluster vs all other clusters
  DARs <- FindAllMarkers(
    object = seurat_obj,
    assay = "peaks",
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    only.pos = only.pos,
    verbose = verbose
  )
  
  # Bonferroni correction
  DARs$p_val_BF <- p.adjust(
    DARs$p_val,
    method = "bonferroni"
  )
  
  # Peak coordinates
  DARs$region <- DARs$gene
  
  # Convert peak strings to GRanges
  regions_DARs_gr <- Signac::StringToGRanges(
    DARs$gene,
    sep = c("-", "-")
  )
  
  # Find closest genes
  closest_genes_DARs <- Signac::ClosestFeature(
    object = seurat_obj,
    regions = regions_DARs_gr
  )
  
  # Add gene annotation
  DARs_anno <- dplyr::left_join(
    DARs,
    closest_genes_DARs,
    by = c("region" = "query_region")
  )
  
  return(DARs_anno)
}

dmso.dar.anno <- find_DARs_and_annotate(dmso.from.combined)
resa.dar.anno <- find_DARs_and_annotate(resa.from.combined)
resb.dar.anno <- find_DARs_and_annotate(resb.from.combined)
enz.dar.anno <- find_DARs_and_annotate(enz.from.combined)

#write.table(dmso.dar.anno,file="Differentially_accessible_regions_all_clusters_DMSO.txt",col.names = T,quote = F, sep="\t", row.names = F)
#write.table(resa.dar.anno,file="Differentially_accessible_regions_all_clusters_RESA.txt",col.names = T,quote = F, sep="\t", row.names = F)
#write.table(resb.dar.anno,file="Differentially_accessible_regions_all_clusters_RESB.txt",col.names = T,quote = F, sep="\t", row.names = F)
#write.table(enz.dar.anno,file="Differentially_accessible_regions_all_clusters_ENZ48.txt",col.names = T,quote = F, sep="\t", row.names = F)

run_dmso_dar_motif_analysis <- function(
    treatment.dar.anno,
    treatment.from.combined,
    output_prefix = "DMSO_motifs"
) {

  library(Seurat)
  library(Signac)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2018)
  
  # Convert region into chr, start, end
  
  treatment.dar.anno[, c("chr", "start", "end")] <- do.call(
    rbind,
    strsplit(treatment.dar.anno$region, "-", fixed = TRUE)
  )
  
  treatment.dar.anno$start <- as.integer(treatment.dar.anno$start)
  treatment.dar.anno$end   <- as.integer(treatment.dar.anno$end)
  
  treatment_dar_sig <- treatment.dar.anno[
    treatment.dar.anno$p_val_BF < 0.05,
  ]

  pfm <- getMatrixSet(
    x = JASPAR2018,
    opts = list(
      collection = "CORE",
      tax_group = "vertebrates"
    )
  )

  tf_keep <- c(
    "Ar",
    "CREB1",
    "CTCF",
    "E2F1",
    "ELF1",
    "ELK4",
    "ETV1",
    "FOXA1",
    "FOXP1",
    "GATA2",
    "GRHL2",
    "HOXB13",
    "JUND",
    "MYC",
    "NKX3-1",
    "POU2F1",
    "RELA",
    "SOX9",
    "TFAP4",
    "ZFX"
  )
  
  
  # Get names of motifs
  pfm_names <- sapply(
    pfm,
    function(x) x@name
  )
  
  
  # Keep only requested TFs
  pfm.use <- pfm[
    pfm_names %in% tf_keep
  ]
  
  
  message(
    "Number of TF motifs selected: ",
    length(pfm.use)
  )
  
  DefaultAssay(treatment.from.combined) <- "peaks"
  
  
  # Calculate GC content for all peaks
  treatment.from.combined <- RegionStats(
    object = treatment.from.combined,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  
  # Add motif information
  treatment.from.combined <- AddMotifs(
    object = treatment.from.combined,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm.use,
    verbose = TRUE
  )

  clusters <- sort(
    unique(treatment_dar_sig$cluster)
  )
  
  message(
    "Clusters to analyse: ",
    paste(clusters, collapse = ", ")
  )
  
  motif_results <- list()

  for (cl in clusters) {
    
    message(
      "Processing cluster ", cl
    )
    
    # Get DARs belonging to this cluster
    peaks_cluster <- treatment_dar_sig$gene[
      treatment_dar_sig$cluster == cl
    ]
    
    # Remove duplicated peaks
    peaks_cluster <- unique(peaks_cluster)
    
    # Skip clusters with no DARs
    if (length(peaks_cluster) == 0) {
      
      message(
        "No DARs found for cluster ",
        cl,
        ". Skipping."
      )
      
      next
    }
    
    
    message(
      "Number of DARs: ",
      length(peaks_cluster)
    )
  
    # Find motif enrichment
    
    motif_result <- FindMotifs(
      object = treatment.from.combined,
      features = peaks_cluster,
      pfm = pfm.use,
      verbose = TRUE
    )
    
    motif_result$p.adjust.FDR <- p.adjust(
      motif_result$pvalue,
      method = "BH"
    )
    
    motif_results[[paste0("cluster_", cl)]] <- motif_result
    
    output_file <- paste0(
      output_prefix,
      "_cluster_",
      cl,
      ".csv"
    )
    
    write.csv(
      motif_result,
      file = output_file,
      row.names = FALSE,
      quote = FALSE
    )
    
  }
  
  return(
    list(
      significant_DARs = treatment_dar_sig,
      motif_results = motif_results,
      seurat_object = treatment.from.combined
    )
  )
}

run_dmso_dar_motif_analysis(dmso.dar.anno, dmso.from.combined, "DMSO_motifs") ## ran okay
run_dmso_dar_motif_analysis(enz.dar.anno, enz.from.combined, "ENZ_motifs") ## ran okay
run_dmso_dar_motif_analysis(resa.dar.anno, resa.from.combined, "RESA_motifs") ## ran okay
run_dmso_dar_motif_analysis(resb.dar.anno, resb.from.combined, "RESB_motifs")

## Compare old and new overlapping regions
lncap.dar.anno <- read.table(file="Differentially_accessible_peaks_all_clusters_DMSO.txt",header = T,sep="\t", stringsAsFactors = F)
resa.dar.anno <- read.table(file="Differentially_accessible_regions_all_clusters_RESA.txt",header = T,sep="\t", stringsAsFactors = F)
resb.dar.anno <- read.table(file="Differentially_accessible_regions_all_clusters_RESB.txt",header = T,sep="\t", stringsAsFactors = F)
enz.dar.anno <- read.table(file="Differentially_accessible_regions_all_clusters_ENZ48.txt",header = T,sep="\t", stringsAsFactors = F)

## overlap paper DAR's and currently discovered DAR's
dar_lncap_paper <- read.csv(file="DAR_for_individual_clusters_from_paper_supp_for_LNCap.csv",header = T,stringsAsFactors = F)
dar_resa_paper <- read.csv(file="DAR_for_individual_clusters_from_paper_supp_for_RESA.csv",header = T,stringsAsFactors = F)
dar_resb_paper <- read.csv(file="DAR_for_individual_clusters_from_paper_supp_for_RESB.csv",header = T,stringsAsFactors = F) ## okay
dar_enz48_paper <- read.csv(file="DAR_for_individual_clusters_from_paper_supp_for_LNCap_ENZ48.csv",header = T,stringsAsFactors = F) ## okay 


dmso.dar.sig <- dmso.dar.anno[dmso.dar.anno$p_val_BF < 0.05, ]
resa.dar.sig <- resa.dar.anno[resa.dar.anno$p_val_BF < 0.05, ]
resb.dar.sig <- resb.dar.anno[resb.dar.anno$p_val_BF < 0.05, ]
enz.dar.sig <- enz.dar.anno[enz.dar.anno$p_val_BF < 0.05, ]

dmso_dar_sig_unique <- dmso.dar.sig[ !duplicated(dmso.dar.sig$region), ]
resa_dar_sig_unique <- resa.dar.sig[ !duplicated(resa.dar.sig$region), ]
resb_dar_sig_unique <- resb.dar.sig[ !duplicated(resb.dar.sig$region), ]
enz_dar_sig_unique <- enz.dar.sig[ !duplicated(enz.dar.sig$region), ]

convert_loci_to_granges <- function(treatment_dar_sig_unique){
  
  coords <- do.call(
    rbind,
    strsplit(unique(treatment_dar_sig_unique$region), "-", fixed = TRUE)
  )
  
  treatment_dar_loci <- data.frame(
    chr = coords[, 1],
    start = as.integer(coords[, 2]),
    end = as.integer(coords[, 3]),
    stringsAsFactors = FALSE
  )
  
  treatment_dar_gr <- GenomicRanges::GRanges(
    seqnames = treatment_dar_loci$chr,
    ranges = IRanges::IRanges(
      start = treatment_dar_loci$start,
      end = treatment_dar_loci$end
    )
  )
  
  return(treatment_dar_gr)
}

dmso_dar_gr <- convert_loci_to_granges(dmso_dar_sig_unique) ## conversion okay
resa_dar_gr <- convert_loci_to_granges(resa_dar_sig_unique) ## conversion okay
resb_dar_gr <- convert_loci_to_granges(resb_dar_sig_unique) ## conversion okay
enz48_dar_gr <- convert_loci_to_granges(enz_dar_sig_unique) ## conversion okay

convert_loci_to_granges_paper_table <- function(dar_treatment_paper){
  
  coords <- do.call(
    rbind,
    strsplit(dar_treatment_paper$Chromatin.region, "[:-]")
  )
  
  dar_treatment_paper$chr <- coords[, 1]
  dar_treatment_paper$start <- as.integer(coords[, 2])
  dar_treatment_paper$end <- as.integer(coords[, 3])
  
  
  dar_treatment_paper_gr <- GenomicRanges::GRanges(
    seqnames = dar_treatment_paper$chr,
    ranges = IRanges::IRanges(
      start = dar_treatment_paper$start,
      end = dar_treatment_paper$end
    )
  )
  
  return(dar_treatment_paper_gr)
}

dar_lncap_paper <- dar_lncap_paper[complete.cases(dar_lncap_paper),]
dar_enz48_paper <- dar_enz48_paper[complete.cases(dar_enz48_paper),]
dar_resa_paper <- dar_resa_paper[complete.cases(dar_resa_paper),]
dar_resb_paper <- dar_resb_paper[complete.cases(dar_resb_paper),]


dar_lncap_paper_gr <- convert_loci_to_granges_paper_table(dar_lncap_paper)
dar_enz48_paper_gr <- convert_loci_to_granges_paper_table(dar_enz48_paper)
dar_resa_paper_gr <- convert_loci_to_granges_paper_table(dar_resa_paper)
dar_resb_paper_gr <- convert_loci_to_granges_paper_table(dar_resb_paper)


pdf("Overlap_of_DARs_paper_reanalysis_LNCAP.pdf")
ol <- findOverlapsOfPeaks(dmso_dar_gr,
                          dar_lncap_paper_gr,
                          ignore.strand = FALSE,
                          connectedPeaks="keepAll")
## venn diagram to show the overlaps
makeVennDiagram(ol, connectedPeaks = "keepAll")
dev.off()


pdf("Overlap_of_DARs_paper_reanalysis_ENZ48.pdf")
ol <- findOverlapsOfPeaks(enz48_dar_gr,
                          dar_enz48_paper_gr,
                          ignore.strand = FALSE,
                          connectedPeaks="keepAll")
## venn diagram to show the overlaps
makeVennDiagram(ol, connectedPeaks = "keepAll")
dev.off()

pdf("Overlap_of_DARs_paper_reanalysis_RESA.pdf")
ol <- findOverlapsOfPeaks(resa_dar_gr,
                          dar_resa_paper_gr,
                          ignore.strand = FALSE,
                          connectedPeaks="keepAll")
## venn diagram to show the overlaps
makeVennDiagram(ol, connectedPeaks = "keepAll")
dev.off()

pdf("Overlap_of_DARs_paper_reanalysis_RESB.pdf")
ol <- findOverlapsOfPeaks(resb_dar_gr,
                          dar_resb_paper_gr,
                          ignore.strand = FALSE,
                          connectedPeaks="keepAll")

## venn diagram to show the overlaps
makeVennDiagram(ol, connectedPeaks = "keepAll")
dev.off()


## get scATAC open space for LNCAP ENZ48 vs LNCAP
head(integrated@meta.data)
table(integrated@meta.data$Sample_Id)

integrated@meta.data$seurat_clusters <- integrated@meta.data$peaks_snn_res.0.6


DefaultAssay(integrated) <- "peaks"

enz48_vs_dmso <- FindMarkers(
  object = integrated,
  ident.1 = "ENZ48",
  ident.2 = "DMSO",
  group.by = "Sample_Id",
  assay = "peaks",
  test.use = "LR",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  only.pos = FALSE
)


DefaultAssay(integrated) <- "peaks"

# Keep only the two conditions
obj <- subset(
  integrated,
  subset = Sample_Id %in% c("DMSO", "ENZ48")
)

clusters <- sort(unique(obj$seurat_clusters))

DE_results <- list()

for (cl in clusters) {
  
  message("Processing cluster: ", cl)
  
  # Cells from this cluster
  cells_cl <- rownames(obj@meta.data)[
    obj@meta.data$seurat_clusters == cl
  ]
  
  obj_cl <- subset(
    obj,
    cells = cells_cl
  )
  
  # DMSO vs ENZ48
  de <- FindMarkers(
    object = obj_cl,
    ident.1 = "ENZ48",
    ident.2 = "DMSO",
    group.by = "Sample_Id",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    only.pos = TRUE
  )
  
  de$gene <- rownames(de)
  de$cluster <- cl
  
  DE_results[[paste0("cluster_", cl)]] <- de
}

write.table(DE_results$cluster_1,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_1.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_2,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_2.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_3,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_3.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_4,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_4.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_5,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_5.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_6,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_6.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_7,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_7.txt",col.names = T,quote = F, sep="\t", row.names = F)
write.table(DE_results$cluster_8,file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_8.txt",col.names = T,quote = F, sep="\t", row.names = F)



run_dar_by_cluster <- function(
    integrated,
    condition_1,
    condition_2,
    output_dir = ".",
    min.pct = 0.15,
    logfc.threshold = 0.25
) {
  
  library(Seurat)
  
  DefaultAssay(integrated) <- "peaks"
  
  obj <- subset(
    integrated,
    subset = Sample_Id %in% c(condition_1, condition_2)
  )
  
  clusters <- sort(
    unique(obj$seurat_clusters)
  )
  
  message(
    "Conditions: ",
    condition_1,
    " vs ",
    condition_2
  )
  
  message(
    "Clusters: ",
    paste(clusters, collapse = ", ")
  )
  
  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  DE_results <- list()
  
  for (cl in clusters) {
    
    message(
      "Processing cluster: ",
      cl
    )
  
    cells_cl <- rownames(integrated@meta.data)[
      integrated@meta.data$seurat_clusters == cl &
        integrated@meta.data$Sample_Id %in% c(condition_1, condition_2)
    ]
    
    
    obj_cl <- subset(
      obj,
      cells = cells_cl
    )
  
    condition_counts <- table(
      obj_cl$Sample_Id
    )
    
    print(condition_counts)
    
    if (
      !condition_1 %in% names(condition_counts) ||
      !condition_2 %in% names(condition_counts)
    ) {
      
      message(
        "Skipping cluster ",
        cl,
        ": one of the conditions is absent."
      )
      
      next
    }
    
    # Positive avg_log2FC = more accessible in condition_1
    
    de <- FindMarkers(
      object = obj_cl,
      ident.1 = condition_1,
      ident.2 = condition_2,
      group.by = "Sample_Id",
      test.use = "wilcox",
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos = TRUE
    )
    
    de$gene <- rownames(de)
    
    de$cluster <- cl
    
    DE_results[[paste0("cluster_", cl)]] <- de
    
    output_file <- file.path(
      output_dir,
      paste0(
        "Differentially_accessible_regions_",
        condition_1,
        "_vs_",
        condition_2,
        "_cluster_",
        cl,
        ".txt"
      )
    )
    
    
    write.table(
      de,
      file = output_file,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
    
    
    message(
      "Written: ",
      output_file
    )
  }

  return(DE_results)
  
}


resa_vs_dmso_dar <- run_dar_by_cluster(integrated, "RESA", "DMSO")
resb_vs_dmso_dar <- run_dar_by_cluster(integrated, "RESB", "DMSO")
resa_vs_enz48_dar <- run_dar_by_cluster(integrated, "RESA", "ENZ48")
resb_vs_enz48_dar <- run_dar_by_cluster(integrated, "RESB", "ENZ48")




