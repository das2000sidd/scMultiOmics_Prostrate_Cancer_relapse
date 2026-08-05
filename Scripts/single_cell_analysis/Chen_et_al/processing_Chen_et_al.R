setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")

library(data.table)
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(ggplot2)

file <- "GSM4203181_data.raw.matrix.txt"
con <- file(file, open = "r")

header <- readLines(con, n = 1)

i_all <- integer()
j_all <- integer()
x_all <- numeric()

gene_names <- character()
row_index <- 0

repeat {
  
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  
  row <- fread(
    text = line,
    header = FALSE,
    data.table = FALSE
  )
  
  gene <- row[[1]]
  vals <- as.numeric(row[-1])
  
  nz <- which(vals != 0)
  
  if (length(nz) > 0) {
    i_all <- c(i_all, rep(row_index + 1, length(nz)))
    j_all <- c(j_all, nz)
    x_all <- c(x_all, vals[nz])
  }
  
  gene_names <- c(gene_names, gene)
  row_index <- row_index + 1
}

close(con)

mat_sparse <- sparseMatrix(
  i = i_all,
  j = j_all,
  x = x_all
)

cell_ids <- read.table(file="columns_ids_of_matrix.txt",header = F,sep="\t", stringsAsFactors = F)
length(unique(cell_ids$V1))

colnames(cell_ids)[1] <- "barcode"


cell_ids$new_column <- sub("-.*", "", cell_ids$barcode)
cell_ids$num <- sub(".*-", "", cell_ids$barcode)
cell_ids$Patient_ID <- paste("Patient",cell_ids$num,sep="_")

cell_ids <- cell_ids[,c(1,4)]


rownames(mat_sparse) <- gene_names
colnames(mat_sparse) <- cell_ids$barcode

seu.obj <- CreateSeuratObject(counts = mat_sparse, min.cells = 3, min.features = 200)
seu.obj


View(seu.obj@meta.data)

length(intersect(cell_ids$barcode, rownames(seu.obj@meta.data)))

seu.obj@meta.data$Cell_id <- rownames(seu.obj@meta.data)

meta.data <- seu.obj@meta.data

meta.data <- left_join(meta.data,cell_ids, by = c("Cell_id" = "barcode"))

rownames(meta.data) <- meta.data$Cell_id

seu.obj@meta.data <- meta.data


housekeeping.genes = read.csv(file="Housekeeping_genes_Chen_at_al.csv",header = T,stringsAsFactors = F)



# Calculate mitochondrial percentage
seu.obj[["percent.mt"]] <- PercentageFeatureSet(
  seu.obj,
  pattern = "^MT-"
)

View(seu.obj@meta.data)

# Visualize QC metrics
VlnPlot(
  seu.obj,
  features = c("nFeature_RNA",
               "nCount_RNA",
               "percent.mt"),
  ncol = 3,
  pt.size = 0
)


p1 <- FeatureScatter(
  seu.obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

p2 <- FeatureScatter(
  seu.obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

p1 + p2

seu.obj <- subset(
  seu.obj,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 5200
)


housekeeping <- housekeeping.genes$Housekeeping_genes

sum(housekeeping %in% rownames(seu.obj))

hk.present <- housekeeping[housekeeping %in% rownames(seu.obj)]

hk.genes <- Matrix::colSums(
  GetAssayData(seu.obj, layer="counts")[hk.present, ] > 0
)

seu.obj$HK <- hk.genes

#seu.obj <- subset(seu.obj, subset = HK >= 56)


library(SingleCellExperiment)
library(scran)
library(batchelor)

samples.use <- lapply(samples.use, NormalizeData)

samples.use <- lapply(samples.use, ScaleData)

sce.list <- lapply(
  samples.use,
  as.SingleCellExperiment
)

sce.mnn <- do.call(
  fastMNN,
  c(
    sce.list,
    list(
      k = 5,
      d = 50,
      approximate = TRUE,
      auto.order = TRUE
    )
  )
)

library(BiocNeighbors)

sce.mnn <- fastMNN(
  sce.list[[1]],
  sce.list[[2]],
  sce.list[[3]],
  sce.list[[4]],
  sce.list[[5]],
  sce.list[[6]],
  k = 5,
  d = 50,
  auto.merge = TRUE
)

reducedDimNames(sce.mnn)

dim(reducedDim(sce.mnn, "corrected"))

# Add the MNN embedding back to Seurat
library(Seurat)

seu.obj.subset[["mnn"]] <- CreateDimReducObject(
  embeddings = reducedDim(sce.mnn, "corrected"),
  key = "MNN_",
  assay = "RNA"
)

Reductions(seu.obj.subset)


# Run UMAP on the corrected space
seu.obj.subset <- RunUMAP(
  seu.obj.subset,
  reduction = "mnn",
  dims = 1:50,
  verbose = TRUE
)

# Before MNN
DimPlot(
  seu.obj.subset,
  reduction = "pca",
  group.by = "Patient_ID"
)

# After MNN
DimPlot(
  seu.obj.subset,
  reduction = "umap",
  group.by = "Patient_ID"
)

# Cluster cells using MNN dimensions
seu.obj.subset <- FindNeighbors(
  seu.obj.subset,
  reduction = "RNA_snn",
  dims = 1:50,
  verbose = TRUE
)

## throwing out of memory
seu.obj.subset <- FindClusters(
  seu.obj.subset,
  graph.name = "RNA_snn",
  resolution = 0.2,
  algorithm = 4,
  random.seed = 123
)

seu.obj.subset <- FindClusters(
  seu.obj.subset,
  graph.name = "RNA_snn",
  resolution = 0.4,
  algorithm = 4,
  random.seed = 123
)

seu.obj.subset <- FindClusters(
  seu.obj.subset,
  graph.name = "RNA_snn",
  resolution = 0.6,
  algorithm = 4,
  random.seed = 123
)

seu.obj.subset <- FindClusters(
  seu.obj.subset,
  graph.name = "RNA_snn",
  resolution = 0.8,
  algorithm = 4,
  random.seed = 123
)

table(
  seu.obj.subset@meta.data$RNA_snn_res.0.2,
  seu.obj.subset$Patient_ID
)


table(
  seu.obj.subset@meta.data$RNA_snn_res.0.8,
  seu.obj.subset$Patient_ID
)

DimPlot(
  seu.obj.subset,
  reduction = "umap",
  group.by = "Patient_ID"
)

Idents(seu.obj.subset) = seu.obj.subset@meta.data$RNA_snn_res.0.2

DimPlot(
  seu.obj.subset,
  reduction = "umap",
  label = TRUE
)

seu.obj.subset <- JoinLayers(
  seu.obj.subset,
  assay = "RNA"
)

DefaultAssay(seu.obj.subset) <- "RNA"

markers <- FindAllMarkers(
  seu.obj.subset,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  verbose = TRUE
)

Reductions(seu.obj.subset)


genes_to_plot <- c("ACTA2","PECAM1","VWF", "ENG","CMA1","MS4A2","TPSAB1","TPSB2",
                   "AR", "KRT19", "KRT18", "KRT8", "TP63", "KRT14", "KRT5", "LYZ",
                   "FCGR3A", "CSF1R", "CD68", "CD163", "CD14", "UCHL1", "HAVCR2",
                   "PDCD1", "CTLA4", "CD8A", "SELL", "PTPRC", "CD4", "BTLA", "IL2RA",
                   "IL7R", "CCR7", "CD28", "CD27", "SLAMF1", "DPP4", "CD7","CD2",
                   "CD3G","CD3E", "CD3D")



pdf("Chen_et_al_marker_genes_plot_by_cluster.pdf")
DotPlot(
  seu.obj.subset,
  features = genes_to_plot
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    )
  ) + coord_flip()

dev.off()


'''
Based on paper markers, following are the cluster identity
Using 0.2 resolution
C6 - T cells
C8 - Monolytic
C11- Basal/Intermediate
C2 & C3 - luminal
C5 - Mast
C4 - Endothelial
C7-fibroblast
'''


seu.obj.subset.identified <- subset(seu.obj.subset, idents = c(2,3,4,5,6,7,8,11))


seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 2  ] <- "Luminal"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 3  ] <- "Luminal"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 4  ] <- "Endothelial"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 5  ] <- "Mast"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 6  ] <- "T cells"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 7  ] <- "Fibroblast"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 8  ] <- "Monolytic"
seu.obj.subset.identified@meta.data$Cell_type[ seu.obj.subset.identified@meta.data$RNA_snn_res.0.2 == 11  ] <- "Basal/Intermediate"


View(seu.obj.subset.identified@meta.data)


pdf("Chen_et_al_UMAP_with_celltype_label.pdf")
DimPlot(
  seu.obj.subset.identified,
  reduction = "umap",
  label = TRUE,
  group.by = "Cell_type"
)
dev.off()



## Plot average expression of PERSIST
persist =read.csv(file="Persist_signature.csv",header = T,stringsAsFactors = F)

grep("Persist", colnames(all_rna_f@meta.data), value = TRUE)


#DefaultAssay(all_rna_f) <- "RNA"

persist_genes <- persist$Persist

persist_genes <- as.character(persist_genes)

persist_genes <- persist_genes[!is.na(persist_genes)]

seu.obj.subset.identified <- AddModuleScore(
  object = seu.obj.subset.identified,
  features = persist_genes,
  name = 'PERSIST'
)

View(seu.obj.subset.identified@meta.data)

metadata = seu.obj.subset.identified@meta.data

average_persist <- apply(metadata[,grep("PERSIST",colnames(metadata))],1,mean)

stopifnot(names(average_persist) == colnames(seu.obj.subset.identified))

seu.obj.subset.identified@meta.data$Persist_avg <- average_persist


Idents(seu.obj.subset.identified) <- seu.obj.subset.identified@meta.data$Cell_type

pdf("Module_Score_for_Persist_Signature_Chen_et_al_UMAP.pdf", width = 15)
FeaturePlot(
  seu.obj.subset.identified,
  features = "Persist_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()


## Run TSNE for this object
seu.obj.subset.identified <- RunTSNE(
  seu.obj.subset.identified,
  reduction = "mnn",
  dims = 1:50
)

pdf("Module_Score_for_Persist_Signature_Chen_et_al_TSNE.pdf", width = 15)
FeaturePlot(
  seu.obj.subset.identified,
  features = "Persist_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "tsne",
  label = TRUE
)
dev.off()




## Plot average expression of PROSGENESIS
prosgenesis =read.csv(file="PROSGenesis_gene_signature.csv",header = T,stringsAsFactors = F)

grep("Prosgenesis", colnames(seu.obj.subset.identified@meta.data), value = TRUE)


#DefaultAssay(all_rna_f) <- "RNA"

prosgenesis_genes <- prosgenesis$PROSGenesis

prosgenesis_genes <- as.character(prosgenesis_genes)

prosgenesis_genes <- prosgenesis_genes[!is.na(prosgenesis_genes)]

seu.obj.subset.identified <- AddModuleScore(
  object = seu.obj.subset.identified,
  features = prosgenesis_genes,
  name = 'PROSGENESIS'
)

View(seu.obj.subset.identified@meta.data)

metadata = seu.obj.subset.identified@meta.data

average_prosgenesis <- apply(metadata[,grep("PROSGENESIS",colnames(metadata))],1,mean)

stopifnot(names(average_prosgenesis) == colnames(seu.obj.subset.identified))

seu.obj.subset.identified@meta.data$Prosgenesis_avg <- average_prosgenesis


#Idents(seu.obj.subset.identified) <- seu.obj.subset.identified@meta.data$Cell_type

pdf("Module_Score_for_Prosgenesis_Signature_Chen_et_al_UMAP.pdf", width = 15)
FeaturePlot(
  seu.obj.subset.identified,
  features = "Prosgenesis_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()


## Run TSNE for this object
#seu.obj.subset.identified <- RunTSNE(
#  seu.obj.subset.identified,
#  reduction = "mnn",
#  dims = 1:50
#)

pdf("Module_Score_for_Prosgenesis_Signature_Chen_et_al_TSNE.pdf", width = 15)
FeaturePlot(
  seu.obj.subset.identified,
  features = "Prosgenesis_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "tsne",
  label = TRUE
)
dev.off()



## histogram of Persist average score

hist(seu.obj.subset.identified@meta.data$Persist_avg)

# Persist Vector to plot
x <- seu.obj.subset.identified@meta.data$Persist_avg

# Calculate 90th percentile
p90 <- quantile(x, 0.90, na.rm = TRUE)

# Plot histogram + 90th percentile line
pdf("Persist_signature_Chen_et_al_histogram.pdf", width = 15)
ggplot(data.frame(x = x), aes(x = x)) +
  geom_histogram(
    bins = 50,
    fill = "grey70",
    color = "black"
  ) +
  geom_vline(
    xintercept = p90,
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = p90,
    y = Inf,
    label = paste0("90th percentile = ", round(p90, 2)),
    vjust = 2,
    hjust = 1.1
  ) +
  theme_classic()

dev.off()




# Prosgenesis Vector to plot

prosgen_avg <- seu.obj.subset.identified@meta.data$Prosgenesis_avg

# Calculate 90th percentile
p90 <- quantile(prosgen_avg, 0.90, na.rm = TRUE)

# Plot histogram + 90th percentile line
pdf("Prosgenesis_signature_Chen_et_al_histogram.pdf", width = 15)
ggplot(data.frame(x = prosgen_avg), aes(x = prosgen_avg)) +
  geom_histogram(
    bins = 50,
    fill = "grey70",
    color = "black"
  ) +
  geom_vline(
    xintercept = p90,
    linetype = "dashed",
    linewidth = 1
  ) +
  annotate(
    "text",
    x = p90,
    y = Inf,
    label = paste0("90th percentile = ", round(p90, 2)),
    vjust = 2,
    hjust = 1.1
  ) +
  theme_classic()

dev.off()


Persist.score <- seu.obj.subset.identified@meta.data[,c("Patient_ID", "Persist_avg")]
Prosgenesis.score <- seu.obj.subset.identified@meta.data[,c("Patient_ID", "Prosgenesis_avg")]


patient.id <- unique(Persist.score$Patient_ID)

## Get fxn of cells above 90% for persist signature
table_Persist_Score_fxn <- matrix(0, nrow = 1, ncol = 2)
colnames(table_Persist_Score_fxn) = c("Above_90_percentile", "Below_90_percentile")

for(patient in patient.id){
  Persist.per.patient <- Persist.score[ Persist.score$Patient_ID == patient, ]
  persist.90.percentile <- quantile(Persist.per.patient$Persist_avg, 0.90, na.rm = TRUE)
  Persist.per.patient$percentile_threshold <- ifelse(Persist.per.patient$Persist_avg > persist.90.percentile, "Above_90_percentile", "Below_90_percentile")
  counts <- table(
    Persist.per.patient$Patient_ID,
    Persist.per.patient$percentile_threshold
  )
  fractions <- counts / rowSums(counts)
  table_Persist_Score_fxn <- rbind(table_Persist_Score_fxn, fractions)
}

table_Persist_Score_fxn = table_Persist_Score_fxn[-c(1),]

table_Persist_Score_fxn = as.data.frame(table_Persist_Score_fxn)

pdf("Persist_signature_percent_of_scores_abive_90th_percentile.pdf", width = 15)
barplot(
  table_Persist_Score_fxn$Above_90_percentile,
  names.arg = rownames(table_Persist_Score_fxn),
  ylab = "Fraction Above 90th percentile",
  xlab = "Patient",
  ylim = c(0, max(table_Persist_Score_fxn$Above_90_percentile)+0.05),
  las = 2
)
dev.off()



## Get fxn of cells above 90% for prosgenesis signature
table_Prosgenesis_Score_fxn <- matrix(0, nrow = 1, ncol = 2)
colnames(table_Prosgenesis_Score_fxn) = c("Above_90_percentile", "Below_90_percentile")

for(patient in patient.id){
  Prosgenesis.per.patient <- Prosgenesis.score[ Prosgenesis.score$Patient_ID == patient, ]
  Prosgenesis.90.percentile <- quantile(Prosgenesis.per.patient$Prosgenesis_avg, 0.90, na.rm = TRUE)
  Prosgenesis.per.patient$percentile_threshold <- ifelse(Prosgenesis.per.patient$Prosgenesis_avg > Prosgenesis.90.percentile, "Above_90_percentile", "Below_90_percentile")
  counts <- table(
    Prosgenesis.per.patient$Patient_ID,
    Prosgenesis.per.patient$percentile_threshold
  )
  fractions <- counts / rowSums(counts)
  table_Prosgenesis_Score_fxn <- rbind(table_Prosgenesis_Score_fxn, fractions)
}

table_Prosgenesis_Score_fxn = table_Prosgenesis_Score_fxn[-c(1),]

table_Prosgenesis_Score_fxn = as.data.frame(table_Prosgenesis_Score_fxn)

pdf("Prosgenesis_signature_percent_of_scores_abive_90th_percentile.pdf", width = 15)
barplot(
  table_Prosgenesis_Score_fxn$Above_90_percentile,
  names.arg = rownames(table_Prosgenesis_Score_fxn),
  ylab = "Fraction Above 90th percentile",
  xlab = "Patient",
  ylim = c(0, max(table_Prosgenesis_Score_fxn$Above_90_percentile)+0.05),
  las = 2
)
dev.off()



## Get GSVA plot
# Average expression
avg.exp <- AverageExpression(
  seu.obj.subset.identified,
  group.by = "Cell_type",
  assay = "RNA",
  layer = "data"
)$RNA


## Pathways to plot for fig 6
pathways.plot <- read.csv(file="Fig6b_gene_expression_signature.csv",header = T,stringsAsFactors = F)


# Run GSVA
library(GSVA)

dim(avg.exp)
head(colnames(avg.exp))
head(rownames(avg.exp))

geneSets <- split(
  pathways.plot$Gene,
  pathways.plot$Signature
)

geneSets <- lapply(
  geneSets,
  trimws
)

gsvapar <- gsvaParam(
  exprData = avg.exp,
  geneSets = geneSets,
  kcdf = "Gaussian"
)


gsva.res <- gsva(gsvapar)


gsva.res.df <- as.data.frame(gsva.res)



library(pheatmap)

pdf("GSVA_Score_Plot_Chen_et_al_cell_types.pdf", width = 15)
pheatmap(
  gsva.res[, c(1,3,4)],
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  breaks = seq(-0.75, 0.75, length.out = 101),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "GSVA Plot, Chen et al"
)
dev.off()




