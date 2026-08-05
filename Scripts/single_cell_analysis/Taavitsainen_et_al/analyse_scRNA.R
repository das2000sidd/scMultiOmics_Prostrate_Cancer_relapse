setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(dplyr)

## First analyse scRNA-seq

dmso_r = Read10X("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168668_RAW_RNA/DMSO/")
enz48_r = Read10X("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168668_RAW_RNA/ENZ48/")
resa_r = Read10X("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168668_RAW_RNA/RESA/")
resb_r = Read10X("/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE168668_RAW_RNA/RESB/")


dmso_s = CreateSeuratObject(dmso_r, assay = "RNA", min.cells = 3, min.features = 200)
enz48_s = CreateSeuratObject(enz48_r, assay = "RNA", min.cells = 3, min.features = 200)
resa_s = CreateSeuratObject(resa_r, assay = "RNA", min.cells = 3, min.features = 200)
resb_s = CreateSeuratObject(resb_r, assay = "RNA", min.cells = 3, min.features = 200)

rm(dmso_r, enz48_r, resa_r, resb_r)

## Standard pre-processing workflow
## QC and selecting cells for further analysis

dmso_s[["percent.mt"]] <- PercentageFeatureSet(dmso_s, pattern = "^MT-")
enz48_s[["percent.mt"]] <- PercentageFeatureSet(enz48_s, pattern = "^MT-")
resa_s[["percent.mt"]] <- PercentageFeatureSet(resa_s, pattern = "^MT-")
resb_s[["percent.mt"]] <- PercentageFeatureSet(resb_s, pattern = "^MT-")

VlnPlot(dmso_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(enz48_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(resa_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(resb_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

feature_Scatture = function(seurat_object){
  plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
}

feature_Scatture(dmso_s)
feature_Scatture(enz48_s)
feature_Scatture(resa_s)
feature_Scatture(resb_s)
## cor of mito and ncount is generally -ve and cor of feature and count is generally +ve
## filtering criteria for number of genes detected
hist(dmso_s@meta.data$nFeature_RNA) ## 2000 to 6000
hist(enz48_s@meta.data$nFeature_RNA) ## 2000 to 5000
hist(resa_s@meta.data$nFeature_RNA) ## 1000 to 5000
hist(resb_s@meta.data$nFeature_RNA) ## 1000 to 5000

## filtering criteria for number of RNA molecules detected
hist(dmso_s@meta.data$nCount_RNA) ## 20000 to 50000 detected
hist(enz48_s@meta.data$nCount_RNA) ## 5000 to 25000 detected
hist(resa_s@meta.data$nCount_RNA) ## 4000 to 24000 detected
hist(resb_s@meta.data$nCount_RNA) ## 5000 to 25000 detected


## filtering criteria for percentage of mitochondrial DNA detected
hist(dmso_s@meta.data$percent.mt) ## < 15%
hist(enz48_s@meta.data$percent.mt) ## < 15%
hist(resa_s@meta.data$percent.mt) ## < 15%
hist(resb_s@meta.data$percent.mt) ## < 20%

## now filtering of cells based on criteria estimated from above histograms 
dmso_s_f <- subset(dmso_s, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 15 & nCount_RNA > 20000 & nCount_RNA < 50000)
enz48_s_f <- subset(enz48_s, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA > 5000 & nCount_RNA < 25000)
resa_s_f <- subset(resa_s, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 15 & nCount_RNA > 4000 & nCount_RNA < 24000)
resb_s_f <- subset(resb_s, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 20 & nCount_RNA > 5000 & nCount_RNA < 25000)


## checking histogram after filtering and we can see much higherand narrower distribution of feature counts
hist(dmso_s_f@meta.data$nFeature_RNA) ## 2000 to 6000, slightly skewed
hist(enz48_s_f@meta.data$nFeature_RNA) ## 2000 to 5000, looks normal
hist(resa_s_f@meta.data$nFeature_RNA) ## 1000 to 5000, looks normal
hist(resb_s_f@meta.data$nFeature_RNA) ## 1000 to 5000, looks normal

## checking histogram after filtering and we can see much higherand narrower distribution of feature counts
hist(dmso_s_f@meta.data$nCount_RNA) ## 20000 to 50000 detected, skewed towards left
hist(enz48_s_f@meta.data$nCount_RNA) ## 5000 to 25000 detected, looks normal
hist(resa_s_f@meta.data$nCount_RNA) ## 4000 to 24000 detected, , looks normal
hist(resb_s_f@meta.data$nCount_RNA) ## 5000 to 25000 detected, looks normal


## filtering criteria for percentage of mitochondrial DNA detected
hist(dmso_s_f@meta.data$percent.mt) ## < 15%
hist(enz48_s_f@meta.data$percent.mt) ## < 15%
hist(resa_s_f@meta.data$percent.mt) ## < 15%
hist(resb_s_f@meta.data$percent.mt) ## < 20%


dmso_s_f@meta.data$Sample_id = "DMSO"
enz48_s_f@meta.data$Sample_id = "ENZ48"
resa_s_f@meta.data$Sample_id = "RESA"
resb_s_f@meta.data$Sample_id = "RESB"


## merge alll samples together
all_rna_f = merge(dmso_s_f, y = c(enz48_s_f, resa_s_f, resb_s_f), add.cellids = c("DMSO", "ENZ48", "RESA", "RESB"), project = "Prostate_Cancer")

View(all_rna_f@meta.data)

all_rna_f[["RNA"]] <- JoinLayers(all_rna_f[["RNA"]])


## do strandard procesisng to check for batch effects
all_rna_f <- NormalizeData(all_rna_f, normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE)
all_rna_f <- FindVariableFeatures(all_rna_f, selection.method = "vst", nfeatures = 2000,verbose = TRUE)
#variable.genes=VariableFeatures(all_rna_f,selection.method = "vst",verbose=TRUE)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(all_rna_f), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(all_rna_f)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all_rna_f <- ScaleData(all_rna_f,verbose = TRUE)
all_rna_f <- RunPCA(all_rna_f,verbose = TRUE,ndims.print = 1:50)


VizDimLoadings(all_rna_f, dims = 1:2, reduction = "pca")

DimPlot(all_rna_f, reduction = "pca") + NoLegend()

DimHeatmap(all_rna_f, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(all_rna_f, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(all_rna_f) ## 10 perhaps

cumsum(all_rna_f@reductions$pca@stdev/sum(all_rna_f@reductions$pca@stdev))


# Cluster the cells
all_rna_f <- FindNeighbors(all_rna_f, dims = 1:45, verbose = TRUE)
all_rna_f <- FindClusters(all_rna_f, dims = 1:45, verbose = TRUE, resolution = 0.5, leiden_method = "igraph",algorithm = 4)


# Look at cluster IDs of the first 5 cells
head(Idents(all_rna_f), 5)



all_rna_f <- RunUMAP(all_rna_f, dims = 1:45, verbose = TRUE)


DimPlot(all_rna_f, reduction = "umap")

DimPlot(all_rna_f, reduction = "umap", group.by = "Sample_id")
## Grouping is currently by treatment instead of celll type after clustering
## Hence batch correction needed

library(SeuratData)
library(patchwork)


all_rna_f[["RNA"]] <- split(all_rna_f[["RNA"]], f = all_rna_f$Sample_id)

all_rna_f <- NormalizeData(all_rna_f)
all_rna_f <- FindVariableFeatures(all_rna_f)
all_rna_f <- ScaleData(all_rna_f)
all_rna_f <- RunPCA(all_rna_f)

all_rna_f <- FindNeighbors(all_rna_f, dims = 1:30, reduction = "pca")
all_rna_f <- FindClusters(all_rna_f, resolution = 2, cluster.name = "unintegrated_clusters")

all_rna_f <- RunUMAP(all_rna_f, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(all_rna_f, reduction = "umap.unintegrated", group.by = c("Sample_id")) ## Batch effect



all_rna_f <- IntegrateLayers(
  object = all_rna_f, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE
)

Reductions(all_rna_f) ## confirmed cca exists


all_rna_f <- RunUMAP(
  all_rna_f,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = TRUE
)

DimPlot(
  all_rna_f,
  reduction = "umap",
  group.by = c("Sample_id"),
  label = FALSE,
  
)


all_rna_f <- FindNeighbors(
  all_rna_f,
  reduction = "integrated.cca",
  dims = 1:30,
  verbose = TRUE
)

all_rna_f <- FindClusters(
  all_rna_f,
  resolution = c(0.2,0.4,0.6,0.8),
  algorithm = 4,
  leiden_method = "igraph"
)

all_rna_f <- FindClusters(
  all_rna_f,
  resolution = c(1,1.2,1.4),
  algorithm = 4,
  leiden_method = "igraph"
)


DimPlot(
  all_rna_f,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.2", "RNA_snn_res.0.4", "RNA_snn_res.0.6", "RNA_snn_res.0.8"),
  label = TRUE
) ## clear separation till res 0.8

DimPlot(
  all_rna_f,
  reduction = "umap",
  group.by = c("RNA_snn_res.1", "RNA_snn_res.1.2", "RNA_snn_res.1.4"),
  label = TRUE
) ## too many



DimPlot(
  all_rna_f,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.8"),
  split.by = c("Sample_id"),
  label = FALSE
)



## plot expression of genes from cells
goi <- c("BIRC5", "KRAS", "PLK1", "SS18", "WDR77", "TGIF1","TOP2A",
         "AURKA", "CCNA2", "STMN1", "BUB1", "BUB1B", "LMNB1", "FOXB1",
         "ZWINT", "WEE1", "BRCA2", "RAD51C", "CTCF", "POLQ", "FANC1",
         "ATAD2", "EZH2", "HNRNPA2B1", "HIST1H1B", "GMNN", "CEBPB",
         "EIF4EBP1","HES6")



Idents(all_rna_f) <- all_rna_f@meta.data$RNA_snn_res.0.8


DotPlot(
  all_rna_f,
  features = goi
) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 2),
    oob = scales::squish
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



## Expression of persist signature genes
goi.persist <- c("CDC20", "MKI67", "CCNB2", "CENPF", "DLGAP5", "HMMR",
                 "PLK1", "CCNB1", "PTTG1", "CENPE", "UBE2S","NCAPD2",
                 "TUBA1C","DYNLL1", "CD81","CKAP5", "ODC1","TMEM54",
                 "BDP1", "PPP1R14B", "PSME4", "MYC")

DotPlot(all_rna_f, features = goi.persist) + RotatedAxis() + scale_color_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0
)+theme(
  axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle("Persist signature") ## looks okay


saveRDS(all_rna_f, "Batch_corrected_clustered_scRNAseq_object.rds")


## Integrating VCaP et al with current data
vcap <- read.table(file="VCAP_DMSO_1_counts.txt",header = T,sep="\t",stringsAsFactors = F)
vcap_enz48 <- read.table(file="VCAP_ENZ48_counts.txt",header = T,sep="\t",stringsAsFactors = F)

head(vcap)

vcap[1:6,1:6]
vcap_enz48[1:6,1:6]

rownames(vcap) <- vcap$GENE

rownames(vcap_enz48) <- vcap_enz48$GENE


vcap <- vcap[,-c(1)]

vcap_enz48 <- vcap_enz48[,-c(1)]


seu.vcap <- CreateSeuratObject(counts = vcap, min.cells = 3, min.features = 200)
seu.vcap


seu.vcap.enz48 <- CreateSeuratObject(counts = vcap_enz48, min.cells = 3, min.features = 200)
seu.vcap.enz48


seu.vcap$Sample_Id = "VCaP"

seu.vcap.enz48$Sample_Id = "VCaP_ENZ48"

## Clustering VCaP
seu.vcap <- NormalizeData(seu.vcap, verbose = TRUE)

seu.vcap <- FindVariableFeatures(seu.vcap, verbose = TRUE)

seu.vcap <- ScaleData(seu.vcap, verbose = TRUE)

seu.vcap <- RunPCA(seu.vcap, verbose = TRUE)

seu.vcap <- RunUMAP(seu.vcap, verbose = TRUE,reduction = "pca", dims = 1:30)

seu.vcap <- FindNeighbors(seu.vcap, dims = 1:30)



seu.vcap <- FindClusters(
  seu.vcap,
  resolution = c(0.2,0.4,0.6,0.8,1),
  algorithm = 4,
  random.seed = 123,
  method = "igraph"
)


anchors.for.vcap <- FindTransferAnchors(
  reference = all_rna_f,
  query = seu.vcap,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:30,
  verbose = TRUE
)



predictions.vcap <- TransferData(
  anchorset = anchors.for.vcap,
  refdata = Idents(all_rna_f),
  dims = 1:30,
  verbose = TRUE
)

seu.vcap@meta.data$Prediction_from_LNCap <- predictions.vcap$predicted.id


Idents(seu.vcap) <- seu.vcap@meta.data$Prediction_from_LNCap


pdf("Label_transfer_from_paper_to_VCaP.pdf")
DimPlot(
  seu.vcap,
  reduction = "umap",
  group.by = c("Prediction_from_LNCap"),
  label = FALSE
)
dev.off()



## Clustering VCaP- ENZ48
seu.vcap.enz48 <- NormalizeData(seu.vcap.enz48, verbose = TRUE)

seu.vcap.enz48 <- FindVariableFeatures(seu.vcap.enz48, verbose = TRUE)

seu.vcap.enz48 <- ScaleData(seu.vcap.enz48, verbose = TRUE)

seu.vcap.enz48 <- RunPCA(seu.vcap.enz48, verbose = TRUE)

seu.vcap.enz48 <- RunUMAP(seu.vcap.enz48, verbose = TRUE,reduction = "pca", dims = 1:30)

seu.vcap.enz48 <- FindNeighbors(seu.vcap.enz48, dims = 1:30)



seu.vcap.enz48 <- FindClusters(
  seu.vcap.enz48,
  resolution = c(0.2,0.4,0.6,0.8,1),
  algorithm = 4,
  random.seed = 123,
  method = "igraph"
)


anchors.for.vcap <- FindTransferAnchors(
  reference = all_rna_f,
  query = seu.vcap.enz48,
  reference.reduction = "pca",
  normalization.method = "LogNormalize",
  dims = 1:30,
  verbose = TRUE
)



predictions.vcap <- TransferData(
  anchorset = anchors.for.vcap,
  refdata = Idents(all_rna_f),
  dims = 1:30,
  verbose = TRUE
)

seu.vcap.enz48@meta.data$Prediction_from_LNCap <- predictions.vcap$predicted.id


Idents(seu.vcap.enz48) <- seu.vcap.enz48@meta.data$Prediction_from_LNCap


pdf("Label_transfer_from_paper_to_VCaP_ENZ48.pdf")
DimPlot(
  seu.vcap.enz48,
  reduction = "umap",
  group.by = c("Prediction_from_LNCap"),
  label = FALSE
)
dev.off()



## plotting expression score of Persist in VCaP
persist = read.csv(file="Persist_signature.csv",header = T,stringsAsFactors = F)

grep("Persist", colnames(all_rna_f@meta.data), value = TRUE)


#DefaultAssay(all_rna_f) <- "RNA"

persist_genes <- persist$Persist

persist_genes <- as.character(persist_genes)

persist_genes <- persist_genes[!is.na(persist_genes)]

seu.vcap <- AddModuleScore(
  object = seu.vcap,
  features = persist_genes,
  name = 'PERSIST'
)

View(seu.vcap@meta.data)

metadata = seu.vcap@meta.data

average_persist <- apply(metadata[,grep("PERSIST",colnames(metadata))],1,mean)

stopifnot(names(average_persist) == colnames(seu.vcap))

seu.vcap@meta.data$Persist_avg <- average_persist


#Idents(seu.vcap) <- seu.vcap@meta.data$Cell_type

pdf("Module_Score_for_Persist_Signature_VCaP_UMAP.pdf", width = 15)
FeaturePlot(
  seu.vcap,
  features = "Persist_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()



## plotting expression score of Persist in VCaP-ENZ48
seu.vcap.enz48 <- AddModuleScore(
  object = seu.vcap.enz48,
  features = persist_genes,
  name = 'PERSIST'
)

View(seu.vcap.enz48@meta.data)

metadata = seu.vcap.enz48@meta.data

average_persist <- apply(metadata[,grep("PERSIST",colnames(metadata))],1,mean)

stopifnot(names(average_persist) == colnames(seu.vcap.enz48))

seu.vcap.enz48@meta.data$Persist_avg <- average_persist


#Idents(seu.vcap.enz48) <- seu.vcap.enz48@meta.data$Cell_type

pdf("Module_Score_for_Persist_Signature_VCaP_ENZ48_UMAP.pdf", width = 15)
FeaturePlot(
  seu.vcap.enz48,
  features = "Persist_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()




## plotting expression score of Prosgenesis in VCaP
prosgenesis = read.csv(file="PROSGenesis_gene_signature.csv",header = T,stringsAsFactors = F)

grep("Prosgenesis", colnames(all_rna_f@meta.data), value = TRUE)

present_genes <- intersect(
  prosgenesis_genes,
  rownames(seu.vcap)
)

length(prosgenesis_genes)
length(present_genes)

present_genes


#DefaultAssay(all_rna_f) <- "RNA"

prosgenesis_genes <- prosgenesis$PROSGenesis

prosgenesis_genes <- as.character(prosgenesis_genes)

prosgenesis_genes <- prosgenesis_genes[!is.na(prosgenesis_genes)]

seu.vcap <- AddModuleScore(
  object = seu.vcap,
  features = list(prosgenesis_genes),
  name = 'PROS'
)

View(seu.vcap@meta.data)

metadata = seu.vcap@meta.data

average_prosgenesis <- metadata[,grep("PROS",colnames(metadata))]

stopifnot(names(average_prosgenesis) == colnames(seu.vcap))

seu.vcap@meta.data$Prosgenesis_avg <- average_prosgenesis


#Idents(seu.vcap) <- seu.vcap@meta.data$Cell_type

pdf("Module_Score_for_Prosgenesis_Signature_VCaP_UMAP.pdf", width = 15)
FeaturePlot(
  seu.vcap,
  features = "Prosgenesis_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()


## Plot progenesis genes in VCAP ENZ48 treated sample
seu.vcap.enz48 <- AddModuleScore(
  object = seu.vcap.enz48,
  features = list(prosgenesis_genes),
  name = 'PROS'
)

View(seu.vcap.enz48@meta.data)

metadata = seu.vcap.enz48@meta.data

average_prosgenesis <- metadata[,grep("PROS",colnames(metadata))]

names(average_prosgenesis) <- rownames(metadata)

stopifnot(names(average_prosgenesis) == colnames(seu.vcap.enz48))

seu.vcap.enz48@meta.data$Prosgenesis_avg <- average_prosgenesis


#Idents(seu.vcap.enz48) <- seu.vcap.enz48@meta.data$Cell_type

pdf("Module_Score_for_Prosgenesis_Signature_VCaP_ENZ48_UMAP.pdf", width = 15)
FeaturePlot(
  seu.vcap.enz48,
  features = "Prosgenesis_avg",
  cols = c("lightgrey", "blue" ,"yellow" , "red"),
  reduction = "umap",
  label = TRUE
)
dev.off()



