setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/GSE130534_RAW_BULK_RNA_SEQ")


files <- list.files(pattern = "\\.txt$")

## make sample sheet
sample_info <- data.frame(
  file = files,
  sample = sub("\\.txt$", "", files),
  stringsAsFactors = FALSE
)

sample_info

counts <- lapply(files, read.delim, header = TRUE)
names(counts) <- sub("\\.txt$", "", basename(files))


count_matrix <- do.call(
  cbind,
  lapply(counts, function(x) x[, 2])
)

rownames(count_matrix) <- counts[[1]][, 1]
colnames(count_matrix) <- sub("\\.txt$", "", basename(files))

#count_matrix <- as.data.frame(count_matrix)
head(sample_info)

sample_info$group <- sapply(
  strsplit(sample_info$sample, "_"),
  function(x) paste(x[1:2], collapse = "_")
)

stopifnot(colnames(count_matrix)==sample_info$sample)


library(DESeq2)
library(edgeR)
#dds <- DESeqDataSetFromTximport(txi, colData=samp_tab, ~ Combined_Id)
dds <- DESeqDataSetFromMatrix(count_matrix, sample_info, design = ~ group)

dds <- estimateSizeFactors(dds)

keep <- rowSums(counts(dds) >= 5) >= (ncol(dds)*0.1)

table(keep)

dds <- dds[keep,]

norm_counts <- counts(dds, normalized = TRUE)

log_norm <- log2(norm_counts + 1)

genesets_plot <- read.csv(file="TCGA_PRAD_GSVA_genesets_to_plot.csv", header = T,stringsAsFactors = F)

genesets_plot$Gene <- trimws(genesets_plot$Gene)

genesets_plot$Entrez <- mapIds(org.Hs.eg.db, genesets_plot$Gene,keytype="SYMBOL", column="ENTREZID")
genesets_plot$Ensembl <- mapIds(org.Hs.eg.db, genesets_plot$Entrez,keytype="ENTREZID", column="ENSEMBL")

genesets_plot$Entrez <- as.character(genesets_plot$Entrez)
genesets_plot$Ensembl <- as.character(genesets_plot$Ensembl)

genesets_plot$Ensembl_from_Symbol <- mapIds(org.Hs.eg.db, genesets_plot$Gene,keytype="SYMBOL", column="ENSEMBL")

genesets_supp <- genesets_plot[genesets_plot$Ensembl != "NULL", ]

genesets_list <- split(genesets_supp$Gene,
                       genesets_supp$Group)


library(GSVA)

param <- gsvaParam(
  expr = log_norm,
  geneSets = c(genesets_list),
  kcdf = "Poisson"
)

gsva_scores <- gsva(param)


gsva_scores_reordered <- gsva_scores[c(11, 3, 13, 12, 7, 5, 14, 8, 9, 4, 6, 1, 10, 2), ]

ph <- pheatmap(
  gsva_scores_reordered,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  breaks = seq(-0.5, 0.5, length.out = 101),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "GSVA Plot Bulk RNA-seq",
  silent = TRUE
)

ph


## get average socre matrix
# Average every 3 consecutive columns
gsva_scores_avg <- sapply(seq(1, ncol(gsva_scores_reordered), by = 3),
                          function(i) {
                            rowMeans(gsva_scores_reordered[, i:(i+2)])
                          })

# Preserve row names
rownames(gsva_scores_avg) <- rownames(gsva_scores_reordered)

# Create treatment names from first replicate of each group
colnames(gsva_scores_avg) <- sub("_R[1-3]$", "",
                                 colnames(gsva_scores_reordered)[seq(1, ncol(gsva_scores_reordered), by = 3)])


pdf("GSVA_Heatmap_in_BulkRNASeq_of_single_cell_derived_clusters_fig5a.pdf", width = 15)

ph.avg <- pheatmap(
  gsva_scores_avg,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  breaks = seq(-0.5, 0.5, length.out = 101),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "GSVA Plot Bulk RNA-seq",
  silent = TRUE
)

print(ph.avg)


dev.off()


