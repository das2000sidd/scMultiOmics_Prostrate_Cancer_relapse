setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(Signac)
library(chromVAR)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(ggpubr)

aucell <- read.csv(
  "aucell.csv",
  check.names = FALSE
)

rownames(aucell) <- aucell$Cell

all_rna_f <- readRDS("Batch_corrected_clustered_scRNAseq_object.rds")

Layers(all_rna_f[["RNA"]])

all_rna_f <- JoinLayers(
  all_rna_f,
  assay = "RNA"
)

meta.seurat <- all_rna_f@meta.data

meta.seurat$Cell <- rownames(meta.seurat)

sum(rownames(aucell) %in% meta.seurat$Cell)

aucell_meta <- cbind(aucell, meta.seurat[rownames(aucell), "Sample_id"])

colnames(aucell_meta)[ncol(aucell_meta)] <- "treatment"

prostate_cancer_tfs <- c(
  "AR(+)",
  "FOXO1(+)",
  "NKX3-1(+)",
  "MYC(+)",
  "E2F1(+)",
  "STAT1(+)",
  "JUN(+)",
  "FOS(+)"
)

prostate_cancer_tfs %in% colnames(aucell_meta)


# SCENIC treatment medians
scenic_medians <- aucell_meta %>%
  group_by(treatment) %>%
  summarise(
    across(
      all_of(prostate_cancer_tfs),
      ~median(.x, na.rm = TRUE)
    )
  )


scenic_medians

# Calculate log2 fold changes relative to DMSO
scenic_fc <- scenic_medians %>%
  pivot_longer(
    cols = -treatment,
    names_to = "TF",
    values_to = "AUC"
  ) %>%
  group_by(TF) %>%
  mutate(
    DMSO = AUC[treatment == "DMSO"][1],
    FC_vs_DMSO = AUC / DMSO,
    log2FC_vs_DMSO = log2(FC_vs_DMSO)
  ) %>%
  ungroup()

scenic_fc


# Statistical testing
kw_results <- lapply(prostate_cancer_tfs, function(tf) {
  
  test <- kruskal.test(
    as.formula(paste0("`", tf, "` ~ treatment")),
    data = aucell_meta
  )
  
  data.frame(
    TF = tf,
    pvalue = test$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(pvalue, method = "BH")
  )

kw_results


# Pairwise comparison
library(rstatix)

dunn_results <- lapply(prostate_cancer_tfs, function(tf) {
  
  aucell_meta %>%
    dunn_test(
      formula = as.formula(
        paste0("`", tf, "` ~ treatment")
      ),
      p.adjust.method = "BH"
    ) %>%
    mutate(TF = tf)
  
}) %>%
  bind_rows()

dunn_results


# Plot SCENIC trajectories
scenic_plot_data <- aucell_meta %>%
  select(treatment, all_of(prostate_cancer_tfs)) %>%
  pivot_longer(
    cols = all_of(prostate_cancer_tfs),
    names_to = "TF",
    values_to = "AUC"
  )

ggplot(
  scenic_plot_data,
  aes(x = treatment, y = AUC)
) +
  geom_boxplot(
    outlier.size = 0.2
  ) +
  facet_wrap(
    ~TF,
    scales = "free_y"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "SCENIC regulon activity"
  )


atac_int <- readRDS("combined_atac_harmony_integrated.rds")

DefaultAssay(atac_int)


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

DefaultAssay(atac_int) <- "peaks"


# Calculate GC content for all peaks
atac_int <- RegionStats(
  object = atac_int,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


# Add motif information
atac_int <- AddMotifs(
  object = atac_int,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm.use,
  verbose = TRUE
)


## Run chromVAR
atac_int <- RunChromVAR(
  object = atac_int,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


motifs <- Motifs(atac_int[["peaks"]])

motifs

tfs <- c("AR", "JUN", "FOS", "E2F1", "MYC", "NKX3-1", "FOXO1", "STAT1")

motifs@motif.names[
  sapply(tfs, function(x) any(grepl(x, motifs@motif.names, ignore.case = TRUE)))
]

motif_matrix <- motifs@data

dim(motif_matrix)

counts <- GetAssayData(
  atac_int,
  assay = "peaks",
  layer = "counts"
)

dim(counts)


# Create the chromVAR SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(
    counts = counts
  )
)

rowRanges(se) <- granges(atac_int[["peaks"]])

colData(se) <- DataFrame(
  cell = colnames(counts)
)

## Identify peaks with any accessibility
peak_totals <- Matrix::rowSums(counts)

summary(peak_totals)

sum(peak_totals == 0)
sum(peak_totals > 0)


## Remove zero-count peaks
keep <- peak_totals > 0

counts_filt <- counts[keep, ]

motif_matrix_filt <- motif_matrix[keep, ]

dim(counts_filt)
dim(motif_matrix_filt)

identical(
  rownames(counts_filt),
  rownames(motif_matrix_filt)
)

# Recreate the SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(
    counts = counts_filt
  )
)

rowRanges(se) <- granges(atac_int[["peaks"]])[keep, ]

colData(se) <- S4Vectors::DataFrame(
  cell = colnames(counts_filt)
)

## figuring chromVar computeDeviations length(bias) == length(fragments_per_peak) is not TRUE
counts <- GetAssayData(
  atac_int,
  assay = "peaks",
  layer = "counts"
)

motifs <- Motifs(atac_int[["peaks"]])
motif_matrix <- motifs@data

cat("Counts:", dim(counts), "\n")
cat("Motifs:", dim(motif_matrix), "\n")
cat("Cells:", ncol(counts), "\n")

cat("Zero-count peaks:", sum(Matrix::rowSums(counts) == 0), "\n")

cat("Peak ranges:", length(GenomicRanges::granges(atac_int[["peaks"]])), "\n")


keep <- Matrix::rowSums(counts) > 0

counts_filt <- counts[keep, ]

motif_matrix_filt <- motif_matrix[keep, ]

peaks_filt <- granges(atac_int[["peaks"]])[keep, ]

cat("Peaks retained:", nrow(counts_filt), "\n")
cat("Motif rows:", nrow(motif_matrix_filt), "\n")
cat("Ranges:", length(peaks_filt), "\n")


se <- SummarizedExperiment(
  assays = list(
    counts = counts_filt
  ),
  rowRanges = peaks_filt
)

colData(se)$cell <- colnames(counts_filt)


# Calculate GC bias
se <- addGCBias(
  se,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

dim(se)

head(rowData(se))

summary(rowData(se)$bias)

#Calculate background peaks
bg <- getBackgroundPeaks(
  se
)

dim(bg)

#run chromVAR
dev <- computeDeviations(
  object = se,
  annotations = motif_matrix_filt,
  background_peaks = bg
)

dev

z.score.motifs <- deviationScores(dev)

dim(z.score.motifs)


# check motif names
rownames(z.score.motifs)


# integrate ATAC results with pySCENIC AUCell results.
chromvar_df <- as.data.frame(t(z.score.motifs))

chromvar_df$Cell <- rownames(chromvar_df)

head(chromvar_df[, 1:5])

chromvar_df <- chromvar_df %>%
  rename(
    E2F1 = MA0024.3,
    MYC  = MA0147.3,
    JUND = MA0491.1,
    ETV1 = MA0761.1,
    HOXB13 = MA0901.1,
    ELK4 = MA0076.2,
    ELF1 = MA0473.2
  )



# Add treatment
meta_atac <- atac_int@meta.data

meta_atac$Cell <- rownames(meta_atac)

table(meta_atac$Sample_Id)

chromvar_df <- left_join(chromvar_df, meta_atac, by=c("Cell"))

table(chromvar_df$Sample_Id, useNA = "ifany")


motif_to_tf <- sapply(
  motifs@motif.names,
  function(x) x[1]
)

motif.id.tf.name.map <- data.frame(
  motif_id = names(motif_to_tf),
  TF = unname(motif_to_tf)
)

rownames(motif.id.tf.name.map) <- motif.id.tf.name.map$motif_id

chromvar.df.order <- colnames(chromvar_df)[ grep("^MA", colnames(chromvar_df)) ]

motif.id.tf.name.map.o <- motif.id.tf.name.map[chromvar.df.order, ]

stopifnot(motif.id.tf.name.map.o$motif_id == colnames(chromvar_df)[ grep("^MA", colnames(chromvar_df)) ])

colnames(chromvar_df)[1:17] <- motif.id.tf.name.map.o$TF


# Plot the ATAC activity
tf.to.plot <- colnames(chromvar_df)[1:17]


## comparisons to add to plot
comparisons <- list(
  c("DMSO", "ENZ48"),
  c("DMSO", "RESA"),
  c("DMSO", "RESB"),
  c("ENZ48", "RESA"),
  c("ENZ48", "RESB")
)

chromvar_df$Sample_Id <- factor(
  chromvar_df$Sample_Id,
  levels = c("DMSO", "ENZ48", "RESA", "RESB")
)


for(tf in tf.to.plot){

  ggplot(
    chromvar_df,
    aes(x = Sample_Id, y = .data[[tf]])
  ) +
    geom_boxplot(outlier.size = 0.3, width=0.7) +
    geom_jitter(
      width = 0.15,
      size = 0.3,
      alpha = 0.3
    ) +
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      p.adjust.method = "BH",
      label = "p.signif",
      hide.ns = TRUE
    )
    theme_classic() +
    labs(
      x = NULL,
      y = paste0(tf, " chromVAR deviation"),
      title = paste0(tf, " motif accessibility")
    )
  
  ggsave(paste0(tf, "_ATAC_acitivity.pdf"), width = 5, height = 5)
}

#ggplot(
#  chromvar_df,
#  aes(x = Sample_Id, y = SOX9)
#) +
#  geom_boxplot(outlier.size = 0.3) +
#  theme_classic() +
#  labs(
#    x = NULL,
#    y = paste0(tf, " chromVAR deviation"),
#    title = paste0(tf, " motif accessibility")
#  )

ggsave(paste0("SOX9", "_ATAC_acitivity_test.pdf"), width = 5, height = 5)







## derive overall how each TF behaves in differences treatments using SCENIC and chromVAR
## calculate SCENIC treament medians
treatments <- c("DMSO", "ENZ48", "RESA", "RESB")

tf.to.plot <- c(
  "AR(+)",
  "NKX3-1(+)",
  "MYC(+)",
  "E2F1(+)",
  "JUN(+)",
  "FOS(+)",
  "FOXO1(+)",
  "STAT1(+)"
)

## first get medians for SCENIC
scenic.medians <- matrix(
  NA,
  nrow = length(tf.to.plot),
  ncol = length(treatments),
  dimnames = list(tf.to.plot, treatments)
)

for (tf in tf.to.plot) {
  
  for (tr in treatments) {
    
    x <- aucell_meta[aucell_meta$treatment == tr, tf]
    
    scenic.medians[tf, tr] <- median(
      x,
      na.rm = TRUE
    )
  }
}

scenic.medians


## calculate SCENIC log2fc
scenic.log2FC <- matrix(
  NA,
  nrow = length(tf.to.plot),
  ncol = 3,
  dimnames = list(
    tf.to.plot,
    c(
      "ENZ48_vs_DMSO",
      "RESA_vs_DMSO",
      "RESB_vs_DMSO"
    )
  )
)

for (tf in tf.to.plot) {
  
  scenic.log2FC[tf, "ENZ48_vs_DMSO"] <-
    log2(
      scenic.medians[tf, "ENZ48"] /
        scenic.medians[tf, "DMSO"]
    )
  
  scenic.log2FC[tf, "RESA_vs_DMSO"] <-
    log2(
      scenic.medians[tf, "RESA"] /
        scenic.medians[tf, "DMSO"]
    )
  
  scenic.log2FC[tf, "RESB_vs_DMSO"] <-
    log2(
      scenic.medians[tf, "RESB"] /
        scenic.medians[tf, "DMSO"]
    )
}

scenic.log2FC


## some TFs have zero median., hence use difference
scenic.diff <- matrix(
  NA,
  nrow = length(tf.to.plot),
  ncol = 3,
  dimnames = list(
    tf.to.plot,
    c(
      "ENZ48_vs_DMSO",
      "RESA_vs_DMSO",
      "RESB_vs_DMSO"
    )
  )
)

for (tf in tf.to.plot) {
  
  scenic.diff[tf, "ENZ48_vs_DMSO"] <-
    scenic.medians[tf, "ENZ48"] -
    scenic.medians[tf, "DMSO"]
  
  scenic.diff[tf, "RESA_vs_DMSO"] <-
    scenic.medians[tf, "RESA"] -
    scenic.medians[tf, "DMSO"]
  
  scenic.diff[tf, "RESB_vs_DMSO"] <-
    scenic.medians[tf, "RESB"] -
    scenic.medians[tf, "DMSO"]
}

scenic.diff



## ChromVAR treatment medians
chromvar.tfs <- c(
  "E2F1",
  "JUND",
  "ELK4",
  "ETV1",
  "HOXB13",
  "RELA",
  "SOX9",
  "ELF1",
  "GRHL2"
)

#chromvar.tfs[which(chromvar.tfs %in% colnames(chromvar_df))]

chromvar.medians <- matrix(
  NA,
  nrow = length(chromvar.tfs),
  ncol = length(treatments),
  dimnames = list(chromvar.tfs, treatments)
)

for (tf in chromvar.tfs) {
  
  for (tr in treatments) {
    
    x <- chromvar_df[
      chromvar_df$Sample_Id == tr,
      tf
    ]
    
    chromvar.medians[tf, tr] <- median(
      x,
      na.rm = TRUE
    )
  }
}

chromvar.medians



# ChromVAR treatment differences
chromvar.diff <- matrix(
  NA,
  nrow = length(chromvar.tfs),
  ncol = 3,
  dimnames = list(
    chromvar.tfs,
    c(
      "ENZ48_vs_DMSO",
      "RESA_vs_DMSO",
      "RESB_vs_DMSO"
    )
  )
)

for (tf in chromvar.tfs) {
  
  chromvar.diff[tf, "ENZ48_vs_DMSO"] <-
    chromvar.medians[tf, "ENZ48"] -
    chromvar.medians[tf, "DMSO"]
  
  chromvar.diff[tf, "RESA_vs_DMSO"] <-
    chromvar.medians[tf, "RESA"] -
    chromvar.medians[tf, "DMSO"]
  
  chromvar.diff[tf, "RESB_vs_DMSO"] <-
    chromvar.medians[tf, "RESB"] -
    chromvar.medians[tf, "DMSO"]
}

chromvar.diff


#Wilcoxon P-values for chromVAR
chromvar.p <- list()

comparisons <- list(
  c("DMSO", "ENZ48"),
  c("DMSO", "RESA"),
  c("DMSO", "RESB"),
  c("ENZ48", "RESA"),
  c("ENZ48", "RESB"),
  c("RESA", "RESB")
)

for (tf in chromvar.tfs) {
  
  for (comp in comparisons) {
    
    g1 <- chromvar_df[
      chromvar_df$Sample_Id == comp[1],
      tf
    ]
    
    g2 <- chromvar_df[
      chromvar_df$Sample_Id == comp[2],
      tf
    ]
    
    test <- wilcox.test(
      g1,
      g2,
      exact = FALSE
    )
    
    chromvar.p[[paste(tf, comp[1], comp[2], sep = "_")]] <-
      test$p.value
  }
}

# Convert chromvar.p to a dataframe
chromvar.p.df <- data.frame(
  TF = character(),
  Group1 = character(),
  Group2 = character(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (tf in chromvar.tfs) {
  
  for (comp in comparisons) {
    
    key <- paste(
      tf,
      comp[1],
      comp[2],
      sep = "_"
    )
    
    chromvar.p.df <- rbind(
      chromvar.p.df,
      data.frame(
        TF = tf,
        Group1 = comp[1],
        Group2 = comp[2],
        Pvalue = chromvar.p[[key]]
      )
    )
  }
}

# Correction for multiple testing
chromvar.p.df$FDR <- p.adjust(
  chromvar.p.df$Pvalue,
  method = "BH"
)

chromvar.p.df <- chromvar.p.df[
  order(chromvar.p.df$FDR),
]

## Wilcoxon P-values for aucell
aucell_meta <- as.data.frame(aucell_meta)

aucell_meta$treatment <- as.character(aucell_meta$treatment)

# Check
table(aucell_meta$treatment)


# Empty results table
scenic.p.df <- data.frame(
  TF = character(),
  Group1 = character(),
  Group2 = character(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (tf in tf.to.plot) {
  
  # Make sure TF column is numeric
  aucell_meta[[tf]] <- as.numeric(aucell_meta[[tf]])
  
  for (comp in comparisons) {
    
    idx1 <- aucell_meta$treatment == comp[1]
    idx2 <- aucell_meta$treatment == comp[2]
    
    g1 <- aucell_meta[[tf]][idx1]
    g2 <- aucell_meta[[tf]][idx2]
    
    # Remove missing values
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    
    cat(
      tf, ":",
      comp[1], "n =", length(g1),
      "|",
      comp[2], "n =", length(g2), "\n"
    )
    
    test <- wilcox.test(
      g1,
      g2,
      exact = FALSE
    )
    
    scenic.p.df <- rbind(
      scenic.p.df,
      data.frame(
        TF = tf,
        Group1 = comp[1],
        Group2 = comp[2],
        Pvalue = test$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
}




# Build final SCENIC + chromVAR integrated table 
sapply(aucell_meta[, tf.to.plot], class)
table(aucell_meta$treatment)


## SCENIC median activity
scenic.medians <- matrix(
  NA,
  nrow = length(tf.to.plot),
  ncol = length(treatments),
  dimnames = list(tf.to.plot, treatments)
)

for (tf in tf.to.plot) {
  
  for (tr in treatments) {
    
    x <- aucell_meta[aucell_meta$treatment == tr, tf]
    
    scenic.medians[tf, tr] <- median(
      x,
      na.rm = TRUE
    )
  }
}

scenic.medians


## exclude FOS, FOXO1, STAT1 as zero scenic score
rows.remove <- which(rownames(scenic.medians) %in% c("FOS(+)", "FOXO1(+)", "STAT1(+)"))

scenic.medians.use <- scenic.medians[ - rows.remove , ] 

tf.to.plot <- rownames(scenic.medians.use)

scenic.fc <- matrix(
  NA,
  nrow = length(tf.to.plot),
  ncol = length(treatments) - 1,
  dimnames = list(
    tf.to.plot,
    c("ENZ48_vs_DMSO", "RESA_vs_DMSO", "RESB_vs_DMSO")
  )
)

for (tf in tf.to.plot) {
  
  dmso <- scenic.medians.use[tf, "DMSO"]
  
  scenic.fc[tf, "ENZ48_vs_DMSO"] <-
    scenic.medians[tf, "ENZ48"] / dmso
  
  scenic.fc[tf, "RESA_vs_DMSO"] <-
    scenic.medians[tf, "RESA"] / dmso
  
  scenic.fc[tf, "RESB_vs_DMSO"] <-
    scenic.medians[tf, "RESB"] / dmso
}

scenic.fc

scenic.log2fc <- log2(scenic.fc)

scenic.log2fc


## Calculate chromVAR fold changes
tf.to.plot <- c(
  "AR(+)",
  "NKX3-1(+)",
  "MYC(+)",
  "E2F1(+)",
  "JUN(+)",
  "FOS(+)",
  "FOXO1(+)",
  "STAT1(+)"
)

treatments <- c("DMSO", "ENZ48", "RESA", "RESB")


chromvar.delta <- matrix(
  NA,
  nrow = nrow(chromvar.medians),
  ncol = 3,
  dimnames = list(
    rownames(chromvar.medians),
    c(
      "ENZ48_vs_DMSO",
      "RESA_vs_DMSO",
      "RESB_vs_DMSO"
    )
  )
)

for (tf in rownames(chromvar.medians)) {
  
  dmso <- chromvar.medians[tf, "DMSO"]
  
  chromvar.delta[tf, "ENZ48_vs_DMSO"] <-
    chromvar.medians[tf, "ENZ48"] - dmso
  
  chromvar.delta[tf, "RESA_vs_DMSO"] <-
    chromvar.medians[tf, "RESA"] - dmso
  
  chromvar.delta[tf, "RESB_vs_DMSO"] <-
    chromvar.medians[tf, "RESB"] - dmso
}

round(chromvar.delta, 3)





 # Calculate SCENIC deltas
scenic.delta <- matrix(
  NA,
  nrow = nrow(scenic.medians),
  ncol = 3,
  dimnames = list(
    rownames(scenic.medians),
    c(
      "ENZ48_vs_DMSO",
      "RESA_vs_DMSO",
      "RESB_vs_DMSO"
    )
  )
)

for (tf in rownames(scenic.medians)) {
  
  scenic.delta[tf, "ENZ48_vs_DMSO"] <-
    scenic.medians[tf, "ENZ48"] -
    scenic.medians[tf, "DMSO"]
  
  scenic.delta[tf, "RESA_vs_DMSO"] <-
    scenic.medians[tf, "RESA"] -
    scenic.medians[tf, "DMSO"]
  
  scenic.delta[tf, "RESB_vs_DMSO"] <-
    scenic.medians[tf, "RESB"] -
    scenic.medians[tf, "DMSO"]
}

round(scenic.delta, 3)



## make integrated table of scenic and chromVar
scenic.delta.df <- as.data.frame(scenic.delta)

scenic.delta.df$TF <- rownames(scenic.delta.df)

rownames(scenic.delta.df) <- NULL

scenic.delta.df <- scenic.delta.df[, c(
  "TF",
  "ENZ48_vs_DMSO",
  "RESA_vs_DMSO",
  "RESB_vs_DMSO"
)]



## Map scenic names to chromVar motif names
tf.mapping <- data.frame(
  SCENIC = c(
    "AR(+)",
    "NKX3-1(+)",
    "MYC(+)",
    "E2F1(+)",
    "JUN(+)",
    "FOS(+)",
    "FOXO1(+)",
    "STAT1(+)"
  ),
  
  chromVAR = c(
    NA,
    NA,
    "MYC",
    "E2F1",
    "JUND",
    NA,
    NA,
    NA
  ),
  
  stringsAsFactors = FALSE
)

tf.mapping



integrated <- tf.mapping

integrated$SCENIC_ENZ48_vs_DMSO <- NA
integrated$SCENIC_RESA_vs_DMSO <- NA
integrated$SCENIC_RESB_vs_DMSO <- NA

integrated$chromVAR_ENZ48_vs_DMSO <- NA
integrated$chromVAR_RESA_vs_DMSO <- NA
integrated$chromVAR_RESB_vs_DMSO <- NA

## fill in scenic
for (i in 1:nrow(integrated)) {
  
  tf <- integrated$SCENIC[i]
  
  integrated$SCENIC_ENZ48_vs_DMSO[i] <-
    scenic.delta[tf, "ENZ48_vs_DMSO"]
  
  integrated$SCENIC_RESA_vs_DMSO[i] <-
    scenic.delta[tf, "RESA_vs_DMSO"]
  
  integrated$SCENIC_RESB_vs_DMSO[i] <-
    scenic.delta[tf, "RESB_vs_DMSO"]
}

## fill in chromvar
for (i in 1:nrow(integrated)) {
  
  tf <- integrated$chromVAR[i]
  
  if (!is.na(tf) && tf %in% rownames(chromvar.delta)) {
    
    integrated$chromVAR_ENZ48_vs_DMSO[i] <-
      chromvar.delta[tf, "ENZ48_vs_DMSO"]
    
    integrated$chromVAR_RESA_vs_DMSO[i] <-
      chromvar.delta[tf, "RESA_vs_DMSO"]
    
    integrated$chromVAR_RESB_vs_DMSO[i] <-
      chromvar.delta[tf, "RESB_vs_DMSO"]
  }
}

integrated




## Add significance stars for scenic findings
scenic.heatmap <- scenic.delta

scenic.stars <- matrix(
  "",
  nrow = nrow(scenic.heatmap),
  ncol = ncol(scenic.heatmap),
  dimnames = dimnames(scenic.heatmap)
)

for (tf in rownames(scenic.heatmap)) {
  
  for (j in 1:ncol(scenic.heatmap)) {
    
    comparison <- colnames(scenic.heatmap)[j]
    
    if (comparison == "ENZ48_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "ENZ48"
    }
    
    if (comparison == "RESA_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "RESA"
    }
    
    if (comparison == "RESB_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "RESB"
    }
    
    idx <- scenic.p.df$TF == tf &
      scenic.p.df$Group1 == g1 &
      scenic.p.df$Group2 == g2
    
    if (any(idx)) {
      
      fdr <- scenic.p.df$FDR[idx][1]
      
      if (fdr < 0.001) {
        scenic.stars[tf, j] <- "***"
      } else if (fdr < 0.01) {
        scenic.stars[tf, j] <- "**"
      } else if (fdr < 0.05) {
        scenic.stars[tf, j] <- "*"
      }
    }
  }
}

## scenic heatmap
pdf("SCENIC_heatmap_of_TF_activity_across_between_treatments.pdf")
pheatmap(
  scenic.heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = scenic.stars,
  number_color = "black",
  fontsize_row = 11,
  fontsize_col = 11,
  main = "SCENIC TF activity change"
)
dev.off()


## Add significance stars for chromvar findings
chromvar.heatmap <- chromvar.delta


chromvar.stars <- matrix(
  "",
  nrow = nrow(chromvar.heatmap),
  ncol = ncol(chromvar.heatmap),
  dimnames = dimnames(chromvar.heatmap)
)

for (tf in rownames(chromvar.heatmap)) {
  
  for (j in 1:ncol(chromvar.heatmap)) {
    
    comparison <- colnames(chromvar.heatmap)[j]
    
    if (comparison == "ENZ48_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "ENZ48"
    }
    
    if (comparison == "RESA_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "RESA"
    }
    
    if (comparison == "RESB_vs_DMSO") {
      g1 <- "DMSO"
      g2 <- "RESB"
    }
    
    idx <- chromvar.p.df$TF == tf &
      chromvar.p.df$Group1 == g1 &
      chromvar.p.df$Group2 == g2
    
    if (any(idx)) {
      
      fdr <- chromvar.p.df$FDR[idx][1]
      
      if (fdr < 0.001) {
        chromvar.stars[tf, j] <- "***"
      } else if (fdr < 0.01) {
        chromvar.stars[tf, j] <- "**"
      } else if (fdr < 0.05) {
        chromvar.stars[tf, j] <- "*"
      }
    }
  }
}

pdf("Chromvar_heatmap_of_TF_activity_across_between_treatments.pdf")
pheatmap(
  chromvar.heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = chromvar.stars,
  number_color = "black",
  fontsize_row = 11,
  fontsize_col = 11,
  main = "chromVAR motif accessibility change"
)
dev.off()



## AR activity progressively increases toward the resistant state, while prostate-lineage/NKX3-1 and MYC activity decline, and AP-1-associated JUN/FOS activity becomes elevated, particularly in RESB.


## files to export
write.csv(
  scenic.medians,
  "SCENIC_TF_activity_medians.csv",
  row.names = TRUE
)

write.csv(
  scenic.delta,
  "SCENIC_TF_delta_vs_DMSO.csv",
  row.names = TRUE
)

write.csv(
  scenic.p.df,
  "SCENIC_pairwise_Wilcoxon_results.csv",
  row.names = FALSE
)

write.csv(
  chromvar.medians,
  "chromVAR_TF_motif_accessibility_medians.csv",
  row.names = TRUE
)

write.csv(
  chromvar.delta,
  "chromVAR_TF_delta_vs_DMSO.csv",
  row.names = TRUE
)

write.csv(
  chromvar.p.df,
  "chromVAR_pairwise_Wilcoxon_results.csv",
  row.names = FALSE
)


write.csv(
  integrated,
  "SCENIC_chromVAR_integrated_results.csv",
  row.names = FALSE
)


write.csv(
  tf.mapping,
  "SCENIC_chromVAR_TF_mapping.csv",
  row.names = FALSE
)


scenic_cell_level <- aucell_meta[ ,c("Cell", "treatment", prostate_cancer_tfs)]

write.csv(
  scenic_cell_level,
  "SCENIC_cell_level_activity.csv",
  row.names = FALSE
)

chromvar_cell_level <- chromvar_df[ ,c("Cell", "Sample_Id", chromvar.tfs)]


write.csv(
  chromvar_cell_level,
  "chromVAR_cell_level_deviation.csv",
  row.names = FALSE
)



