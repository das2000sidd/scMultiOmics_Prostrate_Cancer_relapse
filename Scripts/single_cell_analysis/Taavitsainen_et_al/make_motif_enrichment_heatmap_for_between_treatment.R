setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


library(pheatmap)
library(grid)
library(gridExtra)

# Example values for the barplot
c1 <- read.csv(file="RESB_vs_ENZ48_motifs_cluster_1.csv",header = T,stringsAsFactors = F)
c2 <- read.csv(file="RESB_vs_ENZ48_motifs_cluster_2.csv",header = T,stringsAsFactors = F)
c3 <- read.csv(file="RESB_vs_ENZ48_motifs_cluster_3.csv",header = T,stringsAsFactors = F)
c4 <- read.csv(file="RESB_vs_ENZ48_motifs_cluster_4.csv",header = T,stringsAsFactors = F)
c5 <- read.csv(file="RESB_vs_ENZ48_motifs_cluster_5.csv",header = T,stringsAsFactors = F)

c1$cluster <- 1
c2$cluster <- 2
c3$cluster <- 3
c4$cluster <- 4
c5$cluster <- 5

rownames(c1) <- c1$motif.name
rownames(c2) <- c2$motif.name
rownames(c3) <- c3$motif.name
rownames(c4) <- c4$motif.name
rownames(c5) <- c5$motif.name


all_clusters_combined <- rbind(c1,c2,c3,c4,c5)


sig_motif <- all_clusters_combined[ all_clusters_combined$p.adjust.FDR < 0.05, ]


# Get unique motifs and clusters
motifs <- unique(all_clusters_combined$motif.name)
clusters <- sort(unique(all_clusters_combined$cluster))

# Create empty matrix
motif_fold_enrichment <- matrix(
  NA,
  nrow = length(motifs),
  ncol = 1
  )

rownames(motif_fold_enrichment) <- motifs

motif_fold_enrichment <- cbind(motif_fold_enrichment,c1[rownames(motif_fold_enrichment), "fold.enrichment"],
                               c2[rownames(motif_fold_enrichment), "fold.enrichment"],
                               c3[rownames(motif_fold_enrichment), "fold.enrichment"],
                               c4[rownames(motif_fold_enrichment), "fold.enrichment"],
                               c5[rownames(motif_fold_enrichment), "fold.enrichment"])

colnames(motif_fold_enrichment)[2:6] <- paste("C", c(1,2,3,4,5),sep="") 



motif_fold_enrichment [is.na(motif_fold_enrichment)] <- 0

## Add zeros for two clusters with no enrichment
motif_fold_enrichment <- cbind(
  motif_fold_enrichment,
  C6 = 0,
  C7 = 0
)



heatmap_max_val <- max(motif_fold_enrichment[, -c(1)])

rownames_order <- c("Ar", "CREB1", "CTCF", "E2F1", "ELF1", "ELK4", "ETV1", 
                    "FOXA1", "FOXP1", "GATA2", "GRHL2", "HOXB13", "JUND", 
                    "MYC", "POU2F1", "RELA", "SOX9", "TFAP4")


motif_fold_enrichment.o <- motif_fold_enrichment[rownames_order, ]



resb_vs_enz48.heatmap <- pheatmap(
  motif_fold_enrichment.o[, -c(1)],
  scale = "none",
  color = colorRampPalette(c("white", "red"))(200),
  breaks = seq(0, heatmap_max_val, length.out = 201),
  main = "RES-B vs ENZ48",
  cluster_cols = FALSE,
  cluster_rows = FALSE
)

pdf(
  "RESB_vs_ENZ48_motif_heatmap.pdf",
  width = 6,
  height = 4
  )

resb_vs_enz48.heatmap

dev.off()





make_motif_heatmap <- function(
    condition1,
    condition2,
    n_clusters = 5,
    output_file = NULL,
    output_width = 6,
    output_height = 4
) {
  
  # ============================================================
  # 1. Read motif enrichment files
  # ============================================================
  
  motif_files <- paste0(
    condition1,
    "_vs_",
    condition2,
    "_motifs_cluster_",
    1:n_clusters,
    ".csv"
  )
  
  cluster_data <- lapply(
    motif_files,
    read.csv,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Add cluster number
  for (i in seq_along(cluster_data)) {
    cluster_data[[i]]$cluster <- i
  }
  
  # ============================================================
  # 2. Combine all clusters
  # ============================================================
  
  all_clusters_combined <- do.call(
    rbind,
    cluster_data
  )
  
  # ============================================================
  # 3. Select significant motifs
  # ============================================================
  
  sig_motif <- all_clusters_combined[
    all_clusters_combined$p.adjust.FDR < 0.05,
  ]
  
  # ============================================================
  # 4. Get motifs
  # ============================================================
  
  motifs <- unique(
    all_clusters_combined$motif.name
  )
  
  # ============================================================
  # 5. Create fold-enrichment matrix
  # ============================================================
  
  motif_fold_enrichment <- matrix(
    0,
    nrow = length(motifs),
    ncol = n_clusters
  )
  
  rownames(motif_fold_enrichment) <- motifs
  
  colnames(motif_fold_enrichment) <- paste0(
    "C",
    1:n_clusters
  )
  
  # ============================================================
  # 6. Fill matrix with fold enrichment
  # ============================================================
  
  for (i in seq_along(cluster_data)) {
    
    dat <- cluster_data[[i]]
    
    # Make sure motif names are unique
    dat <- dat[
      !duplicated(dat$motif.name),
    ]
    
    rownames(dat) <- dat$motif.name
    
    common_motifs <- intersect(
      motifs,
      rownames(dat)
    )
    
    motif_fold_enrichment[
      common_motifs,
      i
    ] <- dat[
      common_motifs,
      "fold.enrichment"
    ]
  }
  
  # ============================================================
  # 7. Order motifs
  # ============================================================
  
  rownames_order <- c(
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
    "POU2F1",
    "RELA",
    "SOX9",
    "TFAP4"
  )
  
  rownames_order <- rownames_order[
    rownames_order %in% rownames(motif_fold_enrichment)
  ]
  
  motif_fold_enrichment.o <- motif_fold_enrichment[
    rownames_order,
    ,
    drop = FALSE
  ]
  
  # ============================================================
  # 8. Heatmap colour scale
  # ============================================================
  
  heatmap_max_val <- max(
    motif_fold_enrichment.o,
    na.rm = TRUE
  )
  
  # ============================================================
  # 9. Make heatmap
  # ============================================================
  
  heatmap <- pheatmap(
    motif_fold_enrichment.o,
    scale = "none",
    color = colorRampPalette(
      c("white", "red")
    )(200),
    breaks = seq(
      0,
      heatmap_max_val,
      length.out = 201
    ),
    main = paste(
      condition1,
      "vs",
      condition2
    ),
    cluster_cols = FALSE,
    cluster_rows = FALSE
  )
  
  # ============================================================
  # 10. Export
  # ============================================================
  
  if (!is.null(output_file)) {
    
    pdf(
      output_file,
      width = output_width,
      height = output_height
    )
    
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    
    dev.off()
  }
  
  # ============================================================
  # 11. Return useful objects
  # ============================================================
  
  return(
    list(
      combined_results = all_clusters_combined,
      significant_motifs = sig_motif,
      fold_enrichment = motif_fold_enrichment.o,
      heatmap = heatmap
    )
  )
}

resb_enz48_motif <- make_motif_heatmap(
  condition1 = "RESB",
  condition2 = "ENZ48",
  n_clusters = 5,
  output_file = "RESB_vs_ENZ48_motif_heatmap.pdf"
)

resa_enz48_motif <- make_motif_heatmap(
  condition1 = "RESA",
  condition2 = "ENZ48",
  n_clusters = 5,
  output_file = "RESA_vs_ENZ48_motif_heatmap.pdf"
)

resb_dmso_motif <- make_motif_heatmap(
  condition1 = "RESB",
  condition2 = "DMSO",
  n_clusters = 5,
  output_file = "RESB_vs_DMSO_motif_heatmap.pdf"
)

enz48_dmso_motif <- make_motif_heatmap(
  condition1 = "ENZ48",
  condition2 = "DMSO",
  n_clusters = 2,
  output_file = "ENZ48_vs_DMSO_motif_heatmap.pdf"
)

