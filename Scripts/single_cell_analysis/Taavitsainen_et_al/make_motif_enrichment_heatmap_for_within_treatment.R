setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")

library(pheatmap)
library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)

make_motif_heatmap <- function(
    treatment,
    motif_dir = ".",
    dar_file,
    output_dir = "."
) {
  
  clusters_dmso_enz <- c(1, 2, 4, 5, 6, 7)
  
  clusters_resb <- c(1,2,5,6,7,8)
  
  clusters_resa <- c(1,2,4,5,6,7,8)
  
  if(treatment == "DMSO" | treatment == "ENZ"){
    clusters <- clusters_dmso_enz
  } else if(treatment == "RESB") {
    clusters <- clusters_resb
  } else if(treatment == "RESA"){
    clusters <- clusters_resa
  }
  
  motif_files <- file.path(
    motif_dir,
    paste0(treatment, "_motifs_cluster_", clusters, ".csv")
  )
  
  motif_list <- lapply(
    motif_files,
    read.csv,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  names(motif_list) <- paste0("cluster_", clusters)
  
  motif_list <- lapply(
    motif_list,
    function(x) {
      rownames(x) <- x$motif.name
      x
    }
  )
  
  # Use motifs from the first cluster as reference
  motif_names <- rownames(motif_list[[1]])
  
  motif_fold_enrichment <- do.call(
    cbind,
    lapply(
      motif_list,
      function(x) {
        x[motif_names, "fold.enrichment"]
      }
    )
  )
  
  colnames(motif_fold_enrichment) <- paste0("C", clusters)
  
  rownames(motif_fold_enrichment) <- motif_names
  
  tf_order <- c(
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
  
  # Keep only TFs actually present
  tf_order <- tf_order[
    tf_order %in% rownames(motif_fold_enrichment)
  ]
  
  motif_fold_enrichment_reordered <-
    motif_fold_enrichment[tf_order, , drop = FALSE]
  
  dar <- read.csv(
    file = dar_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  dar_sig <- dar[
    dar$p_val_BF < 0.05,
  ]
  
  # Count significant DARs per cluster
  dar_counts <- sapply(
    clusters,
    function(cl) {
      sum(dar_sig$cluster == cl)
    }
  )
  
  dar_count_df <- data.frame(
    DARs = dar_counts,
    Cluster = paste0("C", clusters)
  )
  
  barplot <- ggplot(
    dar_count_df,
    aes(x = Cluster, y = DARs)
  ) +
    geom_col() +
    theme_classic() +
    labs(
      x = "Cluster",
      y = "Number of DARs"
    )
  
  pdf(
    file.path(
      output_dir,
      paste0(treatment, "_bar_count_by_cluster_repeat.pdf")
    ),
    width = 6,
    height = 4
  )
  
  print(barplot)
  
  dev.off()
  
  heatmap_max_val <- 10
  
  heatmap <- pheatmap(
    motif_fold_enrichment_reordered,
    scale = "none",
    color = colorRampPalette(
      c("white", "red")
    )(100),
    breaks = seq(
      0,
      heatmap_max_val,
      length.out = 101
    ),
    main = "",
    cluster_cols = FALSE,
    cluster_rows = FALSE
  )
  
  pdf(
    file.path(
      output_dir,
      paste0(treatment, "_motif_heatmap_repeat.pdf")
    ),
    width = 6,
    height = 4
  )
  
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  
  dev.off()
  
  return(
    list(
      motif_list = motif_list,
      motif_fold_enrichment =
        motif_fold_enrichment_reordered,
      dar_count_df = dar_count_df,
      dar_sig = dar_sig,
      barplot = barplot,
      heatmap = heatmap
    )
  )
}

dmso_results <- make_motif_heatmap(
  treatment = "DMSO",
  motif_dir = ".",
  dar_file = "Differentially_accessible_peaks_all_clusters_DMSO.csv",
  output_dir = "."
)


enz48_results <- make_motif_heatmap(
  treatment = "ENZ",
  motif_dir = ".",
  dar_file = "Differentially_accessible_regions_all_clusters_ENZ48.csv",
  output_dir = "."
)

resa_results <- make_motif_heatmap(
  treatment = "RESA",
  motif_dir = ".",
  dar_file = "Differentially_accessible_regions_all_clusters_RESA.csv",
  output_dir = "."
)


resb_results <- make_motif_heatmap(
  treatment = "RESB",
  motif_dir = ".",
  dar_file = "Differentially_accessible_regions_all_clusters_RESB.csv",
  output_dir = "."
)

