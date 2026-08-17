setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


motif_files <- c(
  "DMSO_motifs_cluster_1.csv",
  "DMSO_motifs_cluster_2.csv",
  "DMSO_motifs_cluster_4.csv",
  "DMSO_motifs_cluster_5.csv",
  "DMSO_motifs_cluster_6.csv",
  "DMSO_motifs_cluster_7.csv"
)

motif_dmso <- lapply(
  motif_files,
  read.csv,
  header = TRUE,
  stringsAsFactors = FALSE
)

names(motif_dmso) <- c(
  "cluster_1",
  "cluster_2",
  "cluster_4",
  "cluster_5",
  "cluster_6",
  "cluster_7"
)


c1 <- motif_dmso$cluster_1
c2 <- motif_dmso$cluster_2
c4 <- motif_dmso$cluster_4
c5 <- motif_dmso$cluster_5
c6 <- motif_dmso$cluster_6
c7 <- motif_dmso$cluster_7

rownames(c1) <- c1$motif.name
rownames(c2) <- c2$motif.name
rownames(c4) <- c4$motif.name
rownames(c5) <- c5$motif.name
rownames(c6) <- c6$motif.name
rownames(c7) <- c7$motif.name

c1$p.adjust.FDR_log <- -log10(c1$p.adjust.FDR)
c2$p.adjust.FDR_log <- -log10(c2$p.adjust.FDR)
c4$p.adjust.FDR_log <- -log10(c4$p.adjust.FDR)
c5$p.adjust.FDR_log <- -log10(c5$p.adjust.FDR)
c6$p.adjust.FDR_log <- -log10(c6$p.adjust.FDR)
c7$p.adjust.FDR_log <- -log10(c7$p.adjust.FDR)



dmso_motif_fold_enrichment <- cbind(c1[rownames(c1), "fold.enrichment"], c2[rownames(c1), "fold.enrichment"],
                                    c4[rownames(c1), "fold.enrichment"], c5[rownames(c1), "fold.enrichment"],
                                    c6[rownames(c1), "fold.enrichment"], c7[rownames(c1), "fold.enrichment"]);

rownames(dmso_motif_fold_enrichment) <- rownames(c1);


colnames(dmso_motif_fold_enrichment) <- c("C1", "C2", "C4", "C5", "C6", "C7")

dmso_motif_fold_enrichment_reordered <- dmso_motif_fold_enrichment[c("Ar", "CREB1", "CTCF", "E2F1",
                                                                     "ELF1", "ELK4", "ETV1", "FOXA1",
                                                                     "FOXP1", "GATA2", "GRHL2", "HOXB13",
                                                                     "JUND", "MYC", "POU2F1", "RELA",
                                                                     "SOX9", "TFAP4"), ]


library(pheatmap)
library(grid)
library(gridExtra)

# Example values for the barplot
dar_dmso <- read.csv(file="Differentially_accessible_peaks_all_clusters_DMSO.csv",header = T,stringsAsFactors = F)

dar_dmso_sig <- dar_dmso[ dar_dmso$p_val_BF < 0.05, ]

dar_dmso_sig$Significance <- "Yes"

dmso.sig.dar.count <- as.vector(table(dar_dmso_sig$cluster, dar_dmso_sig$Significance))

dmso.sig.dar.count.df <- cbind(dmso.sig.dar.count, c(paste("C", c(1,2,4,5,6,7),sep="")))

dmso.sig.dar.count.df[,1] <- as.numeric(dmso.sig.dar.count.df[,1])

dmso.sig.dar.count.df <- data.frame(
  DARs = dmso.sig.dar.count,
  Cluster = c(paste("C", c(1,2,4,5,6,7),sep=""))
)

pdf(
  "DMSO_bar_count_by_cluster.pdf",
  width = 6,
  height = 4
)

ggplot(
  dmso.sig.dar.count.df,
  aes(x = Cluster, y = DARs)
) +
  geom_col() +
  theme_classic() +
  labs(
    x = "Number of DARs",
    y = "Cluster"
  )

dev.off()


heatmap_max_val <- max(dmso_motif_fold_enrichment_reordered)

dmso.heatmap <- pheatmap(
  dmso_motif_fold_enrichment_reordered,
  scale = "none",
  color = colorRampPalette(c("white", "red"))(200),
  breaks = seq(0, heatmap_max_val, length.out = 201),
  main = "",
  cluster_cols = FALSE,
  cluster_rows = FALSE
)

pdf(
  "DMSO_motif_heatmap.pdf",
  width = 6,
  height = 4
  )

dmso.heatmap

dev.off()






make_motif_heatmap <- function(
    treatment,
    motif_dir = ".",
    dar_file,
    output_dir = "."
) {
  
  # ============================================================
  # Libraries
  # ============================================================
  
  library(pheatmap)
  library(ggplot2)
  
  # ============================================================
  # 1. Read motif files
  # ============================================================
  
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
  
  # ============================================================
  # 2. Set motif names as row names
  # ============================================================
  
  motif_list <- lapply(
    motif_list,
    function(x) {
      rownames(x) <- x$motif.name
      x
    }
  )
  
  # ============================================================
  # 3. Create fold-enrichment matrix
  # ============================================================
  
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
  
  # ============================================================
  # 4. Reorder TFs
  # ============================================================
  
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
  
  # ============================================================
  # 5. Read DARs and calculate significant DAR counts
  # ============================================================
  
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
  
  # ============================================================
  # 6. Barplot
  # ============================================================
  
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
  
  # ============================================================
  # 7. Export barplot
  # ============================================================
  
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
  
  # ============================================================
  # 8. Heatmap
  # ============================================================
  
  #heatmap_max_val <- max(
  #  motif_fold_enrichment_reordered,
  #  na.rm = TRUE
  #)
  
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
  
  # ============================================================
  # 9. Export heatmap
  # ============================================================
  
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
  
  # ============================================================
  # 10. Return useful objects
  # ============================================================
  
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


#make_motif_heatmap("DMSO", ".", "")

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

