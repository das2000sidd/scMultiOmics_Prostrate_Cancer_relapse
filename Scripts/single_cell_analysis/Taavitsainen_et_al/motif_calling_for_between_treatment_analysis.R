run_dmso_dar_motif_analysis <- function(
    treatment.dar.anno,
    treatment.from.combined,
    output_prefix = "motifs"
) {
  
  library(Seurat)
  library(Signac)
  library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2018)
  
  # ============================================================
  # 1. Convert region into chr, start, end
  # ============================================================
  
  treatment.dar.anno$gene <- as.character(treatment.dar.anno$gene)
  
  coords <- do.call(
    rbind,
    strsplit(treatment.dar.anno$gene, "-", fixed = TRUE)
  )
  
  treatment.dar.anno$chr   <- coords[, 1]
  treatment.dar.anno$start <- as.integer(coords[, 2])
  treatment.dar.anno$end   <- as.integer(coords[, 3])
  
  # Check for failed coordinate parsing
  if (any(is.na(treatment.dar.anno$start)) ||
      any(is.na(treatment.dar.anno$end))) {
    
    stop(
      "Some genomic coordinates could not be parsed from the 'gene' column."
    )
  }
  
  # ============================================================
  # 2. Select significant DARs
  # ============================================================
  
  
  
  treatment_dar_sig <- treatment.dar.anno[
    treatment.dar.anno$p_val_BF < 0.05,
  ]
  
  # ============================================================
  # 3. Load JASPAR2018 motifs
  # ============================================================
  
  pfm <- getMatrixSet(
    x = JASPAR2018,
    opts = list(
      collection = "CORE",
      tax_group = "vertebrates"
    )
  )
  
  # ============================================================
  # 4. Select TFs of interest
  # ============================================================
  
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
  
  pfm_names <- sapply(
    pfm,
    function(x) x@name
  )
  
  pfm.use <- pfm[
    pfm_names %in% tf_keep
  ]
  
  message(
    "Number of TF motifs selected: ",
    length(pfm.use)
  )
  
  # ============================================================
  # 5. Prepare ATAC object
  # ============================================================
  
  DefaultAssay(treatment.from.combined) <- "peaks"
  
  treatment.from.combined <- RegionStats(
    object = treatment.from.combined,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  treatment.from.combined <- AddMotifs(
    object = treatment.from.combined,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm.use,
    verbose = TRUE
  )
  
  # ============================================================
  # 6. Identify clusters
  # ============================================================
  
  clusters <- sort(
    unique(treatment_dar_sig$cluster)
  )
  
  message(
    "Clusters to analyse: ",
    paste(clusters, collapse = ", ")
  )
  
  # ============================================================
  # 7. Motif enrichment for each cluster
  # ============================================================
  
  motif_results <- list()
  
  for (cl in clusters) {
    
    message(
      "\n====================================\n",
      "Processing cluster ", cl,
      "\n===================================="
    )
    
    peaks_cluster <- treatment_dar_sig$gene[
      treatment_dar_sig$cluster == cl
    ]
    
    peaks_cluster <- unique(peaks_cluster)
    
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
      "_",
      cl,
      ".csv"
    )
    
    write.csv(
      motif_result,
      file = output_file,
      row.names = FALSE,
      quote = FALSE
    )
    
    message(
      "Written: ",
      output_file
    )
  }
  
  # ============================================================
  # 8. Return results
  # ============================================================
  
  return(
    list(
      significant_DARs = treatment_dar_sig,
      motif_results = motif_results,
      seurat_object = treatment.from.combined
    )
  )
}


integrated <- readRDS("combined_atac_harmony_integrated.rds")

integrated_subset <- subset(
  integrated,
  subset = Sample_Id %in% c("RESB", "ENZ48")
)

resb_enz48_c1 <- read.table(file="Differentially_accessible_regions_RESB_vs_ENZ48_cluster_1.txt",header = T,sep="\t",stringsAsFactors = F)
resb_enz48_c2 <- read.table(file="Differentially_accessible_regions_RESB_vs_ENZ48_cluster_2.txt",header = T,sep="\t",stringsAsFactors = F)
resb_enz48_c3 <- read.table(file="Differentially_accessible_regions_RESB_vs_ENZ48_cluster_3.txt",header = T,sep="\t",stringsAsFactors = F)
resb_enz48_c4 <- read.table(file="Differentially_accessible_regions_RESB_vs_ENZ48_cluster_4.txt",header = T,sep="\t",stringsAsFactors = F)
resb_enz48_c5 <- read.table(file="Differentially_accessible_regions_RESB_vs_ENZ48_cluster_5.txt",header = T,sep="\t",stringsAsFactors = F)

resb_enz48_all_c <- rbind(resb_enz48_c1, resb_enz48_c2, resb_enz48_c3, resb_enz48_c4, resb_enz48_c5)

resb_enz48_motifs <- run_dmso_dar_motif_analysis(
  resb_enz48_all_c,
  integrated_subset,
  "RESB_vs_ENZ48_motifs_cluster"
)



resa_enz48_c1 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_1.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c2 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_2.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c3 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_3.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c4 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_4.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c5 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_5.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c6 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_6.txt",header = T,sep="\t",stringsAsFactors = F)
resa_enz48_c8 <- read.table(file="Differentially_accessible_regions_RESA_vs_ENZ48_cluster_8.txt",header = T,sep="\t",stringsAsFactors = F)

resa_enz48_all_c <- rbind(resa_enz48_c1, resa_enz48_c2, resa_enz48_c3, resa_enz48_c4, resa_enz48_c5,
                          resa_enz48_c6,resa_enz48_c8)

integrated_subset <- subset(
  integrated,
  subset = Sample_Id %in% c("RESA", "ENZ48")
)


resa_enz48_motifs <- run_dmso_dar_motif_analysis(
  resa_enz48_all_c,
  integrated_subset,
  "RESA_vs_ENZ48_motifs_cluster"
) ## Nothing???



resa_dmso_c1 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_1.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c2 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_2.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c3 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_3.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c4 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_4.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c5 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_5.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c6 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_6.txt",header = T,sep="\t",stringsAsFactors = F)
resa_dmso_c8 <- read.table(file="Differentially_accessible_regions_RESA_vs_DMSO_cluster_8.txt",header = T,sep="\t",stringsAsFactors = F)


resa_dmso_c1$p_val_BF <- p.adjust(resa_dmso_c1$p_val, method = "bonferroni")
resa_dmso_c2$p_val_BF <- p.adjust(resa_dmso_c2$p_val, method = "bonferroni")
resa_dmso_c3$p_val_BF <- p.adjust(resa_dmso_c3$p_val, method = "bonferroni")
resa_dmso_c4$p_val_BF <- p.adjust(resa_dmso_c4$p_val, method = "bonferroni")
resa_dmso_c5$p_val_BF <- p.adjust(resa_dmso_c5$p_val, method = "bonferroni")
resa_dmso_c6$p_val_BF <- p.adjust(resa_dmso_c6$p_val, method = "bonferroni")
resa_dmso_c8$p_val_BF <- p.adjust(resa_dmso_c8$p_val, method = "bonferroni")


resa_dmso_all_c <- rbind(resa_dmso_c1, resa_dmso_c2, resa_dmso_c3, resa_dmso_c4, resa_dmso_c5,
                          resa_dmso_c6,resa_dmso_c8)

resa_dmso_sig <- resa_dmso_all_c[resa_dmso_all_c$p_val_BF < 0.05, ]

integrated_subset <- subset(
  integrated,
  subset = Sample_Id %in% c("RESA", "DMSO")
)

resa_dmso_motifs <- run_dmso_dar_motif_analysis(
  resa_dmso_sig,
  integrated_subset,
  "RESA_vs_DMSO_motifs_cluster"
) ## Nothing???



resb_dmso_c1 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_1.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c2 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_2.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c3 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_3.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c4 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_4.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c5 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_5.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c6 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_6.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c7 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_7.txt",header = T,sep="\t",stringsAsFactors = F)
resb_dmso_c8 <- read.table(file="Differentially_accessible_regions_RESB_vs_DMSO_cluster_8.txt",header = T,sep="\t",stringsAsFactors = F)


resb_dmso_c1$p_val_BF <- p.adjust(resb_dmso_c1$p_val, method = "bonferroni")
resb_dmso_c2$p_val_BF <- p.adjust(resb_dmso_c2$p_val, method = "bonferroni")
resb_dmso_c3$p_val_BF <- p.adjust(resb_dmso_c3$p_val, method = "bonferroni")
resb_dmso_c4$p_val_BF <- p.adjust(resb_dmso_c4$p_val, method = "bonferroni")
resb_dmso_c5$p_val_BF <- p.adjust(resb_dmso_c5$p_val, method = "bonferroni")
resb_dmso_c6$p_val_BF <- p.adjust(resb_dmso_c6$p_val, method = "bonferroni")
resb_dmso_c7$p_val_BF <- p.adjust(resb_dmso_c7$p_val, method = "bonferroni")
resb_dmso_c8$p_val_BF <- p.adjust(resb_dmso_c8$p_val, method = "bonferroni")


resb_dmso_all_c <- rbind(resb_dmso_c1, resb_dmso_c2, resb_dmso_c3, resb_dmso_c4, resb_dmso_c5,
                         resb_dmso_c6, resb_dmso_c7,resb_dmso_c8)

resb_dmso_sig <- resb_dmso_all_c[resb_dmso_all_c$p_val_BF < 0.05, ]

integrated_subset <- subset(
  integrated,
  subset = Sample_Id %in% c("RESB", "DMSO")
)

resb_dmso_motifs <- run_dmso_dar_motif_analysis(
  resb_dmso_sig,
  integrated_subset,
  "RESB_vs_DMSO_motifs_cluster"
)



enz48_dmso_c1 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_1.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c2 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_2.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c3 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_3.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c4 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_4.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c5 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_5.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c6 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_6.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c7 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_7.txt",header = T,sep="\t",stringsAsFactors = F)
enz48_dmso_c8 <- read.table(file="Differentially_accessible_regions_ENZ48_vs_DMSO_cluster_8.txt",header = T,sep="\t",stringsAsFactors = F)


enz48_dmso_c1$p_val_BF <- p.adjust(enz48_dmso_c1$p_val, method = "bonferroni")
enz48_dmso_c2$p_val_BF <- p.adjust(enz48_dmso_c2$p_val, method = "bonferroni")
enz48_dmso_c3$p_val_BF <- p.adjust(enz48_dmso_c3$p_val, method = "bonferroni")
enz48_dmso_c4$p_val_BF <- p.adjust(enz48_dmso_c4$p_val, method = "bonferroni")
enz48_dmso_c5$p_val_BF <- p.adjust(enz48_dmso_c5$p_val, method = "bonferroni")
enz48_dmso_c6$p_val_BF <- p.adjust(enz48_dmso_c6$p_val, method = "bonferroni")
enz48_dmso_c7$p_val_BF <- p.adjust(enz48_dmso_c7$p_val, method = "bonferroni")
enz48_dmso_c8$p_val_BF <- p.adjust(enz48_dmso_c8$p_val, method = "bonferroni")


enz48_dmso_all_c <- rbind(enz48_dmso_c1, enz48_dmso_c2, enz48_dmso_c3, enz48_dmso_c4, enz48_dmso_c5,
                          enz48_dmso_c6, enz48_dmso_c7,enz48_dmso_c8)

enz48_dmso_sig <- enz48_dmso_all_c[enz48_dmso_all_c$p_val_BF < 0.05, ]

integrated_subset <- subset(
  integrated,
  subset = Sample_Id %in% c("ENZ48", "DMSO")
)

enz48_dmso_motifs <- run_dmso_dar_motif_analysis(
  enz48_dmso_sig,
  integrated_subset,
  "ENZ48_vs_DMSO_motifs_cluster"
)


