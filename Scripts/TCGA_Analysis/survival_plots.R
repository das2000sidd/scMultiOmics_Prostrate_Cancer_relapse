setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")


gsva_scores <- read.table(file="TCGA_PRAD_RNASeq_GSVA_scores_from_paper.txt",header = T,sep="\t",stringsAsFactors = F)
clinical_data <- read.table(file="TCGA_PRAD_Clinical_data.txt",header = T,sep="\t",stringsAsFactors = F)

gleason_7_higher <- clinical_data[clinical_data$gleason_score >= 7, ]

gleason_score <- gleason_7_higher[,c("submitter_id", "gleason_score")]


query.clinical <- read.csv(file="PCGA_PRAC_OS_PFS.csv",header = T, sep=",")


## rename END indiced clusters column
which(rownames(gsva_scores)=="ENZ-induced clusters")

rownames(gsva_scores)[which(rownames(gsva_scores)=="ENZ-induced clusters")] <- "ENZ_induced_clusters"


gsva_sample_ids <- colnames(gsva_scores)

gsva_sample_ids <- as.data.frame(gsva_sample_ids)

head(gsva_sample_ids)

gsva_sample_ids$patient_id <- sapply(
  strsplit(gsva_sample_ids$gsva_sample_ids, "\\."),
  function(x) paste(x[1:3], collapse = "-")
)


gsva_sample_ids_gleason <- left_join(gsva_sample_ids, gleason_score, by=c("patient_id" = "submitter_id"))


gsva_sample_ids_gleason$Gleason_category <- ifelse(gsva_sample_ids_gleason$gleason_score > 7, "8+", "7")


stopifnot(gsva_sample_ids_gleason$gsva_sample_ids == colnames(gsva_scores))


rownames(gsva_sample_ids_gleason) <- gsva_sample_ids_gleason$gsva_sample_ids


annotation_col <- gsva_sample_ids_gleason[colnames(gsva_scores), , drop = FALSE]




order_of_cols <- c("Initial clusters", "ENZ_induced_clusters", "Persistent clusters",
                   "NEPC upregulated", "NEPC downregulated",
                   "mTORC1 signaling", "PROSGenesis",
                   "Persist", "MYC targets V1", "MYC targets V2",
                   "AR activity", "ARFL", "ARV", "BRCAness",
                   paste0("Cluster",1:11))

gsva_scores_reordered <- gsva_scores[ order_of_cols, ]



ph <- pheatmap(
  gsva_scores_reordered,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col[, c(3:4)],
  breaks = seq(-0.5, 0.5, length.out = 101),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "GSVA Plot, TCGA PRAD",
  silent = TRUE
)

print(ph)

branch <- cutree(
  ph$tree_col,
  k = 2
)

table(branch)


branch1_samples <- names(branch[branch == 1])

branch2_samples <- names(branch[branch == 2])

annotation_col$Branch <- NA

annotation_col[names(branch), "Branch"] <- paste0(
  "Branch ",
  branch
)


branch_df <- data.frame(
  patient_id = names(branch),
  Branch = paste0("Branch_", branch)
)

branch_df$patient_id <- substr(
  branch_df$patient_id,
  1,
  12
)

branch_df$patient_id <- gsub("\\.", "-", branch_df$patient_id)

survival_df <- left_join(
  query.clinical,
  branch_df,
  by = c("patientId"="patient_id")
)

survival_df <- survival_df[ !is.na(survival_df$Branch), ]


survival_df$PFS_event <- ifelse(
  survival_df$PFS_STATUS == "1:PROGRESSION",
  1,
  0
)




## Survival analysis using ENZ induced clusters
gsva_enz_induced <- gsva_scores[ "ENZ_induced_clusters", ]

gsva_enz_induced <- as.data.frame(gsva_enz_induced)

gsva_enz_induced <- t(gsva_enz_induced)

gsva_enz_induced <- as.data.frame(gsva_enz_induced)

gsva_enz_induced$Fill_id <- rownames(gsva_enz_induced)

gsva_enz_induced$Sample_id <- gsub("\\.", "-", gsva_enz_induced$Fill_id)

gsva_enz_induced$Patient_ID <- sapply(
  strsplit(gsva_enz_induced$Sample_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

#gsva_enz_induced <- gsva_enz_induced[,c("gsva_enz_induced", "Patient_ID")]


survival_df_enz_induced <- left_join(survival_df, gsva_enz_induced, by = c("patientId"="Patient_ID"))


#survival_df_enz_induced <- survival_df_enz_induced[!is.na(survival_df_enz_induced$gsva_enz_induced),]


survival_df_enz_induced$enz_induced_threshold <- ifelse(survival_df_enz_induced$ENZ_induced_clusters > median(survival_df_enz_induced$ENZ_induced_clusters), "High", "Low")


survival_df_enz_induced$PFS_event <- ifelse(
  survival_df_enz_induced$enz_induced_threshold == "High",
  1,
  0
)

survival_df_enz_induced$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_enz_induced$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_enz_induced$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_enz_induced$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_enz_induced$PFS_event_status <- as.numeric(
  survival_df_enz_induced$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_enz_induced$PFS_MONTHS,
  event = survival_df_enz_induced$PFS_event_status
)

fit <- survfit(
  surv_object ~ enz_induced_threshold,
  data = survival_df_enz_induced
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_ENZ_induced_clusters.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_enz_induced,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off() 





