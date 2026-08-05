setwd("~/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC")

library(Seurat)
library(pheatmap)

rna_int <- readRDS("Batch_corrected_clustered_scRNAseq_object.rds")

genesets_plot <- read.csv(file="TCGA_PRAD_GSVA_genesets_to_plot.csv", header = T,stringsAsFactors = F)


DimPlot(
  rna_int,
  reduction = "umap",
  group.by = c("RNA_snn_res.0.8"),
  split.by = c("Sample_id"),
  label = FALSE
)
rna_int <- JoinLayers(rna_int, assay="RNA")

#all.markers.rna.int <- FindAllMarkers(rna_int, test.use = "MAST", latent.vars = "nCount_RNA", only.pos = TRUE,min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE)


#write.table(all.markers.rna.int,file = "FindAllMarkers_output_scRNA_seq_clusters.txt",col.names = T,row.names = F,sep="\t",quote = F)


# select the top markers per cluster
library(dplyr)


markers.by.cluster <- read.table(file="FindAllMarkers_output_scRNA_seq_clusters.txt", header = T,sep="\t", stringsAsFactors = F)

markers.by.cluster.sig <- markers.by.cluster[markers.by.cluster$p_val_adj < 0.05, ]


library(org.Hs.eg.db)

markers.by.cluster$Entrez <- mapIds(org.Hs.eg.db, markers.by.cluster$gene,keytype="SYMBOL", column="ENTREZID")
markers.by.cluster$Ensembl <- mapIds(org.Hs.eg.db, markers.by.cluster$Entrez,keytype="ENTREZID", column="ENSEMBL")


genesets_plot$Gene <- trimws(genesets_plot$Gene)

genesets_plot$Entrez <- mapIds(org.Hs.eg.db, genesets_plot$Gene,keytype="SYMBOL", column="ENTREZID")
genesets_plot$Ensembl <- mapIds(org.Hs.eg.db, genesets_plot$Entrez,keytype="ENTREZID", column="ENSEMBL")

genesets_plot$Entrez <- as.character(genesets_plot$Entrez)
genesets_plot$Ensembl <- as.character(genesets_plot$Ensembl)

genesets_plot$Ensembl_from_Symbol <- mapIds(org.Hs.eg.db, genesets_plot$Gene,keytype="SYMBOL", column="ENSEMBL")

genesets_supp <- genesets_plot[genesets_plot$Ensembl != "NULL", ]

genesets_list <- split(genesets_supp$Ensembl,
                       genesets_supp$Group)

markers.by.cluster$Ensembl <- as.character(markers.by.cluster$Ensembl)

gene_sets <- markers.by.cluster %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50)

gene_sets <- gene_sets[gene_sets$Ensembl != "NULL", ]

signatures <- split(gene_sets$Ensembl,
                    gene_sets$cluster)

names(signatures) <- paste("Cluster", names(signatures), sep="")



## download TCGA PRAD data
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

#query.clinical <- GDCquery(project = "TCGA-PRAD",
#                  data.category = "Clinical", data.type = "Clinical Supplement")

query.clinical <- GDCquery_clinic(
  project = "TCGA-PRAD",
  type = "clinical"
)


# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)
#GDCdownload(query.clinical)


# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
PRAD.Rnaseq.SE <- GDCprepare(query)

PRADMatrix <- assay(PRAD.Rnaseq.SE,"unstranded") 


attributes(PRADMatrix)


#write.table(PRADMatrix,file="TCGA_PRAD_Raw_Count_Matrix.txt",col.names = T,row.names = T,sep="\t",quote = F)
#write.table(query.clinical,file="TCGA_PRAD_Clinical_data.txt",col.names = T,row.names = F,sep="\t",quote = F)





gleason_7_higher <- query.clinical[query.clinical$gleason_score >= 7, ]


query.rna.gleason.based <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = gleason_7_higher$submitter_id
)

#GDCdownload(query.rna.gleason.based)



PRAD.Rnaseq.gleason.7.higher <- GDCprepare(query.rna.gleason.based)




library(DSS)
library(DESeq2)

PRADMatrix.gleason.based <- assay(PRAD.Rnaseq.gleason.7.higher,"unstranded") 

write.table(PRADMatrix.gleason.based, "TCGA_PRAD_RNASeq_count_minimum_Gleason_score_7.txt",col.names = T,row.names = T,sep="\t", quote = FALSE)

dds <- DESeqDataSetFromMatrix(
  countData = PRADMatrix.gleason.based,
  colData = data.frame(row.names = colnames(PRADMatrix.gleason.based)),
  design = ~1
)

dds <- estimateSizeFactors(dds)

keep <- rowSums(counts(dds) >= 20) >= (ncol(dds)*0.1)

table(keep)

dds <- dds[keep,]

norm_counts <- counts(dds, normalized = TRUE)

log_norm <- log2(norm_counts + 1)


log_norm[1:5, 1:5]

rownames(log_norm) <- sub("\\..*$", "", rownames(log_norm))


library(GSVA)

param <- gsvaParam(
  expr = log_norm,
  geneSets = c(signatures, genesets_list),
  kcdf = "Poisson"
)

gsva_scores <- gsva(param)

## rename END indiced clusters column
which(rownames(gsva_scores)=="ENZ-induced clusters")

rownames(gsva_scores)[which(rownames(gsva_scores)=="ENZ-induced clusters")] <- "ENZ_induced_clusters"


write.table(gsva_scores, "TCGA_PRAD_RNASeq_GSVA_scores_from_paper.txt",col.names = T,row.names = T,sep="\t", quote = FALSE)


gsva_sample_ids <- colnames(gsva_scores)

gsva_sample_ids <- as.data.frame(gsva_sample_ids)

head(gsva_sample_ids)

gsva_sample_ids$patient_id <- sapply(
  strsplit(gsva_sample_ids$gsva_sample_ids, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gleason_score <- gleason_7_higher[,c("submitter_id", "gleason_score")]


gsva_sample_ids_gleason <- left_join(gsva_sample_ids, gleason_score, by=c("patient_id" = "submitter_id"))


gsva_sample_ids_gleason$Gleason_category <- ifelse(gsva_sample_ids_gleason$gleason_score > 7, "8+", "7")


stopifnot(gsva_sample_ids_gleason$gsva_sample_ids == colnames(gsva_scores))


rownames(gsva_sample_ids_gleason) <- gsva_sample_ids_gleason$gsva_sample_ids


annotation_col <- gsva_sample_ids_gleason[colnames(gsva_scores), , drop = FALSE]


ph <- pheatmap(
  gsva_scores,
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

head(annotation_col)

rownames(gsva_scores)

order_of_cols <- c("Initial clusters", "ENZ_induced_clusters", "Persistent clusters",
                   "NEPC upregulated", "NEPC downregulated",
                   "mTORC1 signaling", "PROSGenesis",
                   "Persist", "MYC targets V1", "MYC targets V2",
                   "AR activity", "ARFL", "ARV", "BRCAness",
                   paste0("Cluster",1:11))

gsva_scores_reordered <- gsva_scores[ order_of_cols, ]

#pdf("TCGA_PRAD_GSVA_scores_plot.pdf", width = 15)

p <- pheatmap(
  gsva_scores_reordered,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col[, c(4:5)],
  breaks = seq(-0.5, 0.5, length.out = 101),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "GSVA Plot, TCGA PRAD"
)

print(p)


dev.off()


#pfs_data <- read.table(file="TCGA_PRAD_data_with_progression_free_survival_data.tsv",header = T,stringsAsFactors = F,sep="\t")

query.clinical <- read.csv(file="PCGA_PRAC_OS_PFS.csv",header = T, sep=",")



## Survival analysis
#clinical_surv <- query.clinical



gsva_patients <- substr(
  colnames(gsva_scores),
  1,
  12
)

#clinical_surv$patient_id <- clinical_surv$submitter_id

branch_df <- data.frame(
  patient_id = names(branch),
  Branch = paste0("Branch_", branch)
)

branch_df$patient_id <- substr(
  branch_df$patient_id,
  1,
  12
)

survival_df <- left_join(
  query.clinical,
  branch_df,
  by = c("patientId"="patient_id")
)

survival_df$PFS_event <- ifelse(
  survival_df$PFS_STATUS == "1:PROGRESSION",
  1,
  0
)


library(survival)
library(survminer)

surv_object <- Surv(
  time = survival_df$PFS_MONTHS,
  event = survival_df$PFS_event
)

fit <- survfit(
  surv_object ~ Branch,
  data = survival_df
)

pdf("Survival_plot_TCGA_PRAD_Overall_GSVA_clustering.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = "GSVA Branch"
)
dev.off()


survival_df$PFI_time <- survival_df$days_to_recurrence

#survival_df$PFI_event <- ifelse(
#  survival_df$progression_or_recurrence == "Yes",
#  1,
#  0
#)


## Survival analysis using Persistent
rownames(gsva_scores)[rownames(gsva_scores) == "Persistent clusters"] <- "Persistent_clusters"

gsva_Persistent_clusters <- gsva_scores[ "Persistent_clusters", ]

gsva_Persistent_clusters <- as.data.frame(gsva_Persistent_clusters)

gsva_Persistent_clusters$Fill_id <- rownames(gsva_Persistent_clusters)

gsva_Persistent_clusters$Patient_ID <- sapply(
  strsplit(gsva_Persistent_clusters$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_Persistent_clusters <- gsva_Persistent_clusters[,c("gsva_Persistent_clusters", "Patient_ID")]

survival_df_Persistent_clusters <- left_join(survival_df, gsva_Persistent_clusters, by = c("patientId"="Patient_ID"))

survival_df_Persistent_clusters <- survival_df_Persistent_clusters[!is.na(survival_df_Persistent_clusters$gsva_Persistent_clusters),]

survival_df_Persistent_clusters$Persistent_clusters_threshold <- ifelse(survival_df_Persistent_clusters$gsva_Persistent_clusters > median(survival_df_Persistent_clusters$gsva_Persistent_clusters), "High", "Low")


survival_df_Persistent_clusters$PFS_event <- ifelse(
  survival_df_Persistent_clusters$Persistent_clusters_threshold == "High",
  1,
  0
)

survival_df_Persistent_clusters$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_Persistent_clusters$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_Persistent_clusters$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_Persistent_clusters$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_Persistent_clusters$PFS_event_status <- as.numeric(
  survival_df_Persistent_clusters$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_Persistent_clusters$PFS_MONTHS,
  event = survival_df_Persistent_clusters$PFS_event_status
)

fit <- survfit(
  surv_object ~ Persistent_clusters_threshold,
  data = survival_df_Persistent_clusters
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_Persistent_clusters.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_Persistent_clusters,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off()



## Survival analysis using AR activity
gsva_AR_activity <- gsva_scores[ "AR activity", ]

gsva_AR_activity <- as.data.frame(gsva_AR_activity)

gsva_AR_activity$Fill_id <- rownames(gsva_AR_activity)

gsva_AR_activity$Patient_ID <- sapply(
  strsplit(gsva_AR_activity$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_AR_activity <- gsva_AR_activity[,c("gsva_AR_activity", "Patient_ID")]


survival_df_AR_activity <- left_join(survival_df, gsva_AR_activity, by = c("patientId"="Patient_ID"))


survival_df_AR_activity <- survival_df_AR_activity[!is.na(survival_df_AR_activity$gsva_AR_activity),]


survival_df_AR_activity$AR_activity_threshold <- ifelse(survival_df_AR_activity$gsva_AR_activity > median(survival_df_AR_activity$gsva_AR_activity), "High", "Low")


survival_df_AR_activity$PFS_event <- ifelse(
  survival_df_AR_activity$AR_activity_threshold == "High",
  1,
  0
)

survival_df_AR_activity$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_AR_activity$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_AR_activity$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_AR_activity$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_AR_activity$PFS_event_status <- as.numeric(
  survival_df_AR_activity$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_AR_activity$PFS_MONTHS,
  event = survival_df_AR_activity$PFS_event_status
)

fit <- survfit(
  surv_object ~ AR_activity_threshold,
  data = survival_df_AR_activity
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_AR_activity.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_AR_activity,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off()



## Survival analysis using persist
gsva_persist <- gsva_scores[ "Persist", ]

gsva_persist <- as.data.frame(gsva_persist)

gsva_persist$Fill_id <- rownames(gsva_persist)

gsva_persist$Patient_ID <- sapply(
  strsplit(gsva_persist$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_persist <- gsva_persist[,c("gsva_persist", "Patient_ID")]


survival_df_persist <- left_join(survival_df, gsva_persist, by = c("patientId"="Patient_ID"))


survival_df_persist <- survival_df_persist[!is.na(survival_df_persist$gsva_persist),]


survival_df_persist$persist_threshold <- ifelse(survival_df_persist$gsva_persist > median(survival_df_persist$gsva_persist), "High", "Low")


survival_df_persist$PFS_event <- ifelse(
  survival_df_persist$persist_threshold == "High",
  1,
  0
)

survival_df_persist$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_persist$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_persist$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_persist$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_persist$PFS_event_status <- as.numeric(
  survival_df_persist$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_persist$PFS_MONTHS,
  event = survival_df_persist$PFS_event_status
)

fit <- survfit(
  surv_object ~ persist_threshold,
  data = survival_df_persist
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_persist.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_persist,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ", ]


## Survival analysis using prosgenesis
gsva_prosgenesis <- gsva_scores[ "PROSGenesis", ]

gsva_prosgenesis <- as.data.frame(gsva_prosgenesis)

gsva_prosgenesis$Fill_id <- rownames(gsva_prosgenesis)

gsva_prosgenesis$Patient_ID <- sapply(
  strsplit(gsva_prosgenesis$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_prosgenesis <- gsva_prosgenesis[,c("gsva_prosgenesis", "Patient_ID")]


survival_df_prosgenesis <- left_join(survival_df, gsva_prosgenesis, by = c("patientId"="Patient_ID"))


survival_df_prosgenesis <- survival_df_prosgenesis[!is.na(survival_df_prosgenesis$gsva_prosgenesis),]


survival_df_prosgenesis$prosgenesis_threshold <- ifelse(survival_df_prosgenesis$gsva_prosgenesis > median(survival_df_prosgenesis$gsva_prosgenesis), "High", "Low")


survival_df_prosgenesis$PFS_event <- ifelse(
  survival_df_prosgenesis$prosgenesis_threshold == "High",
  1,
  0
)

survival_df_prosgenesis$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_prosgenesis$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_prosgenesis$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_prosgenesis$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_prosgenesis$PFS_event_status <- as.numeric(
  survival_df_prosgenesis$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_prosgenesis$PFS_MONTHS,
  event = survival_df_prosgenesis$PFS_event_status
)

fit <- survfit(
  surv_object ~ prosgenesis_threshold,
  data = survival_df_prosgenesis
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_prosgenesis.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_prosgenesis,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off()    




## Survival analysis using arfl
gsva_arfl <- gsva_scores[ "ARFL", ]

gsva_arfl <- as.data.frame(gsva_arfl)

gsva_arfl$Fill_id <- rownames(gsva_arfl)

gsva_arfl$Patient_ID <- sapply(
  strsplit(gsva_arfl$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_arfl <- gsva_arfl[,c("gsva_arfl", "Patient_ID")]


survival_df_arfl <- left_join(survival_df, gsva_arfl, by = c("patientId"="Patient_ID"))


survival_df_arfl <- survival_df_arfl[!is.na(survival_df_arfl$gsva_arfl),]


survival_df_arfl$arfl_threshold <- ifelse(survival_df_arfl$gsva_arfl > median(survival_df_arfl$gsva_arfl), "High", "Low")


survival_df_arfl$PFS_event <- ifelse(
  survival_df_arfl$arfl_threshold == "High",
  1,
  0
)

survival_df_arfl$PFS_event_status <- sapply(
  strsplit(as.character(survival_df_arfl$PFS_STATUS), ":"),
  function(x) x[1]
)

survival_df_arfl$PFS_status_label <- sapply(
  strsplit(as.character(survival_df_arfl$PFS_STATUS), ":"),
  function(x) x[2]
)

survival_df_arfl$PFS_event_status <- as.numeric(
  survival_df_arfl$PFS_event_status
)

surv_object <- Surv(
  time = survival_df_arfl$PFS_MONTHS,
  event = survival_df_arfl$PFS_event_status
)

fit <- survfit(
  surv_object ~ arfl_threshold,
  data = survival_df_arfl
)

pdf("Survival_plot_TCGA_PRAD_using_GSVA_for_ARFL.pdf", width = 15)
ggsurvplot(
  fit,
  data = survival_df_arfl,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  title = "TCGA-PRAD GSVA-derived branches",
  legend.title = ""
)
dev.off() 




## Survival analysis using ENZ induced clusters
gsva_enz_induced <- gsva_scores[ "ENZ_induced_clusters", ]

gsva_enz_induced <- as.data.frame(gsva_enz_induced)

gsva_enz_induced$Fill_id <- rownames(gsva_enz_induced)

gsva_enz_induced$Patient_ID <- sapply(
  strsplit(gsva_enz_induced$Fill_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

gsva_enz_induced <- gsva_enz_induced[,c("gsva_enz_induced", "Patient_ID")]


survival_df_enz_induced <- left_join(survival_df, gsva_enz_induced, by = c("patientId"="Patient_ID"))


survival_df_enz_induced <- survival_df_enz_induced[!is.na(survival_df_enz_induced$gsva_enz_induced),]


survival_df_enz_induced$enz_induced_threshold <- ifelse(survival_df_enz_induced$gsva_enz_induced > median(survival_df_enz_induced$gsva_enz_induced), "High", "Low")


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


