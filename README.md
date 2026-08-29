Reproduction of Single-Cell RNA-seq and ATAC-seq Analysis

Overview
This repository contains my independent reproduction of the computational analyses presented in the following publication:

> **Taavitsainen S, Engedal N, Cao S, et al.**  
> Single-cell ATAC and RNA sequencing reveal pre-existing and persistent cells associated with prostate cancer relapse.  
> *Nature Communications* (2021) 12:5307.  
> DOI: 10.1038/s41467-021-25624-1


The objective of this project was to reproduce the published figures and analytical workflow using the publicly available data for the project from the GEO repository and to gain hands-on experience with modern single-cell transcriptomic and epigenomic analysis pipelines.
This repository is intended as a portfolio project demonstrating proficiency in reproducible computational biology, single-cell data analysis, and scientific programming.

## Project Goals
- Reproduce key analyses from the published study.
- Validate reproducibility of the computational workflow.
- Implement the analysis from raw/publicly available data.

Repository Structure
```
Figures/
	Figure1/
	Figure2/
	Figure3/
	Figure4/
	Figure5/
	Figure6/
    Figure7/


Scripts
	FAIRE Seq Analysis
	TCGA_Analysis
	single_cell_analysis/
		Chen_et_al/
		Taavitsainen_et_al/

README.md

```

Analyses Performed in individual omics
```
The workflow includes following for single cell RNA-seq or scATAC-seq:
•	Data import
•	Quality control
•	Cell filtering
•	Normalization
•	Feature selection
•	Dimensionality reduction (PCA, UMAP)
•	Batch correction
•	Clustering
•	Differential expression analysis for RNA seq / differential accessible region  + motif enrichment for ATAC-seq anaysis
•	Reproduction of published figures
```
Additionally, for TCGA Prostate Cancer (PRAD) data, the data was downloaded using TCGAbiolinks R package and processed using DESeq2 and GSVA package to generate scores for predefined gene sets for subsequent survival analysis and KM curve generation.


```
Integrative analysis of scRNA-seq and scATAC-seq to identify
```

Integration of scRNA and scATAC to investigate transcription factor (TF) regulatory activity and chromatin accessibility changes in prostate cancer cells following androgen receptor (AR) pathway inhibition and in resistant cell states.

The main objective is to determine whether changes in TF activity inferred from gene expression are accompanied by corresponding changes in accessibility of TF binding motifs.

Following was the analysis workflow:

SCENIC / AUCell analysis

```
Single-cell RNA-seq data were analysed using SCENIC-derived AUCell scores to estimate transcription factor regulon activity at the single-cell level.

The analysis included:

TF regulon activity distributions across treatments
Treatment-specific median regulon activity
Changes relative to DMSO
Kruskal-Wallis testing across all treatments
Pairwise Wilcoxon/Dunn testing
Multiple-testing correction using the Benjamini-Hochberg FDR procedure

```

chromVAR analysis

```
Single-cell ATAC-seq data were analysed to quantify TF motif accessibility using chromVAR.

Peak regions were annotated with JASPAR transcription factor motifs and GC bias was calculated before estimating chromVAR deviation scores.

Motifs corresponding to TFs including:

E2F1
MYC
JUND
ELK4
ETV1
HOXB13
RELA
SOX9
ELF1
GRHL2

were examined.

For each TF motif, treatment-specific median chromVAR deviation scores and changes relative to DMSO were calculated.

Pairwise Wilcoxon tests were performed between treatment groups, followed by Benjamini-Hochberg FDR correction.

```

SCENIC–chromVAR integration

```

SCENIC and chromVAR results were integrated to compare two complementary measurements:

SCENIC/AUCell:
TF regulatory activity inferred from transcriptional target-gene expression.

chromVAR:
Accessibility of genomic regions containing TF binding motifs.

This provides a framework for distinguishing TFs showing:

Increased regulon activity and increased motif accessibility
Increased regulon activity without increased motif accessibility
Decreased regulon activity and decreased motif accessibility
Opposing transcriptional and chromatin changes

Because not every SCENIC TF had a directly corresponding motif in the selected chromVAR motif set, TF-to-motif mappings were explicitly defined rather than assuming that every TF could be compared directly.

```

Main outputs

```

The analysis generates:

Treatment-level SCENIC TF activity summaries
SCENIC differential activity statistics
Treatment-level chromVAR motif accessibility summaries
chromVAR differential accessibility statistics
SCENIC TF activity heatmaps
chromVAR motif accessibility heatmaps
Integrated SCENIC–chromVAR comparison tables

The heatmaps display changes relative to DMSO, with statistical significance indicated using FDR-based significance annotations.

```

Interpretation

```

The integrated analysis is intended to identify regulatory programs associated with AR inhibition and resistance.

In particular, TFs such as AR, MYC, E2F1, JUN/FOS and other prostate cancer-associated regulators can be examined for coordinated changes between transcriptional activity and chromatin accessibility.

A concordant increase in SCENIC activity and chromVAR deviation may provide evidence that increased TF regulatory activity is accompanied by increased accessibility of its binding sites. Conversely, discordant patterns may indicate that TF activity changes are not simply explained by changes in global motif accessibility and may involve additional regulatory mechanisms.

```

## Technologies Used

```
R
Seurat v5
Signac
Harmony / FastMNN (where applicable)
ggplot2
patchwork
dplyr
Matrix
GenomicRanges
Bioconductor packages
```


Results

The reproduced figures are available in the Figures/ directory, including:

```
•	Quality-control visualizations
•	Dimensionality reduction (UMAP)
•	Cluster annotations
•	Marker gene expression
•	Figures corresponding to the publication
```
Where differences from the original publication exist, they are likely attributable to software version differences, random initialization, or package updates.

Reproducibility
The analysis scripts are provided in the approximate order they were executed.
The cell ranger processed data was downloaded using accession id's provided in the paper. 

Disclaimer
This repository is an independent educational reproduction of a published study. All scientific credit belongs to the original authors. This repository is intended solely to demonstrate computational and bioinformatics skills.

Skills Demonstrated
```
•	Reproducible computational research
•	Single-cell RNA-seq analysis
•	Single-cell ATAC-seq analysis
•	Statistical analysis in R
•	Data visualization
•	Scientific programming
```
