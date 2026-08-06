Reproduction of Single-Cell RNA-seq and ATAC-seq Analysis

Overview
This repository contains my independent reproduction of the computational analyses presented in the following publication:
Citation: Taavitsainen S, Engedal N, Cao S, Handle F, Erickson A, Prekovic S, Wetterskog D, Tolonen T, Vuorinen EM, Kiviaho A, Nätkin R, Häkkinen T, Devlies W, Henttinen S, Kaarijärvi R, Lahnalampi M, Kaljunen H, Nowakowska K, Syvälä H, Bläuer M, Cremaschi P, Claessens F, Visakorpi T, Tammela TLJ, Murtola T, Granberg KJ, Lamb AD, Ketola K, Mills IG, Attard G, Wang W, Nykter M, Urbanucci A. Single-cell ATAC and RNA sequencing reveal pre-existing and persistent cells associated with prostate cancer relapse. Nat Commun. 2021 Sep 6;12(1):5307. doi: 10.1038/s41467-021-25624-1. PMID: 34489465; PMCID: PMC8421417.

The objective of this project was to reproduce the published figures and analytical workflow using the publicly available data for the project from the GEO repository and to gain hands-on experience with modern single-cell transcriptomic and epigenomic analysis pipelines.
This repository is intended as a portfolio project demonstrating proficiency in reproducible computational biology, single-cell data analysis, and scientific programming.
Project Goals
•	Reproduce key analyses from the published study.
•	Validate reproducibility of the computational workflow.
•	Implement the analysis from raw/publicly available data.
•	Organize the workflow into reusable and well-documented scripts.

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
QC scATAC seq/
QC scRNA seq/

Scripts
	Preprocessing/
	Spatial_RNASeq
	TCGA_Analysis
	single_cell_analysis/
		Chen_et_al/
		Taavitsainen_et_al/

README.md
sessionInfo
```

Analyses Performed
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

Additionally, for TCGA Prostate Cancer (PRAD) data, the data was downloaded using TCGAbiolinks R package and processed using DESeq2 and GSVA package to generate scores for predefined gene sets for subsequent survival analysis and KM curve generation.

Technologies Used
•	R
•	Seurat v5
•	Signac
•	Harmony / FastMNN (where applicable)
•	ggplot2
•	patchwork
•	dplyr
•	Matrix
•	GenomicRanges
•	Bioconductor packages

Results
The reproduced figures are available in the Figures/ directory, including:
•	Quality-control visualizations
•	Dimensionality reduction (UMAP)
•	Cluster annotations
•	Marker gene expression
•	Figures corresponding to the publication
Where differences from the original publication exist, they are likely attributable to software version differences, random initialization, or package updates.

Reproducibility
The analysis scripts are provided in the approximate order they were executed.
The workflow assumes access to the publicly available datasets described in the original publication.

Disclaimer
This repository is an independent educational reproduction of a published study. All scientific credit belongs to the original authors. This repository is intended solely to demonstrate computational and bioinformatics skills.

Skills Demonstrated
•	Reproducible computational research
•	Single-cell RNA-seq analysis
•	Single-cell ATAC-seq analysis
•	Multiomic data integration
•	Statistical analysis in R
•	Data visualization
•	Scientific programming
•	Workflow organization




