# scMultiOmics_Prostrate_Cancer_relapse
Reproducible implementation of a Nature Communications 2021 single-cell multi-omics study integrating scRNA-seq and scATAC-seq to characterize persistent cancer cell states.

In this project, I reanalysed a prostate cancer (PC) dataset to identify pre-existing and treatment persistent cell subpopulations that possess regenerative potential when subjected to treatment focusing on androgen receptor activity.

Briefly, in the paper, LNCaP, a metastatic prostate cancer cell line was used to derive ENZ-resistant cell lines RES-A and RES-B generated via long-term exposure to AR targeting agents. ENZ is an AR antagonist that works by inhibiting AR signaling through various mechanisms. The paper identified the emergence of resistance in the epithelial-derived component of prostate tumors in ENZ-exposed and -resistant PC cell lines at a single-cell level to explore how heterogeneous PCs respond to AR signaling inhibitors. Through enrichment analysis of transcriptional signals from molecular gene classifiers derived in this study, they showed evidence of treatment-persistent and pre-existing PC cells that can predict treatment response in
both primary and advanced patients.


First to study, the molecular consequences of AR signaling suppression and drug-resistance dynamics in PC, they compared open chromatin regions in parental cell line to AR blocking through ENZ treatment and two other resistant cell lines. They used scATAC seq for this assay. 
The found that ATAC seq signal transcription-start sites (TSS) decreased in ENZ-resistant cells compared with the parental genomewide (Figure 1b).



