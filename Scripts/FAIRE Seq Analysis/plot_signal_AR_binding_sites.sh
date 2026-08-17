computeMatrix reference-point \
    --referencePoint center \
    -S ./bigwig/LNCaP-DMSO.CPM.average.bw ./bigwig/LNCaP-R-MDV.CPM.average.bw ./bigwig/LNCaP-DHT.CPM.average.bw   \
    -R /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/AR_binding_sites.txt \
    -a 1000 \
    -b 1000 \
	-bs 25 \
	-o /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_AR_binding_sites.gz


plotProfile \
    -m /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_AR_binding_sites.gz \
	--samplesLabel "DMSO, Avg" "ENZ, Avg" "DHT, Avg" \
	--averageType "mean" \
	--plotType lines \
	--refPointLabel "Distance for MYC binding site" \
	--perGroup \
    -out /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_AR_binding_sites.pdf