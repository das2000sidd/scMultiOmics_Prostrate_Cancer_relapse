bigwigAverage -b ./bigwig/LNCaP-DMSO_1.CPM.bw ./bigwig/LNCaP-DMSO_2.CPM.bw ./bigwig/LNCaP-DMSO_3.CPM.bw -o ./bigwig/LNCaP-DMSO.CPM.average.bw

bigwigAverage -b ./bigwig/LNCaP-R-MDV_1.CPM.bw ./bigwig/LNCaP-R-MDV_2.CPM.bw ./bigwig/LNCaP-R-MDV_3.CPM.bw -o ./bigwig/LNCaP-R-MDV.CPM.average.bw

bigwigAverage -b ./bigwig/LNCap_DHT_1.sorted.bw ./bigwig/LNCap_DHT_2.CPM.bw ./bigwig/LNCap_DHT_3.CPM.bw -o ./bigwig/LNCaP-DHT.CPM.average.bw



computeMatrix reference-point \
    --referencePoint center \
    -S ./bigwig/LNCaP-DMSO.CPM.average.bw ./bigwig/LNCaP-R-MDV.CPM.average.bw ./bigwig/LNCaP-DHT.CPM.average.bw   \
    -R /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/MYC_binding_sites.txt \
    -a 1000 \
    -b 1000 \
	-bs 25 \
	-o /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_MYC_binding_sites.gz


plotProfile \
    -m /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_MYC_binding_sites.gz \
	--samplesLabel "DMSO, Avg" "ENZ, Avg" "DHT, Avg" \
	--averageType "mean" \
	--plotType lines \
	--refPointLabel "Distance for MYC binding site" \
	--perGroup \
    -out /Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/Plot_at_MYC_binding_sites.pdf

