#!/bin/bash

# Input directory
FASTQ_DIR="/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS"

# Reference genome
GENOME="/Users/siddharthadas/Desktop/PROSTRATE_CANCER_RELAPSE_PAIRED_scRNA_AND_scATAC/FAIRE_SEQ_ANALYSIS/hg19.fa"

# Threads
THREADS=8

#cd "$FASTQ_DIR"

for FASTQ in *.fastq
do
    SAMPLE=$(basename "$FASTQ" .fastq)

    echo "Processing ${SAMPLE}"

    # Align reads
    bwa mem -t $THREADS \
        "$GENOME" \
        "$FASTQ" \
        > "${SAMPLE}.sam"

    # Convert SAM to sorted BAM
    samtools view -bS "${SAMPLE}.sam" | \
        samtools sort -@ $THREADS \
        -o "${SAMPLE}.sorted.bam"

    # Index BAM
    samtools index "${SAMPLE}.sorted.bam"

    # Mark duplicates
    java -jar /Users/siddharthadas/tools/picard/picard.jar MarkDuplicates \
        I="${SAMPLE}.sorted.bam" \
        O="${SAMPLE}.markdup.bam" \
        M="${SAMPLE}.metrics.txt" \
        REMOVE_DUPLICATES=false

    # Index duplicate-marked BAM
    samtools index "${SAMPLE}.markdup.bam"

    # Remove intermediate SAM file
    rm "${SAMPLE}.sam"

    echo "${SAMPLE} completed"
done