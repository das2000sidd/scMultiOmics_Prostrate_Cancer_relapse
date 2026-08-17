for bam in *.sorted.bam
do
    sample=$(basename "$bam" .sorted.bam)

    bamCoverage \
        -b "$bam" \
        -o "./bigwig/${sample}.CPM.bw" \
        --normalizeUsing CPM \
        --binSize 10 \
        --ignoreDuplicates \
        -p 8
done