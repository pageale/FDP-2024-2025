#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Use: $0 <bam_folder> <output_folder>"
    exit 1
fi

bam_folder=$1
output_folder="$2/bam"
follder_picard_output="$2/picard"

mkdir -p "$output_folder" "$follder_picard_output"

for bam in "$bam_folder"/*/*.bam; do
    sample=$(basename "$bam" .bam)

    echo "Processing $sample..."

    samtools view -b -F 4 -q 20 "$bam" > "$output_folder/${sample}.filtered.bam"

    samtools index "$output_folder/${sample}.filtered.bam"
    
    picard CollectInsertSizeMetrics I="$output_folder/${sample}.filtered.bam" \
            O="$follder_picard_output/$sample.txt" \
            H="$follder_picard_output/$sample.pdf" M=0.5
    
    echo "Done: $sample."
done

echo "All BAM files have been filtered, indexed and extracted size metrics."

