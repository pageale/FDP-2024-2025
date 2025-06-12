#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Use: $0 <bam_folder> <output_folder>"
    exit 1
fi

bam_folder=$1
output_folder="$2"

mkdir -p "$output_folder"

for bam in "$bam_folder"/*/*.bam; do
    sample=$(basename "$bam" .sorted.aligned.bam)

    echo "Processing $sample..."
    methylartist wgmeth \
        -b "$bam" \
        -m h \
        -p 6 \
        --ref ~/Reference/bwa_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -o "${output_folder}/${sample}_h_wgmeth" \
        --dss \
        -f ~/Reference/bwa_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai \
        --motif C

    echo "Done: $sample."
done




