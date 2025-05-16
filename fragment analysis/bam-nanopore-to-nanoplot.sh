#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Use: $0 <folder_bam> <output_folder>"
    exit 1
fi

bam_folder=$1
output_folder=$2
filtered_bam_folder=$output_folder/bam_filtered
nanoplot_output_folder=$output_folder/nanoplot

mkdir -p "$filtered_bam_folder"
mkdir -p "$nanoplot_output_folder"

echo "Filtering BAMs..."
for bam in "$bam_folder"/*.bam; do
    sample=$(basename "$bam" .bam)

    samtools view -b -F 4 -q 20 "$bam" > "$filtered_bam_folder/${sample}.filtered.bam"
    samtools index "$filtered_bam_folder/${sample}.filtered.bam"
done
echo "Filtering completed."

echo "Nanoplot in parallel..."
ls "$filtered_bam_folder"/*.filtered.bam | parallel -j 4 '
    sample=$(basename {} .filtered.bam)
    NanoPlot --bam {} --raw --tsv_stats --alength -p "${sample}_" -o '"$nanoplot_output_folder"'/${sample}
'

echo "Fragment size extraction using Nanoplot finished."
