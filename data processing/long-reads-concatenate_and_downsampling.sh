#!/bin/bash

GENOME_SIZE=3100000000  # genome size
SEED=42  # Seed
LOG_FILE="downsampling_nanopore_picard_log.txt" # log file name
OUT_DIR="downsampling_nanopore_second_run" # outpur directory
mkdir -p "$OUT_DIR"

echo "Date: $(date)" | tee "$LOG_FILE"

# USE: indicate what folder to pass
if [ $# -ne 1 ]; then
    echo "Use: $0 <directorio_fastq>"
    exit 1
fi

INPUT_DIR="$1"
echo "Input directory: $INPUT_DIR" | tee -a "$LOG_FILE"

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' not found." | tee -a "$LOG_FILE"
    exit 1
fi

MERGED_DIR="${INPUT_DIR}/merged"

if [ ! -d "$MERGED_DIR" ]; then
    echo "Creating directory '$MERGED_DIR'." | tee -a "$LOG_FILE"
    mkdir -p "$MERGED_DIR"
    for barcodeDir in "$INPUT_DIR"/IMN*; do
        barcode=$(basename "$barcodeDir")

        FASTQ_COUNT=$(find "$barcodeDir" -type f -name "*.fastq.gz" | wc -l)

        if [ "$FASTQ_COUNT" -lt 10 ]; then
            echo "Warning: Folder '$barcodeDir' only has $FASTQ_COUNT FASTQ files. Skipping concatenation." | tee -a "$LOG_FILE"
            continue
        fi

        cat "$barcodeDir"/*.fastq.gz > "$MERGED_DIR/${barcode}.all.fastq.gz"
        echo "Concatenating: $barcode with $FASTQ_COUNT files" | tee -a "$LOG_FILE"
    done
else
    echo "Skipping barcoding concatenation step..."
fi

echo "FASTQ concatenation finished. Next step: calculating metrics and downsampling..." | tee -a "$LOG_FILE"

# Downsampling function to call via parallel
downsample_sample() {
    COV=$1
    LC_NUMERIC=C
    COV=$(echo "$COV" | sed "s/,/./g")
    INIT_COV=$(echo "$INIT_COV" | sed "s/,/./g")
    FRAC=$(echo "scale=3; $COV / $INIT_COV" | bc -l)
    if (( $(echo "$FRAC <= 0" | bc -l) )); then
        FRAC=0.01
    fi

    OUT_R1="${OUT_DIR}/${SAMPLE}_${COV}x.fastq.gz"

    echo "  Downsampling ${SAMPLE} a ${COV}x (factor: $FRAC)" | tee -a "$LOG_FILE"
    seqtk sample -s"$SEED" "$r1" "$FRAC" | gzip > "$OUT_R1"

    NEW_READS=$(zcat "$OUT_R1" | wc -l)
    NEW_READS=$((NEW_READS / 4))
    NEW_BASES=$(zcat "$OUT_R1" | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
    NEW_COV=$(echo "scale=5; $NEW_BASES / $GENOME_SIZE" | bc -l)

    echo "  Reads after downsampling: $NEW_READS" | tee -a "$LOG_FILE"
    echo "  FASTQ coverage after downsampling: ${NEW_COV}x" | tee -a "$LOG_FILE"
    echo "-------------------------------------------------" | tee -a "$LOG_FILE"
}
export -f downsample_sample

# Iterate over each sample
for r1 in "$MERGED_DIR"/*.all.fastq.gz; do
    SAMPLE=$(basename "$r1" .all.fastq.gz)

    echo "-------------------------------------------------" | tee -a "$LOG_FILE"
    echo "Processing sample: $SAMPLE" | tee -a "$LOG_FILE"
    echo "-------------------------------------------------" | tee -a "$LOG_FILE"

    echo "Calculating total bases directly from FASTQ..." | tee -a "$LOG_FILE"

    TOTAL_BASES=$(zcat "$r1" | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
    TOTAL_READS=$(zcat "$r1" | wc -l)
    TOTAL_READS=$((TOTAL_READS / 4))
    LC_NUMERIC=C
    INIT_COV=$(echo "scale=5; $TOTAL_BASES / $GENOME_SIZE" | bc -l)

    echo "Total reads: $TOTAL_READS" | tee -a "$LOG_FILE"
    echo "Total bases: $TOTAL_BASES" | tee -a "$LOG_FILE"
    echo "Cobertura inicial: ${INIT_COV}x" | tee -a "$LOG_FILE"

    MIN_COV=0.01
    COVERAGE_LIST=$(awk -v start="$INIT_COV" -v min="$MIN_COV" '
    BEGIN {
        step = 0.1;
        for (i = start; i >= 0.1; i -= step)
            printf "%.5f\n", i;
        step = 0.01;
        for (i = 0.1; i >= min; i -= step)
            printf "%.5f\n", i;
    }')

    LC_NUMERIC=C
    export r1 SEED OUT_DIR GENOME_SIZE LOG_FILE SAMPLE INIT_COV

    if ! command -v seqtk &> /dev/null; then
        echo "Error: seqtk is not installed or not in PATH." | tee -a "$LOG_FILE"
        exit 1
    fi

    echo "$COVERAGE_LIST" | parallel -j 8 downsample_sample {}

    echo "Downsampling completado para $SAMPLE."
done

echo "Proceso finalizado."
