#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Use: $0 <fastq_folder>"
    exit 1
fi

GENOME_SIZE=3100000000
FASTQ_FOLDER="$1"

if [ ! -d "$FASTQ_FOLDER" ]; then
    echo "Error: Folder '$FASTQ_FOLDER' does not exist."
    exit 1
fi

if ! ls "$FASTQ_FOLDER"/*_R1.fastq.gz &>/dev/null; then
    echo "Error: FASTQ files not found in '$FASTQ_FOLDER'"
    exit 1
fi

LOG_FILE="illumina_fastq_downsampling-2.txt"
OUT_DIR="illumina_fastq_downsampling"
mkdir -p "$OUT_DIR"

SEED=42

for r1 in "$FASTQ_FOLDER"/*_R1.fastq.gz; do
    [ -e "$r1" ] || continue  # check empty folder

    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$r1" | sed 's/_R1.fastq.gz//')

    echo "-------------------------------------------------" | tee -a "$LOG_FILE"
    echo "Processing sample: $SAMPLE" | tee -a "$LOG_FILE"
    echo "-------------------------------------------------" | tee -a "$LOG_FILE"

    TOTAL_READS=$(zcat "$r1" | wc -l)
    TOTAL_READS=$((TOTAL_READS / 4))
    READ_LENGTH=$(zcat "$r1" | sed -n '2p' | tr -d '\n' | wc -c)
    TOTAL_BASES=$((TOTAL_READS * READ_LENGTH * 2))
    LC_NUMERIC=C
    INIT_COV=$(echo "scale=5; $TOTAL_BASES / $GENOME_SIZE" | bc -l)
    INIT_COV=$(printf "%.5f" "$INIT_COV")

    echo "Sample: $SAMPLE" | tee -a "$LOG_FILE"
    echo "Total reads: $TOTAL_READS" | tee -a "$LOG_FILE"
    echo "Read length: $READ_LENGTH" | tee -a "$LOG_FILE"
    echo "Total sequenced bases: $TOTAL_BASES" | tee -a "$LOG_FILE"
    echo "Initial coverage: ${INIT_COV}x" | tee -a "$LOG_FILE"
    echo "-------------------------------------------------" | tee -a "$LOG_FILE"

    MIN_COV=0.01
    COVERAGE_LIST=$(awk -v start="$INIT_COV" -v min="$MIN_COV" '
    BEGIN {
        step = 0.1;
        for (i = start; i > 0.1; i -= step) {
            printf "%.5f\n", i;
        }
        step = 0.01;
        for (i = 0.1; i >= min; i -= step) {
            printf "%.5f\n", i;
        }
    }')

    export r1 r2 SEED OUT_DIR GENOME_SIZE READ_LENGTH LOG_FILE SAMPLE INIT_COV

    echo "$COVERAGE_LIST" | parallel -j 7 --env INIT_COV --env r1 --env r2 --env SEED --env OUT_DIR --env GENOME_SIZE --env READ_LENGTH --env LOG_FILE --env SAMPLE '
        COV={}
        LC_NUMERIC=C
        COV=$(echo $COV | sed "s/,/./g")
        INIT_COV=$(echo $INIT_COV | sed "s/,/./g")
        FRAC=$(echo "scale=3; $COV / $INIT_COV" | bc -l)

        if (( $(echo "$FRAC <= 0" | bc -l) )); then
            FRAC=0.01
        fi

        OUT_R1="${OUT_DIR}/${SAMPLE}_${COV}x_R1.fastq.gz"
        OUT_R2="${OUT_DIR}/${SAMPLE}_${COV}x_R2.fastq.gz"

        echo "  Downsampling ${SAMPLE} a ${COV}x (factor: $FRAC)" | tee -a "$LOG_FILE"

        seqtk sample -s$SEED "$r1" $FRAC | gzip > "$OUT_R1"
        seqtk sample -s$SEED "$r2" $FRAC | gzip > "$OUT_R2"

        NEW_READS=$(zcat "$OUT_R1" | wc -l)
        NEW_READS=$((NEW_READS / 4))
        NEW_BASES=$((NEW_READS * READ_LENGTH * 2))
        NEW_COV=$(echo "scale=5; $NEW_BASES / $GENOME_SIZE" | bc -l)

        echo "  Reads after downsampling: $NEW_READS" | tee -a "$LOG_FILE"
        echo "  Coverage after downsampling: ${NEW_COV}x" | tee -a "$LOG_FILE"
        echo "-------------------------------------------------" | tee -a "$LOG_FILE"
    '

    echo "Downsampling completed for $SAMPLE."

done
