#!/bin/bash

calc_coverage() {
  BAM=$1
  echo "Processing '$BAM'..."
  STATS=$(samtools depth -a "$BAM" | awk '{sum+=$3} END { print sum/NR }')
  echo -e "${BAM}\t${STATS}"
  echo "End '$BAM'"

}

export -f calc_coverage

find . -name "*.bam" | parallel --citation -j 10 calc_coverage > bam_coverage_all.txt
