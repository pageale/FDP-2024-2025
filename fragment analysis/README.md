# Fragmentomics Analysis

## Content

Collecting fragment length for:

**NANOPORE**
1. **bam-nanopore-to-nanoplot.sh:** This script collects the metrics
2. **change_format_nanoplot.py:**

**ILLUMINA**
1. **bam-to-picard.sh:**

## Usage

* **bam-nanopore-to-nanoplot.sh:**
requirements:
NanoPlot
samtools
parallel


    ```sh
    ./bam-nanopore-to-nanoplot.sh <folder_bam> <output_folder>
    ```


Installation:

----
* **bam-to-picard.sh:**
samtools
picard 

  **Command:**
  
  ```sh
  ./bam-to-picard.sh <folder_bam> <output_folder>
  ```
  
  **Installation:**
  
  ```sh
  conda install bioconda::picard
  ```

