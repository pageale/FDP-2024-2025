# TISSUE-OF-ORIGIN ANALYSIS
This folder contains workflows and scripts for analyzing tissue-of-origin signals based on cell-free DNA (cfDNA) fragmentation and epigenetic patterns, adapted from the published pipelines by the Shendure and Kircher labs.


## Contents
1. **Modified snakemake workflows**
   * **snakefile_WPS.smk:**: Workflow for calculating the Windowed Protection Score (WPS).
     
     **Input:** A TSV file containing sample IDs and their respective BAM file paths.
   * **snakefile_GE_analysis.smk:** Performs gene expression (GE) analysis. Uses FPKM values for 20,344 Ensembl gene IDs across 44 human cell lines and 32 primary tissues from the Human Protein Atlas (Uhlén et al., 2015). It was downloaded from [http://www.proteinatlas.org/download/rna.csv.zip](http://www.proteinatlas.org/download/rna.csv.zip).  Fast Fourier Transformation (FFT) is applied to cfDNA signals, and the resulting profiles are correlated with gene expression.
   * **GE_unsupervised.smk:**  Unsupervised analysis workflow. Contains utility functions to calculate and visualize similarity relationships between samples.
  
     **NOTE:** More detailed information and original scripts can be found at [Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](https://github.com/kircherlab/cfDNA?tab=readme-ov-file) repository. Original code can be found [Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](https://github.com/shendurelab/cfDNA) repository.
     **PAPER:**
       * [Cell type signatures in cell-free DNA fragmentation profiles reveal disease biology](https://www.nature.com/articles/s41467-024-46435-0)
       * [Cell-free DNA comprises an in vivo nucleosome footprint that informs its tissues-of-origin](https://pmc.ncbi.nlm.nih.gov/articles/PMC4715266/)

3. **Analysis**
   * **ranking-tissue-of-origin.R:**  R script for plotting ranked tissue/cell line correlations over time points to visualize shifts in tissue-of-origin signatures.
     
4. **Workflow (scripts)**
   * **fft_table.py**
   * **snakemake_correlation_rank_table.R**
   * **nakemake_correlation_table.R**

    **NOTE:** These scripts were modified from their original versions in the [Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](https://github.com/kircherlab/cfDNA?tab=readme-ov-file). Other necessary codes are provided there. 
   
6. **Reference**
   * **healthy_Song_100bp_hg38.bed:**  Used as **regions.tsv**. BED file containing regions of interest (e.g., transcription factor binding sites), all of the same length.


## Usage
* **snakefile_WPS.smk:**
Follow instructions from the [Kircher lab cfDNA](https://github.com/kircherlab/cfDNA?tab=readme-ov-file) workflow to install dependencies and activate the conda environment. The different analyses have to be executed separately. To specify the respective workflow use the -s switch followed by the path of the Snakefile (e.g.: ./snakefile_WPS.smk)

  **Input:**
  - configured by the user (samples.tsv):
    - **ID:**  Project or analysis identifier.
    - **sample:** Name of the sample.
    - **path:** Absolute path to the sample's BAM file.
    - **ref_samples:** Comma-separated list of reference samples for comparisons.
    - **genome_build:** Genome assembly used for alignment (e.g., GRCh37 or GRCh38).
      
    
  | ID            | sample        | path                          | ref_samples              | genome_build |
  |:--------------|:--------------|:-------------------------------|:--------------------------|:--------------|
  | experimentID   | testsample1    | /path/to/testsample1.bam        | testsample2,testsample3    | GRCh38        |
  | experimentID   | testsample2    | /path/to/testsample2.bam        | testsample1,testsample3    | GRCh38        |
  | experimentID   | testsample3    | /path/to/testsample3.bam        | testsample1,testsample2    | GRCh38        |

  - configured by the user (regions.tsv](config/regions.tsv):
      - bed file containing regions of interest (e.g. TFBS), all having the same length
        
  | target            | path        |
  |:--------------|:--------------|
  | region_of_interest   | path/to/bed.bed    | 
  | region_of_interest2   | path/to/bed.bed       | 

   **Output:**
  - Table containing bp specific WPS for regions listed in bed
  - Line plot showing normalized WPS of multiple samples

  **Commands:**
  
  Activate the conda environment:
  
  ```bash
  conda activate snakemake
  ```
  
  Test your configuration by performing a dry-run via
  
  ```bash
  snakemake -s snakefile_WPS.smk --use-conda -n
  ```

---
* **snakefile_GE_analysis.smk:**

  **Input**
  
  - included in the [original repository](https://github.com/kircherlab/cfDNA?tab=readme-ov-file):
      - annotations
      - labels
      - RNAtable from Protein Atlas 
          - blood atlas ["Blood"]
          - protein atlas tissues ["Tissue"]
          - protein atlas tissues + cell lines ["Extended"]
  - configured by the user (samples.tsv):
      - ID / patient
      - samples
      - path to sample BAM files
      - genome build per sample
  
  **Output**
  
  - fft_summary tables (results/intermediate/body/fft_summaries)
  - plots showing intensities across tissues (results/plots)
  - table showing correlation with tissues/cell lines
  - table showing correlation rank difference to reference sample

  **Commands:**
  
  Activate the conda environment:
  
  ```bash
  conda activate snakemake
  ```
  
  Test your configuration by performing a dry-run via
  
  ```bash
  snakemake -s snakefile_GE_analysis.smk --use-conda -n
  ```
---

* **ranking-tissue-of-origin.R:** 

  **Input:**
    * *_correlation_results.tsv  files from the **snakefile_GE_analysis.smk** workflow.
    * Sample order list (change variable inside the code).
      
  **Output:** Plots displaying ranked tissue/cell line correlations across multiple time points.

  **Libraries:**


   ```R
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(paletteer)
  ```
  
  **Installation:**

  ```R
  install.packages(c("ggplot2", "dplyr", "readr", "RColorBrewer", "paletteer"))
  ```





