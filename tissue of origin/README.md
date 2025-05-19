# TISSUE-OF-ORIGIN ANALYSIS

## Contents
1. **Modified snakemake workflow**
   * **snakefile_WPS.smk:**: This workflow calculates the Windowed Protection Score (WPS). It receives as input a tsv with the sample ID and its respective path to the sample BAM file. 
   * **snakefile_GE_analysis.smk:** Performs the gene expression (GE) analysis. FPKM gene expression (GE) values measured for 20,344 Ensembl gene identifiers in 44 human cell lines and 32 primary tissues by the Human Protein Atlas (Uhlén et al., 2015) was downloaded from [http://www.proteinatlas.org/download/rna.csv.zip](http://www.proteinatlas.org/download/rna.csv.zip). It does Fast Fourier transformation (FFT), which the  will be correlated with expression. 
   * **GE_unsupervised.smk:** This workflow performs an unsupervised analysis. Contains utility functions to calculate and visualize similarities between samples. 
  
     **NOTE:** More detailed information and original scripts can be found at [Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](https://github.com/kircherlab/cfDNA?tab=readme-ov-file) repository. Original code can be found [Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA](https://github.com/shendurelab/cfDNA) repository.

  ORIGINAL CODE AND PAPER
  Github:   https://github.com/shendurelab/cfDNA 
  Documentation:  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4715266/ 
  STANLEY PAPER AND SNAKEMAKE WORKFLOW
  Github:   https://github.com/kircherlab/cfDNA?tab=readme-ov-file 
  Documentation:  https://www.nature.com/articles/s41467-024-46435-0 

3. **Analysis**
   * **ranking-tissue-of-origin.R:** R script containing functions for plotting the ranking across different time points and see the evolution of tissue cell signatures. 
     
4. **Workflow (scripts)**
   * **fft_table.py**
   * **snakemake_correlation_rank_table.R**
   * **nakemake_correlation_table.R**

    **NOTE:** These scripts were modified from their original version. You can find the remaining in the [Snakemake workflow: Analysis of epigenetic signals captured by fragmentation patterns of cell-free DNA].
   
6. **Reference**
   * **healthy_Song_100bp_hg38.bed:** It was used as **regions.tsv**. Bed file containing regions of interest.


## Usage
* **snakefile_WPS.smk:**

  **Input**
  - configured by the user ([samples.tsv](config/samples.tsv)):
      - analysis ID
      - samples
      - path to sample .bam files
      - reference samples fro plotting
      - genome build per sample
  - configured by the user ([regions.tsv](config/regions.tsv)):
      - bed file containing regions of interest (e.g. TFBS), all having the same length
  **Output**
  - table containing bp specific WPS for regions listed in bed
  - line plot showing normalized WPS of multiple samples

  
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
      - analysis ID
      - samples
      - path to sample BAM files
      - genome build per sample
  
  **Output**
  
  - fft_summary tables (results/intermediate/body/fft_summaries)
  - plots showing intensities across tissues (results/plots)
  - table showing correlation with tissues/cell lines
  - table showing correlation rank difference to reference sample


---

* **ranking-tissue-of-origin.R:** 

  **Input:** *_correlation_results.tsv from the snakefile_GE_analysis.snk workflow and sample order list.

  **Output:** Plot with ranked correlationss for each time point. 

  **Libraries:**


   ```R
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(RColorBrewer)
  library(paletteer)
  
  ```
  
  **Installation:**





