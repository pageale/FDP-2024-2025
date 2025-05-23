# Data Processing
This folder contains all the scripts and workflows used for data preprocessing and other data manipulation.

## Contents

1. **Illumina pre-processing:** In this folder you will find the Nextflow workflow por data pre-preprocessing for illumina paired-end reads fastq (R1 and R2). It does the trimming (CutAdapt), alignment (bwa), indexing (samtools), marking duplicates (Picard), creating the wig files hmm_readCounter), and finally, the copy number analysis (ichorCNA). It receives as input a fastq folder and gives as output each file corresponding to the processes mentioned previously.
    * **illumina_preprocessing-and-ichor.nf:** Workflow with all the steps.
    * **nextflow.config:** Config file where you add the reference genome (.fa), point to input and output folder, and other parameters tuning needed for the workflow.
      
3. **Nanopore pre-processing:** In this folder, you will find the Nextflow workflow for data pre-processing for nanopore long reads (single end). The steps include alignment to the reference genome (using a epi2me-labs workflow), creating the wig files (hmm_readCounter), performing the copy number analysis (ichorCNA), and finally it summarize the tumor fraction (and other parameters) from all samples into a table for easy accesibility. 
    * **nanopore_preprocessing-and-ichor.nf:** This is the nextflow workflow that receives as input a folder with uBAM (direct output from PromethION sequencer) organized into different folders by barcode. As output, it gives the corresponding file from each process. 
    * **nextflow.config:** Config file that contains the input and output path, the name of the final table, path to reference genome, and other parameters tuning for the correct use of the workflow.
      
5. **Other scripts:**
    * **calculate_coverage.sh:** THis sctipt calculates the coverage of the BAM files. It looks for all the .bam extension and then calculate the coverage in parallel for multiples samples. 
    * **illumina-reads_downsampling_for_folder.sh:** This script perfroms downsampling from the fastq (R1 and R2). FIrst it computes sit initial raw coverage (no pre-processing), then it creates a list with differnet cov values (from inital to 0.1x in steps of 0.1 and then from 0.1x to 0.01x in steps of 0.01), and finally, it does the subsampling for both reads using seqtk. 
    * **long-reads-concatenate_and_downsampling.sh:** This performs downsampling for long-read data. It rewceieves as input a folder with one directory for each barcode of fastq (direct output from the sequencer). First, it merges all the fastq into one per barcode, then it calculates the initial raw coverage, create the list with differnet cov values (from inital to 0.1x in steps of 0.1 and then from 0.1x to 0.01x in steps of 0.01), and finally, it does the subsampling in parallel. 


## Usage
* **illumina_preprocessing-and-ichor.nf:**
  
  **Input:**

  **Output:**


  **Requirements:**
  * Edit config file: write input, output, and path to local reference genome (.fa, .fai, .dict).
  * nextflow version 22.10.0.5826
  * Docker version 26.1.3
  * Docker image
    
  **Installation:**
  
  **Command:**
  
  ```bash
  nextflow run illumina_preprocessing-and-ichor.nf
  ```
  
***
* **nanopore_preprocessing-and-ichor.nf:**
  
  **Input:**
  + config file
 
  **Output:**
  * all the files correpsondung to each step (.bam, .bai, .wig, ichor folder)

  **Requirements:**
  * Edit config file: write input, output, name of the table and path to local reference genome (.fa, .fai, .dict).
  * nextflow version 22.10.0.5826 or above. 
  * Install epi2me-labs/wf-alignment, hmm_readCounter, ichorcna

  **Installation:**
  
  **Command:**
  
  ```bash
  nextflow run illumina_nanopore_preprocessing-and-ichor.nf
  ```
  
<fastq_folder>

***
* **calculate_coverage.sh:**
  
  **Input:**
  * BAM files
 
  **Output:**
   * a .txt with all the samples names and its respective coverage

  **Requirements:**
  * samtools
    
  **Installation:**
  
  **Command:**
  
  ```bash
  chmod +x calculate_coverage.sh
  ./calculate_coverage.sh
  ```

***
* **illumina-reads_downsampling_for_folder.sh:**
  
  **Input:**
  * Folder with FASTQ paired-end
 
  **Output:**
   * **illumina_fastq_downsampling-2.txt:** Log file. 
   * **illumina_fastq_downsampling/:** Folder containig all the FASTQ downsmpled at diferent coverages.

  **Requirements:**
  * seqtk
  * GNU parallel
    
  **Installation:**
   ```bash
    conda install bioconda::seqtk
   ```

  **Command:**
  
  ```bash
  chmod +x illumina-reads_downsampling_for_folder.sh
  ./illumina-reads_downsampling_for_folder.sh <fastq_folder>"
  ```

***

* **long-reads-concatenate_and_downsampling.sh:**
  
  **Input:**
  * Folder with barcoded FASTQ divided into their respective folder (output from sequencer). 
 
  **Output:**
   * **downsampling_nanopore_picard_log.txt:** Log file. 
   * **downsampling_nanopore_second_run/:** Output directory containig all the FASTQ downsmpled at diferent coverages.

  **Requirements:**
  * seqtk
  * GNU parallel
   
  **Installation:**
   ```bash
    conda install bioconda::seqtk
   ```

  **Command:**
  
  ```bash
  chmod +x long-reads-concatenate_and_downsampling.sh
  ./long-reads-concatenate_and_downsampling.sh <fastq_folder>"
  ```




