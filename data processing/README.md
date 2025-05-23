# Data Processing
This folder contains all scripts and workflows used for data preprocessing and manipulation for both Illumina and Nanopore sequencing data.

## Contents

1. **Illumina pre-processing:** This section includes a [Nextflow](https://www.nextflow.io/docs/latest/install.html) workflow for preprocessing paired-end Illumina FASTQ files (R1 and R2).
   **Workflow Steps:**
      * Adapter trimming (Cutadapt)
      * Alignment to reference genome (BWA)
      * BAM indexing (Samtools)
      * Marking duplicates (Picard)
      * Wig file creation (hmm_readCounter)
      * Copy number analysis (ichorCNA)


   **Files:**
      * **illumina_preprocessing-and-ichor.nf:**  Main Nextflow workflow.
      * **nextflow.config:** Configuration file to specify paths for reference genomes, input/output directories, and parameter settings.

3. **Nanopore pre-processing:** This section contains a [Nextflow](https://www.nextflow.io/docs/latest/install.html) workflow for preprocessing Nanopore long-read (single-end) data.

   **Workflow Steps:**
     * Alignment to the reference genome (epi2me-labs/wf-alignment)
     * Wig file generation (hmm_readCounter)
     * Copy number analysis (ichorCNA)
     * Summarizing tumor fraction and other ichorCNA outputs into a single results table.

   **Files:**
     * **nanopore_preprocessing-and-ichor.nf:** Main Nextflow workflow.
     * **nextflow.config:** Configuration file for paths to inputs, outputs, reference genomes, and final summary table.
      
5. **Utility Scripts:**
    * **calculate_coverage.sh:** Calculates coverage for all BAM files in a folder in parallel, producing a text file with sample names and respective coverage values.

    * **illumina-reads_downsampling_for_folder.sh:** Downsamples Illumina FASTQ paired-end reads.

        **Steps:**

         * Calculate initial raw coverage.
         * Create a list of target coverages (from initial to 0.1x in 0.1x steps, then from 0.1x to 0.01x in 0.01x steps).
         * Perform parallel downsampling using seqtk.

    * **long-reads-concatenate_and_downsampling.sh:** 

       For Nanopore data:
        * long-reads-concatenate_and_downsampling.sh – For Nanopore data:
        * Concatenates all FASTQ files within each barcode directory.
        * Calculates initial raw coverage.
        * Creates a list of target coverages (as above).
        * Performs parallel downsampling with seqtk.

## Usage
* **illumina_preprocessing-and-ichor.nf:**
  
  **Input:**
   * FASTQ paired-end reads folder
   * Edited nextflow.config specifying input/output paths and reference genome (.fa, .fai, .dict)

  **Output:**
   * Processed trimmed FASTQ, BAM, wig, ichorCNA results.

  **Requirements:**
  * Nextflow v22.10.0.5826
  * Docker v26.1.3
  * Docker image
    
  **Installation:**
  
  **Command:**
  
  ```bash
  nextflow run illumina_preprocessing-and-ichor.nf
  ```

***
* **nanopore_preprocessing-and-ichor.nf:**
  
  **Input:**
  + Folder with uBAM files organized by barcode
  * Edited nextflow.config
 
  **Output:**
  * Processed BAM, wig, ichorCNA results, and final summary table

  **Requirements:**
  * Nextflow v22.10.0.5826 or higher
  * Install epi2me-labs/wf-alignment, hmm_readCounter, ichorCNA.

  **Installation:**
  
  **Command:**
  
  ```bash
  nextflow run illumina_nanopore_preprocessing-and-ichor.nf
  ```

***
* **calculate_coverage.sh:**
  
  **Input:**
  * Folder with BAM files
 
  **Output:**
   * Text file with sample names and coverage values


  **Requirements:**
  * samtools
    
  **Installation:**

  ```bash
  conda install -c bioconda samtools
  ```

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




