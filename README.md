<img align="right" width="100" src="img/esci-logo.png" />

# MULTI-DIMENSIONAL COMPARISON OF CELL-FREE DNA MOLECULAR SIGNATURES USING SHORT AND LONG READ SEQUENCING TECHNOLOGIES IN CLINICAL APPLICATIONS

![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white) 
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) 
![Bash Script](https://img.shields.io/badge/bash_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white) 
![Anaconda](https://img.shields.io/badge/Anaconda-%2344A833.svg?style=for-the-badge&logo=anaconda&logoColor=white) 
![Matplotlib](https://img.shields.io/badge/Matplotlib-%23ffffff.svg?style=for-the-badge&logo=Matplotlib&logoColor=black) 
![Nextflow](https://img.shields.io/badge/Nextflow-%23007ACC.svg?style=for-the-badge&logo=nextflow&logoColor=white) 
![Snakemake](https://img.shields.io/badge/Snakemake-%23E6943B.svg?style=for-the-badge&logo=snakemake&logoColor=white)

---

## Overview

This repository contains the code, workflows, and data processing pipelines developed for my thesis project, which focuses on the **multi-dimensional characterization of cell-free DNA (cfDNA)** using both **short-read and long-read sequencing technologies** in clinical oncology applications.

---

## Abstract

> **Motivation:**  
> Liquid biopsies, leveraging circulating tumor DNA (ctDNA), have emerged as a valuable non-invasive source for cancer biomarkers for detection, diagnosis, and monitoring in patients. These circulating DNA fragments, shed into the bloodstream and other bodily fluids by processes such as apoptosis and necrosis, potentially contain a vast reservoir of molecular information. The plasma cfDNA encapsulates molecular signatures reflecting tumor type, stage, and treatment response, providing insights into fragmentomics, genomics, and epigenetics. Current cfDNA analyses predominantly rely on short-read sequencing technologies, which can introduce biases during library preparation and sequencing. Long-read sequencing, such as
>
> Oxford Nanopore Technology (ONT), offers a promising alternative, enabling rapid analysis and preserving native cfDNA fragment profiles, through its comprehensive in oncology is still emerging.  

> **Results:**  
> This study compared the results of short-read and long-read sequencing from 23 plasma cfDNA samples and demonstrates the feasibility of using shallow whole-genome sequencing (sWGS) with ONT to characterize cfDNA fragment size profiles, detect copy number alterations, and identify tissue-of-origin signals. Furthermore, by leveraging ONT’s real-time methylation calling capabilities, this work demonstrates the potential to additionally extract methylation data without requiring additional library preparation steps or dedicated sequencing protocols. The findings demonstrate that ONT (long-read sequencing) facilitates an enhanced characterization of these cfDNA-derived signals, highlighted by capturing a fuller spectrum of cfDNA fragment lengths for enriched fragmentomic analysis, showing strong positive correlation for tumor fraction (TF) estimates (R = 0.939) and highly concordant CNA profiles with Illumina (short-read sequencing), and increased sensitivity in detecting the tissue-of-origin in the tested samples. We also demonstrate its capability to perform methylation analysis. Thereby, it underscores its potential as a powerful and versatile tool for multifaceted cancer analysis.

---

## Repository Structure

1. **Data processing:** Scripts and pipelines for preprocessing raw data (FASTQ for Illumina, uBAM for Nanopore), downsampling, and other utilities.
3. **Fragment analysis:** Scripts for obtaining fragment size profiles and generating fragmentomics comparison plots.
4. **Copy number analysis:** Code for performing copy number alteration analyses and generating related plots.
5. **Tissue of origin:** Adapted workflows for tissue-of-origin analysis and scripts for figure generation.
6. **Methylation:** Scripts for methylation extraction from BAM files and downstream visualization.
7. **Supplementary Material:** PDF files containing supplementary material referenced in the manuscript.

---

## Contact

- **Author:** Alessandra Bonilla Salon  
- **Email:** [alessandra.bonilla@alum.esci.upf.edu](mailto:alessandra.bonilla@alum.esci.upf.edu)  
- **Institution:** ESCI-UPF  
- **Year:** 2024–2025  




