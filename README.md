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
> Liquid biopsies, leveraging cell-free DNA (cfDNA), have emerged as a valuable non-invasive source for cancer biomarkers for detection, diagnosis, and disease monitoring. These circulating DNA fragments, shed into the bloodstream through biological processes such as apoptosis and necrosis, harbor a rich reservoir of molecular information. Plasma cfDNA reflects tumor type, stage, and treatment response, offering insights into fragmentomics, genomics, and epigenetics. However, most current cfDNA analyses rely on short-read sequencing, which can introduce biases during library preparation and sequencing.  
> 
> Oxford Nanopore Technologies (ONT) presents a promising long-read alternative, capable of preserving native cfDNA profiles and enabling rapid analysis — although its comprehensive evaluation in oncology is still emerging.
>
> **Results:**  
> This study demonstrates the feasibility of using **shallow whole-genome sequencing (sWGS)** with ONT to analyze cfDNA fragment size profiles, detect copy number alterations (CNAs), and infer tissue-of-origin signals. Additionally, by leveraging ONT’s real-time methylation calling capabilities, we show its potential to extract methylation data without additional library preparation or dedicated protocols.  
> 
> Our findings highlight ONT’s enhanced capacity for cfDNA signal characterization — capturing a broader spectrum of fragment lengths for enriched fragmentomics, achieving a strong correlation for tumor fraction estimates (R = 0.939), highly concordant CNA profiles with Illumina, and increased sensitivity in tissue-of-origin detection. Furthermore, its ability to perform direct methylation analysis underscores ONT’s potential as a versatile, multifaceted tool for liquid biopsy-based cancer analysis.

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




