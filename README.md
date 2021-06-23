# Multimodal Analysis of HNSCC Plasma - Scripts used in Analysis
*Author: Justin M. Burgener*

Processed data and code used to generate figures for manuscript titled: 'Tumor-naive multimodal profiling of circulating tumor DNA in localized head and neck squamous cell carcinoma'. Code Ocean capsule including Data and Results (referenced in the various scripts) can be found at DOI: 10.24433/CO.6511868.v1

## *About*

The code within the bratman lab/multimodal-hnscc-manuscript repository is used to demonstrate key findings including: (1) the identification of ctDNA-derived mutations and methylation by CAPP-Seq and cfMeDIP-seq respectively, (2) associated properties of ctDNA fragments, and (3) clinical utility of multimodal detection of ctDNA. The code uses proccessed data from BAM files (FASTQ to SAM by bwa-mem, SAM to BAM by samtools), requiring bash programming and pre-installed software. Access to a high-performance computing cluster is also recommended. Due to privacy concerns, the BAM files used are not available (with the exception of the HNSCC cell-line FaDu which is available upon request).

Analysis is intended to be performed via the workflow.R script, each script may be run independently however some files should be generated from prior scripts into the results folder. Please note that all data used can be located at the Code Ocean capsule referenced above.

## *Requirements:*
  * Computer running a Linux system (≥ 8 GB RAM) **Cluster computing is HIGHLY recommended when working with FASTQ/BAM files**
    * Modules: bwa (version 0.7.15), bowtie2 (version 2.2.6), samtools (version 1.3.1), igenome-human/hg19
  * R/RStudio (version 3.5 or greater)
    * Bioconductor Packages: BSgenome.hsapiens.UCSC.hg19 (v 1.4.3), GenomicRanges (v 1.42.0), AnnotationHub (v 2.22.0), Repitools (v 1.36.0), Biobase (v 2.50.0), DESeq2 (v1.30.0), Homo.sapiens (v 1.3.1), TCGAbiolinks (v 2.18.0), biovizBase (v 1.38.0), ggbio (v 1.38.0), limma (v 3.46.0), minfi (v 1.36.0)
    * R (CRAN) Packages: DT (v 0.17), RColorBrewer (v 1.1-2), cowplot (v 1.1.1), cutpointr (v 1.0.32), data.table (v 1.13.6), dplyer (v 1.0.3), egg (v 0.4.5), ggdendro (v 0.1.22), ggplot2 (v 3.3.3), ggpubr (v 0.4.0), ggridges (v 0.5.3), ggsci (v 2.9), ggthemes (v 4.2.0), glmnet (v 4.1), pROC (v 1.17.0.1), purrr (v 0.3.4), reshape2 (v 1.4.4), survival (v 3.2-7), survminer (v 0.4.8), tidyverse (v 1.3.0)
    
## *FASTQ to BAM File Procedure:*

  1. Align FASTQ files to reference genome (bwa/0.7.15, igenome-human/hg19)
     ```bash
     bwa mem -M -t4 $BWAINDEX $R1.fastq $R2.fastq > $file.sam
     ```
  2. Convert SAM file to BAM file (samtools/1.3.1)
     ```bash 
     samtools view -buS -f 2 -F 4 $file.sam | samtools sort -@4 -o $file.bam
     samtools index $file.bam
     ```
  3. Remove duplicates
     ```bash
     samtools rmdup $file.bam $bam_directory/$rmdup_file.bam
     samtools index $bam_directory/$rmdup_file.bam
     ```
  4. Convert BAM file to BED file
     ```bash
     bedtools bamtobed -i $bam_directory/$rmdup_file.bam
     ```
     
For further clarification of the bioinformatic analysis, please e-mail Justin Burgener @ justin.burgener@uhnresearch.ca
