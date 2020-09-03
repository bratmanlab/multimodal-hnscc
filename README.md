# Multimodal Analysis of HNSCC Plasma - Data for Manuscript
*Author: Justin M. Burgener*

Processed data and code used to generate figures for manuscript titled: 'Tumor-naive multimodal profiling of circulating tumor DNA in localized head and neck squamous cell carcinoma'.

## *About*

The code within the bratman lab/multimodal-hnscc-manuscript repository is used to demonstrate key findings including: (1) the identification of ctDNA-derived mutations and methylation by CAPP-Seq and cfMeDIP-seq respectively, (2) associated properties of ctDNA fragments, and (3) clinical utility of multimodal detection of ctDNA. The code uses proccessed data from BAM files (FASTQ to SAM by bwa-mem, SAM to BAM by samtools), requiring bash programming and pre-installed software. Access to a high-performance computing cluster is also recommended. Due to privacy concerns, the BAM files used are not available (with the exception of the HNSCC cell-line FaDu which is available upon request).

## *Requirements:*
  * Computer running a Linux system (≥ 8 GB RAM) **Cluster computing is HIGHLY recommended when working with FASTQ/BAM files**
    * Modules: bwa (version 0.7.15), bowtie2 (version 2.2.6), samtools (version 1.3.1), igenome-human/hg19
  * R/RStudio (version 3.5 or greater)
    * Packages: BSgenome.hsapiens.UCSC.hg19, GenomicRanges, AnnotationHub, Repitools
    
## *Procedure:*

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
