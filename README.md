# rna_sequencing_workflow_marthlab_rotation
Dr. Gabor Marth Lab Rotation - Proposed RNA Sequencing Analysis Pipeline

This is a proposed analysis pipeline to perform gene-level differential expression (DE) and gene set enrichment analysis (GSEA) on bulk or single-cell RNA sequencing data. 

1) Quality check on fastq files
	a) fastqc
2) Trim/demultiplex fastq files
	a) trimmomatic
3) Align reads to reference genome/transcriptome
	a) hisat2 - genome
	b) salmon - transcriptome (does normalization and quantification) 
4) Normalize and quantify counts
	a) NOISeq - normalize
	b) alpine - gc content normalization
	c) edgeR - rna composition/batch effects
5) DE analysis and GSEA
	a) DESeq2/edgeR/limma+voom/NOISeq
	b) GSVA

Directory Structure: 
$HOME
|-- software # (extracted files)
|	|-- Trimmomatic-0.38
|	|-- salmon-0.11.1-linux_x86_64
|-- scRNA-seq_Expression_Analysis (Set $homeDIR to this directory)
|	|-- fastq_files
|	|	|-- This directory contains all of the raw fastq files with ${sampleName}_1.fq.gz format
|	|-- output
|	|	|-- fastqc_output
|	|	|-- trimmed_fastq_files
|       |       |-- trimmed_fastqc
|       |       |-- hisat2_alignment
|       |       |-- hisat2_alignment_coverage_metrics
|       |       |-- samtools_flagstat
|       |       |-- featureCounts
|       |       |-- salmon
|       |       |-- DESeq2
|	|-- reference
|       |       |-- cdna_fasta
|       |       |-- GRCH37_gtf
|       |       |-- hg19_fasta
|       |       |-- hisat2
|       |       |-- salmon_reference
|	|-- scripts
|       |       |-- adapter_detection.R
|       |       |-- create_batch_file.R
|       |       |-- Rsubread_featureCounts.R
|       |       |-- NOISeq.R
|       |       |-- Rlibs (R library directory with all of the installed packages)









