# rna_sequencing_workflow_marthlab_rotation

Dr. Gabor Marth Lab Rotation - Proposed RNA Sequencing Analysis Pipeline
This is a pipeline that can use raw fastq files form either bulk or single cell rna sequencing. The outputs will include a gene-level differential expression matrix (DESeq2) as well as a list of differentially activated pathways identified through a Gene Set Enrichment Analysis (GSEA). The pipeline also includes steps for removal of adapter features (trimmomatic), demultiplexing in the case of single cell sequencing, genomic/transcriptomic alignment (Hisat2/Salmon respectively), filtering of mapped reads (samtools), and normalization (salmon/NOISeq/edgeR).

## Getting Started

To get started (assuming the user if working in the Marth lab), a few tools need to be installed. The pipeline requires that the tools be installed in a "software" that exists in the home directory. The fastqc, hisat2, samtools, and salmon modules will be used and loaded automatically when running the pipeline.

The first tool that needs to be installed in "trimmomatic". This is needed when analyzing data from single cell RNA sequencing to remove adapter sequences. To install trimmomatic, please use the following command:
 



Required input parameters include:
1) Raw fastq files
2)

Additional options:
1)

1) Quality check on fastq files <br />
	a) fastqc <br />
2) Trim/demultiplex fastq files <br />
	a) trimmomatic <br />
3) Align reads to reference genome/transcriptome <br />
	a) hisat2 - genome <br />
	b) salmon - transcriptome (does normalization and quantification) <br />
4) Normalize and quantify counts <br />
	a) NOISeq - normalize <br />
	b) alpine - gc content normalization <br />
	c) edgeR - rna composition/batch effects <br />
5) DE analysis and GSEA <br />
	a) DESeq2/edgeR/limma+voom/NOISeq <br />
	b) GSVA <br />

Directory Structure: <br />
$HOME <br />
|-- software (extracted files) <br />
|-------|-- Trimmomatic-0.38 <br />
|-------|-- salmon-0.11.1-linux_x86_64 <br />
|-- scRNA-seq_Expression_Analysis (Set $homeDIR to this directory) <br />
|-------|-- fastq_files <br />
|-------|-------|-- This directory contains all of the raw fastq files with ${sampleName}_1.fq.gz format <br />
|-------|-- output <br />
|-------|-------|-- fastqc_output <br />
|-------|-------|-- trimmed_fastq_files <br />
|-------|-------|-- trimmed_fastqc <br />
|-------|-------|-- hisat2_alignment <br />
|-------|-------|-- hisat2_alignment_coverage_metrics <br />
|-------|-------|-- samtools_flagstat <br />
|-------|-------|-- featureCounts <br />
|-------|-------|-- salmon <br />
|-------|-------|-- DESeq2 <br />
|-------|-- reference <br />
|-------|-------|-- cdna_fasta <br />
|-------|-------|-- GRCH37_gtf <br />
|-------|-------|-- hg19_fasta <br />
|-------|-------|-- hisat2 <br />
|-------|-------|-- salmon_reference <br />
|-------|-- scripts <br />
|-------|-------|-- adapter_detection.R <br />
|-------|-------|-- create_batch_file.R <br />
|-------|-------|-- Rsubread_featureCounts.R <br />
|-------|-------|-- NOISeq.R <br />
|-------|-------|-- Rlibs (R library directory with all of the installed packages) <br />
