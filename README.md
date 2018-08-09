# rna_sequencing_workflow_marthlab_rotation

Dr. Gabor Marth Lab Rotation - Proposed RNA Sequencing Analysis Pipeline
This is a pipeline that can use raw fastq files form either bulk or single cell rna sequencing. The outputs will include a gene-level differential expression matrix (`DESeq2`) as well as a list of differentially activated pathways identified through a Gene Set Enrichment Analysis (`GSEA`). The pipeline also includes steps for removal of adapter features ([trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), demultiplexing in the case of single cell sequencing, genomic/transcriptomic alignment ([Hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)/[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) respectively), filtering of mapped reads ([samtools](http://www.htslib.org/doc/samtools.html)), and normalization ([Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)/[NOISeq](https://www.bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)/[edgeR](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)).

## Getting Started
### Installing Necessary Software
To get started (assuming the user works in the Marth lab), a few tools need to be installed. The pipeline requires that the tools be installed in a "software" folder that exists in the home directory. The `fastqc`, `hisat2`, `samtools`, and `salmon` modules will be used and loaded automatically when running the pipeline.

The first tool that needs to be installed is `trimmomatic`. This is needed when analyzing data from single cell RNA sequencing to remove adapter sequences. To install `trimmomatic`, please use the following command:

First make the software directory if it doesn't exist:
```
mkdir ~/software
cd ~/software
```
Next, install `trimmomatic` (v0.38). Use the following commands to download the program from source and unzip the file:
```
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
```

The second tool that needs to be installed is the most up-to-date version of `salmon`. Please refer to https://salmon.readthedocs.io/en/latest/building.html#installation for more specific installation instructions.
```
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.11.2/salmon-0.11.2-linux_x86_64.tar.gz
tar -xvf salmon-0.11.2-linux_x86_64.tar.gz
cd salmon_reference
mkdir build
cd build
make
make install
```

### Establish Directory Structure
The next step is to orgnaize the directories. Within the home directory (`$homeDIR`), there must exist a parent directory that contains a directory of fastq files, an output directory, a directory for the reference files, and a scripts directory that houses the scripts used in the pipeline.
1) Copy the raw fastq files into the `fastq_files` directory (see below)
2) Subdirectories in the `output` folder will be made automatically by the pipeline
3) Files in the `reference` directory include the reference genome, transcriptome, and index files for alignment. Please see commands below to obtain and index these files.
4) `Scripts` directory contains all of the R scripts that the pipeline uses. Get these files from the GitHub repo.

See the commands below to make and supply the necessary files to the directory:
```
mkdir ~/$homeDIR/fastq_files
mkdir ~/$homeDIR/output
mkdir ~/$homeDIR/reference
mkdir ~/$homeDIR/scripts

## Get the reference genome file (build 37)
mkdir ~/$homeDIR/reference/hg19_fasta
cd ~/$homeDIR/reference/hg19_fasta
echo {1..22} | parallel --eta -j+0 'rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr{}.fa.gz ~/$homeDIR/reference/hg19_fasta/'
cat
 ~/$homeDIR/reference/hg19_fasta/*.fa.gz > ~/$homeDIR/reference/hg19_fasta/hg19.fa
rm -rf ~/$homeDIR/reference/hg19_fasta/*.fa.gz

## Index the reference genome file
samtools faidx ~/$homeDIR/reference/hg19_fasta/hg19.fa

## Get the reference transcriptome file (build 37: gtf and cdna fasta)
mkdir ~/$homeDIR/reference/cdna_fasta
cd ~/$homeDIR/reference/cdna_fasta
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-75/fasta/homo_sapiens/cdna/*
gunzip *.gz

cd ~/$homeDIR/reference/GRCH37_gtf
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

## Build the hisat2 index file
module load hisat2
hisat2_extract_splice_sites.py ~/scRNA-seq_Expression_Analysis/reference/GRCh37_gtf/Homo_sapiens.GRCh37.75.gtf > ~/scRNA-seq_Expression_Analysis/reference/hisat2/GRCh37.splicesite
hisat2_extract_exons.py ~/scRNA-seq_Expression_Analysis/reference/GRCh37_gtf/Homo_sapiens.GRCh37.75.gtf > ~/scRNA-seq_Expression_Analysis/reference/hisat2/GRCh37.exon
hisat2-build --ss ~/scRNA-seq_Expression_Analysis/reference/hisat2/GRCh37.splicesite --exon ~/scRNA-seq_Expression_Analysis/reference/hisat2/GRCh37.exon ~/scRNA-seq_Expression_Analysis/reference/hg19_fasta/hg19.fa test

## Build the salmon index file
~/software/salmon-0.11.1-linux_x86_64/bin/salmon index -t $homeDIR/reference/cdna_fasta/Homo_sapiens.GRCh37.75.cdna.all.fa -i $homeDIR/reference/salmon_reference/salmon_transcript_index --type quasi -k 31
```
Directory Structure:
$homeDIR <br />
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

### Required input parameters
```
homeDIR=~/rna_analysis/
trimmomaticFastaPath=~/
```
## Running the Tests
###Example command:
```
bash ~/$homeDIR/scripts/draft1.sh $homeDIR ""
```
