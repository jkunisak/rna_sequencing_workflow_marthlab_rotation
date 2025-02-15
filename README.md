# RNA Sequencing Pipeline for Bulk and Single Cell Data (Marth Lab Rotation)

This is a pipeline that can use raw fastq files form either cohort-level bulk or single cell rna sequencing. The outputs will include a gene-level differential expression matrix ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)) as well as a list of differentially activated pathways identified through a Gene Set Enrichment Analysis ([GSEA](http://software.broadinstitute.org/gsea/)). The pipeline also includes steps for removal of adapter features ([trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), demultiplexing in the case of single cell sequencing, genomic/transcriptomic alignment ([Hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)/[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) respectively), filtering of mapped reads ([samtools](http://www.htslib.org/doc/samtools.html)), and normalization ([Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)/[NOISeq](https://www.bioconductor.org/packages/devel/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf)/[edgeR](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)). Visualizations will be produced to illustrate data quality, verification of proper normalization, differential expression, and GSEA.

The basic workflow of the pipeline is depicted below:
![alt text](https://github.com/jkunisak/rna_sequencing_workflow_marthlab_rotation/blob/master/Fig%201%20-%20Workflow.png)

## Getting Started
### Installing Necessary Software
To get started (assuming the user works in the Marth lab), a few tools need to be installed. The pipeline requires that the tools be installed in a `software` folder that exists in the home directory (`~`). The `fastqc`, `hisat2`, and `samtools`modules will be used and loaded automatically when running the pipeline.

The first tool that needs to be installed is `trimmomatic` to remove adapter sequences from single cell rna sequencing data. To install `trimmomatic`, please use the following command:

First make the software directory if it doesn't exist:
```
mkdir ~/software
```
Next, install `trimmomatic` (v0.38) in the `software` directory. Use the following commands to download the program from source and unzip the file:
```
cd ~/software
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
```

The second tool that needs to be installed is the most up-to-date version of `salmon`. Please refer to https://salmon.readthedocs.io/en/latest/building.html#installation for more specific installation instructions.
```
cd ~/software
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.11.2/salmon-0.11.2-linux_x86_64.tar.gz
tar -xvf salmon-0.11.2-linux_x86_64.tar.gz
cd salmon_reference
mkdir build
cd build
make
make install
```

#### Defining the R library
Follow the instructions provided [here](https://clas.uiowa.edu/linux/help/applications/rpackage) to create a R library folder. Once you create the `Rlibs` folder, copy the contents of the `Rlibs` directory provided on this GitHub page to that newly created directory.
```
cp $github_Rlibs_directory $new_Rlibs_directory
```
A list of the necessary packages include:
1) [data.table](https://cran.r-project.org/web/packages/data.table/data.table.pdf)
2) [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) - differential expression
3) [stringr](https://cran.r-project.org/web/packages/stringr/stringr.pdf)
4) [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) - GSEA
5) [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) - normalization, differential expression
6) [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
7) [GSEABase](https://bioconductor.org/packages/release/bioc/html/GSEABase.html) - GSEA
8) [GSVA](https://bioconductor.org/packages/release/bioc/html/GSVA.html) - GSEA, differential expression
9) [genefilter](https://bioconductor.org/packages/release/bioc/html/genefilter.html) - normalization
10) [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) - GSEA, differential expression
11) [MSigDB](https://github.com/oganm/MSigDB) ==> requires [dplyr](https://cran.r-project.org/web/packages/dplyr/dplyr.pdf) to install - GSEA
12) [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
13) [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
14) [NOISeq](https://bioconductor.org/packages/release/bioc/html/NOISeq.html) - normalization

### Establish Directory Structure
The next step is to orgnaize the directories. Within the home directory (`$homeDIR`), there must exist a parent directory that contains a directory of fastq files, an output directory, a directory for the reference files, and a scripts directory that houses the scripts used in the pipeline.
1) Copy the raw fastq files into the `fastq_files` directory (see below)
2) Subdirectories in the `output` folder will be made automatically by the pipeline
3) Files in the `reference` directory include the reference genome, transcriptome, and index files for alignment. Please see commands below to obtain and index these files.
4) `scripts` directory contains all of the R scripts that the pipeline uses. Get these files from the GitHub repo.

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
mkdir $homeDIR/reference/hisat2
hisat2_extract_splice_sites.py $homeDIR/reference/GRCh37_gtf/Homo_sapiens.GRCh37.75.gtf > $homeDIR/reference/hisat2/GRCh37.splicesite
hisat2_extract_exons.py $homeDIR/reference/GRCh37_gtf/Homo_sapiens.GRCh37.75.gtf > $homeDIR/reference/hisat2/GRCh37.exon
hisat2-build --ss $homeDIR/reference/hisat2/GRCh37.splicesite --exon $homeDIR/reference/hisat2/GRCh37.exon ~/scRNA-seq_Expression_Analysis/reference/hg19_fasta/hg19.fa $homeDIR/reference/hisat2/test

## Build the salmon index file
mkdir $homeDIR/reference/salmon_reference
~/software/salmon-0.11.1-linux_x86_64/bin/salmon index -t $homeDIR/reference/cdna_fasta/Homo_sapiens.GRCh37.75.cdna.all.fa -i $homeDIR/reference/salmon_reference/salmon_transcript_index --type quasi -k 31
```
Directory Structure: <br />
~/ <br />
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

### Required Input Parameters
1) `homeDIR`: defines the main directory that contains the `fastq_files`, `output`, `reference`, and `scripts` folders. The parent directory to the `homeDIR` should contain the `software` folder.
2) `trimmomaticFastaPath`: defines the fasta file that `trimmomatic` uses to remove adapter sequences (if necessary). Defaults to NULL. Options for this parameter include the fasta files provided in the `trimmomatic/adapters/` directory within the `software` folder. The available files (n=6) are listed below:
```
ls ~/software/Trimmomatic-0.38/adapters
NexteraPE-PE.fa  TruSeq2-PE.fa  TruSeq2-SE.fa  TruSeq3-PE-2.fa  TruSeq3-PE.fa  TruSeq3-SE.fa
```

## Running the Tests
### Example command:
```
bash ~/$homeDIR/scripts/draft1.sh $homeDIR NULL
```

## Useful Papers and Resources to Better Understand RNA Sequencing
https://hemberg-lab.github.io/scRNA.seq.course/tabula-muris.html <br />
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190152 <br />
https://github.com/griffithlab/rnaseq_tutorial <br />
https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/ <br /> 
