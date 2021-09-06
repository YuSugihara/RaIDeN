# RaIDeN User Guide
#### version 0.0.2 (beta)


## Table of contents
- [What is RaIDeN?](#What-is-RaIDeN)
  + [Citation](#Citation)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [How to install](#How-to-install)
- [Usage](#Usage)
  + [What RaIDen requires](#What-RaIDen-requires)
  + [Example 1 : run RaIDeN from FASTQ without trimming](#Example-1--run-RaIDeN-from-FASTQ-without-trimming)
  + [Example 2 : run jiji from VCF](#Example-5--run-QTL-plot-from-VCF)
- [Outputs](#Outputs)



## What is RaIDeN?
RaIDeN is a pipeline for <ins>**Ra**</ins>pid <ins>**IDe**</ins>ntification of causal gene with target motif using <ins>**N**</ins>GS technology.

RaIDeN automatically filter the candidate genes based on the patterns of structual variations and mutations.

<img src="https://github.com/YuSugihara/RaIDeN/blob/master/images/Fig.S1.png" width=700>

Supplementary figure 1 in [Shimizu et al. (2021)](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1)

#### Citation
Shimizu M, Hirabuchi A, Sugihara Y, Abe A, Takeda T, Kobayashi M, Hiraka Y, Kanzaki E, Oikawa K, Saitoh H, Langner T, Banfield MJ, Kamoun S, Terauchi R (2021) [A genetically linked pair of NLR immune receptors show contrasting patterns of evolution](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1). _bioRxiv_


## Installation

### Dependencies
- samtools
- bcftools
- hisat2
- bedtools
- gffread
- bamtools
- pigz
- stringtie
- faqcs
- prinseq-lite

### How to install
```
git clone https://github.com/YuSugihara/RaIDeN.git
cd RaIDeN
pip install .
```

### Installation via bioconda

coming soon....


## Usage

```
usage: raiden -r <FASTA>
              -a <FASTQ1>,<FASTQ2> 
              -w <FASTQ1>,<FASTQ2>
              -o <OUT_DIR> 
             [-t <INT>]

RaIDeN version 0.0.2

optional arguments:
  -h, --help            show this help message and exit
  -r , --ref            Reference fasta. You may use De novo assembled contigs.
  -a , --rna-seq        Fastqs of RNA-seq to annotate the genes on the reference genome.
                        This RNA-seq must contain the reads of the causal gene.
                        If you have separated fastq files, You can use this optiion
                        for each pair of files. Please separate paired fastqs by
                        comma (e.g. -a fastq1,fastq2).
  -w , --whole-genome   Fastqs of whole-genome sequences to select causal genes.
                        Those samples must contain the opposite traits to the sample
                        used in the refernce genome. Please separate paired fastqs by
                        comma (e.g. -a fastq1,fastq2). You can use this optiion
                        multiple times for each sample.
  -o , --out            Output directory. Specified name must not exist.
  -t , --threads        Number of threads. If you specify the number below one,
                        then RaIDeN will use the threads as many as possible. [2]
  -d, --disable-RNAseq-trim
                        Disable the trimming of RNA-seq by FaQCs.
  -D, --disable-WGS-trim
                        Disable the trimming of whole-genome sequences
                        by FaQCs and prinseq-lite.
  -s , --strand         Assume a strand library.
                         'fr' : assume a strand library fr-firststrand.
                         'rf' : assume a strand library fr-secondstrand.
                         'None' : Don't assume a strand library.
                        Default is 'None'.
  -m , --minimum-len    Minimum length allowed for the predicted transcripts in
                        stringtie. [200]
  -q , --min-MQ         Minimum mapping quality in mpileup. [40]
  -Q , --min-BQ         Minimum base quality in mpileup. [18]
  -v, --version         show program's version number and exit
```


### What RaIDen requires:

- #### Reference genome (FASTA format)

    RaIDeN requires an assembled reference genome. This reference genome must contain a causal gene. If it is sure that the assembled reference genome contains the causal gene, RaIDeN allows contiguous (not chromosome-scale) referece genome.

- #### RNA-seq for genome annotation (paired FASTQ format)

    RaIDeN requires RNA-seq for genome annotation. This RNA-seq must contain the causal gene because RaIDen analyzes only on the annotated genome region. RaIDeN allows multiple RNA-seq samples.

- #### Whole-genome sequences (paired FASTQ format)

    RaIDen requires whole-genome sequence (WGS). The sample of this WGS should have an opposite phenotype to that of the reference genome. Since RaIDeN expects that the different phenotypes come from structual variation or mutation on the annotated genome region, this WGS must have a structual variation or mutation on the causal gene. RaIDeN allows multiple WGSs or bulked sequences.


