# RaIDeN User Guide
#### version 0.0.3 (beta)


## Table of contents
- [What is RaIDeN?](#What-is-RaIDeN)
  + [Citation](#Citation)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [How to install](#How-to-install)
- [Quick Start Guide](#Quick-Start-Guide)
- [Output files](#Output-files)



## What is RaIDeN?
RaIDeN is a pipeline for <ins>**Ra**</ins>pid <ins>**IDe**</ins>ntification of causal gene with target motif using <ins>**N**</ins>GS technology. RaIDeN automatically filters the candidate genes based on the patterns of structual variations and mutations.

<img src="https://github.com/YuSugihara/RaIDeN/blob/master/images/Fig.S1.png" width=700>

Supplementary figure 1 in [Shimizu et al. (2021)](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1)

#### Citation
Shimizu M, Hirabuchi A, Sugihara Y, Abe A, Takeda T, Kobayashi M, Hiraka Y, Kanzaki E, Oikawa K, Saitoh H, Langner T, Banfield MJ, Kamoun S, Terauchi R (2021) [A genetically linked pair of NLR immune receptors show contrasting patterns of evolution](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1). _bioRxiv_


## Installation

RaIDeN is written in Python3.

### Dependencies
- Python >= 3.5
- samtools >= 1.7
- bcftools >= 1.7
- hisat2
- bedtools
- gffread
- bamtools
- pigz
- stringtie
- faqcs
- prinseq-lite

### How to install

Please run the following command lines to install RaIDeN.

```
git clone https://github.com/YuSugihara/RaIDeN.git
cd RaIDeN
pip install .
```

### Installation via bioconda

Currently, RaIDeN doesn't support the installation via bioconda. However, you can easily install its dependencies because they are distributed via bioconda. You can try the command below:

```
$ conda install -c bioconda samtools bcftools hisat2 bedtools gffread bamtools pigz stringtie faqcs prinseq
```


## Quick Start Guide

```
$ raiden -r reference.fasta \
         -a rnaseq.1.fastq,rnaseq.2.fastq \
         -w wgs.1.fastq,wgs.2.fastq \
         -o test  \
         -t 2
```

- #### `-r` : Reference genome (FASTA format)

    RaIDeN requires an assembled reference genome. This reference genome must contain a causal gene. If it is sure that the assembled reference genome contains the causal gene, RaIDeN allows contiguous (not chromosome-scale) referece genome.

- #### `-a` : RNA-seq for gene annotation (paired FASTQ format)

    RaIDeN requires RNA-seq for gene annotation. This RNA-seq must contain the sequences of causal gene because RaIDen only analyzes the annotated genomic region. Paired FASTQ files have to be separated by commna (eg. fastq1,fastq2). RaIDeN allows multiple RNA-seq samples. FASTQ files can be zipped.

- #### `-w` : Whole-genome sequences (paired FASTQ format)

    RaIDen requires whole-genome sequence (WGS). The sample of this WGS should have an opposite trait to that of the reference genome. Since RaIDeN expects that the different traits come from structual variation or mutation on the annotated genomic region, this WGS must have a structual variation or mutation on the causal gene. Paired FASTQ files have to be separated by commna (eg. fastq1,fastq2). RaIDeN allows multiple WGSs. FASTQ files can be zipped.

- #### `-o` : Name of output directory

    Specified name cannot exist.

- #### `-t` : Number of thread

#### **WARNING**

Selection by target motif has to be run separately by yourself.


## Output files

preparing now...
