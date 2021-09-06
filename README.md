# RaIDeN User Guide
#### version 0.0.2 (beta)


## Table of contents
- [What is RaIDeN?](#What-is-RaIDeN)
  + [Citation](#Citation)
- [Installation](#Installation)
  + [Dependencies](#Dependencies)
  + [How to install](#How-to-install)
- [What RaIDeN requires](#What-RaIDeN-requires)
- [Usage](#Usage)
  + [Example 1 : run RaIDeN from FASTQ with trimming](#Example-1--run-RaIDeN-from-FASTQ-with-trimming)
  + [Example 2 : run jiji from VCF](#Example-5--run-QTL-plot-from-VCF)
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

Currently, RaIDeN doesn't support the installation via bioconda. However, you can easily install its dependencies because they are deposited on bioconda. You can try the command below:

```
$ conda install -c bioconda samtools bcftools hisat2 bedtools gffread bamtools pigz stringtie faqcs prinseq
```


## What RaIDeN requires:

- #### Reference genome (FASTA format; option `-r`)

    RaIDeN requires an assembled reference genome. This reference genome must contain a causal gene. If it is sure that the assembled reference genome contains the causal gene, RaIDeN allows contiguous (not chromosome-scale) referece genome.

- #### RNA-seq for gene annotation (paired FASTQ format; option `-a`)

    RaIDeN requires RNA-seq for gene annotation. This RNA-seq must contain the causal gene because RaIDen only analyzes the annotated genomic region. RaIDeN allows multiple RNA-seq samples.

- #### Whole-genome sequences (paired FASTQ format; option `-w`)

    RaIDen requires whole-genome sequence (WGS). The sample of this WGS should have an opposite trait to that of the reference genome. Since RaIDeN expects that the different traits come from structual variation or mutation on the annotated genomic region, this WGS must have a structual variation or mutation on the causal gene. RaIDeN allows multiple WGSs.



## Usage

```
$ raiden -h

usage: raiden -r <FASTA>
              -a <FASTQ1>,<FASTQ2> 
              -w <FASTQ1>,<FASTQ2>
              -o <OUT_DIR> 
             [-t <INT>]

RaIDeN version 0.0.3

optional arguments:
  -h, --help            show this help message and exit
  -r , --ref            Reference fasta. You may use De novo assembled contigs.
  -a , --rna-seq        Fastqs of RNA-seq for genome annotation on input reference genome.
                        RNA-seq must contain the seqeunce reads of a causal gene.
                        You can use this optiion multiple times for each sample
                        Please separate paired fastqs by comma (e.g. -a fastq1,fastq2).
  -w , --whole-genome   Fastqs of whole-genome sequence to filter causal genes.
                        This sample should have an opposite trait to the sample
                        used for the refernce genome. You can use this optiion multiple
                        times for each sample. Please separate paired fastqs by comma
                        (e.g. -w fastq1,fastq2).
  -o , --out            Output directory. Specified name must not exist.
  -t , --threads        Number of threads. If you specify the number below one,
                        then RaIDeN will use the maximum number of threads. [2]
  -d, --disable-RNAseq-trim
                        Disable the trimming of RNA-seq reads by FaQCs.
  -D, --disable-WGS-trim
                        Disable the trimming of whole-genome sequence
                        by FaQCs and prinseq-lite.
  -s , --strand         Strand library option.
                         'fr' : assume a strand library fr-first strand.
                         'rf' : assume a strand library fr-second strand.
                         'None' : Don't assume a strand library.
                        Default is 'None'.
  -m , --minimum-len    Minimum length of the transcripts predicted by
                        stringtie. [200]
  -q , --min-MQ         Minimum mapping quality during mpileup. [40]
  -Q , --min-BQ         Minimum base quality during mpileup. [18]
  -v, --version         show program's version number and exit
```

RaIDeN can start from different steps.

+ [Example 1 : run RaIDeN from FASTQ with trimming](#Example-1--run-RaIDeN-from-FASTQ-with-trimming)
+ [Example 2 : run jiji from VCF](#Example-5--run-QTL-plot-from-VCF)


### Example 2 : run jiji from VCF

```
$ jiji -h

usage: jiji -a <GFF/GTF>
            -v <VCF>
            -b <BED_DIR>
            -o <OUT_DIR> 
             

jiji version 0.0.3

optional arguments:
  -h, --help       show this help message and exit
  -a , --gff       GFF/GTF containing gene annotation.
  -v , --vcf       VCF generated by RaIDeN.
  -b , --bed       Directory of BED files generated by RaIDeN.
  -o , --out       Output directory. Specified name must not exist.
  --inconsistent   Number of the allowed inconsistent markers for SNPs and indels. [0]
  --miss           Number of the allowed missing markers for SNPs and indels. [0]
  --count-het      Count heterozygotes as causal mutation.
  --coverage       Skip a gene whose coverage more than this parameter for P/A marker. [0.5]
  --version        show program's version number and exit
```


## Output files

preparing now...
