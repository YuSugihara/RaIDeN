# RaIDeN

## What is RaIDeN?
RaIDeN is a pipeline for **Ra**pid **IDe**ntification of causal gene with target motif using **N**GS technology

<img src="https://github.com/YuSugihara/RaIDeN/blob/master/images/Fig.S1.png" width=700>

Supplementary figure 1 in [Shimizu et al. (2021)](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1)

### Citation
Shimizu M, Hirabuchi A, Sugihara Y, Abe A, Takeda T, Kobayashi M, Hiraka Y, Kanzaki E, Oikawa K, Saitoh H, Langner T, Banfield MJ, Kamoun S, Terauchi R (2021) [A genetically linked pair of NLR immune receptors show contrasting patterns of evolution](https://www.biorxiv.org/content/10.1101/2021.09.01.458560v1). _bioRxiv_

### How to install
```
git clone https://github.com/YuSugihara/RaIDeN.git
cd RaIDeN
pip install .
```

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
