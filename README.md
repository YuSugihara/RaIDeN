# RaIDeN

## What is RaIDeN?
RaIDeN is a pipeline for **Ra**pid **Ide**ntification of causal gene with target motif using **N**GS technology

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
