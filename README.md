# Projeto
## Context
Variant Calling Pipeline written in shell for Plasmodium falciparum, using conda interface and packages:
- Bioconda: bwa
- Bioconda: samtools
- Bioconda: picard
- Bioconda: gatk4
- Bioconda: snpeff
- Conda-forge: r
- Bioconda: r-gsalib

To test this Pipeline you can use the files that are in the following links: 

https://uminho365-my.sharepoint.com/:f:/g/personal/pg45971_uminho_pt/EkwFr4fTbBRHqTdxQvL2lP8BEEBgTzDnuoUDGxBld39JJA?e=GkTGSBn (only UMinho) 

https://osf.io/ftkqw/?view_only=9dc1a3496e99420bb42ba5501ac1bb51 

## Abstract
Malaria is still today one of the major global diseases, affecting almost half of the world's population. Caused by several plasmodia parasites, *Plasmodium falciparum* stands out as one of the main causes of the disease. The use of antimalarial drugs is essential to fight the disease, however, due to its wide use, another problem arises, which is the parasites' resistance to antimalarial treatments.
The resistance usually occurs due to mutations (SNPs and Indels) in the Plasmodium's genome, which makes the treatment of the patient very difficult. Therefore, it is necessary to identify the variants causing the resistance. For that, in this study we created a pipeline adapted for *Plasmodium falciparum*, which can be run in any operating system, being only necessary a Conda interface.
The results obtained were VCF files containing both SNPs and Indels, the type of mutations observed and their corresponding location on chromosomes and in genes. This information can then be used to predict which antimalarials drugs the parasite has gained resistance to and its provenance, requiring only the comparison of the data obtained with data from databases such as MalariaGen.
