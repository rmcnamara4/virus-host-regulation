# 0.4-CreateSequenceStatsTables

Here I convert the sequence tables of human, mosquito, and the dengue isolates
into sequence stats tables. The Excel files for each species contain 4 sheets:
Nucleotides, Codons, Amino Acids, and RSCU.

## Scripts

+ **./make_sequence_stats_tables.R**: calculates nucleotide frequency, codon frequency,
amino acid frequency, and RSCU for the CDSs

+ **./run_all.sh**: runs the R script to create the sequence stats tables

## Data

+ **./data/human_hg38_stats.xlsx**: sequence stats tables of *Homo sapien* CDSs

+ **./data/mosquito_AaloF1_stats.xlsx**: sequence stats tables of *Aedes albopictus* CDSs

+ **./data/dengue_1_isolates_stats.xlsx**: sequence stats tables of Dengue 1 isolates CDSs

+ **./data/dengue_2_isolates_stats.xlsx**: sequence stats tables of Dengue 2 isolates CDSs

+ **./data/dengue_3_isolates_stats.xlsx**: sequence stats tables of Dengue 3 isolates CDSs

+ **./data/dengue_4_isolates_stats.xlsx**: sequence stats tables of Dengue 4 isolates CDSs
---
To reproduce the analysis:

```bash
bash run_all.sh
```
