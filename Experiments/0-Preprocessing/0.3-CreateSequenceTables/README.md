# 0.3-CreateSequenceTables

Here I convert the CDS fasta files of human, mosquito, and the dengue isolates
into tables of species, gene_ID, and sequence.

# Scripts

+ **./make_human_mosquito_seq_tables.R**: converts human and mosquito fastas to
sequence tables

+ **./make_isolates_seq_tables.R**: converts Dengue isolates fasta to sequence
tables

+ **./run_all.sh**: Runs the R scripts to produce the sequence tables

## Data

+ **./data/human_hg38_seq.csv**: sequence table of *Homo sapien* CDSs

+ **./data/mosquito_AaloF1_seq.csv**: sequence table of *Aedes albopictus* CDSs

+ **./data/dengue_1_isolates_seq.csv**: sequence table of Dengue 1 isolates CDSs

+ **./data/dengue_2_isolates_seq.csv**: sequence table of Dengue 2 isolates CDSs

+ **./data/dengue_3_isolates_seq.csv**: sequence table of Dengue 3 isolates CDSs

+ **./data/dengue_4_isolates_seq.csv**: sequence table of Dengue 4 isolates CDSs
---
To reproduce the analysis:

```bash
bash run_all.sh
```
