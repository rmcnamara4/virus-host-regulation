# 0.4-CreateSequenceStatsTables

Here I convert the sequence tables of human, mosquito, and the dengue isolates
into sequence stats tables. The Excel files for each species contain 4 sheets:
Nucleotides, Codons, Amino Acids, and RSCU.

The tables are produced in the *data* folder.

## Data

+ **./data/human_hg38_stats.xlsx**: sequence stats tables of *Homo sapien* CDSs

+ **./data/mosquito_AaloF1_stats.xlsx**: sequence stats tables of *Aedes albopictus* CDSs

+ **./data/dengue_1_isolates_stats.xlsx**: sequence stats tables of Dengue 1 isolates CDSs

+ **./data/dengue_2_isolates_stats.xlsx**: sequence stats tables of Dengue 2 isolates CDSs

+ **./data/dengue_3_isolates_stats.xlsx**: sequence stats tables of Dengue 3 isolates CDSs

+ **./data/dengue_4_isolates_stats.xlsx**: sequence stats tables of Dengue 4 isolates CDSs

The script *run_all.sh* will run the R script to create these tables.

To reproduce the analysis:

```bash
bash run_all.sh
```
