# 0.3-CreateSequenceTables

Here I convert the CDS fasta files of human, mosquito, and the dengue isolates
into tables of species, gene_ID, and sequence.

The tables are produced in the *data* folder.

## Data

+ **./data/human_hg38_seq.csv**: sequence table of *Homo sapien* CDSs

+ **./data/mosquito_AaloF1_seq.csv**: sequence table of *Aedes albopictus* CDSs

+ **./data/dengue_1_isolates_seq.csv**: sequence table of Dengue 1 isolates CDSs

+ **./data/dengue_2_isolates_seq.csv**: sequence table of Dengue 2 isolates CDSs

+ **./data/dengue_3_isolates_seq.csv**: sequence table of Dengue 3 isolates CDSs

+ **./data/dengue_4_isolates_seq.csv**: sequence table of Dengue 4 isolates CDSs

The script *run_all.sh* will run all of the R scripts to perform the analysis.

To reproduce the analysis:

```bash
bash run_all.sh
```
