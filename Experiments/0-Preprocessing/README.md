# 0-Preprocessing

In this folder I have subdirectories for individual experiments performed as
preprocessing steps for other experiments.

## Subdirectories

+ **./0.1-GetHostCDSFastas/**: collect the CDS fasta files for *Homo sapiens* (human)
and *Aedes albopictus** (mosquito)

+ **./0.2-DengueIsolatesCDSFastas/**: collect the CDS fasta files for Dengue sequences
isolated from around the world

+ **./0.3-CreateSequenceTables/**: convert fasta files into tables containing species and
sequence for ease of use later

+ **./0.4-CreateSequenceStatsTables/**: calculate the nucleotide composition, codon
composition, amino acid composition, and relative synonymous codon usage (RSCU) for each
of the sequences collected (i.e. human, mosquito, Dengue isolates)

+ **./0.5-GetHumanOptimality/**: collect the codon stability coefficients (CSC) for humans,
which we've already known

+ **./0.6-OtherOptimalities/**: collect the CSCs for other known species
