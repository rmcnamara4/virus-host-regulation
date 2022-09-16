# Figure 3 - Main - Data

In this folder I have the data files that are used to plot the subfigures for
Figure 3 - Main. In the *Scripts/* folder I include the scripts that are used to
collect and create these data files.

## Scripts

+ **./Scripts/collect_opt.sh**: Bash script to copy codon stability coefficient values
(CSC) of various species to this folder

+ **./Scripts/collect_viral_rscu_fc_human.sh**: Bash script to copy the H-relative RSCU
values of all viral genomes downloaded from RefSeq

+ **./Scripts/create_dengue_isolates_m-relative_rscu_fc.R**: R script to calculate the
M-relative RSCU of all of the Dengue isolate sequences

+ **./Scripts/create_viral_rscu_fc_mosquito.R**: R script to calculate M-relative RSCU
values of all viral genomes downloaded from RefSeq; also calculates the Spearman correlation
of these values with mosquito CSC values and the p-value

+ **./Scripts/run_all.sh**: Bash script to run all of the data collection/creation scripts
listed above

## Data

+ **./codon_optimalities.csv**: contains CSC values of various species

+ **./viral_rscu_fc_human.csv**: contains H-relative RSCU values of all viral genomes
downloaded from RefSeq, along with their Spearman correlation with human CSC and the p-value

+ **./dengue_isolates_m-relative_rscu_fc.csv**: contains M-relative RSCU values of all of the
Dengue isolate sequences

+ **./viral_rscu_fc_mosquito.csv**: contains M-relative RSCU values of all viral genomes
downloaded from RefSeq, along with their Spearman correlation with mosquito CSC and the p-value 
