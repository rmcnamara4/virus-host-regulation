# Figure 3 - Supplemental 1 - Data

In this folder I have the data files that are used to plot the subfigures for
Figure 3 - Supplemental 1. In the *Scripts/* folder I include the scripts that are used to
collect and create these data files.

## Scripts

+ **./Scripts/collect_data.sh**: Bash script to copy existing data to this folder

+ **./Scripts/create_viral_aa_fc_mosquito.R**: R script to calculate the amino acid
frequency fold change relative to mosquito for all viral genomes downloaded
from RefSeq

+ **./Scripts/run_all.sh**: Bash script to run all of the data collection/creation
scripts listed above

## Data

+ **./codon_optimalities.csv**: contains CSC values of various species

+ **./dengue_isolates_m-relative_rscu_fc.csv**: contains M-relative RSCU values of all of the
Dengue isolate sequences

+ **./viral_rscu_fc_human.csv**: contains H-relative RSCU values of all viral genomes
downloaded from RefSeq, along with their Spearman correlation with human CSC and the p-value

+ **./viral_rscu_fc_mosquito.csv**: contains M-relative RSCU values of all viral genomes
downloaded from RefSeq, along with their Spearman correlation with mosquito CSC and the p-value

+ **./viral_aa_fc_mosquito.csv**: contains the amino acid frequency fold change relative to
mosquito for all viral genomes downloaded from RefSeq
