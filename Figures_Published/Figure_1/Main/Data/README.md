# Figure 1 - Main - Data

In this folder I have the data files that are used to plot the subfigures for
Figure 1 - Main. In the *Scripts/* folder I include the scripts that are used to
collect and create these data files.

## Scripts

+ **./Scripts/collect_human_opt.sh**: Bash script to copy human codon stability
coefficient (CSC) (optimality) values to this folder

+ **./Scripts/create_dengue_2_iso_and_human_codon_freq.R**: R script to calculate the
codon composition of the Dengue 2 isolates and human gene sequences; also calculates the
optimal to non-optimal ratio for each sequence (relative to human)

+ **./Scripts/create_dengue_isolates_h-relative_rscu_fc.R**: R script to calculate the
RSCU fold change relative to human for all of the Dengue isolate sequences

+ **./Scripts/create_viral_rscu_fc_human.R**: R script to calculate the RSCU fold change
relative to human for all of the viral genomes downloaded from RefSeq; also calculates the
Spearman correlation of H-relative RSCU with human optimality and the p-value

+ **./Scripts/run_all.sh**: Bash script to run all of the data collection/creation scripts
listed above

## Data

+ **./human_endo_csc.csv**: contains human endogenous CSC values previously discovered by
Wu, et al. in *Translation affects mRNA stability in a codon-dependent manner in human cells*

+ **./dengue_2_iso_and_human_codon_freq.csv**: contains codon frequency values for the Dengue
2 isolates and human gene sequences; also contains the optimal to non-optimal ratio for each
sequence (relative to human)

+ **./dengue_isolates_h-relative_rscu_fc.csv**: contains the H-relative RSCU values for all of
the Dengue isolates

+ **./viral_rscu_fc_human.csv**: contains the H-relative RSCU values for all of the viral
genomes downloaded from RefSeq; also contains the Spearman correlation of these values with
human optimality and the p-value
