# Figure 1 - Supplemental 1 - Data

In this folder I have the data files that are used to plot the subfigures for
Figure 1 - Supplemental 1. In the *Scripts/* folder I include the scripts that are used to
collect and create these data files.

## Scripts

+ **./Scripts/collect_data.sh**: Bash script to copy the data that has already been created
to this folder

+ **./Scripts/create_viral_aa_fc_human.R**: R script to calculate the amino acid frequency
fold change relative to human of all of the viral genomes downloaded from RefSeq

+ **./Scripts/run_all.sh**: Bash script to run all of the data collection/creation scripts
listed above

## Data

+ **./human_endo_csc.csv**: contains the human CSC values previously discovered by
Wu, et al. in *Translation affects mRNA stability in a codon-dependent manner in human cells*

+ **./viral_rscu_fc_human.csv**: contains the H-relative RSCU values for all of the viral
genomes downloaded from RefSeq; also contains the Spearman correlation of these values with
human optimality and the p-value

+ **./viral_aa_fc_human.csv**: contains amino acid frequency fold change relative to human of
all of the viral genomes downloaded from RefSeq

+ **./dengue_alignment_nucleotide.txt**: contains the nucleotide similarity percentage between
the reference sequences of all four Dengue serotypes

+ **./dengue_alignment_protein.txt**: contains the amino acid similarity percentage between
the reference sequences of all four Dengue serotypes 
