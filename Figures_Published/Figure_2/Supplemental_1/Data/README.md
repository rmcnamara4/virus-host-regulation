# Figure 2 - Supplemental 1 - Data

In this folder I have the data files that are used to plot the subfigures for
Figure 2 - Supplemental 1. In the *Scripts/* folder I include the scripts that are used to
collect and create these data files.

## Scripts

+ **./Scripts/gather_optimality_data.sh**: Bash script to copy the codon stability coefficient
(CSC) values for various species and *Aedes albopictus* C6/36 cells calculated using
3 different transcription-blocking drugs to this folder

+ **./Scripts/create_egfp_reporters_codon_frequency_ratio.R**: R script to calculate the codon
composition of the EGFP reporters and calculate the optimal to non-optimal ratio of each
sequence (relative to mosquito)

+ **./Scripts/create_mosquito_codon_frequency_ratio.R**: R script to calculate the codon composition
of the genes in the *Aedes albopictus* genome and calculate the optimal to non-optimal
ratio of each sequence (relative to mosquito)

+ **./Scripts/create_nt_out_of_frame_reporters_codon_frequency_ratio.R**: R script to calculate the
codon composition of the 1nt-out of frame reporters and calculate the optimal to non-optimal
ratio of each sequence along with the frequency of neutral codons in each sequence (relative
to mosquito)

## Data

+ **./codon_optimalities.csv**: contains the CSC values for various species

+ **./mosquito_csc.csv**: contains the CSC values calculated for *Aedes albopictus* C6/36
cells using 3 different transcription-blocking drugs

+ **./egfp_reporters_codon_frequency_ratio.csv**: contains the codon frequency values of the EGFP
reporters and the optimal to non-optimal ratio of each sequence (relative to mosquito)

+ **./mosquito_codon_frequency_ratio.csv**: contains the codon frequency values of the genes
in the *Aedes albopictus* genome and the optimal to non-optimal ratio of each sequence (relative
to mosquito)

+ **./nt_out_of_frame_reporters_codon_frequency_ratio.csv**: contains the codon frequency values
of the 1nt-out of frame reporters and the optimal to non-optimal ratio of each sequence along
with the frequency of neutral codons in each sequence (relative to mosquito)
