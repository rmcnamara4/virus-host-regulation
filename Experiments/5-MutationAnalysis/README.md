# 5-MutationAnalysis

Here I analyze the Dengue mutation data that is produced in the paper: *Principles of dengue virus evolvability derived from genotype-fitness maps in human and mosquito cells* by Dolan, et al. The paper can be found [here](https://elifesciences.org/articles/61921#content).

In this folder I investigate the fitness of Dengue virus 2 strain 11681's mutations in both human and mosquito hosts. I seek to understand the way in which Dengue adapts to its different hosts and whether the codon usage of the virus relative to its host plays a role in determining which mutations are beneficial and deleterious.

## Scripts

+ **./download_data.sh**: Bash script to download the mutation data produced in the paper mentioned above

+ **./extract_passage_9.R**: Rscript to filter for only the mutations that were recorded between the reference strain and the 9th passage of the virus and create a separate dataset

+ **./create_rscu_fc_df.R**: Rscript to create a file with the RSCU fold change of Dengue 2 strain 16681 with respect to both mosquito and human (for later use in the analysis)

+ **./create_mutation_dfs.R**: Rscript to filter and aggregate the mutation data in different ways for later plotting; script creates a new file for each of the data.frames created

+ **./helper_functions.R**: Rscript of defined functions used in the analysis for loading and plotting data

+ **./run_all.sh**: Bash script to create the required folders and run the scripts above

## RMarkdown

+ **./dengue_mutation_analysis.Rmd**: RMarkdown file that analyzes the mutation data and produces a number of visualizations

## HTML

+ **./dengue_mutation_analysis.html**: HTML file that is produced when the RMarkdown file is knitted

## Data

+ **./data/point_mutations.csv**: mutation data downloaded from the paper: *Principles of dengue virus evolvability derived from genotype-fitness maps in human and mosquito cells* by Dolan, et al. The dataset is Supplementary File 1 and can be found [here](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjE5MjEvZWxpZmUtNjE5MjEtc3VwcDEtdjIuY3N2LnppcA--/elife-61921-supp1-v2.csv.zip?_hash=nZifKX90cICmBSFUNS40PqGj5KPt2qJy1eZYOy6Iejk%3D) (downloaded by **./download_data.sh** script)

+ **./data/point_mutations_p9.csv**: mutation data filtered for only the 9th passage (produced by **./extract_passage_9.R** script)

+ **./data/rscu_fc.csv**: RSCU fold change data of Dengue 2 strain 16681 with respect to both mosquito and human (produced by **./create_rscu_fc_df.R** script)

+ **./data/condensed_point_mutations_p9.csv**: Mean Relative Fitness of all 9th passage mutations; replicates are aggregated for each species; wtcodon and mutcodon CSC and RSCU fold change are included, along with the delta resulting from the mutation (produced by **./create_mutation_dfs.R** script)

+ **./data/condensed_synonymous_point_mutations_p9.csv**: Mean Relative Fitness of all 9th passage synonymous mutations; all synonymous mutations are aggregated by mutcodon for each species and replicate; one entry per codon per species per replicate (59 * 2 * 2 = 236 rows); fitness class of each mutation is included (produced by **./create_mutation_dfs.R** script)

+ **./data/fitness_class_percentages_per_syn_mutation_p9.csv**: all 9th passage synonymous mutations; aggregated by replicate for each species; fitness class is included; number of instances of each mutation is included; for each mutation the percentage is the fraction of the mutations that are of a particular fitness class (produced by **./create_mutation_dfs.R** script)

+ **./data/syn_mut_mrf_per_species_p9.csv**: Mean Relative Fitness of all 9th passage synonymous mutations; all synonymous mutations are aggregated by mutcodon and replicate for each species; one entry per codon per species (59 * 2 = 118 rows); this dataset is very similar to the **./data/condensed_synonymous_point_mutations_p9.csv** file, except the replicates are also aggregated (produced by **./create_mutation_dfs.R** script)

+ **./data/codon_groups.xlsx**: includes the codons that are in the Frequently Used, Infrequently Used, Preferentially Used, Unpreferentially Used, Denguenized, and Non-denguenized groups for both mosquito and human; Excel file contains two sheets: "Human" and "Mosquito" (produced by **./dengue_mutation_analysis.Rmd** file)

+ **./data/syn_mrf_codon_groups_p9.csv**: Mean Relative Fitness of 9th passage synonymous mutations towards each codon included in the codon groups mentioned above; aggregated by mutcodon and replicate for each species; includes the host and codon group (produced by **./dengue_mutation_analysis.Rmd** file)

## Figures

Note: all figures produced by **./dengue_mutation_analysis.Rmd** file
Note: a neutra mutation is defined as having a Mean Relative Fitness of 1 (as defined in Dolan, et al.)

+ **./figures/syn_mut_all_against_csc.pdf**: plots of Mean Relative Fitness of 9th passage synonymous mutations against CSC of the mutcodon; each point represents the codon that is mutated to and its Mean Relative Fitness/CSC; the four plots correspond to the A & B replicates of human and mosquito; correlations are Spearman

+ **./figures/syn_mut_per_aa_against_csc.pdf**: plots of Mean Relative Fitness of 9th passage synonymous mutations against CSC of the mutcodon; each point represents the codon that is mutated to and its Mean Relative Fitness/CSC; each plot contains the codons for one amino acid; replicates are represented by different point shapes; correlations are Spearman

+ **./figures/rep_corr_all_against_csc.pdf**: plots of the Spearman correlation of replicate B vs the Spearman correlation of replicate A for each amino acid (as plotted in **./figures/syn_mut_per_aa_against_csc.pdf**); plots correspond to human and mosquito hosts; the number of codons in an amino acid is represented by the size of the point; the color of the point corresponds to the CSC of the most optimal codon minus the CSC of the most non-optimal codon

+ **./figures/syn_mut_all_against_rscu_fc.pdf**: plots of Mean Relative Fitness of 9th passage synonymous mutations against RSCU fold change of the mutcodon; each point represents the codon that is mutated to and its Mean Relative Fitness/RSCU fold change; the four plots correspond to the A & B replicates of human and mosquito; correlations are Spearman

+ **./figures/syn_mut_per_aa_against_rscu_fc.pdf**: plots of Mean Relative Fitness of 9th passage synonymous mutations against RSCU fold change of the mutcodon; each point represents the codon that is mutated to and its Mean Relative Fitness/RSCU fold change; each plot contains the codons for one amino acid; replicates are represented by different point shapes; correlations are Spearman

+ **./figures/rep_corr_all_against_rscu_fc.pdf**: plots of the Spearman correlation of replicate B vs the Spearman correlation of replicate A for each amino acid (as plotted in **./figures/syn_mut_per_aa_against_rscu_fc.pdf**); plots correspond to human and mosquito hosts; the number of codons in an amino acid is represented by the size of the point; the color of the point corresponds to the largest RSCU fold change minus the smallest RSCU fold change within the amino acid

+ **./figures/syn_mut_to_aa_codons_with_csc_as_fill.pdf**: barplots showing the Mean Relative Fitness of each 9th passage synonymous mutation, separated by amino acid; standard deviations are shown; x-axis corresponds to the codon that is mutated from; the top of the plots corresponds to the codon mutated to; color of the bar corresponds to the change in CSC caused by the mutation (red = mutation to more optimal codon; blue = mutation to less optimal codon); numbers over bars represent the number of occurences of the mutation; replicates are aggregated together

+ **./figures/syn_mut_to_aa_codons_with_rscu_fc_as_fill.pdf**: barplots showing the Mean Relative Fitness of each 9th passage synonymous mutation, separated by amino acid; standard deviations are shown; x-axis corresponds to the codon that is mutated from; the top of the plots corresponds to the codon mutated to; color of the bar corresponds to the change in RSCU fold change caused by the mutation (green = mutation to a codon with a higher RSCU fold change; purple = mutation to a codon with a lower RSCU fold change); numbers over bars represent the number of occurences of the mutation; replicates are aggregated together

+ **./figures/syn_mut_to_aa_codons_with_fitness_class_as_fill.pdf**: barplots showing the fraction of 9th passage synonymous mutations that are beneficial, neutral, deleterious, and lethal; separated by amino acid; x-axis corresponds to the codon that is mutated from; the top of the plots corresponds to the codon mutated to; color of the bar corresponds to the fitness class of the mutation; numbers over bars represent the number of occurences of the mutation; replicates are aggregated together

+ **./figures/syn_mut_toward_codon_groups.pdf**: boxplot of the Mean Relative Fitness of 9th passage synonymous mutations toward codons present in the codon groups defined in **./data/codon_groups.xlsx**; color of the box represents the codon group; stars indicate significance (unpaired Wilcox test); numbers above the boxes represent the number of codons in each of the codon groups; replicates aggregated together

+ **./figures/mrf_corr_human_mosquito.pdf**: scatter plot of the Mean Relative Fitness of 9th passage synonymous mutations in mosquito vs human; each point represents the codon that is mutated to; correlation is Spearman 
