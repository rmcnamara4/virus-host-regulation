#!/bin/bash

# script to collect data
# collect human optimality
cp ../../../../../Experiments/0-Preprocessing/0.5-GetHumanOptimality/data/human_endo_csc.csv ..

# collect differential expression data
cp ../../../../../Experiments/4-AnalyzeSingleCellSeq/data/edgeR/uninfected-low.csv ..
cp ../../../../../Experiments/4-AnalyzeSingleCellSeq/data/edgeR/uninfected-med.csv ..
cp ../../../../../Experiments/4-AnalyzeSingleCellSeq/data/edgeR/uninfected-high.csv ..

# collect human genome codon frequency
cp ../../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data/human_hg38_stats.xlsx ..

# collect codon groups
cp ../../../../../Experiments/5-MutationAnalysis/data/codon_groups.xlsx ..
