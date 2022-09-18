#!/bin/bash

# script to collect required data
# collect human sequence stats data
cp ../../../../Experiments/0-Preprocessing/0.4-CreateSequenceStatsTables/data/human_hg38_stats.xlsx ..

# collect edgeR data during Dengue infection
cp ../../../../Experiments/4-AnalyzeSingleCellSeq/data/edgeR/uninfected-high.csv ..

# collect GO term data
cp ../../../../Experiments/4-AnalyzeSingleCellSeq/data/go_terms/human_go_terms.csv ..
