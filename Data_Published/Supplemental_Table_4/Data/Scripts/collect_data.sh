#!/bin/bash

# script to collect all required data
# collect human and mosquito optimality
cp ../../../../Experiments/0-Preprocessing/0.5-GetHumanOptimality/data/human_endo_csc.csv ..
cp ../../../../Experiments/1-CalculateMosquitoOptimality/data/mosquito_csc.csv ..

# collect viral RSCU fold change for human and mosquito
cp ../../../../Experiments/3-DownloadAllViruses/data/rscu_ratios/all_rscu_ratio_total_* ..

# collect Dengue isolates RSCU fold change for human and mosquito
cp ../../../../Figures_Published/Figure_1/Main/Data/dengue_isolates_h-relative_rscu_fc.csv ..
cp ../../../../Figures_Published/Figure_3/Main/Data/dengue_isolates_m-relative_rscu_fc.csv ..
