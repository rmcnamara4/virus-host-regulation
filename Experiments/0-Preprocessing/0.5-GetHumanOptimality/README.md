# 0.5-GetHumanOptimality

Here I download the human codon stability coefficients produced by the Bazzini
lab. The dataset is downloaded from the paper: *Translation affects mRNA stability
in a codon-dependent manner in human cells* by Wu, et al (eLife 2019).

The dataset is then filtered for the endogenous CSCs and the mean is taken
across the three cell types.

## Scripts

+ **./get_human_csc.sh**: downloads the whole dataset from the paper

+ **./gather_human_endo_csc.R**: filters for and takes the mean of the endogenous CSCs

+ **./run_all.sh**: Runs both scripts to reproduce the data

## Data

+ **./data/human_csc_all_methods**: whole dataset from the paper, showing all CSC values

+ **./data/human_endo_csc.csv**: filtered dataset, showing endogenous CSCs and average CSC


To reproduce the data:

```bash
bash run_all.sh
```
