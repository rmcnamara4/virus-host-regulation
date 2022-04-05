# 2-ConcatenateOptimalities

Here I simply concatenate the codon optimalities of human (*H. sapiens*),
mosquito (*Ae. albopictus*), and others (zebrafish, frog, yeast, and fly) into
a single data.frame.

## Scripts

+ **./concatenate_optimalities.R**: R script to concatenate the codon optimality
data.frames

+ **./run_all.sh**: Bash script to run the R script and create folder to store the
data

## Data

+ **./data/codon_optimalities.csv**: file containing all of the known codon
optimalities of our lab
---
To reproduce the analysis:

```bash
bash run_all.sh
```
