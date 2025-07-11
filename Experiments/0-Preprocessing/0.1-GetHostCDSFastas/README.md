# 0.1-GetHostCDSFastas

Collect fasta files for the CDSs of Homo sapien (human) and Aedes albopictus (mosquito) genomes.
These fasta files are available through Ensembl.

## Scripts

+ **./get_cds_fastas.sh**: downloads the cds fastas of human and mosquito

## Data

+ **./data/hg38.Ens_102.cds.fa**: [Human_CDS_Fasta](https://webfs//n/analysis/indexes/hg38/annotation/Ens_102/fastas/hg38.Ens_102.cds.fa)

+ **./data/AaloF1.EnsGen_50.cds.fa**: [Mosquito_CDS_Fasta](https://webfs//n/analysis/indexes/AaloF1/annotation/EnsGen_50/fastas/AaloF1.EnsGen_50.cds.fa)
---
To reproduce the data:

```bash
bash get_cds_fastas.sh
```
