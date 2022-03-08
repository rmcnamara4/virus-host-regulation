#!/bin/bash

# script to get the cds fastas of Homo sapien and Aedes albopictus
# make dir to put the fastas
mkdir -p data

# get hg38 (human) cds fasta
wget https://webfs//n/analysis/indexes/hg38/annotation/Ens_102/fastas/hg38.Ens_102.cds.fa
mv hg38.Ens_102.cds.fa data/hg38.Ens_102.cds.fa

# get AaloF1 (mosquito) cds fasta
wget https://webfs//n/analysis/indexes/AaloF1/annotation/EnsGen_50/fastas/AaloF1.EnsGen_50.cds.fa
mv AaloF1.EnsGen_50.cds.fa data/AaloF1.EnsGen_50.cds.fa


