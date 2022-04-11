#!/bin/bash

# script to download the metadata and count files from the paper:
# "Single-cell transcriptional dynamics of flavivirus infection"
# script also downloads the Dengue and Zika sequences used for infection in
# the paper.

# make folders for data
mkdir -p data
mkdir -p data/metadata
mkdir -p data/counts
mkdir -p data/fastas

# download the folder
wget -O temp "https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzI5NDIvZWxpZmUtMzI5NDItc3VwcDctdjIuZ3o-/elife-32942-supp7-v2.gz?_hash=EyCaZA9wZcUhyxTkbnyGtHzZUVxVyyp9Qq4Qdmu4V%2Fc%3D"
tar -xvf temp
rm temp

# change file names
mv cell_metadata_Zika.tsv cell_metadata_zika.tsv
mv counts_Zika.tsv counts_zika.tsv

# move files
mv cell_metadata_* data/metadata/
mv counts_* data/counts/

# download the Dengue virus type 2 (strain 16681) cds sequence
wget -O dengue_2_16681.fa 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_na&id=2155257&from=97&to=10272&extrafeat=null&conwithfeat=on&hide-cdd=on'

# download the Zika virus (strain PRVABC59) cds sequence
wget -O zika_PRVABC59.fa 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_na&id=984874581&from=107&to=10378&extrafeat=null&conwithfeat=on&hide-cdd=on'

# move files
mv *.fa data/fastas
