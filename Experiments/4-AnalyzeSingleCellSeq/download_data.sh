#!/bin/bash

# script to download the metadata and count files from the paper:
# "Single-cell transcriptional dynamics of flavivirus infection"

# make folders for data
mkdir -p data
mkdir -p data/metadata
mkdir -p data/counts

# download the folder
wget https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzI5NDIvZWxpZmUtMzI5NDItc3VwcDctdjIuZ3o-/elife-32942-supp7-v2.gz?_hash=EyCaZA9wZcUhyxTkbnyGtHzZUVxVyyp9Qq4Qdmu4V%2Fc%3D
gunzip elife-32942-supp7-v2

# change file names
cd elife-32942-supp7-v2
mv cell_metadata_Zika.tsv cell_metadata_zika.tsv
mv counts_Zika.tsv counts_zika.tsv

# move files
mv cell_metadata_* ../data/metadata/
mv counts_* ../data/counts/

# delete folder
cd ..
rm -r elife-32942-supp7-v2
