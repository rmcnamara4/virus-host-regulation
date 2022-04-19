#!/bin/bash

# script to download the mutation data from the paper:
# "Principles of dengue virus evolvability derived from genotype-fitness maps in human and mosquito cells"

# create directory to store data
mkdir -p data

# download the file
wget -O temp "https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNjE5MjEvZWxpZmUtNjE5MjEtc3VwcDEtdjIuY3N2LnppcA--/elife-61921-supp1-v2.csv.zip?_hash=nZifKX90cICmBSFUNS40PqGj5KPt2qJy1eZYOy6Iejk%3D"
unzip temp
rm -r temp
rm -r __*

# move the file
mv SupplementaryFile_1.csv data/point_mutations.csv
