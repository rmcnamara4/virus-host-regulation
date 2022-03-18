#!/bin/bash

# script to download all of the cds fastas using the urls in the complete assembly file

# extract the urls and save to file 
cat ../complete_assembly.txt \
    | awk 'BEGIN{FS"\t"}{print $20}' \
    | awk 'BEGIN{OFS=FS="/"}{print $0,$NF"_cds_from_genomic.fna.gz"}' \
    > ../urls.txt

# make directory to put the fastas
cd ..
mkdir fastas
cd fastas

# download the fastas
wget -i ../urls.txt
