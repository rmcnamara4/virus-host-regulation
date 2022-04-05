#!/bin/bash

# script to change the names of the fastas files to the name of the virus

for filename in ./data/fastas/*
do

    clean_name=${filename%_cds_from_genomic.fna.gz}
    clean_name=${clean_name##*/}

    virus_name=$(grep $clean_name ../complete_assembly.txt | cut -f8 | tr -s ' ' | tr ' ' '_' \
    | tr '/' ';')

    mv "$filename" "./data/fastas/$virus_name.fna.gz"

done
