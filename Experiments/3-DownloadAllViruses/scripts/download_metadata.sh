#!/bin/bash

# script to download the metadata required for this analysis
# download virus metadata from NCBI
wget -O "metadata.txt" "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download"
mv metadata.txt ../metadata.txt

# download complete assembly of viruses from RefSeq genome database
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
mv assembly_summary.txt ../complete_assembly.txt
