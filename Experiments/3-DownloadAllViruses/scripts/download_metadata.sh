#!/bin/bash

# script to download the metadata required for this analysis
# download virus metadata from NCBI
wget -O "metadata.txt" "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download"
mv metadata.txt ./metadata.txt

# download complete assembly of viruses from RefSeq genome database
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
mv assembly_summary.txt ./complete_assembly.txt

# download and unzip the gbff files from RefSeq to get other virus information
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.genomic.gbff.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.genomic.gbff.gz

# concatenate the gbff files to one file
cat *.gbff.gz > viruses.gbff.gz

# remove the individual files
rm viral.1.genomic.gbff.gz
rm viral.2.genomic.gbff.gz
rm viral.3.genomic.gbff.gz
rm viral.4.genomic.gbff.gz
