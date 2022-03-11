#!/bin/bash

# script to gather all of the required data for the analysis 
# make dir to put the data 
mkdir -p data

# get Aedes albopictus gff file 
wget https://webfs/n/analysis/indexes/AaloF1/annotation/EnsGen_50/gtfs/AaloF1.EnsGen_50.gff
mv AaloF1.EnsGen_50.gff data/AaloF1.EnsGen_50.gff

# get dmso aligned bam file and bai file 
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_dmso/s_dmso.Aligned.sortedByCoord.out.bam 
mv s_dmso.Aligned.sortedByCoord.out.bam data/s_dmso.Aligned.sortedByCoord.out.bam

wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_dmso/s_dmso.Aligned.sortedByCoord.out.bam.bai
mv s_dmso.Aligned.sortedByCoord.out.bam.bai data/s_dmso.Aligned.sortedByCoord.out.bam.bai

# get drb aligned bam file and bai file 
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_drb_1/s_drb_1.Aligned.sortedByCoord.out.bam
mv s_drb_1.Aligned.sortedByCoord.out.bam data/s_drb_1.Aligned.sortedByCoord.out.bam

wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_drb_1/s_drb_1.Aligned.sortedByCoord.out.bam.bai 
mv s_drb_1.Aligned.sortedByCoord.out.bam.bai data/s_drb_1.Aligned.sortedByCoord.out.bam.bai 

# get flavopiridol bam file and bai file 
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_flavopiridol_2/s_flavopiridol_2.Aligned.sortedByCoord.out.bam 
mv s_flavopiridol_2.Aligned.sortedByCoord.out.bam data/s_flavopiridol_2.Aligned.sortedByCoord.out.bam 

wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_flavopiridol_2/s_flavopiridol_2.Aligned.sortedByCoord.out.bam.bai 
mv s_flavopiridol_2.Aligned.sortedByCoord.out.bam.bai data/s_flavopiridol_2.Aligned.sortedByCoord.out.bam.bai 

# get triptolide bam file and bai file 
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_triptolide_3/s_triptolide_3.Aligned.sortedByCoord.out.bam   
mv s_triptolide_3.Aligned.sortedByCoord.out.bam data/s_triptolide_3.Aligned.sortedByCoord.out.bam   

wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/s_triptolide_3/s_triptolide_3.Aligned.sortedByCoord.out.bam.bai   
mv s_triptolide_3.Aligned.sortedByCoord.out.bam.bai data/s_triptolide_3.Aligned.sortedByCoord.out.bam.bai   

# get exon tpm table
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/secundo/RSEM_TPM_table.csv
mv RSEM_TPM_table.csv data/exon_tpms.csv

# get exon counts table
wget https://webfs/n/core/Bioinformatics/secondary/Bazzini/lc2196/MOLNG-3268/secundo/star_count.csv
mv star_count.csv data/exon_counts.csv
