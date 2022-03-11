#!/bin/bash 

# script to separate exons and introns in the gff file 
# create exon gff file 
cat data/AaloF1.EnsGen_50.gff | grep exon | grep protein_coding > AaloF1.exons.gff

# create intron gff file 
cat data/AaloF1.EnsGen_50.gff | grep intron > AaloF1.introns.gff 