#!/bin/bash 

# script to map intron reads 
# map the dmso sample
htseq-count -m union -f bam -t intron -s no -i Parent data/s_dmso.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/s_dmso.introns.counts.txt

# map the drb sample
htseq-count -m union -f bam -t intron -s no -i Parent data/s_drb_1.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/s_drb_1.introns.counts.txt

# map the flavopiridol sample 
htseq-count -m union -f bam -t intron -s no -i Parent data/s_flavopiridol_2.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/s_flavopiridol_2.introns.counts.txt

# map the triptolide sample 
htseq-count -m union -f bam -t intron -s no -i Parent data/s_triptolide_3.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/s_triptolide_3.introns.counts.txt

