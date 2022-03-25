#!/bin/bash 

# script to count the intron and exon reads  
# map the dmso sample
htseq-count -m intersection-strict -f bam -t exon -s no -i Parent data/bam/s_dmso.Aligned.sortedByCoord.out.bam data/AaloF1.exons.gff > data/counts/s_dmso.exons.counts.txt
htseq-count -m union -f bam -t intron -s no -i Parent data/bam/s_dmso.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/counts/s_dmso.introns.counts.txt

# map the drb sample
htseq-count -m intersection-strict -f bam -t exon -s no -i Parent data/bam/s_drb_1.Aligned.sortedByCoord.out.bam data/AaloF1.exons.gff > data/counts/s_drb_1.exons.counts.txt
htseq-count -m union -f bam -t intron -s no -i Parent data/bam/s_drb_1.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/counts/s_drb_1.introns.counts.txt

# map the flavopiridol sample 
htseq-count -m intersection-strict -f bam -t exon -s no -i Parent data/bam/s_flavopiridol_2.Aligned.sortedByCoord.out.bam data/AaloF1.exons.gff > data/counts/s_flavopiridol_2.exons.counts.txt
htseq-count -m union -f bam -t intron -s no -i Parent data/bam/s_flavopiridol_2.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/counts/s_flavopiridol_2.introns.counts.txt

# map the triptolide sample 
htseq-count -m intersection-strict -f bam -t exon -s no -i Parent data/bam/s_triptolide_3.Aligned.sortedByCoord.out.bam data/AaloF1.exons.gff > data/counts/s_triptolide_3.exons.counts.txt
htseq-count -m union -f bam -t intron -s no -i Parent data/bam/s_triptolide_3.Aligned.sortedByCoord.out.bam data/AaloF1.introns.gff > data/counts/s_triptolide_3.introns.counts.txt

