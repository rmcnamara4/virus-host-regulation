# 1-CalculateMosquitoOptimality

Here I visualize the raw reads from the introns, exons, and ercc
for each of the drugs used to block transcription (DRB, Flavopiridol, Triptolide)
and DMSO.
The intron TPMs are calculated and visualized with the exon TPMs.
Lastly, the codon stability coefficients (CSC) are calculated for C6/36 cells (*Aedes albopictus*)
using each of the samples.

## Scripts

+ **./get_files.sh**: downloads the required files to complete the analysis (i.e. raw counts, tpms, gff, bams)

+ **./split_gff.sh**: splits the *Aedes albopictus* gff into introns and exons

+ **./map_introns.sh**: maps the reads to intronic regions and produces a counts file (note: uses htseq-count)

+ **./run_all.sh**: runs all of the bash scripts in order to complete the preprocessing

## RMarkdown

+ **./calc_viz_intron_tpm.Rmd**: calculates TPMs from the raw intron counts; visualizes the raw reads
from the exons, introns, and ercc; visualizes the distribution of the TPMs of the exons and introns

+ **./calc_mosquito_csc**: calculates the CSC for each of the drug samples

## HTML

+ **./calc_viz_intron_tpm.html**: output of the calc_viz_intron_tpm.Rmd file

+ **./calc_mosquito_csc**: output of the calc_mosquito_csc.Rmd file

## Data

### GFF

+ **./data/AaloF1.EnsGen_50.gff**: Ensembl GFF file of *Aedes albopictus*

+ **./data/AaloF1.exons.gff**: Ensembl GFF file of *Aedes albopictus* filtered for exons

+ **./data/AaloF1.introns.gff**: Ensembl GFF file of *Aedes albopictus* filtered for introns

### CSC

+ **./data/mosquito_csc.csv**: CSC calculated for each drug sample and the average

### BAM

+ **./data/bam/s_dmso.Aligned.sortedByCoord.out.bam/bai**: alignment bam and bai file for DMSO sample

+ **./data/bam/s_drb_1.Aligned.sortedByCoord.out.bam/bai**: alignment bam and bai file for DRB sample

+ **./data/bam/s_flavopiridol_2.Aligned.sortedByCoord.out.bam/bai**: alignment bam and bai file for Flavopiridol sample

+ **./data/bam/s_triptolide_3.Aligned.sortedByCoord.out.bam/bai**: alignment bam and bai file for Triptolide sample

### Counts

+ **./data/counts/s_dmso.introns.counts.txt**: raw intron counts of DMSO sample

+ **./data/counts/s_drb_1.introns.counts.txt**: raw intron counts of DRB sample

+ **./data/counts/s_flavopiridol_2.introns.counts.txt**: raw intron counts of Flavopiridol sample

+ **./data/counts/s_triptolide_3.introns.counts.txt**: raw intron counts of Triptolide sample

+ **./data/counts/ercc_counts.csv**: raw ercc counts for all samples

+ **./data/counts/exon_counts.csv**: raw exon counts for all samples

+ **./data/counts/intron_counts.csv**: raw intron counts for all samples

### TPMs

+ **./data/tpms/intron_tpms.csv**: intron TPMs for all samples

+ **./data/tpms/exon_tpms.csv**: exon TPMs for all samples

## Figures

+ **./figures/raw_exons-1.(png|pdf)**: boxplot of raw exon reads from all samples

+ **./figures/raw_introns-1.(png|pdf)**: boxplot of raw intron reads from all samples

+ **./figures/raw_ercc-1.(png|pdf)**: boxplot of raw ercc reads from all samples

+ **./figures/tpms-1.(png|pdf)**: boxplot of intron and exon TPMs from all samples

+ **./figures/pairs_plot-1.(png|pdf)**: pairs plot of CSC calculated from all samples, showing
correlation and density
---
To reproduce the analysis:

```bash
bash run_all.sh
```
The RMarkdown files need to be knitted to reproduce their calculations. However,
all results are already produced. 
