# 4-AnalyzeSingleCellSeq

Here I analyze the single-cell sequencing data produced in the paper: *Single-cell transcriptional dynamics of flavivirus infection* by Zanini, et al. The paper can be found [here](https://elifesciences.org/articles/32942#s4).

In this folder I download the counts produced in the paper and perform edgeR differential expression analysis. I then visualize the results and determine whether the codon composition of the Dengue virus strain 16681 influences the expression of the endogenous genes during the infection.

## Scripts

+ **./download_data.sh**: Bash script to download all of the counts, metadata, and other files that are needed for the analysis

+ **./calculate_cpms.R**: R script to calculate the Counts Per Million from the raw counts files

+ **./run_all.sh**: Bash script to create the required folders and run the two above scripts

## RMarkdown

+ **./differential_expression_dengue.Rmd**: RMarkdown file that analyzes the single-cell sequencing data and produces a number of visualizations

## HTML

+ **./differential_expression_dengue.html**: HTML file that is produced when the RMarkdown file is knitted

## Data

+ **./data/counts/counts_dengue.tsv**: raw counts of endogenous human genes during Dengue infection

+ **./data/counts/counts_zika.tsv**: raw counts of endogenous human genes during Zika infection

+ **./data/cpms/cpm_dengue.csv**: Counts Per Million (CPMs) of endogenous human genes during Dengue infection

+ **./data/cpms/cpm_zika.csv**: Counts Per Million (CPMs) of endogenous human genes during Zika infection

+ **./data/edgeR/uninfected-low.csv**: table of differential expression of human endogenous genes between uninfected and low infected cells (produced with edgeR)

+ **./data/edgeR/uninfected-med.csv**: table of differential expression of human endogenous genes between uninfected and medium infected cells (produced with edgeR)

+ **./data/edgeR/uninfected-high.csv**: table of differential expression of human endogenous genes between uninfected and high infected cells (produced with edgeR)

+ **./data/fastas/dengue_2_16681.fa**: fasta file of the CDS of Dengue 2 strain 16681 (the strain used in the single-cell sequencing paper)

+ **./data/fastas/zika_PRVABC59.fa**: fasta file of the CDS of Zika strain PRVABC59 (the strain used in the single-cell sequencing paper)

+ **./data/genes/downregulated_genes_dengue.txt**: file with a list of endogenous human genes that are significantly downregulated during Dengue infection

+ **./data/genes/no_change_genes_dengue.txt**: file with a list of endogenous human genes that are not significantly downregulated or upregulated during Dengue infection

+ **./data/genes/upregulated_genes_dengue.txt**: file with a list of endogenous human genes that are significantly upregulated during Dengue infection

+ **./data/gtfs/hg38.Ens_102.gtf**: GTF file of the human genome (hg38 - Ensembl version 102)

+ **./data/metadata/cell_metadata_dengue.tsv**: Dengue metadata file from the single-cell paper showing the cell names, Dengue counts, MOI, time, etc. 

+ **./data/metadata/cell_metadata_zika.tsv**: Zika metadata file from the single-cell paper showing the cell names, Zika counts, MOI, time, etc.
