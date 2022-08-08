# Experiments

In this folder I have multiple subdirectories that include the analysis, data, and
plots associated with each key experiment performed in the creation of this paper.

## Subdirectories

+ **./0-Preprocessing/**: all of the steps taken before performing any of the following
experiments, such as the gathering of data, creation of sequence tables, and counting of
codons

+ **./1-CalculateMosquitoOptimality/**: calculation of codon stability coefficient (CSC)
in *Aedes albopictus* (mosquito) C6/36 cells

+ **./2-ConcatenateOptimalities/**: combine all of the known codon optimality codes into
one convenient file

+ **./3-DownloadAllViruses/**: download all known viruses from the RefSeq database and
analyze their codon composition, relative synonymous codon usage (RSCU), and relative
RSCU to human and mosquito

+ **./4-AnalyzeSingleCellSeq/**: analyze the single cell sequencing data produced in the
paper: *Single-cell transcriptional dynamics of flavivirus infection* by Zanini, et al., in
which human Huh7 cells are infected with a strain of Dengue virus 2

+ **./5-MutationAnalysis/**: analyze the mutation data produced in the paper: *Principles of
dengue virus evolvability derived from genotype-fitness maps in human and mosquito cells* by
Dolan, et al., in which a strain of Dengue virus 2 is passaged in human and mosquito cells and
mutational fitness is calculated

+ **./6-OptimalityReporters/**: contains files with the curated sequences for our reporters
used to verify *Aedes albopictus* optimality in C6/36 cells 
