################################################################################
# Load required libraries 

library(seqinr)
library(Biostrings)

################################################################################
# Load data 

# load Dengue strain CDS 
dengue_strain_char = as.character(read.fasta('../4-AnalyzeSingleCellSeq/data/fastas/dengue_2_16681.fa', as.string = TRUE))

# load human sequences 
human_char = read.table('../0-Preprocessing/0.3-CreateSequenceTables/data/human_hg38_seq.csv', 
                        sep = ',', header = TRUE)$coding

# load mosquito sequences 
mosquito_char = read.table('../0-Preprocessing/0.3-CreateSequenceTables/data/mosquito_AaloF1_seq.csv', 
                           sep = ',', header = TRUE)$coding

################################################################################
# Make DNAString objects for each species 

dengue_strain_seq = DNAString(dengue_strain_char)
human_seq = DNAStringSet(human_char)
mosquito_seq = DNAStringSet(mosquito_char)

################################################################################
# Calculate codon usage for each species 

dengue_strain_codons = trinucleotideFrequency(dengue_strain_seq, step = 3)
human_codons = trinucleotideFrequency(human_seq, step = 3)
mosquito_codons = trinucleotideFrequency(mosquito_seq, step = 3)

################################################################################
# Get the total codon usage for the genomes 

dengue_strain_codons_total = dengue_strain_codons # only one mRNA so the total is the same 
human_codons_total = apply(human_codons, 2, sum)
mosquito_codons_total = apply(mosquito_codons, 2, sum)

################################################################################
# Calculate RSCU for each species 

dengue_strain_rscu = calc_rscu(dengue_strain_codons_total)
human_rscu = calc_rscu(human_codons_total)
mosquito_rscu = calc_rscu(mosquito_codons_total)

################################################################################
# Caclulate RSCU fold change 

human_rscu_fc = log2(dengue_strain_rscu / human_rscu)
mosquito_rscu_fc = log2(dengue_strain_rscu / mosquito_rscu)

################################################################################
# Create df 

rscu_df = data.frame(
  codon = names(human_rscu_fc), 
  aa = codon_to_aa(names(human_rscu_fc)), 
  human_rscu = as.numeric(human_rscu), 
  mosquito_rscu = as.numeric(mosquito_rscu), 
  dengue_16681_rscu = as.numeric(dengue_strain_rscu), 
  human_rscu_fc = as.numeric(human_rscu_fc), 
  mosquito_rscu_fc = as.numeric(mosquito_rscu_fc)
)

################################################################################
# Write file for future use 

write.table(
  rscu_df, 
  './data/rscu_fc.csv', 
  sep = ',', 
  col.names = TRUE, 
  row.names = FALSE, 
  quote = FALSE
)

################################################################################