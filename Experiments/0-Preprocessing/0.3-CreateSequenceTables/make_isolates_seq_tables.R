################################################################################
# Load libraries 

library(dplyr)
library(seqinr)
library(Biostrings)
library(stringr)

################################################################################
# Load Dengue isolates fasta

isolates = read.fasta(
  '../0.2-DengueIsolatesCDSFastas/data/dengue_isolates.fa', 
  as.string = TRUE
)

################################################################################
# Extract annotations and sequence 

annotations = sapply(isolates, function(x) attr(x, 'Annot'))
seq = as.character(isolates)

################################################################################
# Create function for info extraction 

extract_organism = function(annot) {
  str_extract(annot, 'Organism:[\\w\\s]*') %>%
    str_replace('Organism:', '')
}

extract_protein_name = function(annot) {
  str_extract(annot, 'Protein Name:[\\w\\s]*') %>%
    str_replace('Protein Name:', '')
}

################################################################################
# Create table of annotation and seq

df = data.frame(
  species = extract_organism(annotations), 
  gene_ID = extract_protein_name(annotations), 
  coding = seq
)

################################################################################
# Keep only polyprotein sequences

df = df %>%
  filter(gene_ID == 'polyprotein')

################################################################################
# Convert sequences to uppercase

df$coding = str_to_upper(df$coding)

################################################################################
# Remove sequences that contain other letters besides A, C, T, G

iv = grep('[^ACGT]', df$coding)
if (length(iv) > 0) {
  df = df[-iv, ] 
}

################################################################################
# Collect different serotypes in different data.frames

dengue_1 = df %>%
  filter(species == 'Dengue virus 1')
dengue_2 = df %>%
  filter(species == 'Dengue virus 2')
dengue_3 = df %>%
  filter(species == 'Dengue virus 3')
dengue_4 = df %>%
  filter(species == 'Dengue virus 4')

################################################################################
# Get most common length for each serotype and remove sequences that are 
# longer/shorter 

get_mode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

nchar_dengue_1 = get_mode(nchar(dengue_1$coding))
dengue_1 = dengue_1 %>%
  filter(nchar(coding) == nchar_dengue_1)

nchar_dengue_2 = get_mode(nchar(dengue_2$coding))
dengue_2 = dengue_2 %>%
  filter(nchar(coding) == nchar_dengue_2)

nchar_dengue_3 = get_mode(nchar(dengue_3$coding)) 
dengue_3 = dengue_3 %>%
  filter(nchar(coding) == nchar_dengue_3)

nchar_dengue_4 = get_mode(nchar(dengue_4$coding))
dengue_4 = dengue_4 %>%
  filter(nchar(coding) == nchar_dengue_4)

################################################################################
# Write tables

write.table(
  dengue_1,
  './data/dengue_1_isolates_seq.csv', 
  sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE
)

write.table(
  dengue_2,
  './data/dengue_2_isolates_seq.csv', 
  sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE
)

write.table(
  dengue_3,
  './data/dengue_3_isolates_seq.csv', 
  sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE
)

write.table(
  dengue_4,
  './data/dengue_4_isolates_seq.csv', 
  sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE
)

################################################################################