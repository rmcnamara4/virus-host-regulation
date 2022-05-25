################################################################################
# Load libraries

library(dplyr)
library(Biostrings)
library(seqinr)
library(stringr)

################################################################################
# Source helper functions

source('../../Src/create_codon_counts_ind_functions.R')

################################################################################
# Load in fasta file

file = snakemake@input[['file']]

fa = read.fasta(
  file = file,
  as.string = TRUE
)
names(fa) = NULL

################################################################################
# Load in host metadata

host_metadata = read.table(
  './host_metadata.txt',
  sep = '\t',
  header = TRUE,
  quote = "\""
)

################################################################################
# Load in virus composition file

virus_comp = read.table(
    './virus_composition.txt',
    sep = '\t',
    header = TRUE,
    quote = "\""
)

################################################################################
# Create tibble of annotation and coding sequence

annotations = sapply(fa, function(x) {
  attr(x, "Annot")
})

viralseq_df = tibble(
  seq_id = annotations,
  coding = as.character(fa)
)

################################################################################
# Use functions to get required columns

# Get accession number
acc_num = get_accession_num(viralseq_df$seq_id)

# Get protein id
protein_id = get_protein_id(viralseq_df$seq_id)

# Get locus tag
locus_tag = get_locus_tag(viralseq_df$seq_id)

# Get virus name
virus_name = get_virus_name_from_path(file)

################################################################################
# Add data to the table

if (virus_name %in% host_metadata$virus) {

  viralseq_df =
    viralseq_df %>%
    mutate(
      virus = virus_name,
      protein_id = protein_id,
      locus_tag = locus_tag,
      accession = acc_num
    ) %>%
    merge(
      host_metadata,
      by = 'virus'
    ) 

} else {

  viralseq_df =
    viralseq_df %>%
    mutate(
      virus = virus_name,
      protein_id = protein_id,
      locus_tag = locus_tag,
      accession = acc_num,
      host = NA,
      realm = NA
    )

}

################################################################################
# Combine with virus composition data 

if (virus_name %in% virus_comp$virus) {
  
  viralseq_df = viralseq_df %>%
    merge(
      virus_comp, 
      by = 'virus'
    ) %>%
    mutate_all(
      trimws
    ) %>%
    select(virus, realm, locus_tag, host, type, shape, accession, protein_id, coding)
  
} else {
  
  viralseq_df = viralseq_df %>%
    mutate(
      type = NA, 
      shape = NA
    ) %>%
    mutate_all(
      trimws
    ) %>%
    select(virus, realm, locus_tag, host, type, shape, accession, protein_id, coding)
  
}

################################################################################
# Calculate codon composition

codon_counts = get_codon_counts(viralseq_df$coding)

################################################################################
# Column bind with existing table

viralseq_df = cbind(viralseq_df, codon_counts)

################################################################################
# Get total counts and create total table

col_nums = grep('^[ACGT]', names(viralseq_df))

if (any(is.na(viralseq_df$locus_tag))) {

  viralseq_total_df = viralseq_df %>%
    select(virus, realm, host, type, shape, all_of(col_nums)) %>%
    group_by(virus, realm, host, type, shape) %>%
    summarize_all(
      sum
    ) %>%
    mutate(
      locus_tag = viralseq_df$locus_tag[1],
      .after = 2
    )

} else {

  viralseq_total_df = viralseq_df %>%
    select(virus, realm, host, type, shape, all_of(col_nums)) %>%
    group_by(virus, realm, host, type, shape) %>%
    summarize_all(
      sum
    ) %>%
    mutate(
      locus_tag = clean_locus_tag(viralseq_df$locus_tag)[1],
      .after = 2
    )

}

################################################################################
# Write tables

write.table(
  viralseq_df,
  snakemake@output[['codon_counts']],
  sep = ',',
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.table(
  viralseq_total_df,
  snakemake@output[['codon_counts_total']],
  sep = ',',
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

################################################################################
