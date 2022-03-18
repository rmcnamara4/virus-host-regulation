################################################################################
# Get virus name from file name 

get_virus_name_from_path = function(path) {
  
  basename(path) %>%
    str_replace('.fna.gz', '') %>%
    str_replace_all('_', ' ') %>%
    str_replace_all(';', '/')
  
}

################################################################################
# Get the locus tag 

get_locus_tag = function(id) {
  
  str_extract(id, 'locus_tag=\\w+') %>%
    str_replace('locus_tag=', '')
  
}

clean_locus_tag = function(tag) {
  
  substr(tag[1], start = 1, stop = lcprefix(tag[1], tag[length(tag)])) %>%
    str_replace('_.*', '') %>%
    str_replace('s$', '')
  
}

################################################################################
# Get protein id

get_protein_id = function(id) {
  
  str_extract(id, 'cds\\S+') %>%
    str_replace('cds_', '') %>%
    str_replace('\\.\\d', '') 
  
}

################################################################################
# Get accession number

get_accession_num = function(id) {
  
  str_extract(id, 'NC_\\d+[[:punct:]]\\d+') %>%
    str_replace('\\.\\d', '')
  
}

################################################################################
# Calculate codon composition 

get_codon_counts = function(coding) {
  
  seq = DNAStringSet(coding) 
  codons = trinucleotideFrequency(seq, step = 3) %>%
    as.data.frame()
  
  return(codons)
  
}
