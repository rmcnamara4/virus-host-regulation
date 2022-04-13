codon_to_aa = function(codons) {

    aa = vector()

    for (c in codons) {

        aa = append(aa, names(synonymous_codons)[grep(c, synonymous_codons)])

    }

    return(aa)
}
