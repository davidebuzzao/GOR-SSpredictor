#!/usr/bin/env Nextflow

def sample = file('db/uniprot_sprot.fasta')
println sample.countFasta()