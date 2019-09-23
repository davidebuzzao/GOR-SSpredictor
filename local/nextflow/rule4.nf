#!/usr/bin/env Nextflow

genomes = Channel.fromPath('BlindingSet/FASTA/*.fasta')
myFile = file('identifiers.txt')
allines  = myFile.readLines()

process simpleCount {
  input:
  val x from allines
  file "${x}.fasta" from genomes
  
  output:
  stdout result
  
  """
  cat ${x}.fasta | grep '>'
  """
}

result.println {it.trim()}