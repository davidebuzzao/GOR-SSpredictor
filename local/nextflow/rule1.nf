#!/usr/bin/env Nextflow

myDir = file('fasta/prova1/prova2')
result = myDir.mkdirs()
println result ? "OK" : "Cannot create directory: $myDir"

myFile = file('../../local/nextflow/rule3.nf')
myFile.mklink('fasta/prova1/prova2/elfunzia')

params.str = 'prova.txt'

var = 'Hello World me Im Davide Buzzao'


process splitLetters {
    afterScript 'rm *.log*'
    publishDir 'fasta/prova1'
    
    output: 
    file 'chunk_*' into letters mode flatten

    """
    printf '$var' | split -b 6 - chunk_
    """
}


process convertToUpper {
    afterScript 'rm *.log*'

    input:
    file x from letters
    
    output:
    stdout results

    """
    cat $x | tr '[a-z]' '[A-Z]' 
    """
}

results.println {it.trim()}