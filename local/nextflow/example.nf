#!/usr/bin/env Nextflow

params.db = './db/uniprot_sprot.fasta'
database = file(params.db)
params.evalue = 0.01


blindSet = './BlindSet'
trainingSet = './TrainingSet'
psiblast = "$trainingSet/PSIBLAST2"
psiblastPssm = "$trainingSet/PSIBLAST2/pssm2"
psiblastAlign = "$trainingSet/PSIBLAST2/alignments2"
psiblastFasta = "$trainingSet/fasta"
dataset = Channel
                .fromPath("$psiblast/fasta")
                .map { file -> tuple(file.baseName, file)}
/****************************
Training Set data preparation
*////////////////////////////
/*
process psiBlast {
    publishDir "$psiblastAlign"

    input:
    set datasetID, file(datasetFile) from dataset
    
    output:
    set datasetID, file("${datasetID}.blast") into aligned_files
    
    script:
    """
    blastp -query ${datasetFile} -db "$database" -out ${datasetID}.blast -evalue "$params.evalue" -outfmt 6 -num_threads 2
    """
}

    script:
        myReader = jPred4.newReader()
        String line
        while( line = myReader.readLine() ) {
            println line 
            }
        myReader.close()
*/

process clustalw2_align {
    publishDir "$psiblastAlign"

    input:
    set datasetID, file(datasetFile) from dataset

    output:
    set datasetID, file("${datasetID}.aln") into aligned_files

    script:
    """
    clustalw2 -INFILE=${datasetFile} -OUTFILE=${datasetID}.aln
    """
}