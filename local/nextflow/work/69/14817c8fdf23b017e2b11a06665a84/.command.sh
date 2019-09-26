#!/bin/bash -ue
blastp -query "BlindSet.clust.sort.fasta" -db "../../dataset/db/jPred4.fasta" -out "BlindSet.clust.sort.fasta.blast" -evalue "0.01" -outfmt 6 -num_threads 2
