#!/bin/bash -ue
for id in $(cut -d ' ' -f 2 "BlindSet.clust.sort" | awk -F ":" '{print $1}'); do grep -A 1 $id "BlindSet.fasta"; done > "BlindSet.clust.sort.fasta"
