#!/bin/bash -ue
tail -n +2 "BlindSet.csv" | head -n +5549 | sed 's/"//g' | awk -F "," '{print ">"$1"_"$2"\n"$4}' > "BlindSet.fasta"
