#!/bin/bash -ue
sed 's/"//g' BlindSet.csv | awk -F "," '{print $1"_"$2,$3}' | tail -n +2 | head -n +5549 > "BlindSet.res"
