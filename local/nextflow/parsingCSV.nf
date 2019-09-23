#!/usr/bin/env Nextflow

myDir = file('../../results/BlindSet/ok.txt')
blindSet_path = myDir.mkdirs()
println blindSet_path ? "$myDir has been successfully created" : "Cannot create directory: $myDir"

params.coverage = 0.0
params.evalue = 30
blindSet_path = '../../results/BlindSet/'

params.csv = '../../dataset/BlindSet/BlindSet.csv'

csv_ch = Channel.fromPath(params.csv)
csv_ch.into { res_ch; fasta_ch }

process parsingCSV_res {
    publishDir "$blindSet_path"
    
    input:
        file raw_csv from res_ch
    
    output:
        file "${raw_csv.baseName}.res" into clustSort1_ch
    shell:
    '''
    sed 's/"//g' !{raw_csv} | awk -F "," '{print $1"_"$2,$3}' | tail -n +2 | head -n +5549 > "!{raw_csv.baseName}.res"
    '''
}

process parsingCSV_fasta {
    publishDir "$blindSet_path"
    
    input:
        file raw_csv from fasta_ch
    output:
        file "${raw_csv.baseName}.fasta" into clustFasta_ch
    shell:
        '''
        tail -n +2 "!{raw_csv}" | head -n +5549 | sed 's/"//g' | awk -F "," '{print ">"$1"_"$2"\\n"$4}' > "!{raw_csv.baseName}.fasta"
        '''
}

process blastclust_c0_id30 {
    publishDir "$blindSet_path"

    input:
        file blindSet_fasta from clustFasta_ch
    output:
        file "${blindSet_fasta.baseName}.clust" into clustSort2_ch
    shell:
        '''
        blastclust -i "!{blindSet_fasta}" -o "!{blindSet_fasta.baseName}.clust" -L "!{params.coverage}" -S "!{params.evalue}"
        '''
}

process sortClusters {
    publishDir "$blindSet_path"

    input:
        file blindSet_res from clustSort1_ch
        file blindSet_clust from clustSort2_ch
    output:
        file "${blindSet_res.baseName}.sort" into clustSort3_ch
    
    '''
    Sort_cluster "!{blindSet_res}" "!{blindSet_clust}" > "!{blindSet_res}.sort"
    '''
}