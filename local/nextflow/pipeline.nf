#!/usr/bin/env Nextflow

myDir = file('../../results/BlindSet/')
blindSet_path = myDir.mkdirs()
println blindSet_path ? "$myDir has been successfully created" : "Cannot create directory: $myDir"

params.coverage = 0.0
params.evalue = 30
blindSet_path = '../../results/BlindSet/'

params.csv = '../../dataset/BlindSet/BlindSet.csv'

csv_ch = Channel.fromPath(params.csv)
csv_ch.into { res_ch; fasta_ch }

process rawBlindSet_IDRESextraction {
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

process rawBlindSet_FASTAextraction {
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

clustFasta_ch.into { clustFasta1_ch; clustFasta2_ch }

process blastclust_c0_id30 {
    publishDir "$blindSet_path"

    input:
        file blindSet_fasta from clustFasta1_ch
    output:
        file "${blindSet_fasta.baseName}.clust" into clustSort2_ch

        """
        blastclust -i "$blindSet_fasta" -o "${blindSet_fasta.baseName}.clust" -L "$params.coverage" -S "$params.evalue"
        """
}

process sortClusters {
    publishDir "$blindSet_path"

    input:
        file blindSet_res from clustSort1_ch
        file blindSet_clust from clustSort2_ch
    output:
        file "${blindSet_clust}.sort" into clustSort3_ch
    
    """
    Sort_cluster "$blindSet_res" "$blindSet_clust" > "${blindSet_clust}.sort"
    """
}

//println 'This is what comes out from sortClusters process'
//clustSort3_ch.view()

process sortedClust_FASTAextraction {
    publishDir "$blindSet_path"

    input:
        file blindSet_sort from clustSort3_ch
        file blindSet_fasta from clustFasta2_ch

    output:
        file "${blindSet_sort}.fasta" into clustSort4_ch
    shell:
    '''
    for id in $(cut -d ' ' -f 2 "!{blindSet_sort}" | awk -F ":" '{print $1}'); do grep -A 1 $id "!{blindSet_fasta}"; done > "!{blindSet_sort}.fasta"
    '''
}

//println 'This is what comes out from sortedClust_FASTAextraction process'
//clustSort4_ch.view()