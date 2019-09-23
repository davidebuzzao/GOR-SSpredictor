#!/usr/bin/env Nextflow

scripts_ch = Channel
                    .fromPath('../src/*.py')
                    .flatten()
                    
Channel
   .fromPath('../src/*.py', type: 'dir')
   .subscribe { println "Script: $it.baseName" }
/*s
process make_links {

    input:
        file script from scripts_ch
    
    script:
        //println "$script"
        "$script".mklink("../bin/$script")
}
*/