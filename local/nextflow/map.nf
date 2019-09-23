#!/usr/bin/env Nextflow

tuple = Channel.from( [1, 'alpha.txt'], [2, 'beta.txt'], [3, 'delta.txt'] )

process setExample {
    input:
    set val(x), file('latin.txt')  from tuple

    output:
    stdout result

    """
    echo Processing $x
    cat - latin.txt > copy.txt
    """

}

result.println {it}

genomes = Channel.fromPath('fasta/*.fasta')

process foo {
  input:
  file fasta from genomes

  output:
  val x into var_channel
  val '5GN2_A' into str_channel
  val "${fasta.baseName}.out" into exp_channel

  script:
  x = fasta.name
  """
  cat $x > file
  """
}
