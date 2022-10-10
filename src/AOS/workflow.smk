include: config['rules']['peaks']

localrules: lnBams
rule all:
  input:
    expand(
      "sieve/{sample}.bam.bai",
      sample=config['samples']
    ),
    'figures/mitofraction.png'
