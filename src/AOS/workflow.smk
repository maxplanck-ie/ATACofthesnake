include: config['rules']['peaks']

localrules: lnBams
rule all:
  input:
    expand(
      "sieve/{sample}.bam",
      sample=config['samples']
    ),
    'figures/mitofraction.png'
