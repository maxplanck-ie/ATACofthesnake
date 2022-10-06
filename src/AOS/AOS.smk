import os

include: config['rules']['peaks']

rule all:
  input:
    expand(os.path.join(
        config['dirs']['bamdir'],
        '{sample}.bam.bai'
    ), sample=config['samples'])
