from AOS.helper import qcplotter

rule lnBams:
  input: 
    os.path.join(config['dirs']['bamdir'], "{sample}.bam")
  output: 
    'input/{sample}.bam'
  threads: 1
  shell:'''
  ln -t input -s {input}
  '''

rule bamIx:
  input:
    'input/{sample}.bam'
  output:
    'input/{sample}.bam.bai'
  conda: config['envs']['seqtools']
  threads: 5
  shell:'''
  samtools index -@ {threads} {input}
  '''

rule alSieve:
  input:
    bam = 'input/{sample}.bam',
    bai = 'input/{sample}.bam.bai'
  output:
    shortb = 'sieve/{sample}.bam',
    qc = 'qc/{sample}_sieve.txt',
  conda: config['envs']['seqtools']
  threads: 15
  params:
    rar = '--blackListFileName {}'.format(config['files']['readattractingregions']),
    size = '--maxFragmentLength {} --minFragmentLength 0'.format(config['vars']['fragsize'])
  shell:'''
  alignmentSieve --bam {input.bam} --outFile {output.shortb} --filterMetrics {output.qc} -p {threads} {params.rar} {params.size}
  samtools index -@ {threads} {output.shortb}
  '''

rule ixStat:
  input:
    bam = 'input/{sample}.bam',
    bai = 'input/{sample}.bam.bai'
  output:
    ixstat = 'qc/{sample}_ix.tsv'
  conda: config['envs']['seqtools']
  threads: 1
  shell:'''
  set +o pipefail;
  samtools idxstats {input.bam} | cut -f1,3 > {output.ixstat}
  '''

rule mitoPlot:
  input:
    expand('qc/{sample}_ix.tsv', sample=config['samples'])
  output:
    png = 'figures/mitofraction.png'
  params:
    qcdir= 'qc'
  run:
    qcplotter(params.qcdir, output.png)