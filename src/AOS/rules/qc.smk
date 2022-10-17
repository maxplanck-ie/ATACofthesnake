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

rule fragsize:
  input:
    bams = expand('input/{sample}.bam', sample=config['samples']),
    bais = expand('input/{sample}.bam.bai', sample=config['samples'])
  output:
    'qc/fragsize.tsv'
  threads: 10
  conda: config['envs']['seqtools']
  shell:'''
  bamPEFragmentSize -b {input.bams} \
    -p {threads} \
    --outRawFragmentLengths {output}
  '''

rule scalefactors:
  input:
    'peakset/counts.tsv'
  output:
    'peakset/scalefactors.txt'
  params:
    config['rscripts']['scalefactors']
  threads: 1
  conda: config['envs']['seqtools']
  shell:'''
  Rscript {params} {input} {output}
  '''

rule bigwigs:
  input:
    scalefactors = 'peakset/scalefactors.txt',
    bam = 'sieve/{sample}.bam',
    bai = 'sieve/{sample}.bam.bai'
  output:
    bigwig = 'bw/{sample}.scalefac.bw'
  params:
    rar = config['files']['readattractingregions'],
    sample = "{sample}"
  threads: 10
  conda: config['envs']['seqtools']
  shell:'''
  SCALEFAC=$(grep {params.sample} {input.scalefactors} | cut -f2 -d ' ')
  bamCoverage --scaleFactor $SCALEFAC \
    -b {input.bam} -o {output.bigwig} \
    -p {threads} -bs 1 -bl {params.rar}
  '''

rule bigwigs_rpkm:
  input:
    bam = 'sieve/{sample}.bam',
    bai = 'sieve/{sample}.bam.bai'
  output:
    bigwig = 'bw/{sample}.RPKM.bw'
  params:
    rar = config['files']['readattractingregions'],
    sample = "{sample}"
  threads: 10
  conda: config['envs']['seqtools']
  shell:'''
  bamCoverage --normalizeUsing RPKM \
    -b {input.bam} -o {output.bigwig} \
    -p {threads} -bs 1 -bl {params.rar}
  '''